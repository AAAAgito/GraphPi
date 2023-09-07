#pragma once
#include "schedule.h"
#include "vertex_set.h"
#include <map>
#include <utility>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <string>
#include <omp.h>
#include <memory>
#include <unordered_map>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

enum PartitionType {
    METIS,
    RANDOM,
    NAIVE_BFS,
    LDG,
    PAGE
};

enum BlockType {
    RANDOM_BLOCK,
    K_CORE_BLOCK,
    CHAIN,
    SIMPLE_BFS
};

enum LoadType {
    SUB_BLOCK,
    ALL_BLOCK
};

class Graphmpi;
struct VertexTable {
    bool is_intra = false;
    int file_id = 0;
};

struct block_length {
    int v_len = 0;
    int e_len = 0;
    block_length(int i, int ii) {
        v_len = i;
        e_len = ii;
    }
    block_length() {
        v_len = 0;
        e_len = 0;
    }
};

struct patch_data
{
    std::vector<int> insertion;
    std::vector<int> deletion;
};

struct pair
{
    int e1;
    int e2;
};



struct cell {
    int v;
    cell *prev;
    cell *next;
    cell(int vtx) {
        v= vtx;
        prev = NULL;
        next = NULL;
    }
};

class Priority_List {
public:
    cell head = cell(-1);
    cell tail = cell(-1);
    int size;
    int e_size;

    Priority_List() {
        head = cell(-1);
        tail = cell(-1);
        size = 0;
        head.next = &tail;
        tail.prev = &head;
    }
    bool check_list() {
        cell *ptr = &head;
        for (int i=0; i< size+1;i++) ptr = ptr->next;
        return ptr==&tail;
    }

    void insert(const std::vector<int> &scores,cell * &record, std::map<int, cell*> &head_of_score) {
        // O(1)
//        printf("insert\n");
        size += 1;
        int insert_score = scores[record->v];
        if (head.next == &tail) {
            head.next = record;
            record->prev = &head;
            record->next = &tail;
            tail.prev = record;
            head_of_score[insert_score] = record;
//            printf("done 1\n");
            assert(check_list());
            return;
        }
        if (head_of_score.find(insert_score)!=head_of_score.end()) {
            record->next = head_of_score[insert_score]->next;
            record->prev = head_of_score[insert_score];
            head_of_score[insert_score]->next->prev = record;
            head_of_score[insert_score]->next = record;
        }
        else {
            int nearest_lower_score = -1;
            for (auto j : head_of_score) {
                if (j.first < insert_score) nearest_lower_score = j.first;
                else break;
            }
            head_of_score[insert_score] = record;
            if (nearest_lower_score >=0) {
                record->next = head_of_score[nearest_lower_score];
                record->prev = head_of_score[nearest_lower_score]->prev;
            }
            else {
                record->next = &tail;
                record->prev = tail.prev;
            }
            record->prev->next = record;
            record->next->prev = record;
        }

        assert(check_list());
//        printf("done\n");
    }
    void resort(const std::vector<int> &scores, const std::vector<cell*> cells,const std::vector<int> &change_vertex, const std::vector<int> &prev_value, std::map<int, cell*> &head_of_score, bool change_positive) {
        // TODO: if change positive, then start from head, move largest one first
        for (auto i: prev_value) assert(head_of_score.find(i)!=head_of_score.end());
        if (change_positive) {
            // O(Size)
            for (int i=0; i<change_vertex.size();i++) {
                int change_v = change_vertex[i];
                int insert_score = scores[change_v];
                int previous_score = prev_value[i];
                cell *change_v_cell = cells[change_v];
//                printf("idx %d/%d find %d %d\n",i,change_vertex.size(),previous_score,insert_score);
                if (head_of_score.find(previous_score)!=head_of_score.end()) {
                    assert(head_of_score.find(previous_score) != head_of_score.end());
                    if (head_of_score.at(previous_score) == change_v_cell) {
                        cell *next = change_v_cell->next;
                        while ((next != &tail && scores[next->v] > previous_score) ||
                               std::find(change_vertex.begin(), change_vertex.end(), next->v) != change_vertex.end()) {
                            next = next->next;
                        }
                        if (scores[next->v] == previous_score) {
                            head_of_score[previous_score] = next;
//                            printf("head of %d change to v %d\n", previous_score, next->v);
                        } else {
//                            printf("danger operation %d\n", previous_score);
                            head_of_score.erase(previous_score);
                        }
                    }
                }
                assert(check_list());
                change_v_cell->prev->next = change_v_cell->next;
                change_v_cell->next->prev = change_v_cell->prev;
                if (head_of_score.find(insert_score)!=head_of_score.end()) {
                    cell *cur_head = head_of_score[insert_score];
                    head_of_score[insert_score] = change_v_cell;
//                    printf("do if %d %d\n",cur_head->v, change_v_cell->v);
                    change_v_cell->next = cur_head;
                    change_v_cell->prev = cur_head->prev;
                    change_v_cell->prev->next = change_v_cell;
                    change_v_cell->next->prev = change_v_cell;
                }
                else {
//                    printf("do else\n");
                    int nearest_lower_score = -1;
                    for (auto j : head_of_score) {
                        if (j.first < insert_score) nearest_lower_score = j.first;
                        else break;
                    }
                    head_of_score[insert_score] = change_v_cell;
                    if (nearest_lower_score >= 0) {
                        change_v_cell->next = head_of_score[nearest_lower_score];
                        change_v_cell->prev = head_of_score[nearest_lower_score]->prev;
                    }
                    else {
                        change_v_cell->next = &tail;
                        change_v_cell->prev = tail.prev;
                    }
                    change_v_cell->prev->next = change_v_cell;
                    change_v_cell->next->prev = change_v_cell;
                }
                assert(check_list());
            }
        }
        else {

            for (int i=0; i<change_vertex.size();i++) {
                int change_v = change_vertex[i];
                int insert_score = scores[change_v];
                int previous_score = prev_value[i];
                cell *change_v_cell = cells[change_v];
                if (head_of_score.find(previous_score)!=head_of_score.end()) {
                    assert(head_of_score.find(previous_score) != head_of_score.end());
                    if (head_of_score.at(previous_score) == change_v_cell) {
                        cell *next = change_v_cell->next;
                        while (true) {
                            if (std::find(change_vertex.begin(), change_vertex.end(), next->v) !=
                                change_vertex.end())
                                break;
                            next = next->next;
                        }
                        if (scores[next->v] == previous_score)
                            head_of_score[previous_score] = next;
                        else head_of_score.erase(previous_score);
                    }
                }
                change_v_cell->prev->next = change_v_cell->next;
                change_v_cell->next->prev = change_v_cell->prev;
                if (head_of_score.find(insert_score)!=head_of_score.end()) {
                    change_v_cell->next = head_of_score[insert_score]->next;
                    change_v_cell->prev = head_of_score[insert_score];
                    head_of_score[insert_score]->next->prev = change_v_cell;
                    head_of_score[insert_score]->next = change_v_cell;
                }
                else {
                    int nearest_lower_score = -1;
                    for (auto j : head_of_score) {
                        if (j.first < insert_score) nearest_lower_score = j.first;
                        else break;
                    }
                    head_of_score[insert_score] = change_v_cell;
                    if (nearest_lower_score >= 0) {
                        change_v_cell->next = head_of_score[nearest_lower_score];
                        change_v_cell->prev = head_of_score[nearest_lower_score]->prev;
                    }
                    else {
                        change_v_cell->next = &tail;
                        change_v_cell->prev = tail.prev;
                    }
                    change_v_cell->prev->next = change_v_cell;
                    change_v_cell->next->prev = change_v_cell;
                }
            }
        }
    }

    cell *pop_head(const std::vector<int> &scores, std::map<int, cell*> &head_of_score) {
        cell *to_del = head.next;
        assert(to_del!=&tail);
        if (to_del->next != &tail) {
            if (scores[to_del->next->v] == scores[to_del->v]) {
                head_of_score[scores[to_del->v]] = to_del->next;
            }
            else {
                head_of_score.erase(scores[to_del->v]);
            }
        }
        to_del->next->prev = &head;
        head.next = to_del->next;
        size -=1;
        assert(check_list());
        return to_del;
    }

};

static omp_lock_t priority_lock;

class Graph {
public:
    int v_cnt; // number of vertex
    int g_vcnt;
    int extra_v_cnt;
    int insert_cnt;
    int delete_cnt;
    unsigned int e_cnt; // number of edge
    unsigned int g_ecnt;
    long long tri_cnt{}; // number of triangle
    double max_running_time = 60 * 60 * 24; // second

    int *edge; // edges
    unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]

    int *mem;
    int *mmp_edge;
    unsigned int *mmp_vertex;
    double extern_upper_thresold = 0.267;
    long extern_upper_size = 0;
    BlockType blockType = RANDOM_BLOCK;
    LoadType loadType = SUB_BLOCK;
    PartitionType partitionType = RANDOM;
    std:: string raw_data_path;
    std::string file_path;

    int *bitmap;
    int INT_BITS = sizeof(int);
    int SHIFT = 5;
    int MASK = 0x1f;

    std::unordered_map<int, int> intra_vertex_dict; // map v_id to index of vertex in vertex array
    std::vector<patch_data> patch_v;
    std::vector<VertexTable> v_state_map;
    std::map<int, block_length> block_lengths;
    std::vector<block_length> file_len;
    std::vector<std::string> file_str;

    // used for single load
    std::vector<std::shared_ptr<int[]>> physicals_load;
    std::vector<int> physicals_length;

    std::vector<int> bfs_order_queue;

    std::queue<int> pri_Q;

    std::vector<omp_lock_t> lock_vertex;
    std::vector<omp_lock_t> lock_block;

    int query_num = 0;
    int io_num = 0;
    double file_time = 0.0;
    double lock_time = 0.0;
    double blocking_manage_time = 0.0;
    double load_extern_time = 0.0;
    double priority_list_time = 0.0;

    Graph(std:: string path = "/home/yanglaoyuan/AsyncSubGraphStorage/docker_graphPi/bin/") {
        v_cnt = 0;
        e_cnt = 0;
        g_vcnt = 0;
        edge = nullptr;
        vertex = nullptr;
        raw_data_path = std::move(path);
        
    }

    ~Graph() {
        delete[] edge;
        delete[] vertex;
    }

    int memory_map(const std::string &path) {
        int fd = open(path.c_str(),O_RDWR);
        mem = static_cast<int*>(mmap(NULL, (g_vcnt+g_ecnt)*sizeof(int), PROT_READ, MAP_SHARED ,fd,0));
        mmp_vertex = (unsigned int *)mem+2;
        mmp_edge = mem+2+g_vcnt;
        extra_v_cnt = g_vcnt;
        // bitmap = new int[g_vcnt/INT_BITS];
        return fd;
    }

    void insert_edge(int v1, int v2) {
        patch_v[v1].insertion.push_back(v2);
        patch_v[v2].insertion.push_back(v1);
    }

    void delete_edge(int v1, int v2) {
        patch_v[v1].deletion.push_back(v2);
        patch_v[v2].deletion.push_back(v1);
    }

    void update(std::vector<pair> &u, std::vector<pair> &d);

    void dump_v();

    void refine_graph();
    
    void free_map(int &fd) {
        munmap(mem, (g_vcnt+g_ecnt)*sizeof(int));
        close(fd);
        // delete[] bitmap;
    }

    void set(int *p, int i) {
        p[i >> SHIFT] |= 1 << (i & MASK);
    }
    //获取第i位
    int test(int *p, int i) {
        return p[i >> SHIFT] & (1 << (i & MASK));
    }
    //清除第i位
    int clear(int *p, int i) {
        return p[i >> SHIFT] & ~(1 << (i & MASK));
    }

    int intersection_size(int v1,int v2);
    int intersection_size_clique(int v1,int v2);

    //single thread triangle counting
    long long triangle_counting();

    //multi thread triangle counting
    long long triangle_counting_mt(int thread_count);

    //general pattern matching algorithm with multi thread
    long long pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);
    
    long long pattern_matching_oc(const Schedule& schedule, int thread_count, bool clique = false);

    void to_partition_csr(int num, std::vector<int> &side_vertex, std::vector<int> &part_map, const std::string& path);

    int to_partition_csr(int size);

    void dump_csr_vset(std::vector<int> &vs, std::vector<int> &cores, const std::string &p, const std::string &pmap, const std::string &porder);

    void init_extern_storage(int v_num, int e_num);

    void to_global_csr(const std::string& path);

    void gen_out_of_core_component(int part_num, int block_size, const std::string& path, int k_core = 2);

    void to_block_csr(int block_size, std::vector<int> &side_vertex, std::vector<int> &part_map, int part_num, const std::string& path, int k_core = 2);

    void load_partition_graph(int pid, int num, const std::string& path);

    void load_global_graph(const std::string &path);

    int max_degree{};

    void get_edge_index(int v, unsigned int& l, unsigned int& r) const;
    
    void get_mmp_edge_index(int v, unsigned int& l, unsigned int& r) const;

    void load_physical_edge_index(int &v, std::shared_ptr<int[]> &ptr, int &size);

    void gen_bfs_query_order();

    void gen_bipartite_order(std::vector<int> &part_map, int cur_part, int inter_upper);

    void load_order(bool full = false);

    void load_order(int num,int size);

    void get_physical_edge_index(int &v, std::shared_ptr<int[]> &ptr, int& size) {
        query_num +=1;
        omp_set_lock(&lock_vertex[v]);
        ptr = physicals_load[v];
        omp_unset_lock(&lock_vertex[v]);
        if (ptr == NULL) {
            load_physical_edge_index(v,ptr,size);

            omp_set_lock(&priority_lock);

            while (pri_Q.size() > extern_upper_size) {

                int del_v = pri_Q.front();
                pri_Q.pop();
                omp_set_lock(&lock_vertex[del_v]);
                physicals_load[del_v] = NULL;
                omp_unset_lock(&lock_vertex[del_v]);
            }
            omp_set_lock(&lock_vertex[v]);


            physicals_load[v] = ptr;
            physicals_length[v] = size;
            omp_unset_lock(&lock_vertex[v]);
            pri_Q.push(v);

            omp_unset_lock(&priority_lock);
        }
        else {
            size = physicals_length[v];
        }
    }

private:
    friend Graphmpi;
    void tc_mt(long long * global_ans);

    void pattern_matching_aggressive_func_oc(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);
    
    void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    void blocking_data_manage_k_core(int len, int &block_id, int &local_block_size, std::vector<int> &block_vid, std::vector<int> &to_insert_block, const std::string& path);

    double get_wall_time2() {
        struct timeval time;
        if(gettimeofday(&time,NULL)) {
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
    }
};
