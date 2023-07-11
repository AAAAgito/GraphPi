#pragma once
#include "schedule.h"
#include "vertex_set.h"
#include <map>
#include <utility>
#include <vector>
#include <assert.h>
#include <string>
#include <omp.h>
#include <memory>

enum PartitionType {
    METIS,
    RANDOM,
    NAIVE_BFS,
    LDG
};

enum BlockType {
    RANDOM_BLOCK,
    K_CORE_BLOCK,
    CHINK_BFS,
    SIMPLE_BFS
};

enum LoadType {
    SUB_BLOCK
};

class Graphmpi;
struct VertexTable {
    bool is_intra = false;
    bool is_rooted = false; // is used as a root(depth=0) in subgraph matching
    unsigned int k_hop = 0;
    int loading = 0;
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

struct path{
    int depth = 0;
    VertexSet *subtraction_set = nullptr;
    VertexSet *vertex_set = nullptr;
    int *vertex = nullptr;
    unsigned int v_size = 0;
    bool load_queue = false;
    ~path() {

        if (vertex!= nullptr) delete[] vertex;
        if (subtraction_set!= nullptr) delete subtraction_set;
        if (vertex_set!= nullptr) delete[] vertex_set;
    }
};

struct loader{
    int k_hop;
    int vertex;
    int file_id;
};

struct vertex_degree{
    int v = 0;
    int d = 0;
};

struct cell {
    int v;
    cell *prev;
    cell *next;
};

class Priority_List {
public:
    cell head;
    cell tail;
    int size;

    void init_list() {
        head = cell{-1,NULL,NULL};
        tail = cell{-1,NULL,NULL};
        size = 0;
    }
    void insert_head(cell * &record, int v) {
//        assert(record==NULL);
        if (tail.prev == NULL) {
            tail.prev = record;
            record->next = &tail;
        }
        if (head.next == NULL) {
            head.next = record;
            record->prev = &head;
        }
        else {
            record->next = head.next;
            record->prev = &head;
            head.next->prev = record;
            head.next = record;
        }
        size +=1;
    }
    int delete_record(cell *&record) {
        int v = record->v;
        record->prev->next = record->next;
        record->next->prev = record->prev;
        size -=1;
        return v;
    }
    void move_head(cell *&record) {
        record->prev->next = record->next;
        record->next->prev = record->prev;

        record->next = head.next;
        record->prev = &head;
        head.next->prev = record;
        head.next = record;
    }
};

static omp_lock_t load_list_write_lock;
static omp_lock_t load_list_append_lock;
static omp_lock_t init_load_list_append_lock;
static omp_lock_t extern_ve_loading_lock;
static omp_lock_t ready_bin_lock;
static omp_lock_t candidate_bin_lock;
static omp_lock_t v_state_lock;
static omp_lock_t v_state_khop_lock;


class Graph {
public:
    int v_cnt; // number of vertex
    int g_vcnt;
    unsigned int e_cnt; // number of edge
    long long tri_cnt{}; // number of triangle
    double max_running_time = 60 * 60 * 24; // second
    bool all_loaded = false;

    int *edge; // edges
    unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]
    unsigned int external_space =1; // how much can be load from inter-partition edges
    unsigned int external_used =0;
    float extern_upper_thresold = 0.1;
    float extern_lower_thresold = 0.01;
    int *inter_edge; // store inter partition edges
    unsigned int *inter_vertex;
    int inter_vtx_num = 0;
    unsigned int inter_edge_num = 0;
    int extern_v_max_num = 0;
    int extern_e_max_num = 0;
    BlockType blockType = RANDOM_BLOCK;
    LoadType loadType = SUB_BLOCK;
    PartitionType partitionType = RANDOM;
    std:: string raw_data_path;

    std::map<int, int> intra_vertex_dict; // map v_id to index of vertex in vertex array
    std::map<int, int> inter_vertex_dict; // map to_load v_id to index of to_load vertex in
    std::map<int, VertexTable> v_state_map;
    std::map<int, block_length> block_lengths;

    std::vector<path*> candidate_bin;
    std::vector<path*> ready_bin;
    std::map<int,std::vector<loader>> load_list;
    std::vector<loader> init_load_list;

    std::map<int,std::shared_ptr<int[]>> physicals_load;
    std::vector<cell*> physicals_priority;
    std::map<int,int> physicals_length;

    Priority_List L;

    int io_num = 0;
    double file_time = 0.0;
    double extern_init_time = 0.0;
    double loading_manage_time = 0.0;
    double blocking_manage_time = 0.0;
    double external_time = 0.0;
    double load_extern_time = 0.0;
    int load_num = 0;

    Graph(std:: string path = "/home/yanglaoyuan/AsyncSubGraphStorage/docker_graphPi/bin/") {
        v_cnt = 0;
        e_cnt = 0;
        g_vcnt = 0;
        edge = nullptr;
        vertex = nullptr;
        inter_vertex = nullptr;
        inter_edge = nullptr;
        raw_data_path = std::move(path);
    }

    ~Graph() {
        delete[] edge;
        delete[] vertex;
        delete[] inter_vertex;
        delete[] inter_edge;

    }

    int intersection_size(int v1,int v2);
    int intersection_size_clique(int v1,int v2);

    //single thread triangle counting
    long long triangle_counting();

    //multi thread triangle counting
    long long triangle_counting_mt(int thread_count);

    //resume the match stopped by inter-partition
    void resume_matching(const Schedule& schedule, path *p, long long &local_ans);

    //general pattern matching algorithm with multi thread
    long long pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

    long long pattern_matching_asyn(const Schedule& schedule, int thread_count, bool clique = false);

    //this function will be defined at code generation
//    long long unfold_pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

    //general pattern matching algorithm with multi thread ans multi process
    long long pattern_matching_mpi(const Schedule& schedule, int thread_count, bool clique = false);

    void to_partition_csr(int num, std::vector<int> &side_vertex, std::map<int,int> &part_map, const std::string& path);

    void init_extern_storage(int v_num, int e_num);

    void to_global_csr(const std::string& path);

    void gen_out_of_core_component(int part_num, int block_size, const std::string& path, int k_core = 0);

    void to_block_csr(int block_size, std::vector<int> &side_vertex, std::map<int,int> &part_map, int part_num, const std::string& path, int k_core = 0);

    void load_partition_graph(int pid, int num, const std::string& path);

    void load_global_graph(const std::string &path);

    int max_degree{};
    void load_extern_data(int file_id, const std::vector<loader>& load_vertex, std::vector<int> &loaded_vertex);

    void get_edge_index(int v, unsigned int& l, unsigned int& r) const;

    void load_physical_edge_index(int v, bool next_hop);

    void load_physical_edge_index(std::vector<int> &v);

    void get_physical_edge_index(int v, std::shared_ptr<int[]> &ptr, int& size, bool next_hop) {
        ptr = physicals_load[v];
        cell *pri_p = physicals_priority[v];
        if (ptr == NULL) {
            printf("do if\n");
            assert(pri_p==NULL);
            cell *p=L.tail.prev;
            // also O(1) delete tail's vertex and its adj_list
            if (L.size > 1) assert(L.tail.prev != NULL);

            while (L.size > g_vcnt*extern_upper_thresold && p != NULL) {

                if (L.size > 1) assert(L.tail.prev->prev != NULL);
                if (physicals_load[p->v].use_count() <= 1) {
                    int del_v = L.delete_record(p);
                    p = p->prev;
                    delete physicals_priority[del_v];
                    physicals_priority[del_v] = NULL;
                    physicals_load[del_v] = NULL;
                    physicals_load.erase(del_v);
                }
                else {
                    p= p->prev;
                }
            }
            load_physical_edge_index(v, next_hop);
        }
        // O(1) time, switch/ create cell of vertex to head
        else {
            L.move_head(physicals_priority[v]);
        }

        ptr = physicals_load[v];
        size = physicals_length[v];

        if (L.size > 1) {
            assert(L.tail.prev != NULL);
        }
    }

    void drop_physical_edge_index(int v) {
        if (physicals_load[v] == NULL) {
            return;
        }
        int cnt = physicals_load[v].use_count();
        if (cnt <= 1 && physicals_load.size() >= g_vcnt*extern_upper_thresold) {
            physicals_load[v] = NULL;
            physicals_load.erase(v);
        }
    }


    void get_extern_edge_index(int v, unsigned int& l, unsigned int& r) const;

    bool scatter_loading_thread (int thread_num, int id, int to_load_num) {
//        printf("\n%d %d %d\n",thread_num,to_load_num,extern_v_max_num);
        if (to_load_num < extern_v_max_num*extern_lower_thresold) {
            return false;
        }
//        return id <= thread_num * to_load_num / extern_v_max_num;
        return id == 0;
    }


private:
    friend Graphmpi;
    void tc_mt(long long * global_ans);


    // do pre-process to candidate_bin and generate the best load order (considering batch load)
    void extern_load_init(unsigned int available_space, const Schedule& schedule);


    // manage loading order of vertex, and batch load vertices in same file
    void loading_manage();

    void load_list_append(const loader &l, std::vector<int> &loaded_vertex);

    // tracking memory use of extern vertex, drop all if larger than thresold
    void extern_drop_manage();

    // return true if it is still recommend to store.
    // return false if it is recommended to drop.
    bool extern_store_manage(const Schedule& schedule);

    void pattern_matching_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, long long& local_ans, int depth, bool clique = false);

    void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    void blocking_data_manage_k_core(int len, int &block_id, int &local_block_size, std::vector<int> &block_vid, const std::map<int, int> &revert_dict, std::vector<int> &to_insert_block, const std::string& path);

    //this function will be defined at code generation
//    void unfold_pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);


};
