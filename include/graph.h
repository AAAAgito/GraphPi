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

class Graphmpi;

struct pair
{
    int e1;
    int e2;
};


class GraphTask {
public:

};

class Graph {
public:
    std::string raw_data_path;
    std::string backup_data_path;

    std::string patch_insert_path;
    std::string patch_delete_path;
    std::string patch_insert_path_backup;
    std::string patch_delete_path_backup;

    int available_threads=1;

    int v_cnt; // number of vertex
    int g_vcnt;
    int extra_v_cnt;
    int insert_cnt;
    int delete_cnt;
    unsigned int e_cnt; // number of edge
    unsigned int g_ecnt;
    size_t lg_ecnt;
    long long tri_cnt{}; // number of triangle

    int *edge; // edges
    unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]

    int *mem;
    int *mmp_edge;
    unsigned int *mmp_vertex;
    size_t *mmp_l_vertex;
    // with length of extra_v_cnt
    unsigned int *mmp_patch_vertex;
    // with length of insert_cnt
    int *mmp_patch_edge;
    
    // with length of extra_v_cnt
    unsigned int *mmp_patch_delete_vertex;
    // with length of insert_cnt
    int *mmp_patch_delete_edge;

    bool patch_read_open = false;
    int pi_fd;
    int pd_fd;

    bool *bitmap=NULL;
    int *dirty_bit=NULL;
    int INT_BITS = sizeof(int);
    int SHIFT = 5;
    int MASK = 0x1f;

    std::vector<std::vector<int>> patch_vi;
    std::vector<std::vector<int>> patch_vd;

    unsigned int *g_back_up_V=NULL;
    int *g_back_up_E = NULL;
    unsigned int backup_ecnt;
    int back_fd;
    int g_fd;

    Graph(std:: string path = "/home/yanglaoyuan/AsyncSubGraphStorage/docker_graphPi/bin/") {
        v_cnt = 0;
        e_cnt = 0;
        g_vcnt = 0;
        edge = nullptr;
        vertex = nullptr;
        raw_data_path = std::move(path);
        backup_data_path = raw_data_path + std::string("_backup");
        patch_insert_path = raw_data_path + "_patch.insert";
        patch_insert_path_backup = backup_data_path + "_patch.insert";
        patch_delete_path = raw_data_path + "_patch.delete";
        patch_delete_path_backup = backup_data_path + "_patch.delete";
        
    }

    ~Graph() {
        delete[] edge;
        delete[] vertex;
    }

    size_t intersection_size_clique(int v1,int v2);
    
    size_t intersection_size_clique(int v1,int v2, int v3);

    long long triangle_counting_mt(int thread_count);

    void tc_mt(long long *global_ans);

    int memory_map();
    
    int memory_lmap();

    void free_map(int &fd);
    
    void free_lmap(int &fd);

    void insert_edge(int v1, int v2);

    void delete_edge(int v1, int v2);

    void update(std::vector<int> &u, std::vector<int> &d);

    void dump_v(int offset);

    void dump_coarse(int offset);

    void init_patch_ptr();

    void release_patch_ptr();

    void refine_graph(bool preprocess = false);

    //general pattern matching algorithm with multi thread
    long long pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);
    
    long long pattern_matching_oc(const Schedule& schedule, int thread_count, bool clique = false);
    
    long long pattern_matching_oca(const Schedule& schedule, int thread_count, bool clique = false);

    void to_global_csr(const std::string& path);

    void load_global_graph(const std::string &path);

    void get_edge_index(int v, unsigned int& l, unsigned int& r) const;

    void get_large_mmp_edge_index(int v, size_t& l, size_t& r) const;
    
    void get_mmp_edge_index(int v, unsigned int& l, unsigned int& r) const;
    
    void get_mmp_patch_edge_index(int v, unsigned int& l, unsigned int& r) const;
    
    void get_mmp_patch_delete_edge_index(int v, unsigned int& l, unsigned int& r) const;
    
    int max_degree{};

    void generate_request(double ins_rate, double del_rate, std::vector<int> &ins, std::vector<int> &del);

    int get_threads() {return available_threads;}

    void set_threads(int num) {available_threads = num;}

    void needle (int *&p, int v, int &size);

    void bandaid(int *p, int v, int &size);

    int init_backup_ptr(int v, unsigned int e);

    void close_backup_ptr(int fd, int v, unsigned int e);

    void patch_recycle();

    void update_patch(const std::string& piu, const std::string& pdu, int in, int dn);

    void generate_patch_request(double ins_rate, double del_rate, int split);

    void set(int *p, int i) {p[i >> SHIFT] |= 1 << (i & MASK);}
    
    int test(int *p, int i) {return p[i >> SHIFT] & (1 << (i & MASK));}
    
    int clear(int *p, int i) {return p[i >> SHIFT] & ~(1 << (i & MASK));}

private:
    friend Graphmpi;

    void pattern_matching_aggressive_func_oc(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);
    
    void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    void pattern_matching_oca(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    double get_wall_time2() {
        struct timeval time;
        if(gettimeofday(&time,NULL)) {
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
    }
};
