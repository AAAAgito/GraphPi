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

class Graph {
public:
    std::string raw_data_path;
    std::string backup_data_path;

    int available_threads=1;

    int v_cnt; // number of vertex
    int g_vcnt;
    int extra_v_cnt;
    int insert_cnt;
    int delete_cnt;
    unsigned int e_cnt; // number of edge
    unsigned int g_ecnt;
    long long tri_cnt{}; // number of triangle

    int *edge; // edges
    unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]

    int *mem;
    int *mmp_edge;
    unsigned int *mmp_vertex;
    // TODO: following 2 items not initial
    // with length of extra_v_cnt
    unsigned int *mmp_patch_vertex;
    // with length of insert_cnt
    int *mmp_patch_edge;

    int *bitmap;
    int INT_BITS = sizeof(int);
    int SHIFT = 5;
    int MASK = 0x1f;

    // TODO: following 6 items not initial

    bool patch_read_open = false;
    int pi_fd;
    int pic_fd;
    int pdc_fd;

    int *patch_insert_idx;
    int *patch_delete_idx;
    // idx, stop_poi
    // start_poi = idx-1 's stop_poi (0 at 0)
    int *patch_insert_code;
    // |(delete idx_0, delete len_0), (idx1, len1)|
    int *patch_delete_code;
    // update when doing dump_v
    int insert_code_len;
    int delete_code_len;
    std::vector<patch_data> patch_v;

    Graph(std:: string path = "/home/yanglaoyuan/AsyncSubGraphStorage/docker_graphPi/bin/") {
        v_cnt = 0;
        e_cnt = 0;
        g_vcnt = 0;
        edge = nullptr;
        vertex = nullptr;
        raw_data_path = std::move(path);
        backup_data_path = raw_data_path + std::string("_backup");
        
    }

    ~Graph() {
        delete[] edge;
        delete[] vertex;
    }

    int memory_map();

    void free_map(int &fd);

    void insert_edge(int v1, int v2);

    void delete_edge(int v1, int v2);

    void update(std::vector<pair> &u, std::vector<pair> &d);

    void update(std::vector<int> &u, std::vector<int> &d);

    void dump_v();

    void init_patch_ptr();

    void release_patch_ptr();

    void refine_graph();
    
    void set(int *p, int i) {p[i >> SHIFT] |= 1 << (i & MASK);}
    
    int test(int *p, int i) {return p[i >> SHIFT] & (1 << (i & MASK));}
    
    int clear(int *p, int i) {return p[i >> SHIFT] & ~(1 << (i & MASK));}

    //general pattern matching algorithm with multi thread
    long long pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);
    
    long long pattern_matching_oc(const Schedule& schedule, int thread_count, bool clique = false);

    void to_global_csr(const std::string& path);

    void load_global_graph(const std::string &path);

    void get_edge_index(int v, unsigned int& l, unsigned int& r) const;
    
    void get_mmp_edge_index(int v, unsigned int& l, unsigned int& r) const;
    
    void get_mmp_patch_edge_index(int v, unsigned int& l, unsigned int& r) const;

    void get_insert_code_index(int v, int &l, int &r) const;

    void get_delete_code_index(int v, int &l, int &r) const;
    
    int max_degree{};

    void generate_request(double ins_rate, double del_rate, std::vector<int> &ins, std::vector<int> &del);

    int get_threads() {return available_threads;}

    void set_threads(int num) {available_threads = num;}

private:
    friend Graphmpi;

    void pattern_matching_aggressive_func_oc(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);
    
    void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    double get_wall_time2() {
        struct timeval time;
        if(gettimeofday(&time,NULL)) {
            return 0;
        }
        return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
    }
};
