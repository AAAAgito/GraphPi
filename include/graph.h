#pragma once
#include "schedule.h"
#include "vertex_set.h"
#include <map>
#include <utility>
#include <vector>
#include <assert.h>
#include <string>

enum PartitionType {
    METIS,
    RANDOM
};

enum BlockType {
    RANDOM_BLOCK,
    K_CORE_BLOCK
};

class Graphmpi;
struct VertexTable {
    bool is_intra = false;
    bool is_rooted = false; // is used as a root(depth=0) in subgraph matching
    unsigned int k_hop = 0;
    int loading = 0;
    int file_id = 0;
};

struct path{
    int depth;
    VertexSet subtraction_set;
    VertexSet *vertex_set;
    std::vector<int> vertex;
    long long local_ans;
    bool load_queue = false;
};

struct loader{
    int k_hop;
    int vertex;
    int file_id;
};


class Graph {
public:
    int v_cnt; // number of vertex
    int g_vcnt;
    unsigned int e_cnt; // number of edge
    long long tri_cnt{}; // number of triangle
    double max_running_time = 60 * 60 * 24; // second

    int *edge; // edges
    unsigned int *vertex; // v_i's neighbor is in edge[ vertex[i], vertex[i+1]-1]
    unsigned int external_space =1; // how much can be load from inter-partition edges
    unsigned int external_used =0;
    float extern_thresold{};
    int *inter_edge; // store inter partition edges
    unsigned int *inter_vertex;
    int inter_vtx_num = 0;
    unsigned int inter_edge_num = 0;
    BlockType blockType = RANDOM_BLOCK;
    std:: string raw_data_path;

    std::map<int, int> intra_vertex_dict; // map v_id to index of vertex in vertex array
    std::map<int, int> inter_vertex_dict; // map to_load v_id to index of to_load vertex in
    std::map<int, VertexTable> v_state_map;

    std::vector<path*> candidate_bin;
    std::vector<path*> ready_bin;
    std::vector<loader> load_list;

    Graph(std:: string path = "/mnt/d/graph/") {
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
    void resume_matching(const Schedule& schedule, path *p);

    //general pattern matching algorithm with multi thread
    long long pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

    //this function will be defined at code generation
//    long long unfold_pattern_matching(const Schedule& schedule, int thread_count, bool clique = false);

    //general pattern matching algorithm with multi thread ans multi process
    long long pattern_matching_mpi(const Schedule& schedule, int thread_count, bool clique = false);

    // TODO: in csr, vertex should be sort, adj list also
    void to_partition_csr(PartitionType t, int num, const std::string& path);

    void init_extern_storage(int v_num, int e_num);

    void to_global_csr(const std::string& path);

    void to_block_csr(int block_size, const std::string& path, int k_core = 0);

    void load_partition_graph(int pid, int num, const std::string& path);

    void load_global_graph(const std::string &path);

    int max_degree{};
    void load_extern_data(int file_id, const std::vector<loader>& load_vertex);

    void get_edge_index(int v, unsigned int& l, unsigned int& r) const;

    void get_extern_edge_index(int v, unsigned int& l, unsigned int& r) const;
private:
    friend Graphmpi;
    void tc_mt(long long * global_ans);


    // do pre-process to candidate_bin and generate the best load order (considering batch load)
    void extern_load_manage(unsigned int available_space, const Schedule& schedule);


    // manage loading order of vertex, and batch load vertices in same file
    void loading_manage();

    void load_list_append(loader l);

    // tracking memory use of extern vertex, drop all if larger than thresold
    void extern_drop_manage();

    // return true if it is still recommend to store.
    // return false if it is recommended to drop.
    bool extern_store_manage(const Schedule& schedule);

    void pattern_matching_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, long long& local_ans, int depth, bool clique = false);

    void pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

    void blocking_data_manage_k_core(int len, int &block_id, int &local_block_size, std::vector<int> &block_vid, const std::map<int, int> &revert_dict, std::vector<int> &to_insert_block, std::vector<int> &vid, std::vector<int> &file_id, const std::string& path);

    //this function will be defined at code generation
//    void unfold_pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth);

};
