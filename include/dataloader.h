#pragma once
#include "graph.h"
#include <map>
#include <algorithm>
#include <vector>
#include <fstream>
#include <ostream>
#include <istream>

enum DataType {
    Patents,
    Orkut,
    complete8,
    LiveJournal,
    MiCo,
    Twitter,
    CiteSeer,
    Wiki_Vote,
    Invalid
};

const long long Patents_tri_cnt = 7515023LL;
const long long LiveJournal_tri_cnt = 177820130LL;
const long long MiCo_tri_cnt = 12534960LL;
const long long CiteSeer_tri_cnt = 1166LL;
const long long Wiki_Vote_tri_cnt = 608389LL;
const long long Orkut_tri_cnt = 627584181LL;
const long long Twitter_tri_cnt = 34824916864LL;

class DataLoader {
public:

    bool load_data(Graph* &g, DataType type, const char* path, int oriented_type = 0);
        // pattern_diameter means max distance between two vertex in graph
        // oriented_type is used to reorder dataset
        // oriented_type == 0 do nothing
        //               == 1 high degree first
        //               == 2 low degree first
    bool load_complete(Graph* &g, int clique_size);


    static void gen_partition_file(int v_cnt, unsigned int e_cnt, unsigned int *v, int *e, const std::string& path);

    static void load_partition_data(int &v_cnt, unsigned int &e_cnt, unsigned int *v, int *e, const std::string& path);

    static void load_data_size(int &v_cnt, unsigned int &e_cnt, const std::string& path);

    static void gen_map_file(int* key, int* value, int size, const std::string& path);

    static void gen_map_file(unsigned int* key, int* value, int size, const std::string& path);

    static void load_map_data(int* key, int* value, int size, const std::string& path);

    static void load_map_data(unsigned int* key, int* value, int size, const std::string& path);

    static void gen_block_size_file(int* key, int* value, int* value2, int size, const std::string& path);

    static void load_block_size_data(int* key, int* value, int* value2, int size, const std::string& path);

    static void gen_block_file(int* vid, unsigned int* vertex, int* edge, int v_len, unsigned int e_len, const std::string& path);

    static void load_block_data(int* vid, unsigned int* vertex, int* edge, int &v_len, unsigned int &e_len, const std::string& path);

    static void load_block_data_aggregate(int *data, int size, const std::string& path);

private:
    static bool cmp_pair(std::pair<int,int>a, std::pair<int,int>b);
    static bool cmp_degree_gt(std::pair<int,int> a,std::pair<int,int> b);
    static bool cmp_degree_lt(std::pair<int,int> a,std::pair<int,int> b);

    long long comb(int n,int k);
    bool general_load_data(Graph* &g, DataType type, const char* path, int oriented_type = 0);
    bool twitter_load_data(Graph* &g, DataType type, const char* path, int oriented_type = 0);

    std::map<int,int> id;
};

