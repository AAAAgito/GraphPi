//
// Created by yanglaoyuan on 7/7/23.
//
#include "../include/graph.h"
#include "../include/graphmpi.h"
#include "../include/vertex_set.h"
#include "../include/common.h"
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <cstring>
#include <atomic>
#include <queue>
#include <iostream>
#include <random>

void Graph::to_global_csr(const std::string &path) {
    std::vector<int> vid, edges, v_order;
    std::vector<unsigned int> vtx_offset;
    unsigned int cursor = 0;
    for (int i=0; i< v_cnt; i++) {
        std::vector<int> local_edges;
        int len;
        if (i == v_cnt - 1) len = e_cnt - vertex[i];
        else len = vertex[i+1] - vertex[i];
        for (int j = 0; j < len; j++) local_edges.push_back(edge[vertex[i]+j]);
        std::sort(local_edges.begin(), local_edges.end());
        for (auto j : local_edges) edges.push_back(j);
        vtx_offset.push_back(cursor);
        vid.push_back(i);
        v_order.push_back(v_order.size());
        cursor+= len;
    }
    DataLoader::gen_partition_file(v_cnt,e_cnt,vtx_offset.data(),edges.data(),path);
    std::string graph_size(path);
    // printf("%d %d %s\n",v_cnt,e_cnt,path.c_str());
    graph_size.append(".size");
    DataLoader::gen_data_size(v_cnt,e_cnt,graph_size);
}

void Graph::load_global_graph(const std::string &path) {
    std::string graph_size(path);
    graph_size.append(".size");
    DataLoader::load_data_size(v_cnt,e_cnt,graph_size);
    vertex = new unsigned int[v_cnt];
    edge = new int[e_cnt];
    DataLoader::load_partition_data(v_cnt,e_cnt,vertex,edge,path);
    g_vcnt = v_cnt;
    g_ecnt = e_cnt;
}

void Graph::generate_patch_request(double ins_rate, double del_rate, int split) {
    std::vector<int> ins,del;
    generate_request(ins_rate,del_rate,ins,del);
    // printf("total %d %d\n",ins.size(),del.size());
    for (int i=0;i<split;i++){
        std::vector<int> sub_ins(ins.begin()+ins.size()*i/split,ins.begin()+ins.size()*(i+1)/split);
        std::vector<int> sub_del(del.begin()+del.size()*i/split,del.begin()+del.size()*(i+1)/split);
        printf("%d %d\n",i,sub_ins.size());
        update(sub_ins,sub_del);
        dump_coarse(i);
    }
    if (split==0) {
        printf("============================\n");
        std::vector<int> emp;
        update(emp,ins);
        double t1 = get_wall_time2();
        dump_v(0);
        double t2 = get_wall_time2();
        printf("PRE-DUMP %.6lf\n", t2 - t1);
        t1 = get_wall_time2();
        refine_graph();
        t2 = get_wall_time2();
        printf("PRE-REFINE %.6lf\n", t2 - t1);
        free_map(g_fd);
        memory_map();
        update(ins,del);
        t1 = get_wall_time2();
        dump_v(0);
        t2 = get_wall_time2();
        printf("FORMAL-DUMP %.6lf\n", t2 - t1);
        printf("============================\n");
    }
}

void Graph::generate_request(double ins_rate, double del_rate, std::vector<int> &ins, std::vector<int> &del) {
    int FLOAT_MIN = 0;
    int FLOAT_MAX = 1;
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
    ins.clear();
    del.clear();
    for (int i=0; i<g_vcnt; i++) {
        // printf("%d %d\n",i,g_vcnt);
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        for (unsigned int j=l;j<r;j++){
            if (mmp_edge[j]<i) continue;
            float r = distr(eng);
            if (r < ins_rate){
                ins.push_back(i);
                ins.push_back(mmp_edge[j]);
            }
            else if (r>=ins_rate && r<ins_rate+del_rate) {
                del.push_back(i);
                del.push_back(mmp_edge[j]);
            }
        }
    }
    std::string ip(raw_data_path);
    ip.append("_insert.edge");
    DataLoader::gen_data_file(ins.data(),ins.size(),ip);
    
    std::string dp(raw_data_path);
    dp.append("_delete.edge");
    DataLoader::gen_data_file(del.data(),del.size(),dp);
}