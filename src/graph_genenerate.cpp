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
