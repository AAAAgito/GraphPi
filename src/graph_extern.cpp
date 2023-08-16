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
#include <mpi.h>
#include <cmath> 
#include <atomic>
#include <queue>
#include <iostream>

void Graph::load_physical_edge_index(int &v, std::shared_ptr<int[]> &ptr, int &size) {
    io_num+=1;
    int file_id = v_state_map[v].file_id;

    std::string block_path(file_path);
    
    if (blockType == BlockType::K_CORE_BLOCK) block_path+=file_str[file_id];
    else if (blockType == BlockType::RANDOM_BLOCK) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_br");
    else if (blockType == BlockType::CHAIN) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bc");
    else if (blockType == BlockType::SIMPLE_BFS) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bs");

    block_length lens = file_len[file_id];
    int load_v_len=lens.v_len;
    unsigned int load_e_len=lens.e_len;
    int* data = new int[2*load_v_len+load_e_len];

    double t1 = get_wall_time();
    DataLoader::load_block_data_aggregate2(data,2*load_v_len+load_e_len,block_path);
    double t2 = get_wall_time();
    int *vid = data;
    unsigned int *vertex_offset = (unsigned int *)data+load_v_len;
    int *adj_list = data+2*load_v_len;
    int vid_cursor = 0;
    while (vid[vid_cursor] != v) {
//        assert(vid_cursor < load_v_len);
//        assert(vid[vid_cursor] < v);
        vid_cursor++;
    }
//    assert(vid[vid_cursor]==v);
    int l,r;
    l=vertex_offset[vid_cursor];
    if (vid_cursor==load_v_len-1) r=load_e_len;
    else r= vertex_offset[vid_cursor+1];
//    assert(r-l >0);

    std::shared_ptr<int[]> p(new int[r-l]());

    size = r-l;
    for (int i = l; i < r; i++) {
        p[i - l] = adj_list[i];
    }

    ptr.swap(p);
    delete[] data;

    if (omp_get_thread_num()==0)
        blocking_manage_time += t2-t1;
}
