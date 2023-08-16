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
#include <atomic>
#include <queue>
#include <iostream>



void Graph::init_extern_storage(int v_num, int e_num) {
    physicals_length.resize(g_vcnt);
    physicals_load.resize(g_vcnt);
}

void Graph::gen_out_of_core_component(int part_num, int block_size, const std::string &path, int k_core) {
    std::vector<int> side_vtx;
    std::vector<int> part_map;
    for (int i=0; i < v_state_map.size();i++) {
        unsigned int l,r;
        get_edge_index(i,l,r);
    }
    to_partition_csr(part_num,side_vtx,part_map,path);
    to_block_csr(block_size,side_vtx,part_map,part_num,path,k_core);
//    bfs_order_queue.clear();
//    double t1 = get_wall_time2();
//    for (int i=0;i<part_num; i++)
//        gen_bipartite_order(part_map,i,extern_upper_thresold*v_cnt);
//    double t2 = get_wall_time2();
//    assert(bfs_order_queue.size()==v_cnt);
//    for (int i =0; i<v_cnt;i++) {
////        printf(" %d",bfs_order_queue[i]);
//        assert(bfs_order_queue[i]>=0 && bfs_order_queue[i] < v_cnt);
//    }
//    std::string order_path = raw_data_path + "_order_BI";
//    DataLoader::gen_order_file(bfs_order_queue.data(),v_cnt,order_path);
//    printf("gen order time %.6lf\n",t2-t1);
}

void Graph::to_global_csr(const std::string &path) {
    std::vector<int> vid, edges, v_order;
    std::vector<unsigned int> vtx_offset;
    unsigned int cursor = 0;
    if (intra_vertex_dict.empty()) {
        for (int i=0; i< v_cnt; i++) intra_vertex_dict[i] = i;
    }
    for (auto i: intra_vertex_dict) {
        std::vector<int> local_edges;
        int len;
        if (i.second == v_cnt - 1) len = e_cnt - vertex[i.second];
        else len = vertex[i.second+1] - vertex[i.second];
        for (int j = 0; j < len; j++) local_edges.push_back(edge[vertex[i.second]+j]);
        std::sort(local_edges.begin(), local_edges.end());
        for (auto j : local_edges) edges.push_back(j);
        vtx_offset.push_back(cursor);
        vid.push_back(i.first);
        v_order.push_back(v_order.size());
        cursor+= len;
    }

    std::string map_data = path + std::string("_m");
    DataLoader::gen_map_file(vid.data(), v_order.data(), v_cnt, map_data);
    DataLoader::gen_partition_file(v_cnt,e_cnt,vtx_offset.data(),edges.data(),path);
}

void Graph::load_global_graph(const std::string &path) {
    intra_vertex_dict.clear();
    v_state_map.clear();
    DataLoader::load_data_size(v_cnt,e_cnt,path);
    v_state_map.resize(v_cnt);
    vertex = new unsigned int[v_cnt];
    edge = new int[e_cnt];
    DataLoader::load_partition_data(v_cnt,e_cnt,vertex,edge,path);
    std::string map_data = path + std::string("_m");
    int *vkey = new int[v_cnt];
    int *vvalue = new int[v_cnt];
    DataLoader::load_map_data(vkey,vvalue,v_cnt,map_data);
    for (int i=0;i<v_cnt;i++) {
        int key = vkey[i];
        int value = vvalue[i];
        v_state_map[key].is_intra = true;
        intra_vertex_dict.insert(std::make_pair(key,value));
    }
    std::string file_map;
    std::string block_size_path;
    if (blockType == BlockType::K_CORE_BLOCK) {
        file_map = path + std::string("_fk");
        block_size_path = path + std::string("_fk_size");
    }
    if (blockType == BlockType::RANDOM_BLOCK) {
        file_map = path + std::string("_fr");
        block_size_path = path + std::string("_fr_size");
    }
    if (blockType == BlockType::CHAIN) {
        file_map = path + std::string("_fc");
        block_size_path = path + std::string("_fc_size");
    }
    if (blockType == BlockType::SIMPLE_BFS) {
        file_map = path + std::string("_fs");
        block_size_path = path + std::string("_fs_size");
    }
    int *fkey = new int[v_cnt];
    int *fvalue = new int[v_cnt];
    DataLoader::load_map_data(fkey,fvalue,v_cnt,file_map);
    int block_size = 0;
    for(int i=0;i<v_cnt;i++) {
        v_state_map[fkey[i]].file_id = fvalue[i];
        block_size = std::max(block_size,fvalue[i]+1);
    }
    int *key = new int[block_size];
    int *val1 = new int[block_size];
    int *val2 = new int[block_size];
    DataLoader::load_block_size_data(key,val1,val2,block_size,block_size_path);
    for (int i=0;i<block_size;i++) block_lengths[key[i]] = block_length(val1[i],val2[i]);
    g_vcnt = v_cnt;
    g_ecnt = e_cnt;
    delete[] key;
    delete[] val1;
    delete[] val2;
    delete[] vkey;
    delete[] vvalue;
    delete[] fkey;
    delete[] fvalue;
}

void Graph::load_partition_graph(int pid, int num, const std::string &path) {
    intra_vertex_dict.clear();
    v_state_map.clear();
    v_state_map.resize(g_vcnt);
    std::string partition_path,partition_map_data;
    if (partitionType==RANDOM) {
        partition_path = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                         std::string("_pr");
        partition_map_data = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                             std::string("_mr");
    }
    if (partitionType==METIS) {
        partition_path = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                         std::string("_pm");
        partition_map_data = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                             std::string("_mm");
    }
    if (partitionType==NAIVE_BFS) {
        partition_path = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                         std::string("_pb");
        partition_map_data = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                             std::string("_mb");
    }
    if (partitionType==LDG) {
        partition_path = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                         std::string("_pl");
        partition_map_data = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                             std::string("_ml");
    }
    if (partitionType==PAGE) {
        partition_path = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                         std::string("_ppg");
        partition_map_data = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) +
                             std::string("_mpg");

    }
    DataLoader::load_data_size(v_cnt,e_cnt,partition_path);
    vertex = new unsigned int[v_cnt];
    edge = new int[e_cnt];
    DataLoader::load_partition_data(v_cnt,e_cnt,vertex,edge,partition_path);
    int *vkey = new int[v_cnt];
    int *vvalue = new int[v_cnt];
    DataLoader::load_map_data(vkey,vvalue,v_cnt,partition_map_data);

    std::string file_map;
    std::string block_size_path;
    if (blockType == BlockType::K_CORE_BLOCK) {
        file_map = path + std::string("_fk");
        block_size_path = path + std::string("_fk_size");
    }
    if (blockType == BlockType::RANDOM_BLOCK) {
        file_map = path + std::string("_fr");
        block_size_path = path + std::string("_fr_size");
    }
    if (blockType == BlockType::CHAIN) {
        file_map = path + std::string("_fc");
        block_size_path = path + std::string("_fc_size");
    }
    if (blockType == BlockType::SIMPLE_BFS) {
        file_map = path + std::string("_fs");
        block_size_path = path + std::string("_fs_size");
    }
    int *fkey = new int[g_vcnt];
    int *fvalue = new int[g_vcnt];
    int block_size = 0;
    DataLoader::load_map_data(fkey,fvalue,g_vcnt,file_map);
    for (int i=0;i<g_vcnt;i++) {
        int key = fkey[i];
        int value = fvalue[i];
        v_state_map[key].file_id = value;
        block_size = std::max(block_size,fvalue[i]+1);
    }
    for (int i=0;i<v_cnt;i++) {
        int key = vkey[i];
        int value = vvalue[i];
        v_state_map[key].is_intra = true;
        intra_vertex_dict.insert(std::make_pair(key,value));
    }
    // TODO: delete
    for (int i=0; i<v_cnt; i++) v_state_map[vkey[i]].is_intra = true;
    int *key = new int[block_size];
    int *val1 = new int[block_size];
    int *val2 = new int[block_size];
    DataLoader::load_block_size_data(key,val1,val2,block_size,block_size_path);
    for (int i=0;i<block_size;i++) {
        block_lengths[key[i]] = block_length{val1[i], val2[i]};
    }

    delete[] key;
    delete[] val1;
    delete[] val2;
    delete[] vkey;
    delete[] vvalue;
    delete[] fkey;
    delete[] fvalue;
}

void Graph::to_block_csr(int block_size, std::vector<int> &side_vertex, std::vector<int> &part_map, int part_num, const std::string& path, int k_core) {
    int len = sizeof (char) * block_size / sizeof (int);
    if (blockType == BlockType::RANDOM_BLOCK) {
        std::vector<int> vid, edges, map_k, map_v;
        std::vector<unsigned int> vtx_offset;
        std::vector<int> file_len;
        file_len.push_back(0);
        // id should be sorted by keys
        for(auto & id : intra_vertex_dict) {

            int adj_len;
            if (id.second < v_cnt - 1) adj_len = vertex[id.second +1] - vertex[id.second];
            else adj_len = e_cnt - vertex[id.second];
            for(int i=0; i< file_len.size(); i++) {
                if (file_len[i]+ 2 + adj_len < len) {
                    file_len[i] += 2+adj_len;
                    map_k.push_back(id.first);
                    map_v.push_back(i);
                    break;
                }
                if (i == file_len.size()-1) {
                    map_k.push_back(id.first);
                    map_v.push_back(file_len.size());
                    file_len.push_back(2+adj_len);
                    break;
                }
            }
        }
        // dump vertex-block map
        std::string block_map_path = path + std::string("_fr");
        DataLoader::gen_map_file(map_k.data(),map_v.data(),map_v.size(),block_map_path);
        // dump each block
        printf("dump eack block\n");
        for (int block_id = 0; block_id < file_len.size(); block_id++) {
            int vtx_cursor = 0;
            vid.clear();
            edges.clear();
            vtx_offset.clear();
            for (int i = 0; i< map_v.size(); i++) {
                if (map_v[i]==block_id) {
                    vtx_offset.push_back(vtx_cursor);
                    vid.push_back(map_k[i]);
                    int adj_len;
                    int adj_poi = vertex[intra_vertex_dict[map_k[i]]];
                    if (intra_vertex_dict[map_k[i]] == v_cnt-1) adj_len = e_cnt - adj_poi;
                    else adj_len = vertex[intra_vertex_dict[map_k[i]]+1] - vertex[intra_vertex_dict[map_k[i]]];
                    for(int j = 0; j < adj_len; j++) edges.push_back(edge[adj_poi+j]);
                    vtx_cursor += adj_len;
                }
            }
            std::string block_csr_path = path + std::string("_blocks/") + std::to_string(block_id) + std::string("_br");
            DataLoader::gen_block_file(vid.data(),vtx_offset.data(),edges.data(),vid.size(),edges.size(),block_csr_path);
            block_lengths[block_id] = block_length{(int)vid.size(),(int)edges.size()};
        }
        std::string block_size_path = path + std::string("_fr_size");
        std::vector<int> key,vl,el;
        for (auto i : block_lengths) {
            key.push_back(i.first);
            vl.push_back(i.second.v_len);
            el.push_back(i.second.e_len);
        }
        DataLoader::gen_block_size_file(key.data(),vl.data(),el.data(),key.size(),block_size_path);
    }
    if (blockType == BlockType::K_CORE_BLOCK) {
        block_lengths.clear();
        int block_id = 0;
        std::vector<int> degree;
        degree.resize(v_cnt);
        std::vector<int> block_vid, vid, file_id;
        int local_block_size = 0;
        std::vector<int> to_insert_block;
        for (int part=0; part < part_num; part++) {
            for (int i = 0; i < v_cnt; i++) {
                if (part_map[i] != part) {
                    degree[i] = -1;
                }
                else {
                    unsigned int l, r;
                    get_edge_index(i, l, r);
                    int minus = 0;
                    for (int j = l; j < r; j++)
                        if (part_map[edge[j]] != part) minus += 1;
                    degree[i] = r - l - minus;
                }
            }
            for (int core = 1; core <= k_core; core++) {
                bool keep_lop = true;
                while (keep_lop) {
                    keep_lop = false;
                    std::vector<int> local_remove;
                    for (int i=0; i<v_cnt; i++) {
                        if (degree[i] < core && degree[i]!=-1) {
                            unsigned int l,r;
                            get_edge_index(i,l,r);
                            for (int j = l; j < r; j++) {
                                if (part_map[edge[j]]!=part) continue;
                                if (degree[edge[j]] != -1) degree[edge[j]] -= 1;
                            }
                            to_insert_block.push_back(i);
                            local_remove.push_back(i);
                        }
                    }
                    for (auto i: local_remove) degree[i]=-1;
                    for (int i=0; i<v_cnt; i++) keep_lop |= degree[i] < core && degree[i]!=-1 && part_map[edge[i]]==part;
                }
                blocking_data_manage_k_core(len, block_id, local_block_size, block_vid,  to_insert_block,
                                            path);
            }
            for (int i=0; i<v_cnt; i++)
                if (degree[i]!=-1)
                    to_insert_block.push_back(i);
            blocking_data_manage_k_core(len, block_id, local_block_size, block_vid, to_insert_block, path);
        }


        std::string block_map_path = path + std::string("_fk");
        std::vector<int> fkey,fvalue;
        for (int i=0;i< v_state_map.size();i++) {
            fkey.push_back(i);
            fvalue.push_back(v_state_map[i].file_id);
        }
        DataLoader::gen_map_file(fkey.data(),fvalue.data(),fkey.size(),block_map_path);

        std::string block_size_path = path + std::string("_fk_size");
        std::vector<int> key,vl,el;
        for (auto i : block_lengths) {
            key.push_back(i.first);
            vl.push_back(i.second.v_len);
            el.push_back(i.second.e_len);
        }
        DataLoader::gen_block_size_file(key.data(),vl.data(),el.data(),key.size(),block_size_path);
    }

    if (blockType == BlockType::CHAIN) {
        
    }

    if (blockType == BlockType::SIMPLE_BFS) {
        int block_id = 0;
        for (int part = 0; part < part_num; part++) {
            std::vector<int> part_bfs_list;
            bool *choose = new bool[v_cnt];
            int intra_num = 0;
            for (int i = 0; i < v_cnt; i++) {
                choose[i] = part_map[i] != part;
                if (!choose[i]) intra_num++;
            }
            std::queue<int> q;
            for (int i=0; i<v_cnt; i++) {
                if (q.empty()) {
                    if (choose[i]) continue;
                    q.push(i);
                    part_bfs_list.push_back(i);
                    choose[i] = true;
                }
                while (!q.empty()) {
                    int v = q.front();
                    q.pop();
                    unsigned int l,r;
                    get_edge_index(v,l,r);
                    for (int j =l;j<r;j++) {
                        if (choose[edge[j]]) continue;
                        q.push(edge[j]);
                        choose[edge[j]] = true;
                        part_bfs_list.push_back(edge[j]);
                    }
                }
            }
            assert(intra_num==part_bfs_list.size());
            std::vector<int> block_vid;
            int block_vsize = 0;
            for (auto v : part_bfs_list) {
                unsigned int l,r;
                get_edge_index(v,l,r);
                int qlen = 2+r-l;
                if (block_vsize + qlen > len && !block_vid.empty()) {
                    std::sort(block_vid.begin(), block_vid.end());

                    std::vector<int> edges;
                    std::vector<unsigned int> vtx_offset;
                    unsigned int cursor = 0;
                    for (auto j : block_vid) {
                        vtx_offset.push_back(cursor);
                        unsigned int l,r;
                        get_edge_index(j,l,r);
                        for (int k = l; k < r; k++) edges.push_back(edge[k]);
                        cursor += r-l;
                        v_state_map[j].file_id = block_id;
                    }

                    std::string block_csr_path = path + std::string("_blocks/") + std::to_string(block_id) + std::string("_bs");
                    DataLoader::gen_block_file(block_vid.data(),vtx_offset.data(),edges.data(),block_vid.size(),edges.size(),block_csr_path);
                    block_lengths[block_id] = block_length{(int)block_vid.size(),(int)edges.size()};
                    block_vid.clear();
                    block_vsize = 0;
                    block_id+=1;

                }
                block_vid.push_back(v);
                block_vsize += qlen;
            }

            std::sort(block_vid.begin(), block_vid.end());

            std::vector<int> edges;
            std::vector<unsigned int> vtx_offset;
            unsigned int cursor = 0;
            for (auto j : block_vid) {
                vtx_offset.push_back(cursor);
                unsigned int l,r;
                get_edge_index(j,l,r);
                for (int k = l; k < r; k++) edges.push_back(edge[k]);
                cursor += r-l;
                v_state_map[j].file_id = block_id;
            }

            std::string block_csr_path = path + std::string("_blocks/") + std::to_string(block_id) + std::string("_bs");
            DataLoader::gen_block_file(block_vid.data(),vtx_offset.data(),edges.data(),block_vid.size(),edges.size(),block_csr_path);
            block_lengths[block_id] = block_length{(int)block_vid.size(),(int)edges.size()};
            block_vid.clear();
            block_vsize = 0;
            block_id+=1;
            delete[] choose;
        }

        std::string block_map_path = path + std::string("_fs");
        std::vector<int> fkey,fvalue;
        for (int i=0;i< v_state_map.size();i++) {
            fkey.push_back(i);
            fvalue.push_back(v_state_map[i].file_id);
        }
        DataLoader::gen_map_file(fkey.data(),fvalue.data(),fkey.size(),block_map_path);

        std::string block_size_path = path + std::string("_fs_size");
        std::vector<int> key,vl,el;
        for (auto i : block_lengths) {
            key.push_back(i.first);
            vl.push_back(i.second.v_len);
            el.push_back(i.second.e_len);
        }
        DataLoader::gen_block_size_file(key.data(),vl.data(),el.data(),key.size(),block_size_path);
    }
}

void Graph::blocking_data_manage_k_core(int len, int &block_id, int &local_block_size, std::vector<int> &block_vid, std::vector<int> &to_insert_block, const std::string& path) {
    std::queue<int> q;
    bool *pick = new bool[v_cnt];
    for (int i=0;i<v_cnt;i++) pick[i]=true;
    for (auto i : to_insert_block) pick[i] = false;
    to_insert_block.clear();

    for (int v0 = 0; v0<v_cnt;v0++) {
        if (!pick[v0]) {
            q.push(v0);
            pick[v0] = true;
        }
        while (!q.empty()) {
            unsigned int l_front,r_front;
            int v_front=q.front();
            get_edge_index(v_front,l_front,r_front);
            int q_len = r_front-l_front;

            if (local_block_size + 2 + q_len > len && !block_vid.empty()) {
                std::sort(block_vid.begin(), block_vid.end());
                int cursor = 0;
                std::vector<unsigned int> vtx_offset;
                std::vector<int> edges;
                for (auto i : block_vid) {
                    vtx_offset.push_back(cursor);
                    unsigned int l,r;
                    get_edge_index(i,l,r);
                    cursor += r-l;

                    for (int j = l; j < r; j++) edges.push_back(edge[j]);
                }
                block_lengths[block_id] = block_length{(int) block_vid.size(), (int) edges.size()};
                std::string block_csr_path =
                        path + std::string("_blocks/bk_") + std::to_string(block_id);
                DataLoader::gen_block_file(block_vid.data(), vtx_offset.data(), edges.data(), block_vid.size(),
                                           edges.size(), block_csr_path);
                for (auto i: block_vid) {
                    v_state_map[i].file_id = block_id;
                }

                block_id += 1;
                local_block_size = 0;
                block_vid.clear();
            }


            block_vid.push_back(v_front);
            local_block_size += 2 + q_len;
            q.pop();

            for (int j = l_front; j < r_front; j++) {
                int neigh_v = edge[j];
                if (!pick[neigh_v]) {
                    q.push(neigh_v);
                    pick[neigh_v] = true;
                }
            }
        }
    }

    if (!block_vid.empty()) {
        std::sort(block_vid.begin(), block_vid.end());
        int cursor = 0;
        std::vector<unsigned int> vtx_offset;
        std::vector<int> edges;
        for (auto i : block_vid) {
            vtx_offset.push_back(cursor);
            unsigned int l,r;
            get_edge_index(i,l,r);
            cursor += r-l;

            for (int j = l; j < r; j++) edges.push_back(edge[j]);
        }
        if (block_id<0) printf("here\n");
        block_lengths[block_id] = block_length{(int) block_vid.size(), (int) edges.size()};
        std::string block_csr_path =
                path + std::string("_blocks/bk_") + std::to_string(block_id);
        DataLoader::gen_block_file(block_vid.data(), vtx_offset.data(), edges.data(), block_vid.size(),
                                   edges.size(), block_csr_path);

        for (auto i: block_vid) {
            v_state_map[i].file_id = block_id;
        }
        block_id += 1;
        local_block_size = 0;
        block_vid.clear();
    }
    delete[] pick;
}

int get_intersect_num(int *d1, int *d2, int l1, int l2) {
    int i=0,j=0,cnt=1;
    while (i<l1 && j <l2) {
        if (d1[i] < d2[j]) i++;
        else if (d1[i] > d2[j]) j++;
        else {
            j++;
            i++;
            cnt++;
        }
    }
    return cnt;
}

void Graph::to_partition_csr(int num, std::vector<int> &side_vertex, std::vector<int> &part_map, const std::string& path) {
    part_map.resize(v_cnt);
    if (partitionType == PartitionType::RANDOM) {
        for (int p = 0; p<num; p++) {
            int partial_vcnt;
            partial_vcnt = v_cnt /num;
            if (p < v_cnt % num) partial_vcnt += 1;
            auto *part_vertex = new unsigned int[partial_vcnt];
            int *mapk = new int[partial_vcnt];
            int *mapv = new int[partial_vcnt];
            unsigned int len;
            std::vector<int> adj_list;
            auto it = intra_vertex_dict.begin();
            for (int i = 0; i< p; i++) it++;
            int idx = 0;
            while (it != intra_vertex_dict.end()) {
                std::vector<int> local_adj_list;
                int i = intra_vertex_dict[it->first];

                //calculate adj list len
                if (i == v_cnt-1) len = e_cnt-vertex[i];
                else len = vertex[i+1]-vertex[i];

                //collect adj list and sort
                for (unsigned int j = vertex[i]; j < vertex[i]+len; j++) local_adj_list.push_back(edge[j]);
                std::sort(local_adj_list.begin(), local_adj_list.end());

                //insert to global partition adj list
                for (auto &&j : local_adj_list) adj_list.push_back(j);
                mapk[idx] = it->first;
                mapv[idx] = idx;
                if (idx==0) {
                    part_vertex[idx] = 0;
                }
                if (idx < partial_vcnt-1){
                    part_vertex[idx+1] = part_vertex[idx] + len;
                }
                for (int l=0; l< num; l++) {
                    if (it == intra_vertex_dict.end()) break;
                    it++;
                }
                idx++;
            }

            std::string partition_map = path + "_" + std::to_string(p)+ "_" + std::to_string(num) + std::string("_mr");
            DataLoader::gen_map_file(mapk,mapv,partial_vcnt,partition_map);
            std::string partition_path = path+ "_" + std::to_string(p)+ "_" + std::to_string(num)  + std::string("_pr");
            DataLoader::gen_partition_file(partial_vcnt, adj_list.size(), part_vertex, adj_list.data(), partition_path);

            delete[] mapk;
            delete[] mapv;
            delete[] part_vertex;
        }
    }
    if (partitionType == PartitionType::METIS) {
    }
    if (partitionType == PartitionType::NAIVE_BFS) {
        std::map<int,bool> choose;
        std::map<int,bool> is_side;
        std::map<int,int> father;
        std::vector<int> v,d;
        v.resize(v_cnt);
        d.resize(v_cnt);
        for (int i = 0; i < v_cnt; i++) {
            unsigned int l,r;
            get_edge_index(i,l,r);
            d[i] = r-l;
            v[i] = i;
        }
        for (int i =0 ;i< v_cnt-1; i++) {
            for (int j=i;j< v_cnt;j++)
                if (d[i] < d[j]) {
                    int mid = d[i];
                    d[i] = d[j];
                    d[j] = mid;
                    mid = v[i];
                    v[i] = v[j];
                    v[j] = mid;
                }
        }
        for (int i=0; i< v_state_map.size();i++) {
            is_side[i] = false;
            choose[i] = false;
        }
        int cur_part = 0;
        int cur_size = 0;
        std::vector<int> part_vertex;
        std::queue<int> q;
        std::queue<int> prev_q;
        for (int i = 0; i<v_cnt; i++) {
            if (choose[v[i]]) continue;
            father[v[i]] = -1;
            q.push(v[i]);
            choose[v[i]] = true;
            while (!q.empty() || !prev_q.empty()) {
                // dump when edges num are enough
                if (q.empty()) {
                    while (choose[prev_q.front()] && !prev_q.empty()) prev_q.pop();
                    if (prev_q.empty()) continue;
                    q.push(prev_q.front());
                    choose[prev_q.front()] = true;
                    prev_q.pop();
                }
                if (cur_size >= e_cnt / num) {
                    cur_size=0;
                    // dump
                    int partial_vcnt = part_vertex.size();
                    int cursor = 0;
                    std::vector<int> adj_list, vertex_idx;
                    std::vector<unsigned int> vtx_offset;
                    std::sort(part_vertex.begin(), part_vertex.end());
                    for (auto j : part_vertex) {
                        vertex_idx.push_back(vertex_idx.size());
                        vtx_offset.push_back(cursor);
                        unsigned int l,r;
                        get_edge_index(j,l,r);
                        for (int idx= l;idx<r;idx++)
                            adj_list.push_back(edge[idx]);
                        cursor += r-l;
                    }
                    std::string partition_map = path + "_" + std::to_string(cur_part)+ "_" + std::to_string(num) + std::string("_mb");
                    DataLoader::gen_map_file(part_vertex.data(),vertex_idx.data(),partial_vcnt,partition_map);
                    std::string partition_path = path+ "_" + std::to_string(cur_part)+ "_" + std::to_string(num)  + std::string("_pb");
                    DataLoader::gen_partition_file(partial_vcnt, adj_list.size(), vtx_offset.data(), adj_list.data(), partition_path);
                    part_vertex.clear();
                    cur_part+=1;
                    while (!q.empty()) {
                        prev_q.push(q.front());
                        if (!is_side[father[q.front()]]) {
                            side_vertex.push_back(father[q.front()]);
                            is_side[father[q.front()]] = true;
                        }
                        choose[q.front()] = false;
                        q.pop();
                    }
                    continue;
                }
                int vtx = q.front();
                part_map[vtx] = cur_part;
                unsigned int l,r;
                get_edge_index(vtx,l,r);
                cur_size += r-l;
                part_vertex.push_back(vtx);
                q.pop();
                bool has_inter = false;
                for (int k = l;k<r;k++) {
                    has_inter |= (is_side[edge[k]] && part_map[edge[k]]!= cur_part);
                    if (choose[edge[k]]) continue;
                    father[edge[k]] = vtx;
                    q.push(edge[k]);
                    choose[edge[k]] = true;
                }
                if (has_inter && !is_side[vtx]) {
                    side_vertex.push_back(vtx);
                    is_side[vtx] = true;
                }
            }
        }
        // dump remain
        if (!part_vertex.empty()) {

            // dump
            int partial_vcnt = part_vertex.size();
            int cursor = 0;
            std::vector<int> adj_list, vertex_idx;
            std::vector<unsigned int> vtx_offset;
            std::sort(part_vertex.begin(), part_vertex.end());
            for (auto j : part_vertex) {
                vertex_idx.push_back(vertex_idx.size());
                vtx_offset.push_back(cursor);
                unsigned int l,r;
                get_edge_index(j,l,r);
                for (int idx= l;idx<r;idx++) {

                    adj_list.push_back(edge[idx]);
                }
                cursor += r-l;
            }
            std::string partition_map = path + "_" + std::to_string(cur_part)+ "_" + std::to_string(num) + std::string("_mb");
            DataLoader::gen_map_file(part_vertex.data(),vertex_idx.data(),partial_vcnt,partition_map);
            std::string partition_path = path+ "_" + std::to_string(cur_part)+ "_" + std::to_string(num)  + std::string("_pb");
            DataLoader::gen_partition_file(partial_vcnt, adj_list.size(), vtx_offset.data(), adj_list.data(), partition_path);

            part_vertex.clear();
            cur_part+=1;
        }
    }

    if (partitionType == PartitionType::LDG) {

        for (int i=0; i<v_cnt; i++) part_map[i] = -1;
        int capacity = e_cnt / num;
        std::vector<int> part_vertices[num];
        int p_size[num];
        for (int i = 0; i < num; i++) p_size[i] = 0;
        int *d = new int[v_cnt];
        int *v = new int[v_cnt];
        for (int i = 0; i < v_cnt; i++) {
            unsigned int l,r;
//            printf("%d %d\n",i, intra_vertex_dict.empty());
            get_edge_index(i,l,r);
            d[i] = r-l;
            v[i] = i;
        }
        for (int i = 0; i < num - 1; i++) {
            for (int j = i+1; j < num; j++) {
                if (d[i] < d[j]) {
                    auto mid = d[i];
                    d[i] = d[j];
                    d[j] = mid;
                    mid = v[i];
                    v[i] = v[j];
                    v[j] = mid;
                }
            }
        }
        int loading_progress = 0;
        for (int i = 0; i < v_cnt; i++) {
            int vtx= v[i];
            unsigned int l,r;
            get_edge_index(vtx,l,r);
            double score = -1.0 *v_cnt;
            int max = 0;
            for (int part = 0; part < num; part++) {
                int intersect_num = 0;
                double w_i = 1.0 - (double ) p_size[part] / capacity;
                for (int j=l;j<r;j++) {
                    if (part_map[edge[j]]==part) intersect_num++;
                }
//                double new_score = w_i / (double )(r-l) * (double )get_intersect_num(edge+l,part_vertices[part].data(),r-l,part_vertices[part].size());
                double new_score = w_i / (double )(r-l) * (double )intersect_num;
                if (score < new_score) {
                    score = new_score;
                    max = part;
                }
            }
            part_vertices[max].push_back(vtx);
            part_map[vtx] = max;
            p_size[max] += r-l;
        }
        for (int i = 0; i < num; i++) {
            std::vector<int> part_vertex = part_vertices[i];
            int partial_vcnt = part_vertex.size();
            int cursor = 0;
            std::vector<int> adj_list, vertex_idx;
            std::vector<unsigned int> vtx_offset;
//            std::sort(part_vertex.begin(), part_vertex.end());
            for (auto j: part_vertex) {
                vertex_idx.push_back(vertex_idx.size());
                vtx_offset.push_back(cursor);
                unsigned int l, r;
                get_edge_index(j, l, r);
                for (int idx = l; idx < r; idx++)
                    adj_list.push_back(edge[idx]);
                cursor += r - l;
            }
            std::string partition_map =
                    path + "_" + std::to_string(i) + "_" + std::to_string(num) + std::string("_ml");
            DataLoader::gen_map_file(part_vertex.data(), vertex_idx.data(), partial_vcnt, partition_map);
            std::string partition_path =
                    path + "_" + std::to_string(i) + "_" + std::to_string(num) + std::string("_pl");
            DataLoader::gen_partition_file(partial_vcnt, adj_list.size(), vtx_offset.data(), adj_list.data(),
                                           partition_path);
        }

//        for (auto v : intra_vertex_dict) {
//            unsigned int l,r;
//            get_edge_index(v.first,l,r);
//            for (int i = l; i < r; i++) {
//                if (part_map[edge[i]] != part_map[v.first]) {
//                    side_vertex.push_back(v.first);
//                    break;
//                }
//            }
//        }

        delete[] d;
        delete[] v;
    }
}

int Graph::to_partition_csr(int size_byte) {
    int size = size_byte / sizeof(int);
    int process =0;
    int page_num = 0;
    int cur_page_size=0;
    std::vector<bool> chosen, picked;
    std::vector<int> select_per_round, core_per_round;
    std::set<int> candidate;
    chosen.resize(v_cnt);
    picked.resize(v_cnt);
    for (int i=0; i<v_cnt; i++) {
        if (chosen[i]) continue;
        select_per_round.push_back(i);
        picked[i] = true;
        chosen[i] = true;
        if ((process+1)/1000 > process/1000) printf("process %dk\n",(process+1)/1000);
        process++;
        core_per_round.push_back(i);
        unsigned int l,r;
        get_edge_index(i,l,r);
        cur_page_size += 2+r-l;
        for (int j=l;j<r;j++) {
            if (picked[edge[j]]) continue;
            select_per_round.push_back(edge[j]);
            unsigned int l,r;
            get_edge_index(edge[j],l,r);
            cur_page_size += 2+r-l;
            picked[edge[j]] = true;
            if (!chosen[edge[j]])
                candidate.insert(edge[j]);
        }
        while (!candidate.empty()) {
//            printf("%d\n",process);
            if ((process+1)/1000 > process/1000) printf("process %dk\n",(process+1)/1000);
            int max_score = -1;
            int max_vtx = -1;
            for (auto v : candidate) {

                unsigned int l,r;
                get_edge_index(v,l,r);
                int cur_score = 0;
                for (int c =l;c<r;c++) {
                    if (picked[edge[c]]) cur_score++;
                }
                if (cur_score > max_score) {
                    max_score = cur_score;
                    max_vtx = v;
                }
            }
            unsigned int l,r;
            get_edge_index(max_vtx,l,r);
            int need_space =0;

            for (int k = l;k < r; k++) {
                if (picked[edge[k]]) continue;
                unsigned int h,t;
                get_edge_index(edge[k],h,t);
                need_space += 2+t-h;
            }
            if (cur_page_size +need_space <= size) {
                candidate.erase(max_vtx);
                chosen[max_vtx] = true;
                process++;
                core_per_round.push_back(max_vtx);
                for (int k = l;k < r; k++) {
                    if (picked[edge[k]]) continue;
                    picked[edge[k]] = true;
                    select_per_round.push_back(edge[k]);
                    if (!chosen[edge[k]])
                        candidate.insert(edge[k]);
                }
                cur_page_size+= need_space;
            }
            else {
                //drop
                std::string porder = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_opg";
                std::string pmap = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_mpg";
                std::string p = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_ppg";
                std::string pbasic = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_pcore";
                int core_num = core_per_round.size();
                DataLoader::gen_order_file(&core_num,1,pbasic);
                std::sort(select_per_round.begin(), select_per_round.end());
                std::sort(core_per_round.begin(), core_per_round.end());
                dump_csr_vset(select_per_round,core_per_round,p,pmap,porder);
                //drop done
                select_per_round.clear();
//                printf("%d\n",core_per_round.size());
                core_per_round.clear();
                for (int i=0; i<v_cnt;i++) picked[i] = false;
                candidate.clear();
                cur_page_size = 0;
                page_num++;
            }

        }
    }

    //drop
    std::string porder = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_opg";
    std::string pmap = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_mpg";
    std::string p = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_ppg";
    std::string pbasic = raw_data_path +"_" + std::to_string(page_num) + "_" + std::to_string(size_byte) + "_pcore";
    int core_num = core_per_round.size();
    DataLoader::gen_order_file(&core_num,1,pbasic);
    std::sort(select_per_round.begin(), select_per_round.end());
    std::sort(core_per_round.begin(), core_per_round.end());
    dump_csr_vset(select_per_round,core_per_round,p,pmap,porder);
    //drop done
    select_per_round.clear();
    core_per_round.clear();
    for (int i=0; i<v_cnt;i++) picked[i] = false;
    candidate.clear();
    cur_page_size = 0;
    page_num++;
    return page_num;
}
void Graph::dump_csr_vset(std::vector<int> &vs, std::vector<int> &cores, const std::string &p, const std::string &pmap, const std::string &porder) {
    DataLoader::gen_order_file(cores.data(),cores.size(),porder);
    std::vector<int> edges, vertex_idx;
    std::vector<unsigned int> vtx_offset;
    for (int i=0; i<vs.size();i++) {
        vertex_idx.push_back(vertex_idx.size());
        vtx_offset.push_back(edges.size());
        unsigned int l,r;
        get_edge_index(vs[i],l,r);
        for (int j=l;j<r;j++) edges.push_back(edge[j]);
    }
    DataLoader::gen_partition_file(vtx_offset.size(),edges.size(),vtx_offset.data(),edges.data(),p);
    DataLoader::gen_map_file(vs.data(),vertex_idx.data(),vs.size(),pmap);
}