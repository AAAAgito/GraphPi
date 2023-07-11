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
    inter_edge = new int[e_num];
    inter_vertex = new unsigned int[v_num];
    extern_e_max_num = e_num;
    extern_v_max_num = v_num;
}

void Graph::gen_out_of_core_component(int part_num, int block_size, const std::string &path, int k_core) {
    std::vector<int> side_vtx;
    std::map<int,int> part_map;
    for (auto i: v_state_map) {
        unsigned int l,r;
        get_edge_index(i.first,l,r);
    }
    to_partition_csr(part_num,side_vtx,part_map,path);
    to_block_csr(block_size,side_vtx,part_map,part_num,path,k_core);
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
        VertexTable d;
        v_state_map.insert({key, d});
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
    if (blockType == BlockType::CHINK_BFS) {
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
    if (blockType == BlockType::CHINK_BFS) {
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
        VertexTable d;
        v_state_map.insert({key, d});
        v_state_map[key].file_id = value;
        block_size = std::max(block_size,fvalue[i]+1);
    }
    for (int i=0;i<v_cnt;i++) {
        int key = vkey[i];
        int value = vvalue[i];
        v_state_map[key].is_intra = true;
        intra_vertex_dict.insert(std::make_pair(key,value));
    }
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

void Graph::to_block_csr(int block_size, std::vector<int> &side_vertex, std::map<int,int> &part_map, int part_num, const std::string& path, int k_core) {
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
        std::map<int,int> degree;
        std::map<int,int> revert_dict;
        for (auto i : intra_vertex_dict) {
            revert_dict[i.second] = i.first;
        }
        int block_id = 0;
        std::vector<int> block_vid, vid, file_id;
        int local_block_size = 0;
        std::vector<int> to_insert_block;
        for (int i = 0; i < v_cnt; i++) {
            unsigned int l,r;
            get_edge_index(i,l,r);
            degree[i] = r-l;
        }
        for (int core = 1; core <= k_core; core++) {
            bool keep_lop = true;
            while (keep_lop) {
                keep_lop = false;
                std::vector<int> local_remove;
                for (auto i: degree) {
                    if (i.second < core) {
                        for (int j = vertex[i.first]; j < vertex[i.first + 1]; j++) {
                            if (degree.find(edge[j]) != degree.end()) degree[edge[j]] -= 1;
                        }
                        to_insert_block.push_back(i.first);
                        local_remove.push_back(i.first);
                    }
                }
                for (auto i : local_remove) degree.erase(i);
                for (auto i : degree) keep_lop |= i.second < core;
            }
            blocking_data_manage_k_core(len,block_id,local_block_size,block_vid,revert_dict,to_insert_block, path);
        }
        for (auto i : degree) to_insert_block.push_back(i.first);
        blocking_data_manage_k_core(len,block_id,local_block_size,block_vid,revert_dict,to_insert_block, path);

        std::string block_map_path = path + std::string("_fk");
        std::vector<int> fkey,fvalue;
        for (auto i : v_state_map) {
            fkey.push_back(i.first);
            fvalue.push_back(i.second.file_id);
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

    if (blockType == BlockType::CHINK_BFS) {
        unsigned int l,r;
        std::map<int,bool> choose;
        std::queue<int> q_s[part_num];
        std::vector<std::vector<int>> block_vtxs;
        for (auto i: intra_vertex_dict) {
            choose[i.first] = false;
        }
        std::vector<int> same_part_side[part_num];
        for (auto i : side_vertex) {
            int p = part_map[i];
            same_part_side[p].push_back(i);
        }
        //each loop evaluate the connection of each part
        for (int judge_cluster = 0; judge_cluster < part_num; judge_cluster++) {
            int judge_size = same_part_side[judge_cluster].size();
            int *connection_status = new int[judge_size*judge_size];
            for (int i=0;i<judge_size;i++)
                for (int j=0;j<judge_size;j++)
                    connection_status[i*judge_size + j] = 0;
            std::map<int,int> judgeVid_to_idx;
            for (auto judge_vertex : same_part_side[judge_cluster]) judgeVid_to_idx[judge_vertex] = judgeVid_to_idx.size();
            for (int other_part = 0; other_part < part_num; other_part++) {
                if (other_part==judge_cluster) continue;
                for (auto critic_vertex : same_part_side[other_part]) {
                    get_edge_index(critic_vertex,l,r);
                    std::vector<int> judge_idx_find_collecter;
                    for (int i = l; i < r; i++) {
                        if(std::find(judge_idx_find_collecter.begin(), judge_idx_find_collecter.end(),edge[i])!= judge_idx_find_collecter.end())
                            judge_idx_find_collecter.push_back(judgeVid_to_idx[edge[i]]);
                    }

                    std::sort(judge_idx_find_collecter.begin(), judge_idx_find_collecter.end());
                    for (int i = 0; i < judge_idx_find_collecter.size(); i++) {
                        for (int j = i+1; j < judge_idx_find_collecter.size(); j++) {
                            int small_idx = judge_idx_find_collecter[i];
                            int large_idx = judge_idx_find_collecter[j];
                            connection_status[small_idx*judge_size + large_idx] += 1;
                            connection_status[large_idx*judge_size + small_idx] += 1;
                        }
                    }
                }
            }
            // get connect relation in each partition
            // now, get connected component
            std::vector<std::vector<int>> connect_component;
            bool visited[judge_size];
            for (int i = 0; i < judge_size; i++) visited[i] = false;
            std::queue<int> connect_bfs_q;
            for (int i = 0; i < judge_size; i++) {
                if (visited[i]) continue;

                std::vector<int> cur_component;
                connect_bfs_q.push(i);
                visited[i] = true;
                cur_component.push_back(i);

                while (!connect_bfs_q.empty()) {
                    int v0 = connect_bfs_q.front();
                    connect_bfs_q.pop();
                    for (int neigh = 0; neigh < judge_size; neigh++)
                        if (connection_status[v0*judge_size + neigh] > 0)
                            if (!visited[neigh]) {
                                connect_bfs_q.push(neigh);
                                cur_component.push_back(neigh);
                            }
                }
                connect_component.push_back(cur_component);
            }

            // now, in each component, put vertex in block queue in order of importance and connectness
            std::map<int,int> vertex_importance;
            for (int i = 0; i < judge_size; i++) {
                int score = 0;
                for (int j = 0; j < judge_size; j++) {
                    score += connection_status[i*judge_size + j];
                }
                vertex_importance[i] = score;
            }
            bool chosen[judge_size];
            bool neighbored[judge_size];
            for (int i = 0; i < judge_size; i++) {
                chosen[i] = false;
                neighbored[i] = false;
            }
            for (auto component : connect_component) {
                std::vector<int> neighbors;
                // initialize a most importance vertex;
                int max_v = component[0];
                for (auto m : component) {
                    if (vertex_importance[m] > vertex_importance[max_v]) {
                        max_v = m;
                    }
                }
                chosen[max_v] = true;
                //here insert first vertex to queue;
                q_s[judge_cluster].push(same_part_side[judge_cluster][max_v]);
                for (int neigh = 0; neigh < judge_size; neigh++)
                    if (connection_status[max_v*judge_size + neigh] > 0 && !neighbored[neigh]) {
                        neighbors.push_back(neigh);
                        neighbored[neigh] = true;
                    }
                // repeat until all vertex in this component pushed in queue
                while (!neighbors.empty()) {
                    max_v = neighbors[0];

                    for (auto m : neighbors) {

                        if (vertex_importance[m] > vertex_importance[max_v]) {
                            max_v = m;
                        }
                    }
                    chosen[max_v] = true;
                    neighbors.erase(std::find(neighbors.begin(),neighbors.end(),max_v));
                    q_s[judge_cluster].push(same_part_side[judge_cluster][max_v]);
                    for (int neigh = 0; neigh < judge_size; neigh++)
                        if (connection_status[max_v*judge_size + neigh] > 0 && !chosen[max_v] && !neighbored[neigh]) {
                            neighbors.push_back(neigh);
                            neighbored[neigh] = true;
                        }
                }

            }
            delete[] connection_status;
        }
        int cur_block_id = 0;
        int cur_block_size[part_num];
        for (int i=0; i< part_num; i++) cur_block_size[i]=0;
        bool all_queue_empty = false;
        block_vtxs.resize(part_num);
        while (!all_queue_empty) {
            all_queue_empty = true;
            for (int i = 0; i< part_num; i++) {
                if (q_s[i].empty()) continue;
                int vtx = q_s[i].front();
                get_edge_index(vtx,l,r);
                if (cur_block_size[i] + r-l > len && !block_vtxs[i].empty()) {
                    // dump
                    std::sort(block_vtxs[i].begin(), block_vtxs[i].end());
                    std::vector<int> edges;
                    std::vector<unsigned int> vtx_offset;
                    unsigned int cursor = 0;
                    for (auto j : block_vtxs[i]) {
                        vtx_offset.push_back(cursor);
                        get_edge_index(j,l,r);
                        for (int k = l; k < r; k++) edges.push_back(edge[k]);
                        cursor += r-l;
                        v_state_map[j].file_id = cur_block_id;
                    }

                    std::string block_csr_path = path + std::string("_blocks/") + std::to_string(cur_block_id) + std::string("_bc");
                    DataLoader::gen_block_file(block_vtxs[i].data(),vtx_offset.data(),edges.data(),block_vtxs[i].size(),edges.size(),block_csr_path);
                    block_lengths[cur_block_id] = block_length{(int)block_vtxs[i].size(),(int)edges.size()};

                    cur_block_id += 1;
                    cur_block_size[i] = 0;
                    block_vtxs[i].clear();
                    all_queue_empty &= q_s[i].empty();
                    continue;
                }
                if (part_map[vtx] == i) {
                    block_vtxs[i].push_back(vtx);
                    cur_block_size[i] += 2;
                    get_edge_index(vtx, l, r);
                    for (int j = l; j < r; j++) {
                        if (choose[edge[j]]) continue;
                        q_s[i].push(edge[j]);
                        choose[edge[j]] = true;
                    }
                    cur_block_size[i] += r - l;
                }
                q_s[i].pop();
                all_queue_empty &= q_s[i].empty();
            }
        }
        for (int i = 0; i< part_num; i++) {
            if (block_vtxs[i].empty()) continue;
            //dump

            std::sort(block_vtxs[i].begin(), block_vtxs[i].end());
            std::vector<int> edges;
            std::vector<unsigned int> vtx_offset;
            unsigned int cursor = 0;
            for (auto j : block_vtxs[i]) {
                vtx_offset.push_back(cursor);
                get_edge_index(j,l,r);
                for (int k = l; k < r; k++) edges.push_back(edge[k]);
                cursor += r-l;
                v_state_map[j].file_id = cur_block_id;
            }
            std::string block_csr_path = path + std::string("_blocks/") + std::to_string(cur_block_id) + std::string("_bc");
            DataLoader::gen_block_file(block_vtxs[i].data(),vtx_offset.data(),edges.data(),block_vtxs[i].size(),edges.size(),block_csr_path);
            block_lengths[cur_block_id] = block_length{(int)block_vtxs[i].size(),(int)edges.size()};

            cur_block_id += 1;
            cur_block_size[i] = 0;
            block_vtxs[i].clear();
        }
        for (int k = 0; k < part_num; k++) {
            int b_size = 0;
            std::vector<int> b_vtx;
            for (auto i: intra_vertex_dict) {
                if (part_map[i.first] != k) continue;
                if (choose[i.first]) continue;
                choose[i.first] = true;
                int vtx = i.first;
                get_edge_index(vtx, l, r);
                if (b_size + r - l > len && !b_vtx.empty()) {
                    // dump
                    std::sort(b_vtx.begin(), b_vtx.end());
                    std::vector<int> edges;
                    std::vector<unsigned int> vtx_offset;
                    unsigned int cursor = 0;
                    for (auto j: b_vtx) {
                        vtx_offset.push_back(cursor);
                        get_edge_index(j, l, r);
                        for (int k = l; k < r; k++) edges.push_back(edge[k]);
                        cursor += r - l;

                        v_state_map[j].file_id = cur_block_id;
                    }
                    std::string block_csr_path =
                            path + std::string("_blocks/") + std::to_string(cur_block_id) + std::string("_bc");
                    DataLoader::gen_block_file(b_vtx.data(), vtx_offset.data(), edges.data(), b_vtx.size(),
                                               edges.size(), block_csr_path);
                    block_lengths[cur_block_id] = block_length{(int) b_vtx.size(), (int) edges.size()};

                    cur_block_id += 1;
                    b_size = 0;
                    b_vtx.clear();
                    continue;
                }
                b_vtx.push_back(vtx);
                b_size += 1;
            }
        }

        std::string block_map_path = path + std::string("_fc");
        std::vector<int> fkey,fvalue;
        for (auto i : v_state_map) {
            fkey.push_back(i.first);
            fvalue.push_back(i.second.file_id);
        }
        DataLoader::gen_map_file(fkey.data(),fvalue.data(),fkey.size(),block_map_path);

        std::string block_size_path = path + std::string("_fc_size");
        std::vector<int> key,vl,el;
        for (auto i : block_lengths) {
            key.push_back(i.first);
            vl.push_back(i.second.v_len);
            el.push_back(i.second.e_len);
        }
        DataLoader::gen_block_size_file(key.data(),vl.data(),el.data(),key.size(),block_size_path);
    }

    if (blockType == BlockType::SIMPLE_BFS) {
        int block_id = 0;
        bool choose[v_cnt];
        for (int i = 0; i < v_cnt; i++) choose[i] = false;
        for (int part = 0; part < part_num; part++) {
            std::queue<int> q;
            std::vector<int> block_vid;
            int block_vsize = 0;
            for (auto v : v_state_map) {
                if (part_map[v.first] != part || choose[v.first]) continue;
                q.push(v.first);
                choose[v.first] = true;
                while (!q.empty()) {
                    if (block_vsize > len) {
                        block_vsize = 0;
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
                        block_id+=1;
                        continue;
                    }
                    int vtx = q.front();
                    q.pop();
                    block_vid.push_back(vtx);
                    unsigned int l,r;
                    get_edge_index(vtx,l,r);
                    for (int i = l; i < r; i++) {
                        if (choose[edge[i]] || part_map[edge[i]] != part_map[vtx]) continue;
                        q.push(edge[i]);
                        choose[vtx] = true;
                        block_vsize += 2+r-l;
                    }
                }
            }

            block_vsize = 0;
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
            block_id+=1;
        }

        std::string block_map_path = path + std::string("_fs");
        std::vector<int> fkey,fvalue;
        for (auto i : v_state_map) {
            fkey.push_back(i.first);
            fvalue.push_back(i.second.file_id);
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

void Graph::blocking_data_manage_k_core(int len, int &block_id, int &local_block_size, std::vector<int> &block_vid, const std::map<int, int> &revert_dict, std::vector<int> &to_insert_block, const std::string& path) {
    std::queue<int> q;
    while (!to_insert_block.empty() || !q.empty()) {
        if (q.empty()) {
            q.push(to_insert_block.back());
            to_insert_block.pop_back();
        }
        int q_len;
        if (q.front() == v_cnt -1) q_len = e_cnt-vertex[q.front()];
        else q_len = vertex[q.front()+1] - vertex[q.front()];
        if (local_block_size + 2 + q_len > len && !block_vid.empty()) {
            std::sort(block_vid.begin(), block_vid.end());
            int cursor = 0;
            std::vector<unsigned int> vtx_offset;
            std::vector<int> edges;
            vtx_offset.push_back(cursor);
            for (int i=0; i< block_vid.size();i++) {
                int elen;
                if (intra_vertex_dict[block_vid[i]] == v_cnt-1) elen = e_cnt - vertex[intra_vertex_dict[block_vid[i]]];
                else elen =vertex[intra_vertex_dict[block_vid[i]]+1] - vertex[intra_vertex_dict[block_vid[i]]];
                cursor += elen;
                for (int j = 0; j < elen; j++) edges.push_back(edge[j+vertex[intra_vertex_dict[block_vid[i]]]]);
                // the last one do not need to record its ending index.
                if (i< block_vid.size()-1)
                    vtx_offset.push_back(cursor);
            }
            block_lengths[block_id] = block_length{(int) block_vid.size(),(int) edges.size()};
            std::string block_csr_path = path + std::string("_blocks/") + std::to_string(block_id) + std::string("_bk");
            DataLoader::gen_block_file(block_vid.data(),vtx_offset.data(),edges.data(),block_vid.size(),edges.size(),block_csr_path);

            for (auto i : block_vid) {
                v_state_map[i].file_id = block_id;
            }

            block_id+=1;
            local_block_size = 0;
            block_vid.clear();
        }
        int insert_v = q.front();
        block_vid.push_back(revert_dict.at(insert_v));
        local_block_size += 2+ q_len;
        q.pop();
        for (int j = 0; j < q_len; j++) {
            int target_vtx_offset = intra_vertex_dict[vertex[insert_v]+j];
            if (std::find(to_insert_block.begin(), to_insert_block.end(), target_vtx_offset)!= to_insert_block.end()) {
                q.push(target_vtx_offset);
                to_insert_block.erase(std::remove(to_insert_block.begin(),to_insert_block.end(), target_vtx_offset), to_insert_block.end());
            }
        }
        if (q.empty() && to_insert_block.empty()) {
            std::sort(block_vid.begin(), block_vid.end());
            int cursor = 0;
            std::vector<unsigned int> vtx_offset;
            std::vector<int> edges;
            vtx_offset.push_back(cursor);
            for (int i=0; i< block_vid.size();i++) {
                int elen;
                if (intra_vertex_dict[block_vid[i]] == v_cnt-1) elen = e_cnt - vertex[intra_vertex_dict[block_vid[i]]];
                else elen =vertex[intra_vertex_dict[block_vid[i]]+1] - vertex[intra_vertex_dict[block_vid[i]]];
                cursor += elen;
                for (int j = 0; j < elen; j++) edges.push_back(edge[j+vertex[intra_vertex_dict[block_vid[i]]]]);
                // the last one do not need to record its ending index.
                if (i< block_vid.size()-1)
                    vtx_offset.push_back(cursor);
            }
            block_lengths[block_id] = block_length{(int) block_vid.size(),(int) edges.size()};
            std::string block_csr_path = path + std::string("_blocks/") + std::to_string(block_id) + std::string("_bk");
            DataLoader::gen_block_file(block_vid.data(),vtx_offset.data(),edges.data(),block_vid.size(),edges.size(),block_csr_path);
            for (auto i : block_vid) {
                v_state_map[i].file_id = block_id;
            }
            block_id+=1;
            local_block_size = 0;
            block_vid.clear();
        }

    }
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

void Graph::to_partition_csr(int num, std::vector<int> &side_vertex, std::map<int,int> &part_map, const std::string& path) {
    if (partitionType == PartitionType::RANDOM) {
        for (int p = 0; p<num; p++) {
            int partial_vcnt;
            partial_vcnt = v_cnt /num;
            if (p < v_cnt % num) partial_vcnt += 1;
            auto *part_vertex = new unsigned int[partial_vcnt];
            int mapk[partial_vcnt];
            int mapv[partial_vcnt];
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
        for (auto i : v_state_map) {
            is_side[i.first] = false;
            choose[i.first] = false;
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
        int capacity = e_cnt / num;
        std::vector<int> part_vertices[num];
        int p_size[num];
        for (int i = 0; i < num; i++) p_size[i] = 0;
        int d[v_cnt],v[v_cnt];
        for (int i = 0; i < v_cnt; i++) {
            unsigned int l,r;
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
        for (int i = 0; i < v_cnt; i++) {
            int vtx= v[i];
            unsigned int l,r;
            get_edge_index(vtx,l,r);
            double score = 0.0;
            int max = 0;
            for (int part = 0; part < num; part++) {
                double w_i = 1.0 - (double ) p_size[part] / capacity;
                double new_score = w_i / (double )(r-l) * (double )get_intersect_num(edge+l,part_vertices[part].data(),r-l,part_vertices[part].size());
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
            std::sort(part_vertex.begin(), part_vertex.end());
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

        for (auto v : intra_vertex_dict) {
            unsigned int l,r;
            get_edge_index(v.first,l,r);
            for (int i = l; i < r; i++) {
                if (part_map[edge[i]] != part_map[v.first]) {
                    side_vertex.push_back(v.first);
                    break;
                }
            }
        }
    }
}
