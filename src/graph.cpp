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


static omp_lock_t load_list_write_lock;
static omp_lock_t load_list_append_lock;
static omp_lock_t init_load_list_append_lock;
static omp_lock_t extern_ve_loading_lock;
static omp_lock_t ready_bin_lock;
static omp_lock_t candidate_bin_lock;
static omp_lock_t v_state_lock;
static omp_lock_t v_state_khop_lock;

int Graph::intersection_size(int v1,int v2) {
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    unsigned int l2, r2;
    get_edge_index(v2, l2, r2);
    int ans = 0;
    while(l1 < r1 && l2 < r2) {
        if(edge[l1] < edge[l2]) {
            ++l1;
        }
        else {
            if(edge[l2] < edge[l1]) {
                ++l2;
            }
            else {
                ++l1;
                ++l2;
                ++ans;
            }
        }
    }
    return ans;
}

int Graph::intersection_size_clique(int v1,int v2) {
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    unsigned int l2, r2;
    get_edge_index(v2, l2, r2);
    int min_vertex = v2;
    int ans = 0;
    if (edge[l1] >= min_vertex || edge[l2] >= min_vertex)
        return 0;
    while(l1 < r1 && l2 < r2) {
        if(edge[l1] < edge[l2]) {
            if (edge[++l1] >= min_vertex)
                break;
        }
        else {
            if(edge[l2] < edge[l1]) {
                if (edge[++l2] >= min_vertex)
                    break;
            }
            else {
                ++ans;
                if (edge[++l1] >= min_vertex)
                    break;
                if (edge[++l2] >= min_vertex)
                    break;
            }
        }
    }
    return ans;
}

long long Graph::triangle_counting() {
    long long ans = 0;
    for(int v = 0; v < v_cnt; ++v) {
        // for v in G
        unsigned int l, r;
        get_edge_index(v, l, r);
        for(unsigned int v1 = l; v1 < r; ++v1) {
            //for v1 in N(v)
            ans += intersection_size(v,edge[v1]);
        }
    }
    ans /= 6;
    return ans;
}

long long Graph::triangle_counting_mt(int thread_count) {
    long long ans = 0;
#pragma omp parallel num_threads(thread_count)
    {
        tc_mt(&ans);
    }
    return ans;
}

void Graph::tc_mt(long long *global_ans) {
    long long my_ans = 0;
#pragma omp for schedule(dynamic)
    for(int v = 0; v < v_cnt; ++v) {
        // for v in G
        unsigned int l, r;
        get_edge_index(v, l, r);
        for(unsigned int v1 = l; v1 < r; ++v1) {
            if (v <= edge[v1])
                break;
            //for v1 in N(v)
            my_ans += intersection_size_clique(v,edge[v1]);
        }
    }
#pragma omp critical
    {
        *global_ans += my_ans;
    }
}


struct vertex_degree_cmp {
    bool operator() (vertex_degree a, vertex_degree b) {return a.d > b.d;}
} vertexDegreeCmp;


struct loader_cmp {
    bool operator() (loader a, loader b) {return a.vertex< b.vertex;}
} loaderCmp;

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


void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    int vtx = intra_vertex_dict.at(v);
    l = vertex[vtx];
    if (vtx == v_cnt -1) r = e_cnt;
    else r = vertex[vtx + 1];
}

void Graph::get_extern_edge_index(int v, unsigned int& l, unsigned int& r) const
{

    omp_set_lock(&extern_ve_loading_lock);
    int vtx = inter_vertex_dict.at(v);
    omp_unset_lock(&extern_ve_loading_lock);
    l = inter_vertex[vtx];
    if (vtx == inter_vtx_num - 1) r = inter_edge_num;
    else r = inter_vertex[vtx + 1];
}

void Graph::load_extern_data(int file_id, const std::vector<loader>& load_vertex, std::vector<int> &loaded_vertex) {
    io_num+=1;
    load_num += load_vertex.size();
    //read file and get N(v)
    // fread, scanf, getline..
    int *vid = nullptr;
    unsigned int *vertex_offset = nullptr;
    int *adj_list = nullptr;

    std::string block_path;
    if (blockType == BlockType::K_CORE_BLOCK) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bk");
    else if (blockType == BlockType::RANDOM_BLOCK) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_br");
    else if (blockType == BlockType::CHINK_BFS) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bc");
    else if (blockType == BlockType::SIMPLE_BFS) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bs");
    int load_v_len=block_lengths[file_id].v_len;
    unsigned int load_e_len=block_lengths[file_id].e_len;
//    DataLoader::load_data_size(load_v_len,load_e_len,block_path);
    vid = new int[load_v_len];
    vertex_offset = new unsigned int[load_v_len];
    adj_list = new int[load_e_len];
    int* data = new int[2*load_v_len+load_e_len];

    double t1 = get_wall_time();
    // io reading
    DataLoader::load_block_data_aggregate(data,2*load_v_len+load_e_len,block_path);
    memcpy(vid,data,load_v_len*sizeof(int));
    memcpy(vertex_offset,data+load_v_len,load_v_len*sizeof(int));
    memcpy(adj_list,data+2*load_v_len,load_e_len* sizeof(int));
//    DataLoader::load_block_data(vid,vertex_offset,adj_list, load_v_len, load_e_len, block_path);
    delete[] data;
    double t2 = get_wall_time();

    file_time += t2-t1;
    t1 = get_wall_time();
    int vid_cursor = 0;
    if (loadType == SUB_BLOCK){
        for (auto i: load_vertex) {
            // find vertex
            while (vid[vid_cursor] != i.vertex && vid_cursor < load_v_len) {
                vid_cursor++;
                assert(vid[vid_cursor] <= i.vertex);
                if (vid[vid_cursor] == i.vertex) break;
            }
            // add this vertex to list
            // TODO: when extern storage is not enough how to solve it? or avoid loading too much before.

            omp_set_lock(&extern_ve_loading_lock);
            inter_vertex_dict[i.vertex] = inter_vtx_num;
            inter_vertex[inter_vtx_num] = inter_edge_num;
            unsigned int edge_cursor = inter_edge_num;
            unsigned int len;
            if (vid_cursor < load_v_len - 1) len = vertex_offset[vid_cursor + 1] - vertex_offset[vid_cursor];
            else len = load_e_len - vertex_offset[vid_cursor];
#pragma omp atomic
            inter_edge_num += len;
#pragma omp atomic
            inter_vtx_num += 1;
            for (int j = 0; j < len; j++) {
                inter_edge[edge_cursor] = adj_list[vertex_offset[vid_cursor] + j];
                edge_cursor += 1;
            }

            omp_unset_lock(&extern_ve_loading_lock);
            // release lock
            if (i.k_hop == 1) {

                omp_set_lock(&extern_ve_loading_lock);
                omp_set_lock(&v_state_khop_lock);
                v_state_map[i.vertex].k_hop = 1;
                omp_unset_lock(&v_state_khop_lock);
                omp_unset_lock(&extern_ve_loading_lock);
            } else if (i.k_hop > 1) {
                loaded_vertex.push_back(i.vertex);
                for (int j = 0; j < len; j++) {
                    omp_set_lock(&extern_ve_loading_lock);
                    omp_set_lock(&v_state_khop_lock);
                    int origin_hop = v_state_map[i.vertex].k_hop;
                    v_state_map[i.vertex].k_hop = std::max(1,origin_hop);
                    omp_unset_lock(&v_state_khop_lock);
                    omp_unset_lock(&extern_ve_loading_lock);
                    int v = adj_list[vertex_offset[vid_cursor] + j];
                    //here add recursive function (to achieve, just need append load_list, as loading_manage will keep calling this function until load_list becomes empty)
                    if (!v_state_map[v].is_intra) {

                        load_list_append(loader{i.k_hop - 1, v, v_state_map[v].file_id}, loaded_vertex);
                    }
                }
            }
        }
    }
    delete[] vid;
    delete[] vertex_offset;
    delete[] adj_list;

    t2 = get_wall_time();
    blocking_manage_time += t2-t1;

}

void Graph::loading_manage() {
    std::vector<int> loaded_vtx;

    omp_set_lock(&init_load_list_append_lock);
    std::vector<loader> ini_loader;
//    printf("%d %d %d\n", omp_get_thread_num(),load_list.size(), init_load_list.size());
    ini_loader.swap(init_load_list);
    omp_unset_lock(&init_load_list_append_lock);
    if (omp_test_lock(&load_list_append_lock)) {
        if (load_list.empty()) {
            for (auto load_element: ini_loader) {
                load_list_append(load_element, loaded_vtx);
            }
        }
        omp_unset_lock(&load_list_append_lock);
    }
    while (!load_list.empty()) {
        omp_set_lock(&load_list_write_lock);
        std::vector<loader> load_vertex;

        int load_file_size = 0;
        int file_id;
        for (auto i : load_list) {
            if (i.second.size() > load_file_size) {
                load_file_size = i.second.size();
                file_id = i.first;
            }
        }
        load_vertex.swap(load_list[file_id]);
        load_list.erase(file_id);
        omp_unset_lock(&load_list_write_lock);

        //release lock
        std::sort(load_vertex.begin(), load_vertex.end(),loaderCmp);

        double t1 = get_wall_time();
        load_extern_data(file_id, load_vertex, loaded_vtx);
        double t2 = get_wall_time();
        load_extern_time += t2-t1;

    }
    for (auto i : loaded_vtx) {
        if (v_state_map[i].loading != 0) {
            int k_hop = v_state_map[i].loading;
            omp_set_lock(&extern_ve_loading_lock);
            omp_set_lock(&v_state_khop_lock);
            v_state_map[i].k_hop = std::max((unsigned int) v_state_map[i].loading,v_state_map[i].k_hop);
            omp_unset_lock(&v_state_khop_lock);
            omp_unset_lock(&extern_ve_loading_lock);
//            v_state_map[i].loading = 0;
        }
    }
}

void Graph::load_list_append(const loader &l, std::vector<int> &loaded_vertex){
    // if it is not inter-partition vertex, no considering to load k-hop

    if (v_state_map[l.vertex].is_intra) return;
    if (l.k_hop==0) return;
    // when there already exist a loading process that contain higher k-hop, or current k-hop is larger, no need to append to task queue
    if ( v_state_map[l.vertex].loading >= l.k_hop || v_state_map[l.vertex].k_hop >= l.k_hop) return;

    omp_set_lock(&load_list_write_lock);
    v_state_map[l.vertex].loading = std::max((int) l.k_hop,v_state_map[l.vertex].loading);
    omp_unset_lock(&load_list_write_lock);
    int curr_hop = v_state_map[l.vertex].k_hop;
    if (v_state_map[l.vertex].k_hop >= 1) {
        unsigned int i,j;
        get_extern_edge_index(l.vertex,i,j);
        for (int k = i; k<j; k++) {
            int v = inter_edge[k];

            load_list_append(loader{l.k_hop-1, v, v_state_map[v].file_id},loaded_vertex);
        }
        loaded_vertex.push_back(l.vertex);
        return;
    }

    omp_set_lock(&load_list_write_lock);
    if (load_list.find(l.file_id) == load_list.end()) {
        std::vector<loader> vec;
        vec.push_back(l);
        load_list[l.file_id] = vec;
    }
    else {
        load_list[l.file_id].push_back(l);
    }
    omp_unset_lock(&load_list_write_lock);
}

bool Graph::extern_store_manage(const Schedule& schedule) {
    // first calculate the available space
    std::vector<int> k_hop_matrix = schedule.k_hop_matrix;
    float extern_usage = external_used / external_space;
    // second base on the available determine if drop the data
    if (extern_usage > extern_upper_thresold){
        printf("do drop\n");
        // to avoiding load extern info drop even when it is not used, it will only stop loading if ready_bin is not empty
        if (ready_bin.empty()) extern_drop_manage();
        return false;
    }
    // third base on the available space, determine how much vertex will be loaded.
    unsigned int remain_space = external_space - external_used;
    double t1 = get_wall_time();
    extern_load_init(remain_space, schedule);
    double t2 = get_wall_time();
    extern_init_time += t2-t1;
    // only this function is allocate to multi thread
    t1 = get_wall_time();
    loading_manage();
    t2 = get_wall_time();
    loading_manage_time += t2-t1;

    omp_set_lock(&candidate_bin_lock);
    omp_set_lock(&ready_bin_lock);
    for (int i=candidate_bin.size()-1; i >=0 ; i--) {
        int k_hop = k_hop_matrix[candidate_bin[i]->depth];
        bool vector_load = true;
        for (int j = 0; j < candidate_bin[i]->v_size; j++) {
            int v = candidate_bin[i]->vertex[j];

            vector_load &= v_state_map[v].k_hop >= k_hop;
            if (!vector_load) break;
        }
        if (vector_load) {
            ready_bin.push_back(candidate_bin[i]);
            candidate_bin.erase(candidate_bin.begin()+i);
        }
    }
    omp_unset_lock(&ready_bin_lock);
    omp_unset_lock(&candidate_bin_lock);
    return true;
}

// TODO: add todo and space manage
void Graph::extern_drop_manage() {
    omp_set_lock(&extern_ve_loading_lock);
    for (auto it : inter_vertex_dict) {
#pragma omp atomic
        v_state_map[it.first].k_hop *= 0;
#pragma omp atomic
        v_state_map[it.first].loading *= 0;
    }
    inter_vertex_dict.clear();
    printf("do clear\n");
#pragma omp atomic
    inter_edge_num *= 0;

#pragma omp atomic
    inter_vtx_num *= 0;
    omp_unset_lock(&extern_ve_loading_lock);

}

void Graph::extern_load_init(unsigned int available_space, const Schedule& schedule) {
    // manage candidate_bin to to_load set\

    omp_set_lock(&candidate_bin_lock);
    std::vector<path *> candidates(candidate_bin);
    omp_unset_lock(&candidate_bin_lock);
    for (auto & i : candidates) {
        if (i->load_queue) continue;
        i->load_queue = true;
        int k_hop = schedule.k_hop_matrix[i->depth];
        for (int j = 0; j < i->v_size; j++) {
            int v = i->vertex[j];
            omp_set_lock(&init_load_list_append_lock);

            if (!v_state_map[v].is_intra && v_state_map[v].k_hop < k_hop) init_load_list.push_back(loader{k_hop, v, v_state_map[v].file_id});
            omp_unset_lock(&init_load_list_append_lock);
        }
    }
}

void Graph::pattern_matching_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, long long& local_ans, int depth, bool clique)
{
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;
    int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    if (depth == schedule.get_size() - 1)
    {
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (clique)
            local_ans += loop_size;
        else if (loop_size > 0)
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    int last_vertex = subtraction_set.get_last();
    std::vector<int> to_load_vertex;
    for (int i = 0; i < loop_size; ++i)
    {
        if (last_vertex <= loop_data_ptr[i] && clique)
            break;
        int load_v = loop_data_ptr[i];
        if (subtraction_set.has_data(load_v))
            continue;
        unsigned int l, r;

        if (v_state_map[load_v].is_intra){
            get_edge_index(load_v, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, load_v);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) continue;
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(load_v);
            pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@
        }
        else if (v_state_map[load_v].k_hop > 0) {
            get_extern_edge_index(load_v, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &inter_edge[l], (int)r - l, prefix_id, load_v);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) continue;
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(load_v);
            pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@
        }
        else {
            to_load_vertex.push_back(load_v);
            continue;
        }
    }
    if (!to_load_vertex.empty()){
        path *p = new path;
        p->depth = depth;
        p->load_queue = false;
        p->vertex = new int[to_load_vertex.size()];
        p->v_size = to_load_vertex.size();
        for (int i = 0; i< to_load_vertex.size(); i++) p->vertex[i] = to_load_vertex[i];
        // deep copy
        p->subtraction_set = new VertexSet;
        p->subtraction_set->deepcopy(subtraction_set);
        // deep copy vertex_set
        // will it be too large to save?
        p->vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        for (int i=0; i< schedule.get_total_prefix_num(); i++) {
            p->vertex_set[i].deepcopy(vertex_set[i]);
        }
        candidate_bin.push_back(p);
    }
}

void Graph::resume_matching(const Schedule& schedule, path *p, long long &local_ans) {
    VertexSet tmp_set;
    unsigned int l,r;
    int depth = p->depth;
    for (int i=0; i< p->v_size;i++) {
        int v = p->vertex[i];
        get_extern_edge_index(v, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            p->vertex_set[prefix_id].build_vertex_set(schedule, p->vertex_set, &inter_edge[l], (int)r - l, prefix_id, v);
            if( p->vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        p->subtraction_set->push_back(v);
        pattern_matching_aggressive_func(schedule,p->vertex_set,*p->subtraction_set,tmp_set,local_ans,depth+1);
//        pattern_matching_func(schedule,p->vertex_set,*p->subtraction_set,local_ans,depth+1);
        p->subtraction_set->pop_back();
    }
}

long long Graph::pattern_matching(const Schedule& schedule, int thread_count, bool clique)
{
    load_num = 0;
//    io_num = 0;
    for(auto i : v_state_map) v_state_map[i.first].is_rooted = false;
    long long global_ans = 0;
    std::vector<int> k_hop_table = schedule.k_hop_matrix;
    omp_init_lock(&ready_bin_lock);
    omp_init_lock(&candidate_bin_lock);
    omp_init_lock(&extern_ve_loading_lock);
    omp_init_lock(&load_list_write_lock);
    omp_init_lock(&load_list_append_lock);
    omp_init_lock(&init_load_list_append_lock);

    omp_init_lock(&v_state_lock);

#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        long long local_ans = 0;
        int v_num = 0;
        // TODO : try different chunksize
#pragma omp for schedule(dynamic) nowait
        for (int vertex_id = 0; vertex_id < v_state_map.size(); vertex_id++)
        {
            // translate it into vid

            omp_set_lock(&v_state_lock);
            if (!v_state_map[vertex_id].is_intra) {
                omp_unset_lock(&v_state_lock);
                continue;
            }
            if (v_state_map[vertex_id].is_rooted) {
                omp_unset_lock(&v_state_lock);
                continue;
            }

            std::queue<int> bfs_queue;
            bfs_queue.push(vertex_id);
            v_state_map[vertex_id].is_rooted = true;
            omp_unset_lock(&v_state_lock);

            while (!bfs_queue.empty()) {
                // clear remain task in the ready bin is the highest priority
                while (true) {
                    // Highest priority to allocate task from candidate bin to each thread.
                    // do remain work in ready bin
                    omp_set_lock(&ready_bin_lock);
                    if (ready_bin.empty()) {
                        omp_unset_lock(&ready_bin_lock);
                        break;
                    }
                    else {
                        path *p = ready_bin.back();
                        ready_bin.pop_back();
                        omp_unset_lock(&ready_bin_lock);

                        resume_matching(schedule, p, local_ans);
                        delete p;
                        p = nullptr;
                    }
                }
                int vtx = bfs_queue.front();
                bfs_queue.pop();
                v_num++;
                unsigned int l, r;
                get_edge_index(vtx, l, r);
                for (int i=l; i<r;i++) {
                    int enqueue_vtx = edge[i];

                    // for vertex that do not belong to this partition, it will not be selected as a root
                    omp_set_lock(&v_state_lock);
                    if (!v_state_map[enqueue_vtx].is_intra) {
                        omp_unset_lock(&v_state_lock);
                        continue;
                    }
                    if (v_state_map[enqueue_vtx].is_rooted) {
                        omp_unset_lock(&v_state_lock);
                        continue;
                    }
                    bfs_queue.push(enqueue_vtx);
                    v_state_map[enqueue_vtx].is_rooted = true;
                    omp_unset_lock(&v_state_lock);
                }
                for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
                {
                    // get N(v)
                    // As root must be intra-partition, external_edge will not be used
                    vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id);
                }
                subtraction_set.push_back(vtx);
                pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
//                pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, 1);

                subtraction_set.pop_back();
                double t1 = get_wall_time();
                // TODO: add condition for scattering thread for data_managing
                int to_load_num = 0;


                omp_set_lock(&candidate_bin_lock);
                for (auto i : candidate_bin)
                    to_load_num += i->v_size;
                omp_unset_lock(&candidate_bin_lock);

                int id = omp_get_thread_num();
                if (scatter_loading_thread(thread_count,id,to_load_num)) {
                    extern_store_manage(schedule);
                }
                double t2 = get_wall_time();
                external_time += t2-t1;
            }
        }
        double t1 = get_wall_time();
        int id = omp_get_thread_num();
        if (id ==0) {
            while (!candidate_bin.empty()) {
                extern_store_manage(schedule);
            }
        }
        double t2 = get_wall_time();
        external_time += t2-t1;

        while (true) {
            omp_set_lock(&ready_bin_lock);
            if (ready_bin.empty()) {
                omp_unset_lock(&ready_bin_lock);
                break;
            }
            path *p = ready_bin.back();
            ready_bin.pop_back();
            omp_unset_lock(&ready_bin_lock);

            resume_matching(schedule,p,local_ans);
            delete p;
            p = nullptr;
        }
#pragma omp barrier
#pragma omp for schedule(dynamic) nowait
        for (int i=0; i< ready_bin.size();i++) {
            // do remain work in ready bin
            path *p = ready_bin[i];
            resume_matching(schedule,p,local_ans);
            delete p;
            p = nullptr;
        }
//        omp_set_lock(&ready_bin_lock);
//        ready_bin.clear();
//        omp_unset_lock(&ready_bin_lock);
        delete[] vertex_set;
#pragma omp atomic
        global_ans += local_ans;
    }

    omp_destroy_lock(&ready_bin_lock);
    omp_destroy_lock(&candidate_bin_lock);
    omp_destroy_lock(&extern_ve_loading_lock);
    omp_destroy_lock(&load_list_write_lock);
    omp_destroy_lock(&load_list_append_lock);
    omp_destroy_lock(&init_load_list_append_lock);
    omp_destroy_lock(&v_state_lock);
    printf("\n=====log for %d=====\n",blockType);
    printf("\nio num: %d\n", io_num);
    printf("extra vtx num: %d\n", load_num);
    printf("      file read io: %.6lf\n", file_time*100);
    printf("      blocking manage time %.6lf\n",blocking_manage_time*100);
    printf("    load extern io: %.6lf\n", load_extern_time*100);
    printf("  loading manage time %.6lf\n",loading_manage_time*100);
    printf("  extern init time %.6lf\n",extern_init_time*100);
    printf("total load: %.6lf\n", external_time*100);
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

void Graph::pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth) // 3 same # or @ in comment are useful in code generation ###
{
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);// @@@
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;
    int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    //Case: in_exclusion_optimize_num > 1
    if( depth == schedule.get_size() - schedule.get_in_exclusion_optimize_num() ) {
        int in_exclusion_optimize_num = schedule.get_in_exclusion_optimize_num();// @@@
        int loop_set_prefix_ids[ in_exclusion_optimize_num ];
        loop_set_prefix_ids[0] = loop_set_prefix_id;
        for(int i = 1; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id( depth + i );
        for(int optimize_rank = 0; optimize_rank < schedule.in_exclusion_optimize_group.size(); ++optimize_rank) {
            const std::vector< std::vector<int> >& cur_graph = schedule.in_exclusion_optimize_group[optimize_rank];
            long long val = schedule.in_exclusion_optimize_val[optimize_rank];
            for(int cur_graph_rank = 0; cur_graph_rank < cur_graph.size(); ++ cur_graph_rank) {
                //                VertexSet tmp_set;

                //if size == 1 , we will not call intersection(...)
                //so we will not allocate memory for data
                //otherwise, we need to copy the data to do intersection(...)
                if(cur_graph[cur_graph_rank].size() == 1) {
                    int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    val = val * VertexSet::unorderd_subtraction_size(vertex_set[id], subtraction_set);
                }
                else {
                    int id0 = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    int id1 = loop_set_prefix_ids[cur_graph[cur_graph_rank][1]];
                    tmp_set.init(this->max_degree);
                    tmp_set.intersection(vertex_set[id0], vertex_set[id1]);

                    for(int i = 2; i < cur_graph[cur_graph_rank].size(); ++i) {
                        int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][i]];
                        tmp_set.intersection_with(vertex_set[id]);
                    }
                    val = val * VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
                }
                if( val == 0 ) break;

            }
            local_ans += val;
        }
        return;// @@@
    }
    //Case: in_exclusion_optimize_num <= 1
    if (depth == schedule.get_size() - 1)
    {
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0)
        {
            int min_vertex = g_vcnt;
            for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
                if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
            const VertexSet& vset = vertex_set[loop_set_prefix_id];
            int size_after_restrict = std::lower_bound(vset.get_data_ptr(), vset.get_data_ptr() + vset.get_size(), min_vertex) - vset.get_data_ptr();
            if (size_after_restrict > 0)
                local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set, size_after_restrict);
        }
        else
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;// @@@
    }

    int min_vertex = g_vcnt;
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
//    if (depth == 1) Graphmpi::getinstance().get_loop(loop_data_ptr, loop_size);
    int ii = 0;
    std::vector<int> to_load_vertex;
    for (int &i = ii; i < loop_size; ++i)
    {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int load_v = loop_data_ptr[i];
        if (subtraction_set.has_data(load_v))
            continue;
        unsigned int l, r;

        omp_set_lock(&extern_ve_loading_lock);
        omp_set_lock(&v_state_khop_lock);
        int cur_khop = v_state_map[load_v].k_hop;
        omp_unset_lock(&v_state_khop_lock);
        omp_unset_lock(&extern_ve_loading_lock);
        if (v_state_map[load_v].is_intra){
            get_edge_index(load_v, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, load_v);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) continue;
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(load_v);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@

            subtraction_set.pop_back(); // @@@
        }
        else if (cur_khop > 0) {
            get_extern_edge_index(load_v, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &inter_edge[l], (int)r - l, prefix_id, load_v);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) continue;
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(load_v);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@
        }
        else {
            to_load_vertex.push_back(load_v);
            continue;
        }
    }

    if (!to_load_vertex.empty()){
        path *p = new path;
        p->depth = depth;
        p->load_queue = false;
        p->vertex = new int[to_load_vertex.size()];
        p->v_size = to_load_vertex.size();
        for (int i = 0; i< to_load_vertex.size(); i++) p->vertex[i] = to_load_vertex[i];
        // deep copy
        p->subtraction_set = new VertexSet;
        p->subtraction_set->deepcopy(subtraction_set);
        // deep copy vertex_set
        // will it be too large to save?
        p->vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        for (int i=0; i< schedule.get_total_prefix_num(); i++) {
            p->vertex_set[i].deepcopy(vertex_set[i]);
        }
        omp_set_lock(&candidate_bin_lock);
        candidate_bin.push_back(p);
        omp_unset_lock(&candidate_bin_lock);
    }
}
// ###
long long Graph::pattern_matching_mpi(const Schedule& schedule, int thread_count, bool clique)
{
    Graphmpi &gm = Graphmpi::getinstance();
    long long global_ans = 0;
#pragma omp parallel num_threads(thread_count)
    {
#pragma omp master
        {
            gm.init(thread_count, this, schedule);
        }
#pragma omp barrier //mynodel have to be calculated before running other threads
#pragma omp master
        {
            global_ans = gm.runmajor();
        }
        if (omp_get_thread_num()) {
            VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
            long long local_ans = 0;
            VertexSet subtraction_set;
            VertexSet tmp_set;
            subtraction_set.init();
            int last = -1;
            gm.set_loop_flag();
            auto match_edge = [&](int vertex, int *data, int size) {
                if (vertex != last) {
                    if (~last) subtraction_set.pop_back();
                    unsigned int l, r;
                    get_edge_index(vertex, l, r);
                    for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                        vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, edge + l, r - l, prefix_id);
                    }
                    subtraction_set.push_back(vertex);
                    last = vertex;
                }
                gm.set_loop(data, size);
                pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
            };
            for (unsigned int *data; data = gm.get_edge_range();) {
                match_edge(data[1], edge + data[2], data[3] - data[2]);
                /*for (int i = 1; i <= data[4]; i++) {
                    int l, r;
                    get_edge_index(data[1] + i, l, r);
                    match_edge(data[1] + i, edge + l, r - l);
                }*/
            }
            delete[] vertex_set;
            gm.report(local_ans);
        }
    }
    return global_ans;
}
