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

/*int Graph::intersection_size_mpi(int v1, int v2) {
    Graphmpi &gm = Graphmpi::getinstance();
    int ans = 0;
    if (gm.include(v2))
        return intersection_size(v1, v2);
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    int *data = gm.getneighbor(v2);
    for (int l2 = 0; l1 < r1 && ~data[l2];) {
        if(edge[l1] < data[l2]) {
            ++l1;
        }
        else if(edge[l1] > data[l2]) {
            ++l2;
        }
        else {
            ++l1;
            ++l2;
            ++ans;
        }
    }
    return ans;
}
*/

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

void Graph::init_extern_storage(int v_num, int e_num) {
    inter_edge = new int[e_num];
    inter_vertex = new unsigned int[v_num];
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
    if (blockType == BlockType::K_CORE_BLOCK) file_map = path + std::string("_fk");
    if (blockType == BlockType::RANDOM_BLOCK) file_map = path + std::string("_fr");
    int *fkey = new int[v_cnt];
    int *fvalue = new int[v_cnt];
    DataLoader::load_map_data(fkey,fvalue,v_cnt,file_map);
    for(int i=0;i<v_cnt;i++) {
        v_state_map[fkey[i]].file_id = fvalue[i];
    }
    delete[] vkey;
    delete[] vvalue;
    delete[] fkey;
    delete[] fvalue;
}

void Graph::load_partition_graph(int pid, int num, const std::string &path) {
    intra_vertex_dict.clear();
    v_state_map.clear();
    std::string partition_path = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) + std::string("_p");

    DataLoader::load_data_size(v_cnt,e_cnt,partition_path);
    vertex = new unsigned int[v_cnt];
    edge = new int[e_cnt];
    DataLoader::load_partition_data(v_cnt,e_cnt,vertex,edge,partition_path);
    std::string partition_map_data = path + std::string("_") + std::to_string(pid) + std::string("_") + std::to_string(num) + std::string("_m");
    int *vkey = new int[v_cnt];
    int *vvalue = new int[v_cnt];
    DataLoader::load_map_data(vkey,vvalue,v_cnt,partition_map_data);
    for (int i=0; i<v_cnt; i++) intra_vertex_dict[vkey[i]] = vvalue[i];
    std::string file_map;
    if (blockType == BlockType::K_CORE_BLOCK) file_map = path + std::string("_fk");
    if (blockType == BlockType::RANDOM_BLOCK) file_map = path + std::string("_fr");
    int *fkey = new int[g_vcnt];
    int *fvalue = new int[g_vcnt];
    DataLoader::load_map_data(fkey,fvalue,g_vcnt,file_map);
    for (int i=0;i<g_vcnt;i++) {
        v_state_map.insert({fkey[i], VertexTable{}});
        v_state_map[fkey[i]].file_id = fvalue[i];
    }
    for (int i=0; i<v_cnt; i++) v_state_map[vkey[i]].is_intra = true;
    delete[] vkey;
    delete[] vvalue;
    delete[] fkey;
    delete[] fvalue;
}



void Graph::to_block_csr(int block_size, const std::string& path, int k_core) {
    int len = sizeof (char) * block_size / sizeof (int);
    //TODO: is there any blocking algorithm that helps?
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
        }
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
            if (i==v_cnt-1) {
                degree[i] = e_cnt - vertex[i];
                break;
            }
            degree[i] = vertex[i+1] - vertex[i];
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
            blocking_data_manage_k_core(len,block_id,local_block_size,block_vid,revert_dict,to_insert_block, vid, file_id, path);
        }
        for (auto i : degree) to_insert_block.push_back(i.first);
        blocking_data_manage_k_core(len,block_id,local_block_size,block_vid,revert_dict,to_insert_block, vid, file_id, path);

        std::string block_map_path = path + std::string("_fk");
        DataLoader::gen_map_file(block_vid.data(),file_id.data(),block_vid.size(),block_map_path);
    }
}

void Graph::blocking_data_manage_k_core(int len, int &block_id, int &local_block_size, std::vector<int> &block_vid, const std::map<int, int> &revert_dict, std::vector<int> &to_insert_block, std::vector<int> &vid, std::vector<int> &file_id, const std::string& path) {
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
            //TODO: write block_vid to block file
            std::sort(block_vid.begin(), block_vid.end());
            int cursor = 0;
            std::vector<unsigned int> vtx_offset;
            std::vector<int> edges;
            vtx_offset.push_back(cursor);
            for (int i=0; i< block_vid.size()-1;i++) {
                if (intra_vertex_dict[block_vid[i]] == v_cnt-1) cursor += e_cnt - vertex[intra_vertex_dict[block_vid[i]]];
                else cursor +=vertex[intra_vertex_dict[block_vid[i]]+1] - vertex[intra_vertex_dict[block_vid[i]]];
                for (int j = vtx_offset.back(); j < cursor; j++) edges.push_back(edge[j]);
                vtx_offset.push_back(cursor);
            }
            std::string block_csr_path = path + std::string("_blocks/") + std::to_string(block_id) + std::string("_bk");
            DataLoader::gen_block_file(block_vid.data(),vtx_offset.data(),edges.data(),block_vid.size(),edges.size(),block_csr_path);
            for (auto i : block_vid) {
                vid.push_back(i);
                file_id.push_back(block_id);
            }
            block_id+=1;
            printf("%d\n",block_id);
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
    }
}

void Graph::to_partition_csr(PartitionType t, int num, const std::string& path) {
    if (t == PartitionType::RANDOM) {
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
            std::string partition_map = path + "_" + std::to_string(p)+ "_" + std::to_string(num) + std::string("_m");
            DataLoader::gen_map_file(mapk,mapv,partial_vcnt,partition_map);
            std::string partition_path = path+ "_" + std::to_string(p)+ "_" + std::to_string(num)  + std::string("_p");
            DataLoader::gen_partition_file(partial_vcnt,adj_list.size(),part_vertex,adj_list.data(),partition_path);


            delete[] part_vertex;
        }
    }
    if (t == PartitionType::METIS) {
    }
}

struct loader_cmp {
    bool operator() (loader a, loader b) {return a.vertex< b.vertex;}
} loaderCmp;

void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    int vtx = intra_vertex_dict.at(v);
    l = vertex[vtx];
    if (vtx == v_cnt -1) r = e_cnt;
    else r = vertex[vtx + 1];
}
void Graph::get_extern_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    int vtx = inter_vertex_dict.at(v);
    l = inter_vertex[vtx];
    if (vtx == inter_vtx_num - 1) r = inter_edge_num;
    else r = inter_vertex[vtx + 1];
}

void Graph::load_extern_data(int file_id, const std::vector<loader>& load_vertex) {
    //read file and get N(v)
    // fread, scanf, getline..
    int *vid = nullptr;
    unsigned int *vertex_offset = nullptr;
    int *adj_list = nullptr;

    std::string block_path;
    if (blockType == BlockType::K_CORE_BLOCK) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bk");
    else if (blockType == BlockType::RANDOM_BLOCK) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_br");
    int load_v_len=0;
    unsigned int load_e_len=0;
    DataLoader::load_data_size(load_v_len,load_e_len,block_path);
    vid = new int[load_v_len];
    vertex_offset = new unsigned int[load_v_len];
    adj_list = new int[load_e_len];
    DataLoader::load_block_data(vid,vertex_offset,adj_list, load_v_len, load_e_len, block_path);
    int vid_cursor = 0;
    for (auto i : load_vertex) {
        // find vertex
        while (vid[vid_cursor] != i.vertex && vid_cursor < load_v_len) {
            vid_cursor++;
            if (vid_cursor >= load_v_len) printf("expect err");
        }
        // add this vertex to list
        // TODO: add lock
        inter_vertex_dict[i.vertex] = inter_vtx_num;
        inter_vertex[inter_vtx_num] = inter_edge_num;
        inter_vtx_num +=1;
        unsigned int edge_cursor = inter_edge_num;
        unsigned int len;
        if (vid_cursor <load_v_len-1 ) len = vertex_offset[vid_cursor+1] - vertex_offset[vid_cursor];
        else len = load_e_len - vertex_offset[vid_cursor];
        inter_edge_num += len;
        // release lock
        if (i.k_hop==1) {
            for (int j = 0; j < len; j++) {
                inter_edge[edge_cursor] = adj_list[vertex_offset[vid_cursor]+j];
                edge_cursor += 1;
            }
            v_state_map[i.vertex].k_hop = 1;
        } else if (i.k_hop > 1) {
            for (int j = 0; j < len; j++) {
                inter_edge[edge_cursor] = adj_list[vertex_offset[vid_cursor]+j];
                edge_cursor += 1;
                //here add recursive function (to achieve, just need append load_list, as loading_manage will keep calling this function until load_list becomes empty)
                load_list_append(loader{i.k_hop-1, adj_list[vertex_offset[vid_cursor]+j], v_state_map[adj_list[vertex_offset[vid_cursor]+j]].file_id});
            }
        }
    }
    delete[] vid;
    delete[] vertex_offset;
    delete[] adj_list;


}

void Graph::loading_manage() {
    std::set<int> load_vtx;
    while (!load_list.empty()) {
        std::vector<loader> load_vertex;
        load_vertex.push_back(load_list.back());
        int file_id = load_list.back().file_id;
        load_list.pop_back();
        //TODO: add a mutex lock;
        while (load_list.back().file_id == file_id) {
            if (load_vtx.find(load_list.back().vertex)==load_vtx.end()) load_vtx.insert(load_list.back().vertex);
            load_vertex.push_back(load_list.back());
            load_list.pop_back();
        }
        //release lock
        std::sort(load_vertex.begin(), load_vertex.end(),loaderCmp);
        load_extern_data(file_id, load_vertex);
    }
    for (auto i : load_vtx) {
        v_state_map[i].k_hop = v_state_map[i].loading;
        v_state_map[i].loading = 0;
    }
}

void Graph::load_list_append(loader l){
    int idx = 0;
    // if it is not inter-partition vertex, no considering to load k-hop
    if (v_state_map[l.vertex].is_intra) return;
    // when there already exist a loading process that contain higher k-hop, or current k-hop is larger, no need to append to task queue
    if ( v_state_map[l.vertex].loading >= l.k_hop || v_state_map[l.vertex].k_hop >= l.k_hop) return;
    if (v_state_map[l.vertex].k_hop >= 1) {
        unsigned int i,j;
        get_extern_edge_index(l.vertex,i,j);
        for (int k = 0; k<j-i; k++) {
            load_list_append(loader{l.k_hop-1, inter_edge[i+k], l.file_id});
        }
        return;
    }
    for (auto it : load_list) {
        if (it.vertex == l.vertex) {
            it.k_hop = std::max(it.k_hop,l.k_hop);
            v_state_map[it.vertex].loading = it.k_hop;
            return;
        }
        if (it.file_id > l.file_id) {
            load_list.insert(load_list.begin()+idx,l);
            v_state_map[it.vertex].loading = it.k_hop;
            return;
        }
        idx++;
    }
    load_list.push_back(l);
}

bool Graph::extern_store_manage(const Schedule& schedule) {
    // first calculate the available space
    std::vector<int> k_hop_matrix = schedule.k_hop_matrix;
    float extern_usage = external_used / external_space;
    // second base on the available determine if drop the data
    if (extern_usage > extern_thresold){
        // to avoiding load extern info drop even when it is not used, it will only stop loading if ready_bin is not empty
        if (ready_bin.empty()) extern_drop_manage();
        return false;
    }
    // third base on the available space, determine how much vertex will be loaded.
    unsigned int remain_space = external_space - external_used;
    extern_load_manage(remain_space, schedule);
    // only this function is allocate to multi thread
    loading_manage();
    for (int i=candidate_bin.size()-1; i >=0 ; i--) {
        int k_hop = k_hop_matrix[candidate_bin[i]->depth];
        bool vector_load = true;
        for (int v : candidate_bin[i]->vertex) {
            // TODO: change it into per vertex, use smart pointer to vertex_set
            // two version can be formulated, one is one by one. another is satisfying only when all in vectors are satisfied.
            // satisfy conditions
            vector_load &= v_state_map[v].k_hop >= k_hop;
            if (!vector_load) break;
        }
        if (vector_load) {
            ready_bin.push_back(candidate_bin[i]);
            candidate_bin.erase(candidate_bin.begin()+i);
        }
    }
    return true;
}

void Graph::extern_drop_manage() {
    for (auto it : inter_vertex_dict) {
        v_state_map[it.first].k_hop = 0;
        v_state_map[it.first].loading = 0;
    }
    inter_vertex_dict.clear();
    inter_edge_num = 0;
    inter_vtx_num = 0;
}

void Graph::extern_load_manage(unsigned int available_space, const Schedule& schedule) {
    // manage candidate_bin to to_load set
    std::vector<int> k_hop_matrix = schedule.k_hop_matrix;
    for (auto & i : candidate_bin) {
        if (i->load_queue) continue;
        i->load_queue = true;
        int k_hop = k_hop_matrix[i->depth];
        for (auto j : i->vertex) {
            if (v_state_map[j].k_hop < k_hop) load_list_append(loader{k_hop, j, v_state_map[j].file_id});
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
    /*if (clique == true)
      {
      int last_vertex = subtraction_set.get_last();
    // The number of this vertex must be greater than the number of last vertex.
    loop_start = std::upper_bound(loop_data_ptr, loop_data_ptr + loop_size, last_vertex) - loop_data_ptr;
    }*/
    if (depth == schedule.get_size() - 1)
    {
        // TODO : try more kinds of calculation.
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (clique)
            local_ans += loop_size;
        else if (loop_size > 0)
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    int last_vertex = subtraction_set.get_last();
    for (int i = 0; i < loop_size; ++i)
    {
        if (last_vertex <= loop_data_ptr[i] && clique)
            break;
        int v = loop_data_ptr[i];
        if (!clique)
            if (subtraction_set.has_data(v))
                continue;
        unsigned int l, r;
        get_edge_index(v, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, v, clique);
            if( vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(v);
        pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, depth + 1, clique);
        subtraction_set.pop_back();
    }
}

void Graph::resume_matching(const Schedule& schedule, path *p) {

    VertexSet tmp_set;
    unsigned int l,r;
    int depth = p->depth;
    for (int i=0; i< p->vertex.size();i++) {
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
        p->subtraction_set.push_back(v);
        pattern_matching_aggressive_func(schedule,p->vertex_set,p->subtraction_set,tmp_set,p->local_ans,depth+1);
    }
    // after done release memory
    delete[] p->vertex_set;
    delete[] p->subtraction_set.get_data_ptr();
    delete p;
}

long long Graph::pattern_matching(const Schedule& schedule, int thread_count, bool clique)
{
    long long global_ans = 0;
    std::vector<int> k_hop_table = schedule.k_hop_matrix;
#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        double start_time = get_wall_time();
        double current_time;
        VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        long long local_ans = 0;
        // TODO : try different chunksize
#pragma omp for schedule(dynamic) nowait
        for (int v = 0; v < v_cnt; ++v)
        {
            if (v_state_map[v].is_rooted) continue;
            std::queue<int> bfs_queue;
            bfs_queue.push(v);
            while (!bfs_queue.empty()) {
                // clear remain task in the ready bin is the highest priority
                while (!ready_bin.empty()) {
                    // Highest priority to allocate task from candidate bin to each thread.
                    // do remain work in ready bin
                    path *p = ready_bin.back();
                    ready_bin.pop_back();
                    resume_matching(schedule,p);
                }
                int vtx = bfs_queue.front();
                bfs_queue.pop();
                if (v_state_map[vtx].is_rooted) continue;
                v_state_map[vtx].is_rooted = true;
                unsigned int l, r;
                get_edge_index(vtx, l, r);
                for (int i=0; i<=r-l;i++) {
                    int enqueue_vtx = edge[l+i];
                    // for vertex that do not belong to this partition, it will not be selected as a root
                    if (!v_state_map[enqueue_vtx].is_intra) continue;
                    if (v_state_map[enqueue_vtx].is_rooted) continue;
                    bfs_queue.push(enqueue_vtx);
                }
                for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
                {
                    // get N(v)
                    // As root must be intra-partition, external_edge will not be used
                    vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id);
                }
                //subtraction_set.insert_ans_sort(vertex);
                subtraction_set.push_back(vtx);
                //if (schedule.get_total_restrict_num() > 0 && clique == false)
                if(true) {
                    pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
                }
                else
                    pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, 1, clique);
                subtraction_set.pop_back();
                if (!extern_store_manage(schedule)) break;
            }
        }
        // TODO: add condition for scattering thread for data_managing
        while (!candidate_bin.empty()) extern_store_manage(schedule);
#pragma omp for schedule(dynamic) nowait
        for (int i=0; i< ready_bin.size();i++) {
            // do remain work in ready bin
            resume_matching(schedule,ready_bin[i]);
        }
        ready_bin.clear();
        delete[] vertex_set;
        // TODO : Computing multiplicty for a pattern
        global_ans += local_ans;
        
    }
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

void Graph::pattern_matching_aggressive_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth) // 3 same # or @ in comment are useful in code generation ###
{
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);// @@@
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;

    int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
/* @@@ 
    //Case: in_exclusion_optimize_num = 2
    if (depth == schedule.get_size() - 2 && schedule.get_in_exclusion_optimize_num() == 2) { 
        int loop_set_prefix_id_nxt = schedule.get_loop_set_prefix_id( depth + 1);
        int loop_size_nxt = vertex_set[loop_set_prefix_id_nxt].get_size();
        int size1 = VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        int size2 = VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id_nxt], subtraction_set);
        VertexSet tmp_set;
        tmp_set.init();
        tmp_set.intersection(vertex_set[loop_set_prefix_id], vertex_set[loop_set_prefix_id_nxt]);
        int size3 = VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
        local_ans += 1ll * size1 * size2 - size3;
        return;
    }
*/
/*
    //Case: in_exclusion_optimize_num = 3
    if( depth == schedule.get_size() - 3 && schedule.get_in_exclusion_optimize_num() == 3) { 
        int in_exclusion_optimize_num = 3;
        int loop_set_prefix_ids[ in_exclusion_optimize_num];
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id( depth + i );
        
        int loop_sizes[ in_exclusion_optimize_num ];
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            loop_sizes[i] = VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_ids[i]], subtraction_set);
        
        local_ans += 1ll * loop_sizes[0] * loop_sizes[1] * loop_sizes[2];

        for(int i = 1; i < 3; ++i) 
            for(int j = 0; j < i; ++j){
                VertexSet tmp_set;
                tmp_set.init();
                tmp_set.intersection(vertex_set[loop_set_prefix_ids[i]], vertex_set[loop_set_prefix_ids[j]]);
                int tmp_size = VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
                int size2;
                for(int k = 0; k < 3; ++k)
                    if( i != k && j != k) size2 = loop_sizes[k];
                local_ans -= 1ll * tmp_size * size2;
            }
        VertexSet tmp1;
        tmp1.init();
        tmp1.intersection(vertex_set[loop_set_prefix_ids[0]], vertex_set[loop_set_prefix_ids[1]]);
        VertexSet tmp2;
        tmp2.init();
        tmp2.intersection(vertex_set[loop_set_prefix_ids[2]], tmp1);
        local_ans += 1ll * 2 * VertexSet::unorderd_subtraction_size(tmp2, subtraction_set);
        return;
    }
*/
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
        // TODO : try more kinds of calculation. @@@
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0)
        {
            int min_vertex = v_cnt;
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
  
    // TODO : min_vertex is also a loop invariant @@@
    int min_vertex = v_cnt;
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
    if (depth == 1) Graphmpi::getinstance().get_loop(loop_data_ptr, loop_size);
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
        if (v_state_map[load_v].is_intra){
            get_edge_index(load_v, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(1); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
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
        else if(v_state_map[load_v].loading == 0) {
            to_load_vertex.push_back(load_v);
            continue;
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
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@
        }
    }
    if (!to_load_vertex.empty()){
        path *path_ptr = new path;
        path_ptr->depth = depth;
        path_ptr->load_queue = false;
        path_ptr->vertex = to_load_vertex;
        path_ptr->local_ans = local_ans;
        // deep copy
        path_ptr->subtraction_set.deepcopy(subtraction_set);
        // deep copy vertex_set
        // will it be too large to save?
        path_ptr->vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        for (int i=0; i< schedule.get_total_prefix_num(); i++) path_ptr->vertex_set[i].deepcopy(vertex_set[i]);
        candidate_bin.push_back(path_ptr);
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
