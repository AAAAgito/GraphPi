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


struct vertex_degree_cmp {
    bool operator() (vertex_degree a, vertex_degree b) {return a.d > b.d;}
} vertexDegreeCmp;


struct loader_cmp {
    bool operator() (loader a, loader b) {return a.vertex< b.vertex;}
} loaderCmp;

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
            int cur_hop = v_state_map[i.vertex].k_hop;
//            assert(cur_hop==0);
            while (vid[vid_cursor] != i.vertex && vid_cursor < load_v_len) {
                vid_cursor++;
//                assert(vid[vid_cursor] <= i.vertex);
                if (vid[vid_cursor] == i.vertex) break;
            }
            // add this vertex to list
            // TODO: when extern storage is not enough how to solve it? or avoid loading too much before.

            omp_set_lock(&extern_ve_loading_lock);
//            assert(inter_vertex_dict.find(i.vertex) == inter_vertex_dict.end());
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

            memcpy(inter_edge+edge_cursor,adj_list+vertex_offset[vid_cursor],len* sizeof(int));
//            for (int j = 0; j < len; j++) {
//
//                inter_edge[edge_cursor] = adj_list[vertex_offset[vid_cursor] + j];
//                edge_cursor += 1;
//            }
//            assert(edge_cursor < extern_e_max_num);

            omp_unset_lock(&extern_ve_loading_lock);
            // release lock
            if (i.k_hop == 1) {

                omp_set_lock(&extern_ve_loading_lock);
                v_state_map[i.vertex].k_hop = 1;
                omp_unset_lock(&extern_ve_loading_lock);
            } else if (i.k_hop > 1) {
                loaded_vertex.push_back(i.vertex);
                omp_set_lock(&extern_ve_loading_lock);
                int origin_hop = v_state_map[i.vertex].k_hop;
                v_state_map[i.vertex].k_hop = std::max(1,origin_hop);
                omp_unset_lock(&extern_ve_loading_lock);
                for (int j = 0; j < len; j++) {
                    int v = adj_list[vertex_offset[vid_cursor] + j];
                    //here add recursive function (to achieve, just need append load_list, as loading_manage will keep calling this function until load_list becomes empty)
                    load_list_append(loader{i.k_hop - 1, v, v_state_map[v].file_id}, loaded_vertex);
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

    if (omp_get_thread_num()==0) {
        if (load_list.empty()) {
            omp_set_lock(&init_load_list_append_lock);
            std::vector<loader> ini_loader;
            ini_loader.swap(init_load_list);
            omp_unset_lock(&init_load_list_append_lock);

            for (auto load_element: ini_loader) {
                load_list_append(load_element, loaded_vtx);
            }
        }
        while (!load_list.empty()) {
            omp_set_lock(&load_list_write_lock);
            std::vector<loader> load_vertex;

            //        printf("=g %d s %d\n",omp_get_thread_num(),load_list.size());
            int load_file_size = 0;
            int file_id = 0;
            //        printf("=f %d s %d\n",omp_get_thread_num(),load_list.size());
            for (auto i : load_list) {
                if (i.second.size()==0) load_list.erase(i.first);
                //            printf("==%d==%d==%d\n",i.first,i.second.size(),load_file_size);
                if (i.second.size() > load_file_size) {
                    load_file_size = i.second.size();
                    file_id = i.first;
                }
            }


            //        printf("=a %d %d s %d\n",file_id,omp_get_thread_num(),load_list.size());
            //        assert(load_list.find(file_id) != load_list.end());
            //        for (auto i : load_list.at(file_id)) {
            //            assert(v_state_map[i.vertex].k_hop==0);
            //        }


            load_vertex.swap(load_list.at(file_id));
            //        printf("=aa %d s %d\n",omp_get_thread_num(),load_list.size());
            load_list.erase(file_id);
            //        printf("=b %d s %d\n",omp_get_thread_num(),load_list.size());
            omp_unset_lock(&load_list_write_lock);
            for (auto i : load_vertex) {
                loaded_vtx.push_back(i.vertex);
                assert(i.vertex < g_vcnt && i.vertex >= 0);
            }

            //        printf("=c %d s %d\n",omp_get_thread_num(),load_list.size());
            //release lock
            std::sort(load_vertex.begin(), load_vertex.end(),loaderCmp);

            double t1 = get_wall_time();
            load_extern_data(file_id, load_vertex, loaded_vtx);
            double t2 = get_wall_time();
            load_extern_time += t2-t1;
            //        printf("=d %d s %d\n",omp_get_thread_num(),load_list.size());

        }
        //    printf("=e %d s %d\n",omp_get_thread_num(),load_list.size());
        for (auto i : loaded_vtx) {
            if (v_state_map[i].loading != 0) {
                omp_set_lock(&extern_ve_loading_lock);
                //            omp_set_lock(&v_state_khop_lock);
                v_state_map[i].k_hop = std::max((unsigned int) v_state_map[i].loading,v_state_map[i].k_hop);
                //            omp_unset_lock(&v_state_khop_lock);
                omp_unset_lock(&extern_ve_loading_lock);
                //            v_state_map[i].loading = 0;
            }
        }
    }
    else {

        omp_set_lock(&load_list_write_lock);
        std::vector<loader> load_vertex;
        int load_file_size = 0;
        int file_id = 0;
        for (auto i : load_list) {
            if (i.second.size()==0) load_list.erase(i.first);
            if (i.second.size() > load_file_size) {
                load_file_size = i.second.size();
                file_id = i.first;
            }
        }
        load_vertex.swap(load_list.at(file_id));
        load_list.erase(file_id);
        omp_unset_lock(&load_list_write_lock);
        for (auto i : load_vertex) {
            loaded_vtx.push_back(i.vertex);
            assert(i.vertex < g_vcnt && i.vertex >= 0);
        }

        std::sort(load_vertex.begin(), load_vertex.end(),loaderCmp);
        double t1 = get_wall_time();
        load_extern_data(file_id, load_vertex, loaded_vtx);
        double t2 = get_wall_time();
        load_extern_time += t2-t1;
    }
}

void Graph::load_list_append(const loader &l, std::vector<int> &loaded_vertex){
    // if it is not inter-partition vertex, no considering to load k-hop

//    printf("insertion %d\n",omp_get_thread_num());
    if (v_state_map[l.vertex].is_intra) return;
    if (l.k_hop==0) return;
    // when there already exist a loading process that contain higher k-hop, or current k-hop is larger, no need to append to task queue
    if ( v_state_map[l.vertex].loading >= l.k_hop || v_state_map[l.vertex].k_hop >= l.k_hop) return;

    if (v_state_map[l.vertex].k_hop >= 1) {
        omp_set_lock(&extern_ve_loading_lock);
        v_state_map[l.vertex].loading = std::max((int) l.k_hop,v_state_map[l.vertex].loading);
        omp_unset_lock(&extern_ve_loading_lock);
        loaded_vertex.push_back(l.vertex);
        int curr_hop = v_state_map[l.vertex].loading;
        unsigned int i,j;
        get_extern_edge_index(l.vertex,i,j);
        for (int k = i; k<j; k++) {
            int v = inter_edge[k];
            assert(l.k_hop-1 > 0 && l.k_hop < 5);
            load_list_append(loader{l.k_hop-1, v, v_state_map[v].file_id},loaded_vertex);
        }
        return;
    }

    omp_set_lock(&load_list_write_lock);
    if (v_state_map[l.vertex].k_hop == 0 && v_state_map[l.vertex].loading ==0) {
        omp_set_lock(&extern_ve_loading_lock);
        v_state_map[l.vertex].loading = std::max((int) l.k_hop,v_state_map[l.vertex].loading);
        omp_unset_lock(&extern_ve_loading_lock);
        loaded_vertex.push_back(l.vertex);
        int curr_hop = v_state_map[l.vertex].loading;
        if (load_list.find(l.file_id) == load_list.end()) {
            std::vector<loader> vec;
            vec.push_back(l);
            load_list[l.file_id] = vec;
        } else {
            for (auto &i: load_list[l.file_id])
                if (i.vertex == l.vertex) {
                    i.k_hop = std::max(i.k_hop, l.k_hop);
                    omp_unset_lock(&load_list_write_lock);
                    return;
                }
            load_list[l.file_id].push_back(l);
        }
    }
    if (v_state_map[l.vertex].k_hop == 0 && v_state_map[l.vertex].loading > 0) {

        omp_set_lock(&extern_ve_loading_lock);
        v_state_map[l.vertex].loading = std::max((int) l.k_hop,v_state_map[l.vertex].loading);
        omp_unset_lock(&extern_ve_loading_lock);
        loaded_vertex.push_back(l.vertex);
        int curr_hop = v_state_map[l.vertex].loading;
        for (auto &i: load_list[l.file_id])
            if (i.vertex == l.vertex) {
                i.k_hop = std::max(i.k_hop, l.k_hop);
                omp_unset_lock(&load_list_write_lock);
                return;
            }
    }
    omp_unset_lock(&load_list_write_lock);
}

bool Graph::extern_store_manage(const Schedule& schedule) {
    // first calculate the available space
    std::vector<int> k_hop_matrix = schedule.k_hop_matrix;
    double extern_usage = (double)inter_edge_num / (double)extern_e_max_num;
    // second base on the available determine if drop the data
    if (extern_usage > extern_upper_thresold){
//        printf("%d\n",inter_edge_num);
        printf("do drop\n");
        // to avoiding load extern info drop even when it is not used, it will only stop loading if ready_bin is not empty

        omp_set_lock(&ready_bin_lock);
        omp_set_lock(&v_state_lock);
        omp_set_lock(&extern_ve_loading_lock);
        if (ready_bin.empty())
            extern_drop_manage();
        omp_unset_lock(&extern_ve_loading_lock);
        omp_unset_lock(&v_state_lock);
        omp_unset_lock(&ready_bin_lock);
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
//    omp_set_lock(&extern_ve_loading_lock);
    int max_keep_idx = inter_vtx_num / 2;
    int start_v_idx = max_keep_idx+1;
    int start_e_idx = inter_vertex[start_v_idx];
    for (auto it : inter_vertex_dict) {
        if (! it.second >= start_v_idx) {
#pragma omp atomic
            v_state_map[it.first].k_hop *= 0;
#pragma omp atomic
            v_state_map[it.first].loading *= 0;
            inter_vertex_dict.erase(it.first);
        }
    }
//    inter_vertex_dict.clear();
    printf("do clear\n");
//#pragma omp atomic
    inter_edge_num -= start_e_idx;
    assert(inter_edge_num >0);
    memcpy(inter_edge,inter_edge+start_e_idx,inter_edge_num);

//#pragma omp atomic
    inter_vtx_num -= start_v_idx;
    assert(inter_vtx_num > 0);
//    omp_unset_lock(&extern_ve_loading_lock);

}

void Graph::extern_load_init(unsigned int available_space, const Schedule& schedule) {
    // manage candidate_bin to to_load set\

    std::vector<path *> candidates(candidate_bin);
    omp_set_lock(&candidate_bin_lock);
    for (int i = 0; i < candidate_bin.size(); i++) {
        if (candidate_bin[i]->load_queue) continue;
        candidate_bin[i]->load_queue = true;
        int k_hop = schedule.k_hop_matrix[candidate_bin[i]->depth];
        for (int j = 0; j < candidate_bin[i]->v_size; j++) {
            int v = candidate_bin[i]->vertex[j];
            omp_set_lock(&init_load_list_append_lock);
            assert(k_hop >= 0 && k_hop < 5);
            init_load_list.push_back(loader{k_hop, v, v_state_map[v].file_id});
            omp_unset_lock(&init_load_list_append_lock);
        }
    }
    omp_unset_lock(&candidate_bin_lock);
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


long long Graph::pattern_matching_asyn(const Schedule& schedule, int thread_count, bool clique)
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
        bool choose[v_state_map.size()];
        for (int vertex_id = 0; vertex_id < v_state_map.size(); vertex_id++) choose[vertex_id] = false;
#pragma omp for schedule(dynamic) nowait
        for (int vertex_id = 0; vertex_id < v_state_map.size(); vertex_id++)
        {
            // translate it into vid

            omp_set_lock(&v_state_lock);
            if (!v_state_map[vertex_id].is_intra) {
                omp_unset_lock(&v_state_lock);
                continue;
            }
            if (v_state_map[vertex_id].is_rooted || choose[vertex_id]) {
                omp_unset_lock(&v_state_lock);
                continue;
            }

            std::queue<int> bfs_queue;
            choose[vertex_id] = true;
            bfs_queue.push(vertex_id);
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
                if (!bfs_queue.empty())
                    bfs_queue.pop();

                omp_set_lock(&v_state_lock);
                if (v_state_map[vtx].is_rooted) {

                    omp_unset_lock(&v_state_lock);
                    continue;
                }
                v_state_map[vtx].is_rooted = true;
                omp_unset_lock(&v_state_lock);
//                printf("%d %d %d\n",vtx,bfs_queue.size(), omp_get_thread_num());
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
                    if (v_state_map[enqueue_vtx].is_rooted || choose[enqueue_vtx]) {
                        omp_unset_lock(&v_state_lock);
                        continue;
                    }
                    omp_unset_lock(&v_state_lock);
                    choose[enqueue_vtx] = true;
                    bfs_queue.push(enqueue_vtx);
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

        int to_load_num = 0;
        omp_set_lock(&candidate_bin_lock);
        for (auto i : candidate_bin)
            to_load_num += i->v_size;
        omp_unset_lock(&candidate_bin_lock);

        if (scatter_loading_thread(thread_count,id,to_load_num)) {
            while (!candidate_bin.empty()) {
//                printf("can size %d\n",candidate_bin.size());
//                for (auto i : candidate_bin) {
//                    for (int j = 0; j < i->v_size; j++) {
//                        printf("%d %d %d %d %d\n",i->vertex[j], i->depth,v_state_map[i->vertex[j]].k_hop,v_state_map[i->vertex[j]].loading,schedule.k_hop_matrix[i->depth]);
//                    }
//                }
                extern_store_manage(schedule);
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
        omp_set_lock(&ready_bin_lock);
        ready_bin.clear();
        omp_unset_lock(&ready_bin_lock);
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




void Graph::load_physical_edge_index(int v, bool next_hop) {
    io_num+=1;


    //read file and get N(v)
    // fread, scanf, getline..
    int file_id = v_state_map[v].file_id;
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
    if (loadType == SUB_BLOCK){
        int vid_cursor = 0;
        while (vid[vid_cursor] != v) vid_cursor++;
        int l,r;
        l=vertex_offset[vid_cursor];
        if (vid_cursor==load_v_len-1) r=load_e_len;
        else r= vertex_offset[vid_cursor+1];
        assert(r-l >0);
        std::shared_ptr<int[]> ptr(new int[r-l]());
//        int *ptr = new int[r-l];

        assert(physicals_load[v] == NULL);
        assert(physicals_priority[v] == NULL);
        physicals_priority[v] = new cell{v,NULL,NULL};
        L.insert_head(physicals_priority[v],v);

        physicals_length[v] = r-l;
        physicals_load[v] = ptr;


        std::vector<int> local_v;
        std::map<int,std::vector<int>> file_v;
//        printf("loop out side\n");
        for (int i = l; i < r; i++) {
            ptr[i - l] = edge[i];
            if (physicals_load.find(adj_list[i])!=physicals_load.end()) {
                if (v_state_map[adj_list[i]].is_intra) continue;
                if (physicals_load[adj_list[i]] != NULL) {
                    assert(physicals_priority[adj_list[i]] != NULL);
                    continue;
                }
            }
//            printf("==%d %d %d==\n",i,adj_list[i],physicals_priority[adj_list[i]]==NULL);
            assert(physicals_priority[adj_list[i]]==NULL);
            int fid = v_state_map[adj_list[i]].file_id;
            if (fid == file_id) {
                local_v.push_back(adj_list[i]);
            } else {
                if (file_v.find(fid)==file_v.end()) {
                    std::vector<int> vec;
                    vec.push_back(adj_list[i]);
                    file_v[fid] = vec;
                } else {
                    file_v[fid].push_back(adj_list[i]);
                }
            }
        }
        std::sort(local_v.begin(), local_v.end());
        if (next_hop) {
            int vid_cursor = 0;
            for (auto i : local_v) {
                while (vid[vid_cursor] != i) {
                    assert(vid[vid_cursor] < i);
                    vid_cursor++;
                }
                int l,r;
                l=vertex_offset[vid_cursor];
                if (vid_cursor==load_v_len-1) r=load_e_len;
                else r= vertex_offset[vid_cursor+1];
                assert(r-l >0);
                std::shared_ptr<int[]> ptr(new int[r-l]());
                for (int j = l; j < r; j++)
                    ptr[j - l] = adj_list[j];
//                int *ptr = new int[r-l];
                physicals_length[i] = r-l;
                physicals_load[i] = ptr;
                assert(physicals_priority[i] == NULL);
                physicals_priority[i] = new cell{i,NULL,NULL};
                L.insert_head(physicals_priority[i],i);
            }
            for (auto pair : file_v) {
                load_physical_edge_index(pair.second);
            }
            for (int i = l; i < r; i ++) {
                assert(physicals_load.find(adj_list[i])!= physicals_load.end());
                assert(physicals_load[adj_list[i]] != NULL);
//                printf("loaded %d %d\n",edge[i],physicals_load[edge[i]][0]);
            }
        }
    }
    delete[] vid;
    delete[] vertex_offset;
    delete[] adj_list;

    t2 = get_wall_time();
    blocking_manage_time += t2-t1;

}



void Graph::load_physical_edge_index(std::vector<int> &v) {
    io_num+=1;

    for (auto vtx : v) {

    }

    //read file and get N(v)
    // fread, scanf, getline..
    int file_id = v_state_map[v[0]].file_id;
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
    if (loadType == SUB_BLOCK){
        for (auto vtx : v) {
            int vid_cursor = 0;
            while (vid[vid_cursor] != vtx) vid_cursor++;
            int l, r;
            l = vertex_offset[vid_cursor];
            if (vid_cursor == load_v_len - 1) r = load_e_len;
            else r = vertex_offset[vid_cursor + 1];
            std::shared_ptr<int[]> ptr(new int[r - l]());
//            int *ptr = new int[r-l];
            physicals_length[vtx] = r - l;
            assert(physicals_load[vtx] == NULL);
            physicals_load[vtx] = ptr;


//            printf("--%d %d--\n",vtx,physicals_priority[vtx]==NULL);
            assert(physicals_priority[vtx] == NULL);
            physicals_priority[vtx] = new cell{vtx,NULL,NULL};
            L.insert_head(physicals_priority[vtx],vtx);

            std::map<int, std::vector<int>> file_v;
            for (int i = l; i < r; i++) {
                ptr[i - l] = adj_list[i];
            }
        }
    }
    delete[] vid;
    delete[] vertex_offset;
    delete[] adj_list;

    t2 = get_wall_time();
    blocking_manage_time += t2-t1;

}
