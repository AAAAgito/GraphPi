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

void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    v = intra_vertex_dict.at(v);
    l = vertex[v];
    r = vertex[v + 1];
}
void Graph::get_extern_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    v = inter_vertex_dict.at(v);
    l = inter_vertex[v];
    r = inter_vertex[v + 1];
}

void Graph::load_extern_data(int file_id) {
    //read file and get N(v)

    // append further neighbor to load_list
    return;
}

void Graph::load_order_manage() {
    while (!load_list.empty()) {
        std::vector<loader> load_vertex;
        load_vertex.push_back(load_list.back());
        int file_id = load_list.back().file_id;
        load_vertex.pop_back();
        //add a mutex lock;
        while (load_vertex.back().file_id == file_id) {
            load_vertex.push_back(load_list.back());
            load_list.pop_back();
        }
        //release lock
        load_extern_data(file_id);
    }
}

void Graph::load_list_append(loader l){
    int idx = 0;
    for (auto it : load_list) {
        if (v_state_map[it.vertex].is_intra) continue;
        if ((v_state_map[it.vertex].is_loaded || v_state_map[it.vertex].is_loading) && v_state_map[it.vertex].k_hop >= l.k_hop) continue;
        if (it.vertex == l.vertex) {
            it.k_hop = std::max(it.k_hop,l.k_hop);
            return;
        }
        if (it.file_id > l.file_id) {
            load_list.insert(load_list.begin()+idx,l);
            return;
        }
        idx++;
    }
    load_list.push_back(l);
    return;
}

bool Graph::extern_store_manage(const Schedule& schedule) {
    // first calculate the available space
    std::vector<int> k_hop_matrix = schedule.k_hop_matrix;
    float extern_usage = external_used / external_space;
    // second base on the available determine if drop the data
    if (extern_usage > extern_thresold){
        extern_drop_manage();
        return false;
    }
    // third base on the available space, determine how much vertex will be loaded.
    unsigned int remain_space = external_space - external_used;
    extern_load_manage(remain_space, schedule);
    for (int i=candidate_bin.size()-1; i >=0 ; i--) {
        int k_hop = k_hop_matrix[candidate_bin[i]->depth];
        bool vector_load = true;
        for (int j=0; j< candidate_bin[i]->vertex.size();j++) {
            int vertex = candidate_bin[i]->vertex.at(j);
            // two version can be formulated, one is one by one. another is satisfying only when all in vectors are satisfied.
            // satisfy conditions
            vector_load &= v_state_map[vertex].is_loaded && v_state_map[vertex].k_hop >= k_hop;
            if (!vector_load) break;
        }
        if (vector_load) {
            ready_bin.push_back(candidate_bin[i]);
            candidate_bin.erase(candidate_bin.begin()+i);
        }
    }
}

void Graph::extern_drop_manage() {
    for (auto it : inter_vertex_dict) {
        v_state_map[it.first].is_loaded = false;
        v_state_map[it.first].k_hop = 0;
        v_state_map[it.first].is_loading = false;
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
            // case 1: to_load vertex is not loaded
            if (!v_state_map[j].is_loaded && !v_state_map[j].is_loading) load_list_append(loader{k_hop, j, v_state_map[j].file_id});
            // case 2: to_load vertex is loaded and loaded hop is larger
            if (v_state_map[j].is_loaded && v_state_map[j].k_hop >= k_hop) continue;
            // case 3: to_load vertex is loaded but loaded hop is smaller
            if (v_state_map[j].is_loaded && v_state_map[j].k_hop < k_hop) {
                v_state_map[j].is_loaded = false;
                v_state_map[j].is_loading = true;
                v_state_map[j].k_hop = k_hop;
                load_list_append(loader{k_hop, j, v_state_map[j].file_id});
            }
            // case 4: to_load vertex is loading and loading hop is larger
            if (v_state_map[j].is_loading && v_state_map[j].k_hop >= k_hop) continue;
            // case 5: to_load vertex is loading but loading hop is smaller
            if (v_state_map[j].is_loading && v_state_map[j].k_hop < k_hop) {
                v_state_map[j].k_hop = k_hop;
                load_list_append(loader{k_hop, j, v_state_map[j].file_id});
            }
        }
    }
    // a loading function from class DataLoader
    load_order_manage();
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
        int vertex = loop_data_ptr[i];
        if (!clique)
            if (subtraction_set.has_data(vertex))
                continue;
        unsigned int l, r;
        get_edge_index(vertex, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, vertex, clique);
            if( vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex);
        pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, depth + 1, clique);
        subtraction_set.pop_back();
    }
}

void Graph::resume_matching(const Schedule& schedule, path *p) {

    VertexSet tmp_set;
    unsigned int l,r;
    int depth = p->depth;
    for (int i=0; i< p->vertex.size();i++) {
        int vertex = p->vertex[i];
        get_extern_edge_index(vertex, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            p->vertex_set[prefix_id].build_vertex_set(schedule, p->vertex_set, &inter_edge[l], (int)r - l, prefix_id, vertex);
            if( p->vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        p->subtraction_set.push_back(vertex);
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
                // TODO: add a manager to allocate task from candidate bin to each thread.
                // clear remain task in the ready bin is the highest priority
                while (!ready_bin.empty()) {
                    // do remain work in ready bin
                    path *p = ready_bin.back();
                    ready_bin.pop_back();
                    resume_matching(schedule,p);
                }
                int vtx = bfs_queue.front();
                bfs_queue.pop();
                if (v_state_map[vtx].is_rooted == true) continue;
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
                if(true)
                    pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
                else
                    pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, 1, clique);
                subtraction_set.pop_back();
                // TODO: add break condition
                if (!extern_store_manage(schedule)) break;
            }
        }
        // TODO: parallel I/O
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
    std::vector<int> inter_vertex;
    for (int &i = ii; i < loop_size; ++i)
    {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int vertex = loop_data_ptr[i];
        if (subtraction_set.has_data(vertex))
            continue;
        unsigned int l, r;
        VertexTable t = v_state_map[vertex];
        if(!t.is_intra && !t.is_loaded) {
            inter_vertex.push_back(vertex);
            continue;
        }
        else if (t.is_intra){
            get_edge_index(vertex, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id, vertex);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) continue;
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(vertex);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@
        }
        else if (t.is_loaded) {
            get_extern_edge_index(vertex, l, r);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &inter_edge[l], (int)r - l, prefix_id, vertex);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) continue;
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(vertex);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@
        }
    }
    if (!inter_vertex.empty()){

        path *path_ptr = new path;
        path_ptr->depth = depth;
        path_ptr->load_queue = false;
        path_ptr->vertex = inter_vertex;
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
