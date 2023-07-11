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

long long Graph::pattern_matching(const Schedule& schedule, int thread_count, bool clique)
{
    load_num = 0;
    io_num = 0;
    physicals_priority.resize(v_state_map.size());
    for(auto i : v_state_map) physicals_priority[i.first] = NULL;
    L.init_list();
    long long global_ans = 0;

#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        long long local_ans = 0;
        // TODO : try different chunksize
#pragma omp for schedule(dynamic) nowait
        for (int vertex_id = 0; vertex_id < v_state_map.size(); vertex_id++)
        {
            unsigned int l,r;
            int vtx = vertex_id;
            if (!v_state_map[vtx].is_intra) continue;
            printf("%d %d\n",vtx,omp_get_thread_num());
            get_edge_index(vtx,l,r);
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
        }
        global_ans += local_ans;
    }
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
    for (int &i = ii; i < loop_size; ++i)
    {
        if (min_vertex <= loop_data_ptr[i])
            break;
        int load_v = loop_data_ptr[i];
        if (subtraction_set.has_data(load_v))
            continue;

        if (v_state_map[load_v].is_intra){
            unsigned int l, r;
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
        else {
            int adj_size = 0;
            std::shared_ptr<int[]> ptr;
            bool next_hop = false;
            if (schedule.k_hop_matrix[depth] > 1) next_hop = true;
            get_physical_edge_index(load_v, ptr, adj_size, next_hop);
            if (ptr==NULL) {
                printf("is null\n");
            }
//            printf("%d %d\n",adj_size,ptr[0]);
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, ptr.get(), (int)adj_size, prefix_id, load_v);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) {
                ptr = NULL;
//                drop_physical_edge_index(load_v);
                continue;
            }
            //subtraction_set.insert_ans_sort(vertex);
            subtraction_set.push_back(load_v);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@

            ptr = NULL;
//            drop_physical_edge_index(load_v);

        }
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
