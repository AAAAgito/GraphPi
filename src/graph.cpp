#include "../include/graph.h"
#include "../include/graphmpi.h"
#include "../include/vertex_set.h"
#include "../include/common.h"
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <omp.h>
#include <cstring>
#include <atomic>
#include <queue>
#include <stack>
#include <iostream>


void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = vertex[v];
    if (v == v_cnt -1) r = e_cnt;
    else r = vertex[v + 1];
}

void Graph::get_mmp_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = mmp_vertex[v];
    if (v == g_vcnt -1) r = g_ecnt;
    else r = mmp_vertex[v + 1];
}

void Graph::get_mmp_patch_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = mmp_patch_vertex[v];
    if (v == extra_v_cnt -1) r = insert_cnt;
    else r = mmp_patch_vertex[v + 1];
}

void Graph::get_insert_code_index(int v, int &l, int &r) const {
    l = patch_insert_idx[v];
    if (v == extra_v_cnt - 1) r = insert_code_len;
    else r = patch_insert_idx[v+1];
}

void Graph::get_delete_code_index(int v, int &l, int &r) const {
    l = patch_delete_idx[v];
    if (v == extra_v_cnt - 1) r = delete_code_len;
    else r = patch_delete_idx[v+1];
}

long long Graph::pattern_matching(const Schedule& schedule, int thread_count, bool clique)
{
    long long global_ans = 0;

#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        long long local_ans = 0;
#pragma omp for schedule(dynamic) nowait
        for (int vertex_id = 0; vertex_id < v_cnt; vertex_id++)
        {
            unsigned int l,r;
            int vtx = vertex_id;
            get_edge_index(vtx,l,r);
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int)r - l, prefix_id);
            }
            subtraction_set.push_back(vtx);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);

            subtraction_set.pop_back();
        }
        global_ans += local_ans;
    }
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}


long long Graph::pattern_matching_oc(const Schedule& schedule, int thread_count, bool clique)
{
    long long global_ans = 0;

#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        long long local_ans = 0;
#pragma omp for schedule(dynamic) nowait
        for (int vertex_id = 0; vertex_id < g_vcnt; vertex_id++)
        {
            unsigned int l,r;
            int vtx = vertex_id;
            get_mmp_edge_index(vtx,l,r);
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &mmp_edge[l], (int)r - l, prefix_id);
            }
            subtraction_set.push_back(vtx);
            pattern_matching_aggressive_func_oc(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);

            subtraction_set.pop_back();
        }
        global_ans += local_ans;
    }
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

void Graph::pattern_matching_aggressive_func_oc(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth) // 3 same # or @ in comment are useful in code generation ###
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

        unsigned int l, r;
        get_mmp_edge_index(load_v, l, r);
        int *patch_merge_data = NULL;
        int *input_data;
        // if adj list not changed, directly used, else regenerate

        if (false) {
        // if (load_v > g_vcnt || test(bitmap,load_v)) {
            int del_l,del_r,ins_l,ins_r;
            get_insert_code_index(load_v,ins_l,ins_r);
            get_delete_code_index(load_v,del_l,del_r);
            int origin_cursor = 0;
            int ins_cursor = ins_l;
            int del_cursor = del_l;
            patch_merge_data = new int[r-l+ins_r-ins_l-del_r+del_l];
            int merge_cursor = 0;

            while (true)
            {
                if (origin_cursor < r-l)
                {
                    // who first come who first served
                    if (del_cursor < del_r && ins_cursor < ins_r) {
                        // round for doing deletion
                        if (patch_delete_code[del_cursor] < patch_insert_code[ins_cursor]) {
                            // we are now at:
                            int operating_cursor = patch_delete_code[del_cursor];
                            // reserve all data from last saved
                            memcpy(patch_merge_data+merge_cursor,mmp_edge+l+origin_cursor,operating_cursor-origin_cursor);
                            // next storing will move forward as this step has store such length
                            merge_cursor += operating_cursor-origin_cursor;
                            // next reserving start index move to current index
                            origin_cursor = operating_cursor;
                            // number of deletions(ignore)
                            int del_len = patch_delete_code[del_cursor+1];
                            origin_cursor += del_len;
                            // goto next del operation
                            del_cursor += 2;
                        }
                        // round for doing insertion
                        else {
                            // we are now at:
                            int operating_cursor = patch_insert_code[ins_cursor];
                            // reserve all data from last saved
                            memcpy(patch_merge_data+merge_cursor,mmp_edge+l+origin_cursor,operating_cursor-origin_cursor);
                            // next storing will move forward as this step has store such length
                            merge_cursor += operating_cursor-origin_cursor;
                            // next reserving start index move to current index
                            origin_cursor = operating_cursor;
                            // number of insertion
                            int ins_len;
                            if (ins_cursor == 0) ins_len = patch_insert_code[ins_cursor+1]-0;
                            else ins_len = patch_insert_code[ins_cursor+1] - patch_insert_code[ins_cursor-1];
                            // TODO:get data from .insert
                            unsigned int patch_l,patch_r;
                            get_mmp_patch_edge_index(load_v,patch_l,patch_r);
                            int start_poi = patch_insert_code[ins_cursor+1] - ins_len;
                            memcpy(patch_merge_data+merge_cursor,mmp_patch_edge+patch_l+start_poi,ins_len);

                            // insertion has store such lens to merge array. origin would not change.
                            merge_cursor += ins_len;
                            // goto next ins operation
                            ins_cursor += 2;

                        }
                    }
                    // In this case, only deletion after
                    else if (del_cursor < del_r) {
                        // we are now at:
                        int operating_cursor = patch_delete_code[del_cursor];
                        // reserve all data from last saved
                        memcpy(patch_merge_data+merge_cursor,mmp_edge+l+origin_cursor,operating_cursor-origin_cursor);
                        // next storing will move forward as this step has store such length
                        merge_cursor += operating_cursor-origin_cursor;
                        // next reserving start index move to current index
                        origin_cursor = operating_cursor;
                        // number of deletions(ignore)
                        int del_len = patch_delete_code[del_cursor+1];
                        origin_cursor += del_len;
                        // goto next del operation
                        del_cursor += 2;
                    }
                    // In this case, only insertion after
                    else if (ins_cursor < ins_r) {
                        // we are now at:
                        int operating_cursor = patch_insert_code[ins_cursor];
                        // reserve all data from last saved
                        memcpy(patch_merge_data+merge_cursor,mmp_edge+l+origin_cursor,operating_cursor-origin_cursor);
                        // next storing will move forward as this step has store such length
                        merge_cursor += operating_cursor-origin_cursor;
                        // next reserving start index move to current index
                        origin_cursor = operating_cursor;
                        // number of insertion
                        int ins_len;
                        if (ins_cursor == 0) ins_len = patch_insert_code[ins_cursor+1]-0;
                        else ins_len = patch_insert_code[ins_cursor+1] - patch_insert_code[ins_cursor-1];
                        // TODO:get data from .insert
                        unsigned int patch_l,patch_r;
                        get_mmp_patch_edge_index(load_v,patch_l,patch_r);
                        int start_poi = patch_insert_code[ins_cursor+1] - ins_len;
                        memcpy(patch_merge_data+merge_cursor,mmp_patch_edge+patch_l+start_poi,ins_len);

                        // insertion has store such lens to merge array. origin would not change.
                        merge_cursor += ins_len;
                        // goto next ins operation
                        ins_cursor += 2;
                    }
                    // In this case no more change operation
                    else {
                        memcpy(patch_merge_data+merge_cursor,mmp_edge+l+origin_cursor,r-origin_cursor);
                        origin_cursor = r-l;
                        merge_cursor += r-origin_cursor;
                        // assure it exactly reach the end
                        assert(merge_cursor == r-l+ins_r-ins_l-del_r+del_l);
                    }
                }
                // reach end but still has insertion
                else if (ins_cursor < ins_r) {
                    assert(ins_cursor+2 > ins_r);
                    
                    // we are now at:
                    int operating_cursor = patch_insert_code[ins_cursor];
                    // reserve all data from last saved
                    memcpy(patch_merge_data+merge_cursor,mmp_edge+l+origin_cursor,operating_cursor-origin_cursor);
                    // next storing will move forward as this step has store such length
                    merge_cursor += operating_cursor-origin_cursor;
                    // next reserving start index move to current index
                    origin_cursor = operating_cursor;
                    // number of insertion
                    int ins_len;
                    if (ins_cursor == 0) ins_len = patch_insert_code[ins_cursor+1]-0;
                    else ins_len = patch_insert_code[ins_cursor+1] - patch_insert_code[ins_cursor-1];
                    // TODO:get data from .insert
                    unsigned int patch_l,patch_r;
                    get_mmp_patch_edge_index(load_v,patch_l,patch_r);
                    int start_poi = patch_insert_code[ins_cursor+1] - ins_len;
                    memcpy(patch_merge_data+merge_cursor,mmp_patch_edge+patch_l+start_poi,ins_len);

                    // insertion has store such lens to merge array. origin would not change.
                    merge_cursor += ins_len;

                    // assure it exactly reach the end
                    assert(merge_cursor == r-l+ins_r-ins_l-del_r+del_l);
                }
                else
                    break;
            }
            input_data = patch_merge_data;

        }
        else {
            input_data = mmp_edge+l;
        }

        
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, input_data, (int)r - l, prefix_id, load_v);
            if( vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(load_v);
        pattern_matching_aggressive_func_oc(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@

        subtraction_set.pop_back(); // @@@
        if (patch_merge_data != NULL) delete[] patch_merge_data;
        
    }

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

}

