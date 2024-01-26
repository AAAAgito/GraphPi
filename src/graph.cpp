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

size_t Graph::intersection_size_clique(int v1,int v2, int v3) {
    size_t l1, r1;
    get_large_mmp_edge_index(v1, l1, r1);
    size_t l2, r2;
    get_large_mmp_edge_index(v2, l2, r2);
    size_t l3, r3;
    get_large_mmp_edge_index(v3, l3, r3);
    int min_vertex = v3;
    size_t ans = 0;
    if (mmp_edge[l1] >= min_vertex || mmp_edge[l2] >= min_vertex || mmp_edge[l3] >= min_vertex)
        return 0;
    while(l1 < r1 && l2 < r2 && l3 < r3) {
        if(mmp_edge[l1] < std::max(mmp_edge[l3],mmp_edge[l2])) {
            if (mmp_edge[++l1] >= min_vertex)
                break;
        }
        if(mmp_edge[l2] < std::max(mmp_edge[l3],mmp_edge[l1])) {
            if (mmp_edge[++l2] >= min_vertex)
                break;
        }
        else {
            if(mmp_edge[l3] < std::max(mmp_edge[l1],mmp_edge[l2])) {
                if (mmp_edge[++l3] >= min_vertex)
                    break;
            }
            else {
                ++ans;
                if (mmp_edge[++l1] >= min_vertex)
                    break;
                if (mmp_edge[++l2] >= min_vertex)
                    break;
                if (mmp_edge[++l3] >= min_vertex)
                    break;
            }
        }
    }
    return ans;
}

size_t Graph::intersection_size_clique(int v1,int v2) {
    size_t l1, r1;
    assert(v1>=0);
    get_large_mmp_edge_index(v1, l1, r1);
    size_t l2, r2;
    assert(v2>=0);
    get_large_mmp_edge_index(v2, l2, r2);
    int min_vertex = v2;
    size_t ans = 0;
    if (mmp_edge[l1] >= min_vertex || mmp_edge[l2] >= min_vertex)
        return 0;
    while(l1 < r1 && l2 < r2) {
        if(mmp_edge[l1] < mmp_edge[l2]) {
            if (mmp_edge[++l1] >= min_vertex)
                break;
        }
        else {
            if(mmp_edge[l2] < mmp_edge[l1]) {
                if (mmp_edge[++l2] >= min_vertex)
                    break;
            }
            else {
                ++ans;
                if (mmp_edge[++l1] >= min_vertex)
                    break;
                if (mmp_edge[++l2] >= min_vertex)
                    break;
            }
        }
    }
    return ans;
}

long long Graph::triangle_counting_mt(int thread_count) {
    long long ans = 0;
#pragma omp parallel num_threads(thread_count)
    {
        tc_mt(&ans);
    }
    printf("done\n");
    return ans;
}

void Graph::tc_mt(long long *global_ans) {
    long long my_ans = 0;
    int period = g_vcnt/100;
    #pragma omp for schedule(dynamic)
    for(int v = 0; v < g_vcnt; ++v) {
        // for v in G
        if (v%period==0) printf("progress %d\n",v/period);
        size_t l, r;
        get_large_mmp_edge_index(v, l, r);
        for(size_t v1 = l; v1 < r; ++v1) {
            if (v <= mmp_edge[v1])
                break;
            //for v1 in N(v)

            my_ans += intersection_size_clique(v,mmp_edge[v1]);
        }
    }
    #pragma omp critical
    {
        *global_ans += my_ans;
    }
}

void Graph::get_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = vertex[v];
    if (v == v_cnt -1) r = e_cnt;
    else r = vertex[v + 1];
}

void Graph::get_large_mmp_edge_index(int v, size_t& l, size_t& r) const
{
    l = mmp_l_vertex[v];
    if (v == g_vcnt -1) r = lg_ecnt;
    else r = mmp_l_vertex[v + 1];
}

void Graph::get_mmp_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = mmp_vertex[v];
    // if (v == g_vcnt -1) r = g_ecnt;
    // else r = mmp_vertex[v + 1];
    r = mmp_vertex[v+1];

    // // prefetch
    // if (KMD[2*v]!=-1) {
    //     madvise(mmp_edge+KMD[2*v],KMD[2*v+1],MADV_WILLNEED);
    // }
}

void Graph::get_mmp_patch_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = mmp_patch_vertex[v];
    if (v == extra_v_cnt -1) r = insert_cnt;
    else r = mmp_patch_vertex[v + 1];
}

void Graph::get_mmp_patch_delete_edge_index(int v, unsigned int& l, unsigned int& r) const
{
    l = mmp_patch_delete_vertex[v];
    if (v == extra_v_cnt -1) r = delete_cnt;
    else r = mmp_patch_delete_vertex[v + 1];
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
            // do prefetch
            if (do_prefetch){
                int pref_v=vertex_id+thread_count;
                madvise(mmp_edge+KMD[2*pref_v],KMD[2*pref_v+1],MADV_WILLNEED);
            }
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
            if (do_prefetch){
                int pref_v=vertex_id+prefetch_interval;
                if (pref_v < g_vcnt)
                    madvise(mmp_edge+KMD[2*pref_v],KMD[2*pref_v+1],MADV_WILLNEED);
            }
            unsigned int l,r;
            int vtx = vertex_id;
            get_mmp_edge_index(vtx,l,r);
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &mmp_edge[l], (int)r - l, prefix_id);
            }
            subtraction_set.push_back(vtx);
            pattern_matching_mmap_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);

            subtraction_set.pop_back();
        }
        global_ans += local_ans;
    }
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

long long Graph::pattern_matching_oca(const Schedule& schedule, int thread_count, bool clique)
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
        for (int vtx = 0; vtx < g_vcnt; vtx++)
        {
            unsigned int l,r;
            // printf("v=%d %d\n",vtx,!test(bitmap,vtx));
            // for (int a=l;a<r;a++) {
            //     printf("%d ",mmp_edge[a]);
            // }
            // printf("\n");
            // printf("needle \n",);
            int *input_data=NULL;
            int size;
            if(!bitmap[vtx]){
                needle(input_data,vtx,size);
                bitmap[vtx]=true;
            }
            else{
                input_data = g_back_up_E+g_back_up_V[vtx];
                if (vtx==extra_v_cnt-1) size = backup_ecnt - g_back_up_V[vtx];
                else size = g_back_up_V[vtx+1] - g_back_up_V[vtx];
                // printf("hit\n");
            }
            // for (int a=0;a<size;a++) {
            //     printf("%d ",input_data[a]);
            // }
            // printf("\n");
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, input_data, (int)size, prefix_id);
            }
            subtraction_set.push_back(vtx);
            pattern_matching_oca(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);

            subtraction_set.pop_back();
        }
        global_ans += local_ans;
    }
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

void Graph::pattern_matching_mmap_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth)
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
        int size = r-l;
        // if adj list not changed, directly used, else regenerate

        
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, mmp_edge+l, (int)size, prefix_id, load_v);
            if( vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(load_v);
        pattern_matching_mmap_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@

        subtraction_set.pop_back(); // @@@
        
    }

}


void Graph::pattern_matching_oca(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, VertexSet& tmp_set, long long& local_ans, int depth)
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

        int *input_data;
        int size;

        if (!bitmap[load_v]) {
            needle(input_data,load_v,size);
            
            bitmap[load_v]=true;
            // printf("call needle\n");
        }
        else {
            // printf("hit out\n");
            input_data = g_back_up_E+g_back_up_V[load_v];
            if (load_v==extra_v_cnt-1) size = backup_ecnt - g_back_up_V[load_v];
            else size = g_back_up_V[load_v+1] - g_back_up_V[load_v];
        }

        
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, input_data, (int)size, prefix_id, load_v);
            if( vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if( is_zero ) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(load_v);
        pattern_matching_mmap_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@

        subtraction_set.pop_back(); // @@@
        
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

void Graph::pattern_matching_func(const Schedule& schedule, VertexSet* vertex_set, VertexSet& subtraction_set, long long& local_ans, int depth, bool clique)
{
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;
    int* loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    if (depth == schedule.get_size() - 1)
    {
        if (clique == true)
            local_ans += loop_size;
        else if (loop_size > 0)
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    int last_vertex = subtraction_set.get_last();
    for (int i = 0; i < loop_size; ++i)
    {
        if (last_vertex <= loop_data_ptr[i] && clique == true)
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

void Graph::Sampling_exploring(std::vector<std::vector<int>> &edge_list, std::map<int,int> &cc, std::unordered_set<int> & bin, double rate, int depth) {
    double t1 = get_wall_time2();
    std::random_device e; 
    int v=0;
    srand(time(NULL));
    printf("rate %.6lf\n",rate);
    std::vector<int> sample_vertex_degree(g_vcnt,0);
    int sample_ecnt;
    for(unsigned int i=0; i<g_ecnt;i++) {
        if (v!= g_vcnt && i >= mmp_vertex[v+1]) {
            v+=1;
        }
        if (mmp_edge[i] > v) continue;
        double roll = (double)(rand()/(double)RAND_MAX);
        // printf("roll %.6lf %.6lf\n",roll,rate);
        if (roll < rate) {
            edge_list[v].push_back(mmp_edge[i]);
            edge_list[mmp_edge[i]].push_back(v);
            sample_vertex_degree[v]+=1;
            sample_vertex_degree[mmp_edge[i]]+=1;
        }
    }
    printf("collected done\n");
    // printf("sample edge %d\n",sample_edges.size());
    std::vector<bool> visited(g_vcnt,false);
    for (int v=0; v<g_vcnt; v++){
        if (sample_vertex_degree[v] > 0 && !visited[v]) {
            bin.insert(v);
            visited[v]=true;
        }
    }
    
        
    double t2 = get_wall_time2();
    printf("Function Sampling Times: %.6lf Bin %d\n", t2 - t1,bin.size());
}

void Graph::Gathering(int repeat, double rate, int depth, int threshold) {
    double t1 = get_wall_time2();
    std::unordered_map<int,int> counter;
    std::map<int,int> cc;
    std::vector<std::vector<int>> edge_list;
    edge_list.resize(g_vcnt);
    for(int i=0; i< repeat; i++) {
        std::unordered_set<int> bin;
        Sampling_exploring(edge_list, cc, bin, rate, depth);
        
        for (auto j : bin) {
            if (counter.find(j)==counter.end()) counter[j]=0;
            else counter[j]+=1;
        }
    }
    std::unordered_set<int> set;
    for (int i=0; i<repeat; i++) {
        for (auto j: counter) {
            if(j.second < threshold) continue;
            set.insert(j.first);
        }
    }
    
    ConnectedComponent(edge_list,set,cc);
    printf("collected done %d\n",set.size());
    
    for (auto& element : Bins) {
        element = -1;
    }
    Disk_benchmark();

    int component=0;
    for (auto a:cc) {
        component=std::max(a.second,component);
    }
    for(auto c : set) {
        Bins[c]=cc[c];
    }
    Advise_Translation(component);
    double t2 = get_wall_time2();
    printf("Function Gathering Times: %.6lf\n", t2 - t1);
}

int Graph::ConnectedComponent(std::vector<std::vector<int>> &edge_list, const std::unordered_set<int> &set, std::map<int,int> &cc) {
    double t1 = get_wall_time2();
    std::unordered_map<int,bool> visited;
    // for (auto v: edge_list) {
    //     if (v.size()>10)
    //         printf("vsize %d\n",v.size());
    // }
    for(auto j:set) {
        visited[j]=false;
    }
    int component=0;
    // printf("v0 %d\n",edge_list[0].size());
    for (auto j : set) {
        // printf("visiting\n");
        if(!visited[j]) {
            std::stack<int> dfs;
            visited[j]=true;
            dfs.push(j);
            cc[j]=component;
            while (!dfs.empty()){
                int vtx = dfs.top();
                dfs.pop();
                for (auto e : edge_list[vtx]){
                    if (visited[e] || set.find(e)==set.end()){
                        continue;
                    }
                    visited[e]=true;
                    dfs.push(e);
                    cc[e]=cc[j];
                }

            }
            component+=1;
        }
    }
    double t2 = get_wall_time2();
    
    printf("Function CC Times: %.6lf\n", t2 - t1);
    return component;
}

void Graph::Advise_Translation(int component) {
    double t1 = get_wall_time2();
    reading.resize(component);
    printf("component %d\n",component);
    for (int c=0; c<component; c++) {
        int cc_num=0;
        unsigned int prev_start=0,prev_end=0;
        printf("cur com %d\n",c);
        for (int v=0; v< g_vcnt;v++) {
            if (Bins[v]!=c) continue;
            // printf("v %d\n",v);
            cc_num+=1;
            // initial case
            if (prev_end==0) {
                get_mmp_edge_index(v,prev_start,prev_end);
                reading[c].push_back(interval{prev_start,prev_end});
                continue;
            }
            unsigned int present_start,present_end;
            get_mmp_edge_index(v,present_start,present_end);
            double sequential_speed = (present_start-prev_start)*4/system_disk.sequential_read;
            double split= 1/system_disk.iops + (prev_end-prev_start)*4/ system_disk.sequential_read;
            if (split < sequential_speed) {
                prev_start=present_start;
                prev_end=present_end;
                // reading[c].push_back(interval{prev_start,prev_end});
            }
            else {
                prev_end = present_end;
                // reading[c].back().end = prev_end;
            }
        }
        printf("Advise length %d cc num:%d\n",reading[c].size(),cc_num);
    }
    double t2 = get_wall_time2();
    
    printf("Function Advise Times: %.6lf\n", t2 - t1);
}

void Graph::Disk_benchmark() {
    system_disk.iops=98,000;
    system_disk.sequential_read=550,000,000;
}

void Graph::KMeasureDecompose(double density_threshold, double threshold) {
    // density = degree of useful / total space
    // available_degree_threshold
    KMD = new unsigned int[2*g_vcnt];
    printf("gvcnt %d\n",g_vcnt);
    
#pragma omp parallel for
    for (int v=0; v < g_vcnt; v++) {
        // printf("%d\n",v);
        KMD[2*v] = -1;
        unsigned int l,r;
        l=mmp_vertex[v];
        r=mmp_vertex[v+1];
        unsigned int adt = std::floor(threshold*(r-l));
        if (adt==0) continue;
        assert(r>l);
        unsigned int density_index=l;
        double density=0.0;
        unsigned int useful_space=0;
        // printf("%lld\n",r-adt);
        for (unsigned int w=l; w < r-adt; w++) {
            int sum=0;
            for (unsigned int m=0;m<adt;m++) {
                // contained vertex ids
                int neighbor = mmp_edge[w+m];
                sum += mmp_vertex[neighbor+1]-mmp_vertex[neighbor];
            }
            int start_v = mmp_edge[w];
            int end_v = mmp_edge[w+adt-1];
            assert(w+adt-1 < r);
            assert(w+adt -1 >= l);
            assert(start_v <= end_v);
            unsigned int radius = mmp_vertex[end_v+1]-mmp_vertex[start_v];
            double cal = (double)sum/(double)radius;
            assert(cal <= 1.0);
            if (density < cal) {
                density = cal;
                density_index=w;
                useful_space = radius;
            }
        }
        if (density < density_threshold) continue;
        // while (density_index>l && density_index < r-radius-1)
        // {
        //     int start_v = mmp_edge[density_index];
        //     int left_explore_v = mmp_edge[density_index-1];
        //     int right_explore_v = mmp_edge[density_index+radius];

        //     double cal_l = (useful_space + mmp_vertex[left_explore_v+1] - mmp_vertex[left_explore_v]) / (mmp_vertex[start_v+radius]-mmp_vertex[left_explore_v]);
        //     double cal_r = (useful_space + mmp_vertex[right_explore_v+1] - mmp_vertex[right_explore_v]) / (mmp_vertex[start_v+radius+1]-mmp_vertex[start_v]);
        //     if (cal_l < density_threshold && cal_r < density_threshold) break;
        //     radius+=1;
        //     if (cal_r < cal_l) {
        //         density_index-=1;
        //         useful_space += mmp_vertex[left_explore_v+1] - mmp_vertex[left_explore_v];
        //     } else {
        //         useful_space += mmp_vertex[right_explore_v+1] - mmp_vertex[right_explore_v];
        //     }
        // }
        KMD[2*v] = density_index;
        KMD[2*v+1] = useful_space;
    }
}