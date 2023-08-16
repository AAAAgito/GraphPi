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
#include <mpi.h>
#include <atomic>
#include <queue>
#include <stack>
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
//    printf("get v %d\n",intra_vertex_dict.size());
    // int vtx = intra_vertex_dict.at(v);
    l = vertex[v];
    if (v == v_cnt -1) r = e_cnt;
    else r = vertex[v + 1];
}

void Graph::get_mmp_edge_index(int v, unsigned int& l, unsigned int& r) const
{
//    printf("get v %d\n",intra_vertex_dict.size());
    l = mmp_vertex[v];
    if (v == g_vcnt -1) r = g_ecnt;
    else r = mmp_vertex[v + 1];
}

long long Graph::pattern_matching(const Schedule& schedule, int thread_count, bool clique)
{
    io_num = 0;
    lock_vertex.resize(v_state_map.size());
    // file_len.resize(block_lengths.size());
    // file_str.resize(block_lengths.size());
    // for (auto i : block_lengths) {
    //     file_len[i.first]=i.second;
    //     file_str[i.first] = std::to_string(i.first);
    // }
    // block_lengths.clear();
    extern_upper_size = extern_upper_thresold * g_vcnt;
    file_path = raw_data_path.append("_blocks/bk_");
    // // file_path = "/home/yanglaoyuan/patents_blocks/bk_";
    
    omp_init_lock(&priority_lock);
    for (int i=0; i< v_state_map.size(); i++) omp_init_lock(&lock_vertex[i]);
    long long global_ans = 0;

#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        VertexSet* vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        long long local_ans = 0;
#pragma omp for schedule(dynamic) nowait
        for (int vertex_id = 0; vertex_id < v_state_map.size(); vertex_id++)
        {
            unsigned int l,r;
            int vtx = vertex_id;
            if (!v_state_map[vtx].is_intra) continue;
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
    omp_destroy_lock(&priority_lock);
    for (int i=0; i< v_state_map.size(); i++) omp_destroy_lock(&lock_vertex[i]);
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
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
        {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &mmp_edge[l], (int)r - l, prefix_id, load_v);
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
            // int adj_size = 0;
            // std::shared_ptr<int[]> ptr;

            double t1 = get_wall_time2();
            // get_physical_edge_index(load_v, ptr, adj_size);
            unsigned int l,r;
            get_mmp_edge_index(load_v,l,r);
            double t2 = get_wall_time2();
            if (omp_get_thread_num()==0) load_extern_time += t2-t1;
            bool is_zero = false;
            for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id))
            {
                // vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, ptr.get(), (int)adj_size, prefix_id, load_v);
                vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &mmp_edge[l], (int)r-l, prefix_id, load_v);
                if( vertex_set[prefix_id].get_size() == 0) {
                    is_zero = true;
                    break;
                }
            }
            if( is_zero ) {
                // ptr = NULL;
                continue;
            }
            subtraction_set.push_back(load_v);
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, depth + 1);// @@@
            subtraction_set.pop_back(); // @@@

            // ptr = NULL;

        }
    }

}

// void Graph::gen_bfs_query_order() {
//     bfs_order_queue.clear();
//     bool *picked = new bool[g_vcnt];
//     for (int i=0;i<g_vcnt;i++)
//         picked[i] = !v_state_map[i].is_intra;
//     std::queue<int> q;
//     for (int i=0; i<g_vcnt; i++) {
//         if (picked[i]) continue;
//         if (q.empty()) {
//             q.push(i);
//             bfs_order_queue.push_back(i);
//             picked[i] = true;
//         }
//         while (!q.empty()) {
//             int v = q.front();
//             q.pop();
//             unsigned int l,r;
//             get_edge_index(v,l,r);
//             std::vector<int> sort_bfs;
//             std::map<int,int> degree;
//             for (int j =l;j<r;j++) {
//                 if (picked[edge[j]]) continue;
//                 q.push(edge[j]);
//                 picked[edge[j]] = true;
//                 unsigned int s,e;
//                 get_edge_index(edge[j],s,e);
//                 degree[edge[j]] = e-s;
//                 sort_bfs.push_back(edge[j]);
//             }
//             for (int l=0; l<sort_bfs.size();l++) {
//                 for (int m=l+1;m<sort_bfs.size();m++) {
//                     if (degree[sort_bfs[l]] < degree[sort_bfs[m]])
//                     {
//                         std::swap(sort_bfs[l],sort_bfs[m]);
//                     }
//                 }
//             }
//             for (auto a:sort_bfs)
//                 bfs_order_queue.push_back(a);
//         }
// //        bfs_order_queue.push_back(i);
//     }
//     assert(bfs_order_queue.size()==v_cnt);
// }

// void Graph::gen_bipartite_order(std::vector<int> &part_map, int cur_part, int inter_upper) {
// //    bfs_order_queue.clear();
//     printf("part %d\n",cur_part);
//     std::vector<int> cross_num;
//     cross_num.resize(v_cnt);
//     for (int i=0;i<v_cnt;i++) {
//         cross_num[i] = 0;
//         unsigned int l,r;
//         if (part_map[i]!=cur_part) continue;
//         get_edge_index(i,l,r);
//         for (int j=l;j<r;j++) {
//             if (part_map[edge[j]] !=cur_part) cross_num[i] +=1;
//         }
//     }
//     std::vector<int> score;
//     std::queue<int> inter_queue;
//     std::set<int> s;
//     std::vector<bool> picked_visited;
//     picked_visited.resize(v_cnt);
//     score.resize(v_cnt);
//     int process_size = 0;
//     for (int i=0;i<v_cnt;i++) score[i] =0;
//     // O(V)
//     for (int intra_v = 0; intra_v < v_cnt; intra_v++) {
//         if (part_map[intra_v]!=cur_part) continue;
//         if (picked_visited[intra_v]) continue;
// //        printf("%d\n",intra_v);
//         int max_ptr = intra_v;
//         double max_score = 0;
//         unsigned int l,r;
//         bfs_order_queue.push_back(intra_v);
//         picked_visited[intra_v] = true;
//         while (max_ptr>=0) {
//             get_edge_index(max_ptr,l,r);
//             max_score = -1;
//             max_ptr = -1;
//             if ((process_size+1)/10000 > process_size/10000) printf("process %d w\n",process_size/10000);
//             // O(Degree)
//             for (auto i: s) {
//                 if (part_map[i] != cur_part) continue;
//                 if (picked_visited[i]) continue;
//                 if (max_score < (double )score[i]/ cross_num[i]) {
//                     max_ptr = i;
//                     max_score = (double )score[i]/ cross_num[i];
//                 }
//             }
//             for (int i=l;i<r;i++) {
//                 if (part_map[edge[i]]==cur_part) continue;
//                 if(picked_visited[edge[i]]) continue;
//                 inter_queue.push(edge[i]);
//                 picked_visited[edge[i]] = true;
//                 unsigned int l0,r0;
//                 get_edge_index(edge[i],l0,r0);
//                 // O(Degree)
//                 for (int j=l0;j<r0;j++) {
//                     if (part_map[edge[j]]!=cur_part) continue;
//                     if (picked_visited[edge[j]]) continue;
//                     if (score[edge[j]]==0) s.insert(edge[j]);
//                     score[edge[j]] += 1;
//                     if ((double )score[edge[j]] / cross_num[edge[j]] > max_score) {
//                         max_score = (double )score[edge[j]] / cross_num[edge[j]];
//                         max_ptr = edge[j];
//                     }
//                 }
//             }
//             if (max_ptr <0) break;
//             if (cross_num[max_ptr]==0)
//             assert(max_ptr>=0 && max_ptr < v_cnt);
// //            printf("max score %f\n",max_score);
//             bfs_order_queue.push_back(max_ptr);
//             picked_visited[max_ptr] = true;
//             process_size+=1;
//             s.erase(max_ptr);
//             while (inter_queue.size() > inter_upper) {
//                 int quit = inter_queue.front();
//                 assert(part_map[quit] != cur_part);
//                 inter_queue.pop();
//                 picked_visited[quit] = false;
//                 get_edge_index(quit,l,r);
//                 for (int i=l; i<r;i++) {
//                     if (part_map[edge[i]]!=cur_part) continue;
//                     score[edge[i]] -=1;
//                 }
//             }
//         }
//     }
// }


// void Graph::load_order(int num, int size) {
//     std::string porder = raw_data_path +"_" + std::to_string(num) + "_" + std::to_string(size) + "_opg";
//     std::string pbasic = raw_data_path +"_" + std::to_string(num) + "_" + std::to_string(size) + "_pcore";
//     int core_num=0;
//     DataLoader::load_order_data(&core_num,1,pbasic);
//     bfs_order_queue.clear();
//     bfs_order_queue.resize(core_num);
//     DataLoader::load_order_data(bfs_order_queue.data(),core_num,porder);
//     for (auto i: bfs_order_queue) v_state_map[i].is_rooted = true;
// }