#include "gtest/gtest.h"
#include "dataloader.h"
#include "common.h"
//single function test

//TEST(testDataLoad, global_csr) {
//    std::string raw_path = "/mnt/d/graph/patents_input";
//    std::string graph_path = "/mnt/d/graph/asgs_bin/patents";
//    Graph *g;
//    DataLoader D;
//    double t1 = get_wall_time();
//    D.load_data(g,DataType::Patents,raw_path.c_str());
//    double t2 = get_wall_time();
//    printf("raw load time %.6lf\n", t2-t1);
//    g->to_global_csr(graph_path);
//    Graph g2;
//    t1 = get_wall_time();
//    g2.load_global_graph(graph_path);
//    t2 = get_wall_time();
//    printf("csr binary load time %.6lf\n", t2-t1);
//    EXPECT_EQ(g2.v_cnt,g->v_cnt);
//    EXPECT_EQ(g2.e_cnt,g->e_cnt);
//    EXPECT_EQ(g2.intra_vertex_dict.size(),g->intra_vertex_dict.size());
//}


//TEST(testDataLoad, sorted_csr) {
//
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
//    Graph g2;
//    g2.load_global_graph(graph_path);
//    for (int i = 0; i<g2.v_cnt-1; i++) {
//        unsigned int l,r;
//        g2.get_edge_index(i,l,r);
//        for (int j = l ; j< r-1; j++) EXPECT_LE(g2.edge[j], g2.edge[j+1]);
//    }
//}

//TEST(testDataLoad, block_csr_random) {
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
//    Graph g2;
//    std::vector<std::vector<int>> v;
//    g2.blockType = BlockType::RANDOM_BLOCK;
//    g2.load_global_graph(graph_path);
//    g2.to_block_csr(4*1024,v,graph_path,5);
//}


//TEST(testDataload, read_block_csr) {
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote_blocks/";
//    for (int i=0 ;i<300 ;i++) {
//        std::string p = graph_path + std::to_string(i) + std::string("_bk");
////        printf("%s",p.c_str());
//        int v=20;
//        unsigned int e=100;
//        int a[100];
//        double t1 = get_wall_time();
////        DataLoader::load_data_size(v,e,p);
//        DataLoader::load_block_data_aggregate(a,e,p);
//        double t2 = get_wall_time();
//        printf("\nread time: %.6lf\n", t2 - t1);
//
////        printf("\n%d %d %d\n",i,v,e);
////        int vid[v], edge[e];
////        unsigned int vertex[v];
////        DataLoader::load_block_data(vid,vertex,edge,v,e,p);
////        for (int i = 0; i< v-1; i++) {
////            printf("%d\n", vid[i]);
////            EXPECT_LE(vid[i],vid[i+1]);
////        }
//
//    }
//}

//TEST(testDataLoad, block_csr_correctness) {
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
//    Graph g2(graph_path);
//    g2.blockType = BlockType::K_CORE_BLOCK;
//    g2.load_global_graph(graph_path);
//    // TODO: debug K_CORE_BLOCK
//    g2.init_extern_storage(g2.v_cnt,g2.e_cnt);
//    for (auto i : g2.intra_vertex_dict) {
////        if (i.first < 6800) continue;
//        std::vector<loader> v;
//        v.push_back({1,i.first,g2.v_state_map[i.first].file_id});
//
//        unsigned int l0,r0;
//        g2.get_edge_index(i.first,l0,r0);
//
////        printf("vid %d fid %d\n",i.first,g2.v_state_map[i.first].file_id);
//        g2.load_extern_data(g2.v_state_map[i.first].file_id,v);
//        EXPECT_NE(g2.inter_vertex_dict.find(i.first),g2.inter_vertex_dict.end());
//        unsigned int l,r;
//        g2.get_extern_edge_index(i.first,l,r);
//
////        printf("len %d %d\n",r-l,r0-l0);
//        ASSERT_EQ(r-l,r0-l0);
//        for (int j =0; j< r-l; j++) ASSERT_EQ(g2.edge[j+l0],g2.inter_edge[j+l]);
//    }
//}


//TEST(testDataLoad, partition_csr) {
//
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
//    Graph g2(graph_path);
//    g2.load_global_graph(graph_path);
//    g2.partitionType = NAIVE_BFS;
//    int part = 4;
//    std::vector<std::vector<int>> vec;
//    g2.to_partition_csr(part,vec,graph_path);
//    int v=0;
//    int e=0;
//    Graph gs[part];
//    for (int p=0; p< part; p++) {
//        gs[p] = Graph(graph_path);
//        gs[p].g_vcnt = g2.v_cnt;
//        gs[p].partitionType = NAIVE_BFS;
//        gs[p].load_partition_graph(p,part, graph_path);
//        v+=gs[p].v_cnt;
//        e+=gs[p].e_cnt;
//    }
//    EXPECT_EQ(v,g2.v_cnt);
//    EXPECT_EQ(e,g2.e_cnt);
//
//    for (auto i: g2.intra_vertex_dict) {
//        int vid = i.first;
//        int cnt = 0;
//        unsigned int l0,r0;
//        g2.get_edge_index(vid,l0,r0);
//        for (int p=0;p<part;p++) {
//            if (gs[p].v_state_map[vid].is_intra) {
//                unsigned int l,r;
//                gs[p].get_edge_index(vid,l,r);
//                EXPECT_EQ(r-l,r0-l0);
////                printf("%d %d %d\n",vid,l,r);
//                for (int j =0; j< r-l; j++) ASSERT_EQ(g2.edge[j+l0],gs[p].edge[j+l]);
//                cnt++;
//            }
//        }
//        ASSERT_EQ(cnt,1);
//    }
//}

//TEST(testMatching, match_in_global) {
//    // patterns:
//    Pattern tc_pattern(3);
//    tc_pattern.add_edge(0, 1);
//    tc_pattern.add_edge(1, 2);
//    tc_pattern.add_edge(2, 0);
//    Pattern r(Rectangle);
//    // setting
//    bool is_pattern_valid;
//    bool use_in_exclusion_optimize = false;
//    int performance_type = 2;
//    int restricts_type = 2;
//    int thread_num = 1;
//
//    //global matching
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
//    Graph g(graph_path);
//    g.load_global_graph(graph_path);
//    g.tri_cnt = 608389;
//    g.g_vcnt = g.v_cnt;
//    for (int i =0 ; i < g.v_cnt -1; i++) {
//        int len = g.vertex[i+1]-g.vertex[i];
//        VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, len);
//    }
//    int len = g.e_cnt - g.vertex[g.v_cnt-1];
//    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, len);
//
//
//    Schedule tc_schedule(tc_pattern, is_pattern_valid, performance_type, restricts_type, use_in_exclusion_optimize, g.v_cnt, g.e_cnt);
//    double t1 = get_wall_time();
//    long long global_result = g.pattern_matching(tc_schedule,1);
//    double t2 = get_wall_time();
//    printf("candidate_bin %d\n",candidate_insert);
//    printf("%lld\n",global_result);
//    printf("thread %d time: %.6lf\n", thread_num, t2 - t1);
////    ASSERT_EQ(g.tri_cnt,global_result);
//
//}


//TEST(testMatching, gen_k_hops) {
//
//    // patterns:
//    Pattern tc_pattern(3);
//    tc_pattern.add_edge(0, 1);
//    tc_pattern.add_edge(1, 2);
//    tc_pattern.add_edge(2, 0);
//    Pattern r(Rectangle);
//    Pattern d(House);
//    Pattern c(Clique_7_Minus);
//    // setting
//    bool is_pattern_valid;
//    bool use_in_exclusion_optimize = false;
//    int performance_type = 2;
//    int restricts_type = 2;
//    int thread_num = 1;
//    int part =4;
//    //global matching
//
//    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
//    Graph g(graph_path);
//    g.load_global_graph(graph_path);
//    Schedule tc_schedule(c, is_pattern_valid, performance_type, restricts_type, use_in_exclusion_optimize, g.v_cnt, g.e_cnt);
//    tc_schedule.gen_k_hop_matrix();
//    tc_schedule.print_schedule();
//    for (int i : tc_schedule.k_hop_matrix) printf("%d ",i);
//    printf("\n");
//}

TEST(testMatching, match_in_partition_r) {

    // patterns:
    Pattern tc_pattern(3);
    tc_pattern.add_edge(0, 1);
    tc_pattern.add_edge(1, 2);
    tc_pattern.add_edge(2, 0);

    // setting
    bool is_pattern_valid;
    bool use_in_exclusion_optimize = false;
    int performance_type = 2;
    int restricts_type = 2;
    int thread_num = 1;
    int part = 4;
    LoadType loadType = SUB_BLOCK;
    BlockType blockType = RANDOM_BLOCK;
    PartitionType partitionType = RANDOM;
    BlockType block_types[3] = {RANDOM_BLOCK, K_CORE_BLOCK, CHINK_BFS};
    //global matching
    blockType = block_types[0];
    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
    Graph g(graph_path);
    g.load_global_graph(graph_path);
    g.tri_cnt = 608389;
    g.g_vcnt = g.v_cnt;
    for (int i = 0; i < g.v_cnt - 1; i++) {
        unsigned int l, r;
        g.get_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
    g.max_degree = VertexSet::max_intersection_size;
    g.loadType = loadType;
    g.blockType = blockType;
    g.partitionType = partitionType;
    g.gen_out_of_core_component(part, 4 * 1024, graph_path);

    Schedule tc_schedule(tc_pattern, is_pattern_valid, performance_type, restricts_type, use_in_exclusion_optimize,
                         g.v_cnt, g.e_cnt);
    tc_schedule.gen_k_hop_matrix();
    long long ground_truth_result = 0;

    double t1 = get_wall_time();
    ground_truth_result = g.pattern_matching(tc_schedule, 1);
    double t2 = get_wall_time();

    long long global_result = 0;
    double total_time = 0.0;
    for (int i = 0; i < part; i++) {
        Graph g2(graph_path);
        g2.set_static_indicator_zero();
        g2.blockType = blockType;
        g2.loadType = loadType;
        g2.partitionType = partitionType;
        g2.g_vcnt = g.v_cnt;
        g2.load_partition_graph(i, part, graph_path);
        printf("vertex num %d",g2.v_cnt);
        g2.init_extern_storage(g.v_cnt, g.e_cnt);
        g2.tri_cnt = g.tri_cnt / part;
        g2.max_degree = VertexSet::max_intersection_size;
        double t1 = get_wall_time();
        global_result += g2.pattern_matching(tc_schedule, 1);
        double t2 = get_wall_time();
        total_time += t2 - t1;
    }

    printf("\nglobal time: %.6lf\n", t2 - t1);
    printf("total time: %.6lf\n", total_time);
    ASSERT_EQ(ground_truth_result, global_result);
}

TEST(testMatching, match_in_partition_c) {

    // patterns:
    Pattern tc_pattern(3);
    tc_pattern.add_edge(0, 1);
    tc_pattern.add_edge(1, 2);
    tc_pattern.add_edge(2, 0);

    // setting
    bool is_pattern_valid;
    bool use_in_exclusion_optimize = false;
    int performance_type = 2;
    int restricts_type = 2;
    int thread_num = 1;
    int part = 4;
    LoadType loadType = SUB_BLOCK;
    BlockType blockType = RANDOM_BLOCK;
    PartitionType partitionType = NAIVE_BFS;
    BlockType block_types[3] = {RANDOM_BLOCK, K_CORE_BLOCK, CHINK_BFS};
    //global matching
    blockType = block_types[2];
    std::string graph_path = "/mnt/d/graph/asgs_bin/WikiVote";
    Graph g(graph_path);
    g.load_global_graph(graph_path);
    g.tri_cnt = 608389;
    g.g_vcnt = g.v_cnt;
    for (int i = 0; i < g.v_cnt - 1; i++) {
        unsigned int l, r;
        g.get_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
    g.max_degree = VertexSet::max_intersection_size;
    g.loadType = loadType;
    g.blockType = blockType;
    g.partitionType = partitionType;
    g.gen_out_of_core_component(part, 4 * 1024, graph_path);

    Schedule tc_schedule(tc_pattern, is_pattern_valid, performance_type, restricts_type, use_in_exclusion_optimize,
                         g.v_cnt, g.e_cnt);
    tc_schedule.gen_k_hop_matrix();
    long long ground_truth_result = 0;

    double t1 = get_wall_time();
    ground_truth_result = g.pattern_matching(tc_schedule, 1);
    double t2 = get_wall_time();

    long long global_result = 0;
    double total_time = 0.0;
    for (int i = 0; i < part; i++) {
        Graph g2(graph_path);
        g2.set_static_indicator_zero();
        g2.blockType = blockType;
        g2.loadType = loadType;
        g2.partitionType = partitionType;
        g2.g_vcnt = g.v_cnt;
        g2.load_partition_graph(i, part, graph_path);
        g2.init_extern_storage(g.v_cnt, g.e_cnt);
        g2.tri_cnt = g.tri_cnt / part;
        g2.max_degree = VertexSet::max_intersection_size;
        double t1 = get_wall_time();
        global_result += g2.pattern_matching(tc_schedule, 1);
        double t2 = get_wall_time();
        total_time += t2 - t1;
    }

    printf("\nglobal time: %.6lf\n", t2 - t1);
    printf("total time: %.6lf\n", total_time);
    ASSERT_EQ(ground_truth_result, global_result);
}