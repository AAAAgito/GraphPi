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
//    std::string graph_path = "/mnt/d/graph/asgs_bin/patents";
//    Graph g2;
//    g2.load_global_graph(graph_path);
//    for (int i = 0; i<g2.v_cnt-1; i++) {
//        int l = g2.vertex[i];
//        int r = g2.vertex[i+1];
//        for (int j = l ; j< r-1; j++) EXPECT_LE(g2.edge[j], g2.edge[j+1]);
//    }
//}

//TEST(testDataLoad, block_csr) {
//    std::string graph_path = "/mnt/d/graph/asgs_bin/patents";
//    Graph g2;
//    g2.load_global_graph(graph_path);
//    g2.blockType = BlockType::RANDOM_BLOCK;
//    g2.to_block_csr(32*1024,graph_path,5);
//}
//
//TEST(testDataload, read_block_csr) {
//    std::string graph_path = "/mnt/d/graph/asgs_bin/patents_blocks/";
//    for (int i=0 ;i<2 ;i++) {
//        std::string p = graph_path + std::to_string(i) + std::string("_br");
//        printf("%s",graph_path.c_str());
//        int v=0;
//        unsigned int e=0;
//        DataLoader::load_data_size(v,e,p);
//        printf("%d %d %d\n",i,v,e);
//    }
//}

TEST(testDataLoad, block_csr_correctness) {
    std::string graph_path = "/mnt/d/graph/asgs_bin/patents";
    Graph g2(graph_path);
    g2.load_global_graph(graph_path);
    g2.blockType = BlockType::RANDOM_BLOCK;
    g2.init_extern_storage(g2.v_cnt,g2.e_cnt);
    for (auto i : g2.intra_vertex_dict) {
        std::vector<loader> v;
        v.push_back({1,i.first,g2.v_state_map[i.first].file_id});
        printf("vector %zu\n",v.size());
        printf("%d\n",g2.v_state_map[i.first].file_id);
        g2.load_extern_data(g2.v_state_map[i.first].file_id,v);
        ASSERT_NE(g2.inter_vertex_dict.find(i.first),g2.inter_vertex_dict.end());
        int l = g2.inter_vertex[g2.inter_vertex_dict[i.first]];
        int r;
        printf("%d\n",g2.inter_vertex_dict[i.first]);
        if (g2.inter_vertex_dict[i.first] == g2.inter_vtx_num-1) r = g2.inter_edge_num;
        else r = g2.inter_vertex[g2.inter_vertex_dict[i.first]+1];

        int l0 = g2.vertex[g2.intra_vertex_dict[i.first]];
        int r0;
        if (g2.intra_vertex_dict[i.first] == g2.v_cnt-1) r0 = g2.e_cnt;
        else r0 = g2.vertex[g2.intra_vertex_dict[i.first]+1];
//        printf("%d %d\n",l,r);
        EXPECT_EQ(r-l,r0-l0);
//        for (int j =0; j< r-l; j++) EXPECT_EQ(g2.edge[j+l0],g2.inter_edge[j+l]);
    }
}
//TEST(testDataLoad, partition_csr) {
//
//    std::string graph_path = "/mnt/d/graph/asgs_bin/patents";
//    Graph g2;
//    g2.load_global_graph(graph_path);
//    int part = 4;
//    g2.to_partition_csr(PartitionType::RANDOM,part,graph_path);
//    int v=0;
//    int e=0;
//    for (int p=0; p< part; p++) {
//    Graph g;
//    g.g_vcnt = g2.v_cnt;
//    g.load_partition_graph(p,part, graph_path);
//    v+=g.v_cnt;
//    e+=g.e_cnt;
//    }
//    EXPECT_EQ(v,g2.v_cnt);
//    EXPECT_EQ(e,g2.e_cnt);
//}