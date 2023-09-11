#include "gtest/gtest.h"
#include "dataloader.h"
#include "common.h"

// std::string test_graph("friendster");
std::string test_graph("Patents");
//single function test


TEST(testrunning, hello) {
    printf("hello world\n");
    std::cout << __cplusplus << std::endl;
}

TEST(unit_test, generate_update_request) {
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
    graph_path.append(test_graph);
    Graph g(graph_path);
    g.memory_map();
    std::vector<int> insertion, deletion;
    g.generate_request(0.1,0.1,insertion,deletion);
    
    printf("%d %d\n",insertion.size(),deletion.size());
    for (int i=0; i< insertion.size()/2;i++) {
        ASSERT_LE(insertion[2*i],insertion[2*i+1]);
    }
    for (int i=0; i< deletion.size()/2;i++) {
        ASSERT_LE(deletion[2*i],deletion[2*i+1]);
    }
}

TEST(unit_test, verify_request) {
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
    graph_path.append(test_graph);
    Graph g(graph_path);
    g.memory_map();
    std::string insert_p(graph_path);
    insert_p.append("_insert.edge");
    std::string delete_p(graph_path);
    delete_p.append("_delete.edge");
    std::vector<int> data;
    int size = DataLoader::load_data_aggregate(data,insert_p)/2;
    for (int i=0;i<size;i++) {
        int a = data[i*2];
        int b = data[i*2+1];
        unsigned int l,r;
        bool found = false;
        g.get_mmp_edge_index(a,l,r);
        for (int j=l;j<r;j++) {
            if (g.mmp_edge[j]==b) found = true;
        }
        ASSERT_EQ(found,true);
        found = false;
        
        g.get_mmp_edge_index(b,l,r);
        for (int j=l;j<r;j++) {
            if (g.mmp_edge[j]==a) found = true;
        }
        ASSERT_EQ(found,true);
    }
    
    data.clear();
    size = DataLoader::load_data_aggregate(data,delete_p)/2;
    for (int i=0;i<size;i++) {
        int a = data[i*2];
        int b = data[i*2+1];
        unsigned int l,r;
        bool found = false;
        g.get_mmp_edge_index(a,l,r);
        for (int j=l;j<r;j++) {
            if (g.mmp_edge[j]==b) found = true;
        }
        ASSERT_EQ(found,true);
        found = false;
        
        g.get_mmp_edge_index(b,l,r);
        for (int j=l;j<r;j++) {
            if (g.mmp_edge[j]==a) found = true;
        }
        ASSERT_EQ(found,true);
    }
}

TEST(unit_test, insert_edge) {
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
    graph_path.append(test_graph);
    Graph g(graph_path);
    g.memory_map();
    int insert_size = 100;
    for (int i=0; i<insert_size; i++) {
        g.insert_edge(1,100000+i);
    }
    int increment = 0;
    for (auto i : g.patch_v) increment += i.insertion.size();
    ASSERT_EQ(increment,2*insert_size);

}

TEST(unit_test, update) {
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
    graph_path.append(test_graph);
    Graph g(graph_path);
    g.memory_map();
    std::string insert_p(graph_path);
    insert_p.append("_insert.edge");
    std::string delete_p(graph_path);
    delete_p.append("_delete.edge");
    std::vector<int> ins_set, del_set;
    
    int size = DataLoader::load_data_aggregate(del_set,insert_p);
    g.set_threads(20);
    g.update(ins_set,del_set);
    int del_sum=0;
    for (int i=0; i<g.g_vcnt;i++) {
        auto dels = g.patch_v[i].deletion;
        del_sum += dels.size();
    }
    EXPECT_EQ(size, del_sum);

    double tg1 = get_wall_time();
    g.dump_v();
    double tg2 = get_wall_time();
    printf("\nupdate time: %.6lf\n", tg2 - tg1);
    int prev_ecnt = g.g_ecnt;

    for (int i=0; i<g.g_vcnt;i++) {
        auto dels = g.patch_v[i].deletion;
        int l,r;
        g.get_delete_code_index(i,l,r);
        int il,ir;
        g.get_insert_code_index(i,il,ir);
        ASSERT_EQ(ir-il,0);
        int del_degree = 0;
        for (int j=l+1;j<r;j+=2) {
            del_degree += g.patch_delete_code[j];
        }
        for (auto del_v : dels) {
            bool exist = false;
            unsigned int vl,vr;
            g.get_mmp_edge_index(i,vl,vr);
            for (int j=l;j<r;j+=2) {
                int last_len = g.patch_delete_code[j+1];
                for (int last=0;last<last_len;last++)
                    if (g.mmp_edge[g.patch_delete_code[j]+vl+last]==del_v) exist = true;
            }
            ASSERT_TRUE(exist);
        }
        ASSERT_EQ(del_degree,dels.size());
    }
    ASSERT_EQ(g.insert_cnt,0);
    ASSERT_EQ(g.delete_cnt,del_set.size());

    double t1 = get_wall_time();
    g.refine_graph();
    double t2 = get_wall_time();
    printf("\nrefine time: %.6lf\n", t2 - t1);

    g.memory_map();
    int now_ecnt = g.g_ecnt;
    ASSERT_EQ(del_set.size(),prev_ecnt-now_ecnt);

    // test insertion

}

// TEST(unit_test, dump_v) {
//     printf("hello world\n");
//     std::cout << __cplusplus << std::endl;
// }

// TEST(unit_test, init_patch_ptr) {
//     printf("hello world\n");
//     std::cout << __cplusplus << std::endl;
// }

// TEST(unit_test, init_patch_ptr) {
//     printf("hello world\n");
//     std::cout << __cplusplus << std::endl;
// }

// TEST(unit_test, release_patch_ptr) {
//     printf("hello world\n");
//     std::cout << __cplusplus << std::endl;
// }

// TEST(unit_test, refine_graph) {
//     printf("hello world\n");
//     std::cout << __cplusplus << std::endl;
// }

// TEST(unit_test, bit_map) {
//     printf("hello world\n");
//     std::cout << __cplusplus << std::endl;
// }

// TEST(testDataLoad, global_csr) {
//     std::string raw_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/friendster_input";
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/friendster";
//     Graph *g;
//     DataLoader D;
//     double t1 = get_wall_time();
//     D.load_data(g,DataType::Patents,raw_path.c_str());
//     double t2 = get_wall_time();
//     printf("raw load time %.6lf\n", t2-t1);
//     g->to_global_csr(graph_path);
//     Graph g2;
//     t1 = get_wall_time();
//     g2.load_global_graph(graph_path);
//     t2 = get_wall_time();
//     printf("csr binary load time %.6lf\n", t2-t1);
//     EXPECT_EQ(g2.v_cnt,g->v_cnt);
//     EXPECT_EQ(g2.e_cnt,g->e_cnt);
//     for (int v = 0; v < g->v_cnt; v++) {
//         unsigned int l1,r1;
//         unsigned int l2,r2;
//         g->get_edge_index(v,l1,r1);
//         g2.get_edge_index(v,l2,r2);
//         ASSERT_EQ(l1,l2);
//         ASSERT_EQ(r1,r2);
//         for (int i=l1 ; i<r1; i++) {
//             ASSERT_EQ(g->edge[i],g2.edge[i]);
//         }
//     }
// }


// TEST(test_counting, oc) {
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
//     graph_path.append(test_graph);
//     Graph g(graph_path);
//     Graph g2(graph_path);
//     double t1= get_wall_time();
//     g.load_global_graph(graph_path);
    
//     g2.memory_map();
//     for (int i=0;i<g.v_cnt;i++) {
//         unsigned int l1,l2,r1,r2;
//         g.get_edge_index(i,l1,r1);
//         g2.get_mmp_edge_index(i,l2,r2);
//         ASSERT_EQ(l1,l2);
//         ASSERT_EQ(r1,r2);
//         for (unsigned int j=l1;j<r1;j++) {
//             ASSERT_EQ(g.edge[j],g2.mmp_edge[j]);
//         }
//     }
//     for (int i=0; i < 100;i++)
//         g2.insert_edge(1,100000+i);
//     g2.dump_v();
//     double t2= get_wall_time();
//     printf("time %.6lf\n",t2-t1);
//     printf("loaded %dv %llde\n",g.v_cnt,g.e_cnt);
// }

TEST(testMatching, match_in_partition_c) {

    // patterns:
    Pattern tc_pattern(3);
    tc_pattern.add_edge(0, 1);
    tc_pattern.add_edge(1, 2);
    tc_pattern.add_edge(2, 0);
    Pattern house(House);

    // setting
    bool is_pattern_valid;
    bool use_in_exclusion_optimize = false;
    int performance_type = 2;
    int restricts_type = 2;
    int thread_num = 32;
    //global matching
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
    graph_path.append(test_graph);
    Graph g(graph_path);
    g.load_global_graph(graph_path);
    g.g_vcnt = g.v_cnt;
    for (int i = 0; i < g.v_cnt - 1; i++) {
        unsigned int l, r;
        g.get_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
    g.max_degree = VertexSet::max_intersection_size;
    Schedule tc_schedule(tc_pattern, is_pattern_valid, performance_type, restricts_type, use_in_exclusion_optimize, g.v_cnt, g.e_cnt);
    
    unsigned long long ground_truth_result = 0;
    unsigned long long oc_result = 0;

    double t1 = get_wall_time();
    ground_truth_result = g.pattern_matching(tc_schedule, thread_num);
    double t2 = get_wall_time();
    printf("\nglobal time: %.6lf\n", t2 - t1);
    printf("count answer %lld\n",ground_truth_result);

    Graph g2(graph_path);
    int fd = g2.memory_map();
    t1 = get_wall_time();
    oc_result = g2.pattern_matching_oc(tc_schedule, thread_num);
    t2 = get_wall_time();
    printf("out of core time: %.6lf\n", t2 - t1);

    printf("oc count %lld\n",oc_result);

    g2.free_map(fd);

}
