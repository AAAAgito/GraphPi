#include "gtest/gtest.h"
#include "dataloader.h"
#include "common.h"

std::string test_graph("orkut");
// std::string test_graph("friendster");
// std::string test_graph("patents");
//single function test


TEST(testrunning, hello) {
    printf("hello world\n");
    std::cout << __cplusplus << std::endl;
    int test_size = 3*1000*1000;
    int *p = new int[test_size];
    double t1 = get_wall_time();
    for (int i=0; i<test_size;i++) {
        p[i] = i;
    }
    double t2 = get_wall_time();
    printf("%.6lf\n",t2-t1);
    std::vector<int> r;
    t1 = get_wall_time();
    r.reserve(test_size/2);
    for (int i=0; i<test_size;i++) {
        r.push_back(i);
    }
    t2 = get_wall_time();
    printf("%d\n",r[0]);
    printf("%.6lf\n",t2-t1);

    
    std::vector<std::vector<int>> pr;
    pr.resize(3774768);
    t1 = get_wall_time();
    for (int i=0;i<3000000;i++) {
        pr[i].reserve(3);
    }
    // pr[0].reserve(test_size/2);
    for (int i=0; i<test_size;i++) {
        pr[i%3000000].push_back(i);
        pr[i%3000000].push_back(i);
    }
    t2 = get_wall_time();
    printf("%d\n",pr[0][0]);
    printf("%.6lf\n",t2-t1);
}


// TEST(testDataLoad, global_csr) {
//     std::string raw_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/"+test_graph+"_input";
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/"+test_graph;
//     Graph *g;
//     DataLoader D;
//     double t1 = get_wall_time();
//     D.load_data(g,DataType::Patents,raw_path.c_str());
//     // g = new Graph;
//     // g->load_global_graph(graph_path);
//     double t2 = get_wall_time();
//     printf("raw load time %.6lf\n", t2-t1);
//     t1 = get_wall_time();
//     g->to_global_csr(graph_path);
//     t2 = get_wall_time();
//     Graph g2;
//     g2.load_global_graph(graph_path);
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
//         for (unsigned int i=l1 ; i<r1; i++) {
//             ASSERT_EQ(g->edge[i],g2.edge[i]);
//         }
//     }
//     delete g;
// }

TEST(test_counting, oc) {
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
    graph_path.append(test_graph);
    Graph g(graph_path);
    Graph g2(graph_path);
    double t1= get_wall_time();
    g.load_global_graph(graph_path);
    
    g2.memory_map();
    double t2= get_wall_time();
    printf("time %.6lf\n",t2-t1);
    for (int i=0;i<g.v_cnt;i++) {
        unsigned int l1,l2,r1,r2;
        g.get_edge_index(i,l1,r1);
        g2.get_mmp_edge_index(i,l2,r2);
        ASSERT_EQ(l1,l2);
        ASSERT_EQ(r1,r2);
        for (unsigned int j=l1;j<r1;j++) {
            ASSERT_EQ(g.edge[j],g2.mmp_edge[j]);
        }
    }
    printf("loaded %dv %llde\n",g.v_cnt,g.e_cnt);
}

// TEST(unit_test, generate_update_request) {
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
//     graph_path.append(test_graph);
//     Graph g(graph_path);
//     g.memory_map();
//     std::vector<int> insertion, deletion;
//     g.generate_request(0.1,0.1,insertion,deletion);
    
//     printf("%ld %ld\n",insertion.size(),deletion.size());
//     for (int i=0; i< (int)insertion.size()/2;i++) {
//         ASSERT_LE(insertion[2*i],insertion[2*i+1]);
//     }
//     for (int i=0; i< (int)deletion.size()/2;i++) {
//         ASSERT_LE(deletion[2*i],deletion[2*i+1]);
//     }
// }

// TEST(unit_test, verify_request) {
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
//     graph_path.append(test_graph);
//     Graph g(graph_path);
//     g.memory_map();
//     std::string insert_p(graph_path);
//     insert_p.append("_insert.edge");
//     std::string delete_p(graph_path);
//     delete_p.append("_delete.edge");
//     std::vector<int> data;
//     int size = DataLoader::load_data_aggregate(data,insert_p)/2;
//     for (int i=0;i<size;i++) {
//         int a = data[i*2];
//         int b = data[i*2+1];
//         unsigned int l,r;
//         bool found = false;
//         g.get_mmp_edge_index(a,l,r);
//         for (unsigned int j=l;j<r;j++) {
//             if (g.mmp_edge[j]==b) found = true;
//         }
//         ASSERT_EQ(found,true);
//         found = false;
        
//         g.get_mmp_edge_index(b,l,r);
//         for (unsigned int j=l;j<r;j++) {
//             if (g.mmp_edge[j]==a) found = true;
//         }
//         ASSERT_EQ(found,true);
//     }
    
//     data.clear();
//     size = DataLoader::load_data_aggregate(data,delete_p)/2;
//     for (int i=0;i<size;i++) {
//         int a = data[i*2];
//         int b = data[i*2+1];
//         unsigned int l,r;
//         bool found = false;
//         g.get_mmp_edge_index(a,l,r);
//         for (unsigned int j=l;j<r;j++) {
//             if (g.mmp_edge[j]==b) found = true;
//         }
//         ASSERT_EQ(found,true);
//         found = false;
        
//         g.get_mmp_edge_index(b,l,r);
//         for (unsigned int j=l;j<r;j++) {
//             if (g.mmp_edge[j]==a) found = true;
//         }
//         ASSERT_EQ(found,true);
//     }
// }

// TEST(unit_test, insert_edge) {
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
//     graph_path.append(test_graph);
//     Graph g(graph_path);
//     g.memory_map();
//     int insert_size = 100;
//     for (int i=0; i<insert_size; i++) {
//         g.insert_edge(1,100000+i);
//     }
//     int increment = 0;
//     for (auto i : g.patch_v) increment += i.insertion.size();
//     ASSERT_EQ(increment,2*insert_size);

// }

// TEST(unit_test, update) {
//     std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/";
//     graph_path.append(test_graph);
//     Graph g(graph_path);
//     Graph g2(graph_path);
//     g2.load_global_graph(graph_path);
//     g.memory_map();
//     printf("gsize: %d\n",g.g_ecnt);
//     std::string insert_p(graph_path);
//     insert_p.append("_insert.edge");
//     std::string delete_p(graph_path);
//     delete_p.append("_delete.edge");
//     std::vector<int> ins_set, del_set;
    
//     int size = DataLoader::load_data_aggregate(del_set,insert_p);
//     g.set_threads(1);
//     printf("ingest size %d %d\n",ins_set.size(),del_set.size());
//     double tg1 = get_wall_time();
//     g.update(ins_set,del_set);
//     double tg2 = get_wall_time();
//     printf("\n update time: %.6lf\n", tg2 - tg1);
    
//     tg1 = get_wall_time();
//     g.dump_v(0);
//     tg2 = get_wall_time();
//     printf("\n dump time: %.6lf\n", tg2 - tg1);
//     int prev_ecnt = g.g_ecnt;
//     for (int i=0; i<g.extra_v_cnt;i++) {
//         int li,ld,ri,rd;
//         li = g.mmp_patch_vertex[i];
//         ld = g.mmp_patch_delete_vertex[i];
//         if (i==g.extra_v_cnt-1) {
//             ri = g.insert_cnt;
//             rd = g.delete_cnt;
//         }
//         else {
//             ri = g.mmp_patch_vertex[i+1];
//             rd = g.mmp_patch_delete_vertex[i+1];
//         }
//         ASSERT_EQ(g.patch_vd[i].size(),rd-ld);
//         for (int j=0;j<ri-li;j++) {
//             ASSERT_EQ(g.patch_vi[i][j],g.mmp_patch_edge[j+li]);
//         }
//         for (int j=0;j<rd-ld;j++) {
//             ASSERT_EQ(g.patch_vd[i][j],g.mmp_patch_delete_edge[j+ld]);
//         }
//     }

//     ASSERT_EQ(g.insert_cnt,0);
//     ASSERT_EQ(g.delete_cnt,del_set.size());
//     // printf("aaaa\n");
//     double t1 = get_wall_time();
//     g.refine_graph();
//     double t2 = get_wall_time();
//     printf("\nrefine time: %.6lf\n", t2 - t1);

//     g.free_map(g.g_fd);
//     g.memory_map();
//     printf("gsize: %d\n",g.g_ecnt);
//     int now_ecnt = g.g_ecnt;
//     ASSERT_EQ(del_set.size(),prev_ecnt-now_ecnt);
//     for (int i=0; i<g.extra_v_cnt;i++) {
//         unsigned int l,r,ll,rr;
//         g2.get_edge_index(i,ll,rr);
//         g.get_mmp_edge_index(i,l,r);
//     }

//     printf("update\n");
    
//     int is = DataLoader::load_data_aggregate(ins_set,insert_p);
//     int ds = DataLoader::load_data_aggregate(del_set,delete_p);
//     g.update(ins_set,del_set);
    
//     int ins_sum=0,del_sum=0;
//     for (int i=0; i<g.g_vcnt;i++) {
//         auto in = g.patch_vi[i];
//         ins_sum += in.size();
//         auto de = g.patch_vd[i];
//         del_sum += de.size();
//     }
//     EXPECT_EQ(is, ins_sum);
//     EXPECT_EQ(ds, del_sum);


//     g.dump_v(0);
//     ASSERT_EQ(g.delete_cnt,del_set.size());
//     ASSERT_EQ(g.insert_cnt,ins_set.size());
//     t1 = get_wall_time();
//     g.refine_graph();
//     g.free_map(g.g_fd);
//     t2 = get_wall_time();
//     printf("\nrefine time: %.6lf\n", t2 - t1);
//     g.memory_map();
//     printf("gsize: %d\n",g.g_ecnt);
//     ASSERT_EQ((int)g.g_ecnt-now_ecnt,(int)ins_set.size()-del_set.size());

//     // test insertion

// }

TEST(testMatching, match_in_patch) {
    Pattern tc_pattern(3);
    tc_pattern.add_edge(0, 1);
    tc_pattern.add_edge(1, 2);
    tc_pattern.add_edge(2, 0);
    bool is_pattern_valid;
    bool use_in_exclusion_optimize = false;
    int performance_type = 2;
    int restricts_type = 2;
    int thread_num = 40;
    std::string graph_path = "/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/orkut";

    std::string insert_p(graph_path);
    insert_p.append("_insert.edge");
    std::string delete_p(graph_path);
    delete_p.append("_delete.edge");

    graph_path.append("_backup");
    
    Schedule tc_schedule(tc_pattern, is_pattern_valid, performance_type, restricts_type, use_in_exclusion_optimize, 0, 0);
    
    Graph g2(graph_path);
    int fd = g2.memory_map();
    for (int i = 0; i < g2.g_vcnt - 1; i++) {
        unsigned int l, r;
        g2.get_mmp_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
    std::vector<int> in,de;
    int inum = DataLoader::load_data_aggregate(in,insert_p);
    int dnum = DataLoader::load_data_aggregate(de,delete_p);
    printf("updates %d %d\n",inum,dnum);
    std::string p1(graph_path);
    p1.append("_patch.insert");
    std::string p2(graph_path);
    p2.append("_patch.delete");
    g2.update_patch(p1,p2,inum,dnum);
    double t1 = get_wall_time();
    long oc_result = g2.pattern_matching_oca(tc_schedule, thread_num);
    double t2 = get_wall_time();
    printf("out of core time: %.6lf\n", t2 - t1);

    printf("oc count %lld\n",oc_result);
}

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
    int thread_num = 40;
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
