#include "graph.h"
#include "dataloader.h"
#include "common.h"
#include <iostream>
#include "pattern.h"
#include <fstream>
void init(Graph &g) {
    for (int i = 0; i < g.g_vcnt - 1; i++) {
        unsigned int l, r;
        g.get_mmp_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
}

void initl(Graph &g) {
    for (int i = 0; i < g.g_vcnt - 1; i++) {
        size_t l, r;
        g.get_large_mmp_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
}
void init_in(Graph &g) {
    for (int i = 0; i < g.v_cnt - 1; i++) {
        unsigned int l, r;
        g.get_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
}

int main (int argc, char **argv) {
    if (argc < 2) {
        printf("param err\n");
        return 0;
    }
    const char *graph_p = argv[1];
    int enable_thread = std::atoi(argv[2]);
    // 1: matching 2: update 3: needle
    int experiment = std::atoi(argv[3]);
    // 1: 0=in-memory 1=out
    // 2: 0=wait rebuild 1=Draft
    int setting = std::atoi(argv[4]);
    int query = std::atoi(argv[5]);
    int split, tri_cnt=0;
    std::string log_path("/home/yanglaoyuan/AsyncSubGraphStorage/GraphPi/");
    log_path+=std::string(graph_p);
    log_path+=argv[3];
    log_path+=argv[4];
    log_path+=argv[5];
    log_path+=".log";
    std::ofstream MyFile(log_path.c_str());
    if (experiment>=1) {
        split = std::atoi(argv[6]);
    }
    else {
        tri_cnt= std::atoi(argv[6]);
    }
    double rate = (double)query/10;
    PatternType t[5];
    t[0]=P1;
    t[1]=P2;
    t[2]=P3;
    t[3]=P4;
    t[4]=P5;
    Pattern p(t[query-1]);
    Pattern tri(P1);
    bool valid;
    std::string graph_path(graph_p);
    Graph g(graph_path);
    if (experiment==0) {
        if (setting==0) {
            g.load_global_graph(graph_path);
            g.available_threads = enable_thread;
            init_in(g);
            Schedule s(p, valid, 1, 1, false, g.v_cnt, g.e_cnt,tri_cnt);
            
            double t1 = get_wall_time();
            unsigned long long res = g.pattern_matching(s,enable_thread);
            double t2 = get_wall_time();
            printf("%.6lf\n", t2 - t1);
            MyFile << t2-t1 << std::endl;
            // std::cout<<t2-t1<<std::endl;
        }
        else if (setting==1) {
            int fd = g.memory_map();
            init(g);
            g.available_threads = enable_thread;
            
            Schedule s(p, valid, 2, 2, false, g.g_vcnt, g.g_ecnt,tri_cnt);
            double t1 = get_wall_time();
            unsigned long long res = g.pattern_matching_oc(s,enable_thread);
            double t2 = get_wall_time();
            printf("%.6lf\n", t2 - t1);
            MyFile << t2-t1 << std::endl;
            g.free_map(fd);
        }
        else if (setting==2) {
            int fd = g.memory_lmap();
            initl(g);
            double t1 = get_wall_time();
            long long result = g.triangle_counting_mt(enable_thread);
            double t2 = get_wall_time();
            printf("GSH: %.6lf\n", t2 - t1);

            MyFile << "GSH: " << t2-t1 << std::endl;
            g.free_lmap(fd);
        }
    }
    if (experiment==1) {
        // printf("rate %.6lf\n",rate);
        int fd = g.memory_map();
        init(g);
        g.available_threads = enable_thread;
        
        Schedule s(tri, valid, 2, 2, false, g.g_vcnt, g.g_ecnt);
        g.generate_patch_request(rate,rate,split);
        // Draft
        if (setting==0) {
            double t1 = get_wall_time();
            g.refine_graph();
            g.memory_map();
            double t3 = get_wall_time();
            g.pattern_matching_oc(s,enable_thread);
            double t4 = get_wall_time();
            double t2 = get_wall_time();
            printf("%.6lf\n", t2 - t1);
            MyFile << t2-t1 << std::endl;
            printf("base: %.6lf\n", t4 - t3);
            MyFile << "Base:" << t4-t3 << std::endl;
        }
        if (setting==1) {
            int in,dn;
            std::string psize(g.raw_data_path);
            psize+="_patch.size";
            // printf("load size %s\n",psize.c_str());
            // printf("pi %s\n",g.patch_insert_path.c_str());
            // printf("pd %s\n",g.patch_delete_path.c_str());
            DataLoader::load_data_size(in,dn,psize);
            double t1 = get_wall_time();
            g.update_patch(g.patch_insert_path,g.patch_delete_path,in,dn);
            size_t result = g.pattern_matching_oca(s,enable_thread);
            double t2 = get_wall_time();
            printf("%.6lf\n", t2 - t1);
            MyFile << t2-t1 << std::endl;
        }
        g.free_map(fd);
    }
    MyFile.close();
}
