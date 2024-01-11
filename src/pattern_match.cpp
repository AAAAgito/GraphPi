#include "graph.h"
#include "dataloader.h"
#include "common.h"
#include <iostream>
#include "pattern.h"
#include <fstream>
#include <sys/resource.h>
#include <unistd.h>

void init(Graph &g) {
    for (int i = 0; i < g.g_vcnt - 1; i++) {
        unsigned int l, r;
        g.get_mmp_edge_index(i, l, r);
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
}
int main(int argc, char **argv){
    
    if (argc < 2) {
        printf("param err\n");
        return 0;
    }
    // struct rlimit rl;
    // // Set the limit, in bytes
    // rl.rlim_cur = 1024 * 1024 * 10; // 50 MB
    // rl.rlim_max = 1024 * 1024 * 20; // 100 MB

    // // Apply the limit
    // if (setrlimit(RLIMIT_RSS, &rl) == -1) {
    //     perror("setrlimit");
    //     return 1;
    // }
    const char *graph_p = argv[1];
    int query = std::atoi(argv[2]);
    int setting = std::atoi(argv[3]);
    int enable_thread = std::atoi(argv[4]);
    
    PatternType t[5];
    t[0]=P1;
    t[1]=P2;
    t[2]=P3;
    t[3]=P4;
    t[4]=P5;
    Pattern p(t[query-1]);
    
    std::string graph_path(graph_p);
    Graph g(graph_path);

    
    int fd = g.memory_map();
    init(g);
    // g.load_global_graph(graph_path);
    g.available_threads = enable_thread;
    bool valid;
    Schedule s(p, valid, 2, 2, false, g.g_vcnt, g.g_ecnt,0);
    double t1 = get_wall_time();
    unsigned long long res = g.pattern_matching_oc(s,enable_thread);
    double t2 = get_wall_time();
    // res = g.pattern_matching(s,enable_thread);
    printf("%lld, %.6lf\n",res, t2 - t1);
    g.free_map(fd);
}