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
        l=g.mmp_vertex[i];
        r=g.mmp_vertex[i+1];
        int len = r - l;
        VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size, len);
    }
}
int main(int argc, char **argv){
    
    if (argc < 8) {
        printf("param err\n");
        return 0;
    }
    // int a;
    // std::cin >> a;
    // std::cout << a;
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
    // 0: mmap(), 1: mmap() with prefetch, 2: load into memory
    int do_prefetch = std::atoi(argv[3]);
    int enable_thread = std::atoi(argv[4]);
    double threshold = std::atof(argv[5]);
    double density = std::atof(argv[6]);
    int interval = std::atoi(argv[7]);
    
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
    g.do_prefetch = do_prefetch==1;
    g.prefetch_interval = interval;
    g.KMeasureDecompose(density,threshold);
    int count=0;
    int sum=0;
    unsigned int max=0;
    for (int i=0;i<g.g_vcnt;i++) {
        if (g.KMD[2*i]!=-1) {
            count+=1;
            sum+=g.KMD[2*i+1];
            max = std::max(max,g.KMD[2*i+1]);
        }
    }
    printf("KMD count %d %d %lld\n",count, sum/count, max);
    // g.load_global_graph(graph_path);
    g.available_threads = enable_thread;
    bool valid;
    Schedule s(p, valid, 2, 2, false, g.g_vcnt, g.g_ecnt,0);
    double t1 = get_wall_time();
    unsigned long long res = g.pattern_matching_oc(s,enable_thread);
    double t2 = get_wall_time();
    // res = g.pattern_matching(s,enable_thread);
    printf("Result: %lld, Times: %.6lf\n",res, t2 - t1);
    g.free_map(fd);
}