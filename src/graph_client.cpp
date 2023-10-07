#include "graph.h"
#include "dataloader.h"
#include "common.h"
#include <iostream>

void init(Graph &g) {
    for (int i = 0; i < g.g_vcnt - 1; i++) {
        unsigned int l, r;
        g.get_mmp_edge_index(i, l, r);
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
    // 0: needle, 1: reconstruct first
    int mode = std::atoi(argv[3]);
    std::string graph_path(graph_p);
    Graph g(graph_path);
    g.memory_map();
    g.available_threads = enable_thread;
    init(g);
    
    if (mode==0)
        for (int i=4; i<argc; i+=3)
        {
            std::string ins(argv[i]);
            std::string del(argv[i+1]);
            std::string pattern(argv[i+2]);
            g.update_patch(ins,del,0,0);
            // match
            Pattern p(House);
            bool valid;
            Schedule s(p, valid, 2, 2, false, g.v_cnt, g.e_cnt);
            g.pattern_matching_oc(s,enable_thread);

            g.wait_refine();
        }
    else
        for (int i=4; i<argc; i+=3)
        {
            std::string ins(argv[i]);
            std::string del(argv[i+1]);
            std::string pattern(argv[i+2]);
            std::vector<int> ins_set, del_set;
            int isize = DataLoader::load_data_aggregate(ins_set,ins);
            int dsize = DataLoader::load_data_aggregate(del_set,del);
            g.update(ins_set,del_set);
            g.dump_v(0);
            g.refine_graph();
            g.memory_map();
            // match
            Pattern p(House);
            bool valid;
            Schedule s(p, valid, 2, 2, false, g.v_cnt, g.e_cnt);
            g.pattern_matching_oc(s,enable_thread);


        }
}
