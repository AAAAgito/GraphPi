#include "dataloader.h"
#include "graph.h"


int main(int argc, char** argv) {
    if(argc<5) {
        printf("error! Invalid path provided");
        return 0;
    }
    char *input_path = argv[1];
    std::string graph_path(input_path);
    double insert_rate = std::stod(argv[2]);
    double delete_rate = std::stod(argv[3]);
    int split = std::stoi(argv[4]);
    int threads = 1;
    if (argc == 6) threads = std::stoi(argv[5]);
    Graph g(graph_path);
    g.memory_map();
    g.available_threads = threads;
    std::vector<int> insertion, deletion;
    g.generate_request(insert_rate,delete_rate,insertion,deletion);
    g.generate_patch_request(insert_rate,delete_rate,split);
    return 0;
}