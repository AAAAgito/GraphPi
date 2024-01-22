#include "dataloader.h"
#include "graph.h"

int main(int argc, char** argv) {
    if(argc!=2) {
        printf("error! Invalid path provided");
        return 0;
    }
    char *input_path = argv[1];
    std::string graph_path(input_path);
    Graph g(graph_path);
    g.load_global_graph(graph_path);
    printf("export\n");
    g.to_global_csr(graph_path);
    return 0;
}