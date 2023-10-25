#include "dataloader.h"
#include "graph.h"

int main(int argc, char** argv) {
    // printf("a\n");
    if(argc!=5) {
        printf("error! Invalid path provided");
        return 0;
    }
    char *input_path = argv[1];
    char *export_path = argv[2];
    int v = std::atoi(argv[3]);
    // v = DataLoader::large_graph_v(input_path);
    size_t mem = std::atoll(argv[4]);
    printf("input: %s\nexport: %s\n",input_path,export_path);
    DataLoader::load_large_graph(input_path,export_path,v,mem);
    return 0;
}