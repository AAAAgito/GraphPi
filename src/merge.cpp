#include "dataloader.h"
#include "graph.h"

int main(int argc, char** argv) {
    if(argc!=3) {
        printf("error! Invalid path provided");
        return 0;
    }
    char *bin = argv[1];
    int packs = std::atoi(argv[2]);
    DataLoader::merge_graph(std::string(bin),packs);
    DataLoader::merge_size(std::string(bin),packs);
    return 0;
}