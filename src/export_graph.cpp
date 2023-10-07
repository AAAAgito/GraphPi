#include "dataloader.h"
#include "graph.h"

int main(int argc, char** argv) {
    if(argc!=3) {
        printf("error! Invalid path provided");
        return 0;
    }
    char *input_path = argv[1];
    char *export_path = argv[2];
    printf("input: %s\nexport: %s\n",input_path,export_path);
    Graph *g;
    DataLoader D;
    std::string graph_path(export_path);
    D.load_data(g,DataType::Patents,input_path);
    g->to_global_csr(graph_path);
    delete g;
    return 0;
}