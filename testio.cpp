
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <ostream>
#include <sys/time.h>
#include <istream>
#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <cassert>

double get_wall_time() {
    struct timeval time;
    if(gettimeofday(&time,NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

void gen_block_file(int *v,int size,const std::string& path) {

    std::ofstream binaryIo(path.c_str(), std::ios::binary);

    binaryIo.write((char*)v, size * sizeof(v[0]));
    binaryIo.close();
}


void load_block_data_aggregate(int *data, int size, const std::string& path) {
    int fd = open(path.c_str(), O_RDWR);
    int ret = read(fd, data, sizeof(int)* size);
    close(fd);
}

void load_block_data_aggregate2(int *data, int size, const std::string& path) {
    std::ifstream binaryIo(path.c_str(), std::ios::binary);
    binaryIo.read((char*)data, size * sizeof(int));
    binaryIo.close();
}
int main() {
    std::cout << __cplusplus << std::endl;
    int len = 1000;
    // srand((unsigned) time(NULL));
    std::string root_path("/home/yanglaoyuan/test/");
    // std::string root_path("/home/yanglaoyuan/AsyncSubGraphStorage/dataset/test/");
    int *v = new int[len];
    int *data = new int[len];
    for (int i=0;i<len;i++) v[i]=i;
    for (int i=0;i<30000;i++) {
        std::string path(root_path);
        int a = i;
        std::stringstream ss;
        ss << a;
        std::string str = ss.str();
        path += str;
        gen_block_file(v,len,path);
        // load_block_data_aggregate(data,len,path);
    }
    double total_time = 0;
    for (int i=0; i < 400*10000; i++) {
        std::string path(root_path);
        int random = rand();
        int a = (random%30000);
        std::stringstream ss;
        ss << a;
        std::string str = ss.str();
        path += str;
        const char *p = path.c_str();
        double t1= get_wall_time();
        load_block_data_aggregate(data,len,path);
        // load_block_data_aggregate3(data,len,p);
        double t2= get_wall_time();
        total_time+=t2-t1;
    }
    printf("time %.6lf\n",total_time);
    delete[] data;
}