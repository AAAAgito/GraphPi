#include "../include/dataloader.h"
#include "../include/graph.h"
#include "../include/vertex_set.h"
#include "../include/common.h"
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <random>
#include <iostream>

bool DataLoader::load_data(Graph* &g, DataType type, const char* path, int oriented_type) {
    if(type == Patents || type == Orkut || type == complete8 || type == LiveJournal || type == MiCo || type == CiteSeer || type == Wiki_Vote) {
        return general_load_data(g, type, path, oriented_type);
    }

    printf("invalid DataType!\n");
    return false;
}

bool DataLoader::general_load_data(Graph *&g, DataType type, const char* path, int oriented_type) {
    if (freopen(path, "r", stdin) == NULL)
    {
        printf("File not found. %s\n", path);
        return false;
    }
    // printf("Load begin in %s\n",path);
    g = new Graph();

    //load triangle counting information
    switch(type) {
        case DataType::Patents : {
            g->tri_cnt = Patents_tri_cnt;
            break;
        }
        case DataType::LiveJournal : {
            g->tri_cnt = LiveJournal_tri_cnt;
            break;
        }
        case DataType::MiCo : {
            g->tri_cnt = MiCo_tri_cnt;
            break;
        }
        case DataType::CiteSeer : {
            g->tri_cnt = CiteSeer_tri_cnt;
            break;
        }
        case DataType::Wiki_Vote : {
            g->tri_cnt = Wiki_Vote_tri_cnt;
            break;
        }
        case DataType::Orkut : {
            g->tri_cnt = Orkut_tri_cnt;
            break;
        }
        default : {
            g->tri_cnt = -1;
            break;
        }
    }

    scanf("%d%u",&g->v_cnt,&g->e_cnt);
    int* degree = new int[g->v_cnt];
    memset(degree, 0, g->v_cnt * sizeof(int));
    g->e_cnt *= 2;
    std::pair<int,int> *e = new std::pair<int,int>[g->e_cnt];
    id.clear();
    int x,y;
    int tmp_v;
    unsigned int tmp_e;
    tmp_v = 0;
    tmp_e = 0;
    while(scanf("%d%d",&x,&y)!=EOF) {
        if(x == y) {
            // printf("find self circle\n");
            g->e_cnt -=2;
            continue;
            //return false;
        }
        if(!id.count(x)) id[x] = tmp_v ++;
        if(!id.count(y)) id[y] = tmp_v ++;
        x = id[x];
        y = id[y];
        e[tmp_e++] = std::make_pair(x,y);
        e[tmp_e++] = std::make_pair(y,x);
        ++degree[x];
        ++degree[y];
        //if(tmp_e % 1000000u == 0u) {
        //    printf("load %u edges\n",tmp_e);
        //    fflush(stdout);
        //}
    }


    // oriented_type == 0 do nothing
    //               == 1 high degree first
    //               == 2 low degree first
    if ( oriented_type != 0 ) {
        std::pair<int,int> *rank = new std::pair<int,int>[g->v_cnt];
        int *new_id = new int[g->v_cnt];
        for(int i = 0; i < g->v_cnt; ++i) rank[i] = std::make_pair(i,degree[i]);
        if( oriented_type == 1) std::sort(rank, rank + g->v_cnt, cmp_degree_gt);
        if( oriented_type == 2) std::sort(rank, rank + g->v_cnt, cmp_degree_lt);
        for(int i = 0; i < g->v_cnt; ++i) new_id[rank[i].first] = i;
        for(unsigned int i = 0; i < g->e_cnt; ++i) {
            e[i].first = new_id[e[i].first];
            e[i].second = new_id[e[i].second];
        }
        delete[] rank;
        delete[] new_id;
    }
    std::sort(degree, degree + g->v_cnt);

    // The max size of intersections is the second largest degree.
    //TODO VertexSet::max_intersection_size has different value with different dataset, but we use a static variable now.
    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, degree[g->v_cnt - 2]);
    g->max_degree = degree[g->v_cnt - 1];
    delete[] degree;
    if(tmp_v != g->v_cnt) {
        printf("vertex number error!\n");
    }
    if(tmp_e != g->e_cnt) {
        printf("edge number error!\n");
    }
    if(tmp_v != g->v_cnt || tmp_e != g->e_cnt) {
        fclose(stdin);
        delete g;
        delete[] e;
        return false;
    }
    std::sort(e,e+tmp_e,cmp_pair);
    g->e_cnt = unique(e,e+tmp_e) - e;
    for(unsigned int i = 0; i < g->e_cnt - 1; ++i)
        if(e[i] == e[i+1]) {
            printf("have same edge\n");
            fclose(stdin);
            delete g;
            delete[] e;
            return false;
        }
    g->edge = new int[g->e_cnt];
    g->vertex = new unsigned int[g->v_cnt + 1];
    bool* have_edge = new bool[g->v_cnt];
    int lst_v = -1;
    for(int i = 0; i < g->v_cnt; ++i) have_edge[i] = false;
    for(unsigned int i = 0; i < g->e_cnt; ++i) {
        if(e[i].first != lst_v) {
            have_edge[e[i].first] = true;
            g->vertex[e[i].first] = i;
        }
        lst_v = e[i].first;
        g->edge[i] = e[i].second;
    }
    delete[] e;
    // printf("Success! There are %d nodes and %u edges.\n",g->v_cnt,g->e_cnt);
    fflush(stdout);
    g->vertex[g->v_cnt] = g->e_cnt;
    for(int i = g->v_cnt - 1; i >= 0; --i)
        if(!have_edge[i]) {
            g->vertex[i] = g->vertex[i+1];
        }
    delete[] have_edge;

    return true;
}



bool DataLoader::load_complete(Graph* &g, int clique_size) {
    g = new Graph();

    g->v_cnt = clique_size;
    g->e_cnt = clique_size * (clique_size - 1) / 2;

    int* degree = new int[g->v_cnt];
    memset(degree, 0, g->v_cnt * sizeof(int));
    g->e_cnt *= 2;
    std::pair<int,int> *e = new std::pair<int,int>[g->e_cnt];
    id.clear();
    int tmp_v;
    unsigned int tmp_e;
    tmp_v = 0;
    tmp_e = 0;
    for(int i = 0; i < clique_size; ++i)
        for(int j = 0; j < i; ++j) {
            int x = i, y = j;
            if(!id.count(x)) id[x] = tmp_v ++;
            if(!id.count(y)) id[y] = tmp_v ++;
            x = id[x];
            y = id[y];
            e[tmp_e++] = std::make_pair(x,y);
            e[tmp_e++] = std::make_pair(y,x);
            ++degree[x];
            ++degree[y];
        }
    std::sort(degree, degree + g->v_cnt);

    // The max size of intersections is the second largest degree.
    //TODO VertexSet::max_intersection_size has different value with different dataset, but we use a static variable now.
    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, degree[g->v_cnt - 2]);
    g->max_degree = degree[g->v_cnt - 1];
    delete[] degree;
    if(tmp_v != g->v_cnt) {
        printf("vertex number error!\n");
    }
    if(tmp_e != g->e_cnt) {
        printf("edge number error!\n");
    }
    if(tmp_v != g->v_cnt || tmp_e != g->e_cnt) {
        fclose(stdin);
        delete g;
        delete[] e;
        return false;
    }
    std::sort(e,e+tmp_e,cmp_pair);
    g->e_cnt = unique(e,e+tmp_e) - e;
    g->edge = new int[g->e_cnt];
    g->vertex = new unsigned int[g->v_cnt + 1];
    bool* have_edge = new bool[g->v_cnt];
    int lst_v = -1;
    for(int i = 0; i < g->v_cnt; ++i) have_edge[i] = false;
    for(unsigned int i = 0; i < g->e_cnt; ++i) {
        if(e[i].first != lst_v) {
            have_edge[e[i].first] = true;
            g->vertex[e[i].first] = i;
        }
        lst_v = e[i].first;
        g->edge[i] = e[i].second;
    }
    delete[] e;
    g->vertex[g->v_cnt] = g->e_cnt;
    for(int i = g->v_cnt - 1; i >= 0; --i)
        if(!have_edge[i]) {
            g->vertex[i] = g->vertex[i+1];
        }
    delete[] have_edge;
    return true;
}

bool DataLoader::cmp_pair(std::pair<int,int>a, std::pair<int,int>b) {
    return a.first < b.first || (a.first == b.first && a.second < b.second);
}

bool DataLoader::cmp_degree_gt(std::pair<int,int> a,std::pair<int,int> b) {
    return a.second > b.second;
}

bool DataLoader::cmp_degree_lt(std::pair<int,int> a,std::pair<int,int> b) {
    return a.second < b.second;
}

long long DataLoader::comb(int n, int k) {
    long long ans = 1;
    for(int i = n; i > n - k; --i)
        ans = ans * i;
    for(int i = 1; i <= k; ++i)
        ans = ans / k;
    return ans;
}

void DataLoader::gen_partition_file(int v_cnt, unsigned int e_cnt, unsigned int *v, int *e, const std::string &path) {
    std::fstream binaryIo;
    binaryIo.open(path, std::ios::out| std::ios::binary | std::ios::trunc);
    binaryIo.seekp(0);

    binaryIo.write((char*)v, v_cnt * sizeof(v[0]));
    binaryIo.write((char*)e, e_cnt * sizeof(e[0]));
    binaryIo.close();

}

void DataLoader::load_partition_data(int &v_cnt, unsigned int &e_cnt, unsigned int *v, int *e, const std::string& path) {
    std::fstream binaryIo;

    binaryIo.open(path, std::ios::in | std::ios::binary);
    binaryIo.read((char *)v, v_cnt * sizeof(int));
    binaryIo.read((char *)e, e_cnt * sizeof(int));
    binaryIo.close();
}

void DataLoader::load_data_size(int &v_cnt, unsigned int &e_cnt, const std::string& path) {
    std::fstream binaryIo;

    binaryIo.open(path, std::ios::in | std::ios::binary);
    binaryIo.read(reinterpret_cast<char *>(&v_cnt), sizeof(int)); // read the number of elements
    binaryIo.read(reinterpret_cast<char *>(&e_cnt), sizeof(int));
    binaryIo.close();
}

void DataLoader::load_data_size(int &v_cnt, size_t &e_cnt, const std::string& path) {
    std::fstream binaryIo;

    binaryIo.open(path, std::ios::in | std::ios::binary);
    binaryIo.read(reinterpret_cast<char *>(&v_cnt), sizeof(int)); // read the number of elements
    binaryIo.read(reinterpret_cast<char *>(&e_cnt), sizeof(size_t));
    binaryIo.close();
}

void DataLoader::load_data_size(int &v_cnt, int &e_cnt, const std::string& path) {
    std::fstream binaryIo;

    binaryIo.open(path, std::ios::in | std::ios::binary);
    binaryIo.read(reinterpret_cast<char *>(&v_cnt), sizeof(int)); // read the number of elements
    binaryIo.read(reinterpret_cast<char *>(&e_cnt), sizeof(int));
    binaryIo.close();
}

void DataLoader::gen_data_size(int &v_cnt, unsigned int &e_cnt, const std::string& path) {
    std::fstream binaryIo;
    binaryIo.open(path, std::ios::out| std::ios::binary | std::ios::trunc);
    binaryIo.seekp(0);
    binaryIo.write((char *)&v_cnt, sizeof(int));
    binaryIo.write((char *)&e_cnt, sizeof(int));
    binaryIo.close();
}
void DataLoader::gen_data_size(int &v_cnt, int &e_cnt, const std::string& path) {
    std::fstream binaryIo;
    binaryIo.open(path, std::ios::out| std::ios::binary | std::ios::trunc);
    binaryIo.seekp(0);
    binaryIo.write((char *)&v_cnt, sizeof(int));
    binaryIo.write((char *)&e_cnt, sizeof(int));
    binaryIo.close();
}
void DataLoader::gen_data_size(int &v_cnt, size_t &e_cnt, const std::string& path) {
    std::fstream binaryIo;
    binaryIo.open(path, std::ios::out| std::ios::binary | std::ios::trunc);
    binaryIo.seekp(0);
    binaryIo.write((char *)&v_cnt, sizeof(int));
    binaryIo.write((char *)&e_cnt, sizeof(size_t));
    binaryIo.close();
}

int DataLoader::load_data_aggregate(std::vector<int> &data, const std::string& path) {
    std::ifstream binaryIo;
    binaryIo.open(path, std::ios::binary);
    binaryIo.seekg(0, std::ios_base::end);

	std::streampos size = binaryIo.tellg()/sizeof(int);

	binaryIo.seekg(0, std::ios_base::beg);
    data.resize(size);
    binaryIo.read((char*)data.data(), size*sizeof(int));
    binaryIo.close();
    return size;
}


void DataLoader::gen_data_file(int *v,int size,const std::string& path) {

    std::ofstream binaryIo(path.c_str(), std::ios::binary);

    binaryIo.write((char*)v, size * sizeof(v[0]));
    binaryIo.close();
}

int DataLoader::open_mmp_r(const std::string& path, size_t size_int, int *ptr) {
    int fd = open(path.c_str(), O_RDWR);
    lseek (fd, size_int*sizeof(int)-1, SEEK_SET);
    write (fd, "", 1);
    ptr = (int *)mmap(NULL, size_int*sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
    return fd;
}

int DataLoader::open_mmp_w(const std::string& path, size_t size_int, int *ptr) {
    int fd = open(path.c_str(), O_RDWR| O_CREAT,0777);
    lseek (fd, size_int*sizeof(int)-1, SEEK_SET);
    write (fd, "", 1);
    ptr = (int *)mmap(NULL, size_int*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    return fd;
}

void DataLoader::close_mmp(int fd, int size, int *ptr) {
    munmap(ptr, size*sizeof(int));
    close(fd);
}

int DataLoader::large_graph_v(const std::string& raw) {
    
    if (freopen(raw.c_str(), "r", stdin) == NULL)
    {
        printf("File not found. %s\n", raw.c_str());
        return -1;
    }
    int x;
    int l=0;
    while(scanf("%d",&x)!=EOF) {
        l=std::max(l,x);
    }
    return l;
}

void DataLoader::load_large_graph(const std::string& raw, const std::string& bin, int v, size_t mem_limit) {
    if (freopen(raw.c_str(), "r", stdin) == NULL)
    {
        printf("File not found. %s\n", raw.c_str());
        return ;
    }
    int x,y;
    std::vector<std::vector<int>> adj_list;
    adj_list.resize(v);
    size_t global_count=0;
    size_t count=0;
    int sub_file=0;
    printf("start scanning\n");
    std::cout << "start scanning\n";
    double t1 = get_wall_time();
    mem_limit /= sizeof(int)*2;
    while(scanf("%d%d",&x,&y)!=EOF) {
        if (x==y) continue;
        adj_list[x].push_back(y);
        adj_list[y].push_back(x);
        count+=1;
        global_count += 1;
        if (global_count % (1000*1000*1000) ==0) {
            // printf("scan time: %.6lf\n", get_wall_time() - t1);
            std::cout << "scan time: " << get_wall_time() - t1 << std::endl;
            t1=get_wall_time();
            // printf("process %dG\n",global_count/(1000*1000*1000));
            std::cout << "process " << global_count/(1000*1000*1000) << "G" << std::endl;
        }
        if (count>mem_limit) {
            printf("dump e %lld\n",count);
            std::string s(bin);
            s+="_";
            s.append(std::to_string(sub_file));
            size_t offset=0;
            std::vector<size_t> vs;
            for (auto i:adj_list) {
                vs.push_back(offset);
                offset+=i.size();
            }
            std::ofstream binaryIo(s.c_str(), std::ios::binary | std::ios::trunc);

            binaryIo.write((char*)vs.data(), vs.size() * sizeof(vs[0]));
            for (auto i:adj_list) {
                for (auto j:i) assert(j>0);
                binaryIo.write((char*)i.data(), i.size()*sizeof(i[0]));
            }
            binaryIo.close();
            std::string ssize(s);
            ssize+=".size";
            std::ofstream binaryIo2(ssize.c_str(), std::ios::binary | std::ios::trunc);
            
            binaryIo2.write((char*)&v, sizeof(v));
            binaryIo2.write((char*)&offset, sizeof(offset));
            binaryIo2.close();
            adj_list.clear();
            adj_list.resize(v);
            sub_file+=1;
            count=0;
        }
    }
    printf("%d of %lld reach the end\n",sub_file,global_count);
    std::string s(bin);
    s+="_";
    s.append(std::to_string(sub_file));
    size_t offset=0;
    std::vector<size_t> vs;
    for (auto i:adj_list) {
        vs.push_back(offset);
        offset+=i.size();
    }
    std::ofstream binaryIo(s.c_str(), std::ios::binary | std::ios::trunc);

    binaryIo.write((char*)vs.data(), vs.size() * sizeof(vs[0]));
    for (auto i:adj_list) {
        binaryIo.write((char*)i.data(), i.size()*sizeof(i[0]));
    }
    binaryIo.close();
    adj_list.clear();
    std::string ssize(s);
    ssize+=".size";
    std::ofstream binaryIo2(ssize.c_str(), std::ios::binary | std::ios::trunc);
    
    binaryIo2.write((char*)&v, sizeof(v));
    binaryIo2.write((char*)&offset, sizeof(offset));
    binaryIo2.close();
    sub_file+=1;
    merge_graph(bin,sub_file);
    merge_size(bin,sub_file);
}

void DataLoader::merge_size(const std::string &bin, int packs) {
    int v;
    size_t ecnt=0;
    for (int i=0;i<packs;i++) {
        std::string s(bin);
        s+="_";
        s+=std::to_string(i);
        std::string ss(s);
        ss+=".size";
        
        std::fstream sizeIo;
        size_t e;
        sizeIo.open(ss, std::ios::in | std::ios::binary);
        sizeIo.read(reinterpret_cast<char *>(&v), sizeof(int)); // read the number of elements
        sizeIo.read(reinterpret_cast<char *>(&e), sizeof(size_t));

        ecnt +=e;
    }
    std::string ss(bin);
    ss+=".size";
    gen_data_size(v,ecnt,ss);
}

void DataLoader::merge_graph(const std::string &bin, int packs) {
    std::vector<size_t*> vmap;
    std::vector<int*> emap;
    std::vector<int> fdmap;
    std::vector<size_t> es;
    size_t ecnt = 0;
    int v;
    printf("start merge\n");
    for (int i=0;i<packs;i++) {
        std::string s(bin);
        s+="_";
        s+=std::to_string(i);
        std::string ss(s);
        ss+=".size";
        
        std::fstream sizeIo;
        size_t e;
        sizeIo.open(ss, std::ios::in | std::ios::binary);
        sizeIo.read(reinterpret_cast<char *>(&v), sizeof(int)); // read the number of elements
        sizeIo.read(reinterpret_cast<char *>(&e), sizeof(size_t));
        
        
        int fd = open(s.c_str(), O_RDWR);
        fdmap.push_back(fd);
        printf("prev fd %d\n",fd);
        char *ptr = (char *)mmap(NULL, v*sizeof(size_t)+e*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        size_t *vptr = (size_t *)ptr;
        int *eptr = (int *)(vptr+v);
        // for (size_t j=0;j<e;j++){
        //     // if (eptr[j]<0)
        //     assert(eptr[j]>=0 &&eptr[j]<v);
        // }
        vmap.push_back(vptr);
        emap.push_back(eptr);
        ecnt +=e;
        printf("read e %lld\n",e);
        es.push_back(e);
        printf("v %d\n",v);
    }

    
    int fd = open(bin.c_str(), O_RDWR| O_CREAT,0777);
    lseek (fd, v*sizeof(size_t)+ecnt*sizeof(int)-1, SEEK_SET);
    write (fd, "", 1);
    char *ptr = (char *)mmap(NULL, v*sizeof(size_t)+ecnt*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    size_t *vptr = (size_t *)ptr;
    int *eptr = (int *)(vptr+v);
    size_t offset=0;
    printf("fd %d\n",fd);
    printf("ecnt %lld\n",ecnt);
    for (int vtx=0;vtx<v;vtx++) {
        vptr[vtx]=offset;
        if(vtx%(1000*1000)==0) printf("process %dM V\n",vtx/(1000*1000));
        // ptr[vtx] = offset;
        if (vtx != v-1)
            for (auto i: vmap) offset += i[vtx+1]-i[vtx];
        else for (int i=0;i<packs;i++) offset += es[i] - vmap[i][vtx];
        std::vector<int> edges;
        for (int i=0;i<packs;i++) {
            size_t l,r;
            if (vtx==v-1) r=es[i];
            else r=vmap[i][vtx+1];
            l=vmap[i][vtx];
            assert(r<=es[i]);
            // printf("L R %lld %lld\n",l,r);
            for (size_t a=l; a<r; a++) {
                // assert(emap[i][a]>=0);
                // assert(emap[i][a]<v);
                edges.push_back(emap[i][a]);
            }
        }
        std::sort(edges.begin(),edges.end());
        for (auto j:edges) assert(j>=0);
        // printf("%d %d\n",vptr[vtx]+edges.size(),offset);
        assert(vptr[vtx]+edges.size()==offset);
        memcpy(eptr+vptr[vtx],edges.data(),edges.size());
    }
    munmap(ptr,v*sizeof(size_t)+ecnt*sizeof(int));
    close(fd);
}