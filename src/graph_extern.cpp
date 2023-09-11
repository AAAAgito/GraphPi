//
// Created by yanglaoyuan on 7/7/23.
//
#include "../include/graph.h"
#include "../include/graphmpi.h"
#include "../include/vertex_set.h"
#include "../include/common.h"
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <cstring>
#include <cmath> 
#include <atomic>
#include <queue>
#include <iostream>
#include <random>


void Graph::generate_request(double ins_rate, double del_rate, std::vector<int> &ins, std::vector<int> &del) {
    int FLOAT_MIN = 0;
    int FLOAT_MAX = 1;
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
    ins.clear();
    del.clear();
    for (int i=0; i<g_vcnt; i++) {
        // printf("%d %d\n",i,g_vcnt);
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        for (int j=l;j<r;j++){
            if (mmp_edge[j]<i) continue;
            float r = distr(eng);
            if (r < ins_rate){
                ins.push_back(i);
                ins.push_back(mmp_edge[j]);
            }
            else if (r>=ins_rate && r<ins_rate+del_rate) {
                del.push_back(i);
                del.push_back(mmp_edge[j]);
            }
        }
    }
    std::string ip(raw_data_path);
    ip.append("_insert.edge");
    DataLoader::gen_data_file(ins.data(),ins.size(),ip);
    
    std::string dp(raw_data_path);
    dp.append("_delete.edge");
    DataLoader::gen_data_file(del.data(),del.size(),dp);
}

int Graph::memory_map() {
    
    std::string graph_size(raw_data_path);
    graph_size.append(".size");
    DataLoader::load_data_size(g_vcnt,g_ecnt,graph_size);
    
    int fd = open(raw_data_path.c_str(),O_RDWR);
    mem = static_cast<int*>(mmap(NULL, (g_vcnt+g_ecnt)*sizeof(int), PROT_READ, MAP_SHARED ,fd,0));
    mmp_vertex = (unsigned int *)mem;
    mmp_edge = mem+g_vcnt;
    extra_v_cnt = g_vcnt;
    patch_v.resize(g_vcnt);
    bitmap = new int[g_vcnt/INT_BITS];
    return fd;
}

void Graph::free_map(int &fd) {
    munmap(mem, (g_vcnt+g_ecnt)*sizeof(int));
    close(fd);
    delete[] bitmap;
}

void Graph::insert_edge(int v1, int v2) {
    if (patch_v.size() <= v2) {
        patch_v.resize(v2+1);
        extra_v_cnt = v2+1;
    }
    if (patch_v.size() <= v1) {
        patch_v.resize(v1+1);
        extra_v_cnt = v1+1;
    }
    patch_v[v1].insertion.push_back(v2);
    patch_v[v2].insertion.push_back(v1);
}

void Graph::delete_edge(int v1, int v2) {
    patch_v[v1].deletion.push_back(v2);
    patch_v[v2].deletion.push_back(v1);
}

void Graph::update(std::vector<pair> &u, std::vector<pair> &d) {
    // TODO: openMP
#pragma omp parallel num_threads(available_threads)
    {
#pragma omp for schedule(dynamic) nowait
        for (int i=0; i<u.size();i++) {
        insert_edge(u[i].e1,u[i].e2);
        set(bitmap,u[i].e1);
        set(bitmap,u[i].e2);
        }
#pragma omp for schedule(dynamic) nowait
        for (int i=0; i<d.size();i++) {
            delete_edge(d[i].e1,d[i].e2);
            set(bitmap,d[i].e1);
            set(bitmap,d[i].e2);
        }
    }
}

void Graph::update(std::vector<int> &u, std::vector<int> &d) {
    // TODO: openMP
    for (int i=0; i<u.size()/2;i++) {
        insert_edge(u[2*i],u[2*i+1]);
        set(bitmap,u[2*i]);
        set(bitmap,u[2*i+1]);
    }
    for (int i=0; i<d.size()/2;i++) {
        delete_edge(d[2*i],d[2*i+1]);
        set(bitmap,d[2*i]);
        set(bitmap,d[2*i+1]);
    }
}

void Graph::dump_v() {
    release_patch_ptr();
    std::string pi(raw_data_path);
    pi.append("_patch.insert");
    std::string pd(raw_data_path);
    pd.append("_patch.delete");
    std::string picode(raw_data_path);
    picode.append("_patch_code.insert");
    std::string pdcode(raw_data_path);
    pdcode.append("_patch_code.delete");

    // sequential
    // TODO: add reading prev patch info (and merge)
    std::vector<int> ins_idx,del_idx;
    int ins=0,del=0;
    for (int i=0; i<extra_v_cnt; i++) {
        patch_data *p = &patch_v[i];
        ins_idx.push_back(ins);
        ins+=p->insertion.size();
        del_idx.push_back(del);
        del+=p->deletion.size();
    }
    insert_cnt = ins;
    delete_cnt = del;
    // Do Open mmap() patch file as write
    int insFile = open(pi.c_str(), O_RDWR| O_CREAT,0777);
    lseek (insFile, (ins+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (insFile, "", 1);
    int *insBuffer = (int *)mmap(NULL, (ins+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, insFile, 0);
    
    int delFile = open(pd.c_str(), O_RDWR| O_CREAT,0777);
    lseek (delFile, (del+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (delFile, "", 1);
    int *delBuffer = (int *)mmap(NULL, (del+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, delFile, 0);

    std::vector<int> Icode_idx, Dcode_idx, Icode, Dcode;
    Icode_idx.resize(extra_v_cnt);
    Dcode_idx.resize(extra_v_cnt);
    // TODO: openMP
    double t1 = get_wall_time2();
#pragma omp parallel num_threads(available_threads)
{
#pragma omp for
    for (int i=0; i<extra_v_cnt; i++) {
        patch_data *p = &patch_v[i];
        std::sort(p->insertion.begin(),p->insertion.end());
        std::sort(p->deletion.begin(),p->deletion.end());
        // // Do AIO(Normal I/O first) Dump to patch file
        memcpy(insBuffer+i,ins_idx.data()+i,sizeof(int));
        memcpy(insBuffer+ins_idx[i]+extra_v_cnt,p->insertion.data(),p->insertion.size()*sizeof(int));
        memcpy(delBuffer+i,del_idx.data()+i,sizeof(int));
        memcpy(delBuffer+del_idx[i]+extra_v_cnt,p->deletion.data(),p->deletion.size()*sizeof(int));
    }
}
    double t2 = get_wall_time2();
    printf("\n code time: %.6lf\n", t2 - t1);
    // Generate Patch Index Code
    for (int i=0; i<extra_v_cnt; i++) {
        patch_data *p = &patch_v[i];
        Icode_idx[i]=(Icode.size());
        Dcode_idx[i]=(Dcode.size());
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        int Icur=0, Dcur=0;
        int Dlen=0,D_sum=0;
        bool D_next=false;
        for (unsigned int j=0;j<r-l;j++) {
            unsigned int start_idx = j+l;
            if (Dcur < p->deletion.size()) {
                if (p->deletion[Dcur]==mmp_edge[start_idx]) {
                    if (!D_next) Dcode.push_back(j);
                    Dcur++;
                    Dlen++;
                    D_next = true;
                    if (Dcur == p->deletion.size()) {
                        Dcode.push_back(Dlen);
                        D_sum+=Dlen;
                        Dlen = 0;
                        D_next = false;
                    }
                    continue;
                }
                else {
                    if (D_next == true) {
                        Dcode.push_back(Dlen);
                        D_sum+=Dlen;
                    }
                    Dlen = 0;
                    D_next = false;
                }
            }
            if (Icur < p->insertion.size()) {
                if (p->insertion[Icur]<mmp_edge[start_idx]) {
                    Icode.push_back(start_idx);
                    while (Icur < p->insertion.size() && p->insertion[Icur] < mmp_edge[start_idx])
                    {
                        Icur++;
                    }
                    Icode.push_back(Icur);
                }
            }
        }
        // assert(Dcur==p->deletion.size());
        // assert(D_sum==p->deletion.size());
        if (Icur < p->insertion.size()) Icode.push_back(p->insertion.size());
    }
    assert(Dcode_idx.size()==extra_v_cnt);
    assert(Icode_idx.size()==extra_v_cnt);

    
    // close patch file
    munmap(insBuffer, (ins+extra_v_cnt)*sizeof(int));
    close(insFile);

    // // Do Open mmap() patch Index file as write
    // // Sequentially dump Patch Index Code file

    // update code_len used for reading
    insert_code_len = Icode.size();
    delete_code_len = Dcode.size();

    int insCodeFile = open(picode.c_str(), O_RDWR| O_CREAT,0777);
    int icode_size = Icode_idx.size()+Icode.size();
    lseek (insCodeFile, icode_size*sizeof(int)-1, SEEK_SET);
    write (insCodeFile, "", 1);
    int *insCodeBuffer = (int *)mmap(NULL, icode_size*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, insCodeFile, 0);
    memcpy(insCodeBuffer,Icode_idx.data(),Icode_idx.size()*sizeof(int));
    memcpy(insCodeBuffer+Icode_idx.size(),Icode.data(),Icode.size()*sizeof(int));
    munmap(insCodeBuffer, icode_size*sizeof(int));
    close(insCodeFile);

    int delCodeFile = open(pdcode.c_str(), O_RDWR| O_CREAT,0777);
    int dcode_size = Dcode_idx.size()+Dcode.size();
    lseek (delCodeFile, dcode_size*sizeof(int)-1, SEEK_SET);
    write (delCodeFile, "", 1);
    int *delCodeBuffer = (int *)mmap(NULL, dcode_size*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, delCodeFile, 0);
    memcpy(delCodeBuffer,Dcode_idx.data(),Dcode_idx.size()*sizeof(int));
    memcpy(delCodeBuffer+Dcode_idx.size(),Dcode.data(),Dcode.size()*sizeof(int));
    munmap(delCodeBuffer, dcode_size*sizeof(int));
    close(delCodeFile);

    // Update Done
    init_patch_ptr();
    return;
}

void Graph::init_patch_ptr() {
    patch_read_open = true;
    std::string pi(raw_data_path);
    pi.append("_patch.insert");
    std::string picode(raw_data_path);
    picode.append("_patch_code.insert");
    std::string pdcode(raw_data_path);
    pdcode.append("_patch_code.delete");


    pi_fd = open(pi.c_str(), O_RDWR);
    int *insert_info = static_cast<int*>(mmap(NULL, (insert_cnt+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,pi_fd,0));
    mmp_patch_vertex = (unsigned int*)insert_info;
    mmp_patch_edge = insert_info+extra_v_cnt;
    
    pic_fd = open(picode.c_str(), O_RDWR);
    int *idx_insert_info = static_cast<int*>(mmap(NULL, (insert_code_len+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,pic_fd,0));
    patch_insert_idx = idx_insert_info;
    patch_insert_code = idx_insert_info+extra_v_cnt;
    
    pdc_fd = open(pdcode.c_str(), O_RDWR);
    int *idx_delete_info = static_cast<int*>(mmap(NULL, (delete_code_len+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,pdc_fd,0));
    patch_delete_idx = idx_delete_info;
    patch_delete_code = idx_delete_info+extra_v_cnt;
}

void Graph::release_patch_ptr() {
    if (!patch_read_open) return;
    munmap(mmp_vertex, (insert_cnt+extra_v_cnt)*sizeof(int));
    close(pi_fd);

    munmap(patch_insert_idx, (insert_code_len+extra_v_cnt)*sizeof(int));
    close(pic_fd);

    munmap(patch_delete_idx, (delete_code_len+extra_v_cnt)*sizeof(int));
    close(pdc_fd);
}

void Graph::refine_graph() {
    // set backup graph path

    // Open mmap() backup graph file as write
    unsigned int total_cnt = g_ecnt + insert_cnt - delete_cnt;
    // TODO: path problem
    std::string pi(raw_data_path);
    pi.append("_patch.insert");
    std::string pd(raw_data_path);
    pd.append("_patch.delete");
    std::string p_back(backup_data_path);

    std::string graph_size(p_back);
    graph_size.append(".size");
    DataLoader::gen_data_size(extra_v_cnt,total_cnt,graph_size);

    int gfd = open(p_back.c_str(), O_RDWR| O_CREAT,0777);
    int gsize = extra_v_cnt+total_cnt; 
    lseek(gfd, gsize*sizeof(int)-1, SEEK_SET);
    write (gfd, "", 1);
    int *g_ptr = (int *)mmap(NULL, gsize*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, gfd, 0);
    
    // Open mmap() patch File as read
    int idx_ifd = open(pi.c_str(), O_RDWR);
    int *idx_iinfo = static_cast<int*>(mmap(NULL, (insert_cnt+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,idx_ifd,0));
    int *insert_content = idx_iinfo+extra_v_cnt;

    int idx_dfd = open(pd.c_str(), O_RDWR);
    int *idx_dinfo = static_cast<int*>(mmap(NULL, (delete_cnt+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,idx_dfd,0));
    int *delete_content = idx_dinfo+extra_v_cnt;

    // info preparation of insert offset
    int *adj_idx = new int[extra_v_cnt];
    int poi = 0;
    for (int i=0;i<extra_v_cnt; i++) {
        adj_idx[i] = poi;
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        // TODO: calculate degree of each vertex;
        int degree=r-l;
        if (i<extra_v_cnt-1) degree += idx_iinfo[i+1]-idx_iinfo[i];
        else degree += insert_cnt-idx_iinfo[i];
        if (i<extra_v_cnt-1) degree -= idx_dinfo[i+1]-idx_dinfo[i];
        else degree -= delete_cnt-idx_dinfo[i];

        g_ptr[i] = poi;
        poi += degree;
    }
    
    // TODO: This can be done in parallel.
    int *fill_edge_start = g_ptr+extra_v_cnt;
    for (int i=0; i<extra_v_cnt; i++) {
        if ((i<g_vcnt && test(bitmap,i) || i>=g_vcnt)) {
            // TODO: asyn I/O
            int fill_cur = adj_idx[i];
            int fill_end;
            if (i<extra_v_cnt-1) {
                fill_end = adj_idx[i+1];
            }
            else {
                fill_end = total_cnt;
            }
            unsigned int l,r;
            get_mmp_edge_index(i,l,r);
            int ld=idx_dinfo[i],li=idx_iinfo[i];
            int rd,ri;
            if (i<extra_v_cnt-1) {
                rd = idx_dinfo[i+1];
                ri = idx_iinfo[i+1];
            }
            else {
                rd = delete_cnt;
                ri = insert_cnt;
            }
            while (fill_cur < fill_end)
            {
                if (ld < rd && delete_content[ld]==mmp_edge[l]) {
                    ld++;
                    l++;
                }
                else if (insert_content[li] < mmp_edge[l]) {
                    fill_edge_start[fill_cur] = insert_content[li];
                    li++;
                    fill_cur++;
                }
                else {
                    fill_edge_start[fill_cur] = mmp_edge[l];
                    l++;
                    fill_cur++;
                }
            }
            // assertion
            
        }
        else {
            // TODO: asyn I/O
            unsigned int l,r;
            get_mmp_edge_index(i,l,r);
            // printf("a: %d\n",adj_idx[i]);
            memcpy(fill_edge_start+adj_idx[i],mmp_edge+l,(r-l)*sizeof(int));
        }
    }
    delete []adj_idx;
    raw_data_path.swap(backup_data_path);
}