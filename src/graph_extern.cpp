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

void Graph::load_physical_edge_index(int &v, std::shared_ptr<int[]> &ptr, int &size) {
    io_num+=1;
    int file_id = v_state_map[v].file_id;

    std::string block_path(file_path);
    
    if (blockType == BlockType::K_CORE_BLOCK) block_path+=file_str[file_id];
    else if (blockType == BlockType::RANDOM_BLOCK) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_br");
    else if (blockType == BlockType::CHAIN) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bc");
    else if (blockType == BlockType::SIMPLE_BFS) block_path = raw_data_path + std::string("_blocks/") + std::to_string(file_id) + std::string("_bs");

    block_length lens = file_len[file_id];
    int load_v_len=lens.v_len;
    unsigned int load_e_len=lens.e_len;
    int* data = new int[2*load_v_len+load_e_len];

    double t1 = get_wall_time();
    DataLoader::load_block_data_aggregate2(data,2*load_v_len+load_e_len,block_path);
    double t2 = get_wall_time();
    int *vid = data;
    unsigned int *vertex_offset = (unsigned int *)data+load_v_len;
    int *adj_list = data+2*load_v_len;
    int vid_cursor = 0;
    while (vid[vid_cursor] != v) {
//        assert(vid_cursor < load_v_len);
//        assert(vid[vid_cursor] < v);
        vid_cursor++;
    }
//    assert(vid[vid_cursor]==v);
    int l,r;
    l=vertex_offset[vid_cursor];
    if (vid_cursor==load_v_len-1) r=load_e_len;
    else r= vertex_offset[vid_cursor+1];
//    assert(r-l >0);

    std::shared_ptr<int[]> p(new int[r-l]());

    size = r-l;
    for (int i = l; i < r; i++) {
        // p[i - l] = adj_list[i];
    }

    ptr.swap(p);
    delete[] data;

    if (omp_get_thread_num()==0)
        blocking_manage_time += t2-t1;
}

void Graph::update(std::vector<pair> &u, std::vector<pair> &d) {
    // TODO: openMP
    for (int i=0; i<u.size();i++) {
        insert_edge(u[i].e1,u[i].e2);
        set(bitmap,u[i].e1);
        set(bitmap,u[i].e2);
    }
    for (int i=0; i<d.size();i++) {
        delete_edge(d[i].e1,d[i].e2);
        set(bitmap,d[i].e1);
        set(bitmap,d[i].e2);
    }
}

void Graph::dump_v() {
    std::string pi(raw_data_path);
    pi.append("_patch.insert");
    std::string pd(raw_data_path);
    pd.append("_patch.delete");
    std::string picode(raw_data_path);
    picode.append("_patch_code.insert");
    std::string pdcode(raw_data_path);
    pdcode.append("_patch_code.delete");

    // sequential
    // TODO: add reading prev patch info
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
    int insFile = open(pi.c_str(), O_RDWR);
    lseek (insFile, ins+extra_v_cnt-1, SEEK_SET);
    write (insFile, "", 1);
    int *insBuffer = (int *)mmap(NULL, ins+extra_v_cnt, PROT_READ | PROT_WRITE, MAP_SHARED, insFile, 0);
    
    int delFile = open(pd.c_str(), O_RDWR);
    lseek (delFile, del+extra_v_cnt-1, SEEK_SET);
    write (delFile, "", 1);
    int *delBuffer = (int *)mmap(NULL, del+extra_v_cnt, PROT_READ | PROT_WRITE, MAP_SHARED, delFile, 0);

    std::vector<int> Icode_idx, Dcode_idx, Icode, Dcode;
    // TODO: openMP
    for (int i=0; i<extra_v_cnt; i++) {
        patch_data *p = &patch_v[i];
        std::sort(p->insertion.begin(),p->insertion.end());
        std::sort(p->deletion.begin(),p->deletion.end());
        // Do AIO(Normal I/O first) Dump to patch file
        insBuffer[i] = ins_idx[i];
        memcpy(insBuffer+ins_idx[i]+extra_v_cnt,p->insertion.data(),p->insertion.size());
        delBuffer[i] = del_idx[i];
        memcpy(delBuffer+del_idx[i]+extra_v_cnt,p->deletion.data(),p->deletion.size());
        // Gen Patch Index Code
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        int Icur=0, Dcur=0;
        bool D_next=false;
        for (unsigned int j=0;j<r-l;j++) {
            unsigned int start_idx = j+l;
            if (Dcur < p->deletion.size()) {
                if (p->deletion[Dcur]==mmp_edge[start_idx]) {
                    if (!D_next) Dcode.push_back(start_idx);
                    Dcur++;
                    D_next = true;
                    continue;
                }
                else {
                    if (D_next == true) {
                        Dcode.push_back(Dcur);
                    }
                    D_next = false;
                }
            }
            if (Icur < p->insertion.size()) {
                if (p->insertion[Icur]<mmp_edge[start_idx]) {
                    Icode.push_back(start_idx);
                    while (p->insertion[Icur] < mmp_edge[start_idx] && Icur < p->insertion.size())
                    {
                        Icur++;
                    }
                    Icode.push_back(Icur);
                }
            }
        }
        if (Icur < p->insertion.size()) Icode.push_back(p->insertion.size());
        Icode_idx.push_back(Icode.size());
        Dcode_idx.push_back(Dcode.size());
    }
    
    // close patch file
    munmap(insBuffer, ins);
    close(insFile);
    munmap(delBuffer, del);
    close(delFile);

    // Do Open mmap() patch Index file as write
    // Sequentially dump Patch Index Code file
    int insCodeFile = open(picode.c_str(), O_RDWR);
    int icode_size = Icode_idx.size()+Icode.size();
    lseek (insCodeFile, icode_size-1, SEEK_SET);
    write (insCodeFile, "", 1);
    int *insCodeBuffer = (int *)mmap(NULL, icode_size, PROT_READ | PROT_WRITE, MAP_SHARED, insCodeFile, 0);
    memcpy(insCodeBuffer,Icode_idx.data(),Icode_idx.size());
    memcpy(insCodeBuffer+Icode_idx.size(),Icode.data(),Icode.size());
    munmap(insCodeBuffer, icode_size);
    close(insCodeFile);

    int delCodeFile = open(pdcode.c_str(), O_RDWR);
    int dcode_size = Dcode_idx.size()+Dcode.size();
    lseek (delCodeFile, dcode_size-1, SEEK_SET);
    write (delCodeFile, "", 1);
    int *delCodeBuffer = (int *)mmap(NULL, dcode_size, PROT_READ | PROT_WRITE, MAP_SHARED, delCodeFile, 0);
    memcpy(delCodeBuffer,Dcode_idx.data(),Dcode_idx.size());
    memcpy(delCodeBuffer+Dcode_idx.size(),Dcode.data(),Icode.size());
    munmap(delCodeBuffer, dcode_size);
    close(delCodeFile);

    // Update Done
}

void Graph::refine_graph() {
    // set backup graph path

    // Open mmap() backup graph file as write
    int total_cnt = g_ecnt + insert_cnt - delete_cnt;
    std::string p_back,pi,pd;
    int gfd = open(p_back.c_str(), O_RDWR);
    int gsize = extra_v_cnt+total_cnt; 
    lseek(gfd, gsize-1, SEEK_SET);
    write (gfd, "", 1);
    int *g_ptr = (int *)mmap(NULL, gsize, PROT_READ | PROT_WRITE, MAP_SHARED, gfd, 0);
    
    // Open mmap() patch File as read
    int idx_ifd = open(pi.c_str(), O_RDWR);
    int i_size;
    int *idx_iinfo = static_cast<int*>(mmap(NULL, (insert_cnt+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,idx_ifd,0));
    int *insert_content = idx_iinfo+extra_v_cnt;

    int idx_dfd = open(pd.c_str(), O_RDWR);
    int d_size;
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
        if (i<extra_v_cnt) degree += idx_iinfo[i+1]-idx_iinfo[i];
        else degree += insert_cnt-idx_iinfo[i];
        if (i<extra_v_cnt) degree -= idx_dinfo[i+1]-idx_dinfo[i];
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
            if (i<extra_v_cnt) {
                fill_end = adj_idx[i+1];
            }
            else {
                fill_end = total_cnt;
            }
            unsigned int l,r;
            get_mmp_edge_index(i,l,r);
            int ld=idx_dinfo[i],li=idx_iinfo[i];
            int rd,ri;
            if (i<extra_v_cnt) {
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
            memcpy(fill_edge_start+adj_idx[i],mmp_edge+l,r-l);
        }
    }
}