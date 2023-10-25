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

int Graph::memory_map() {
    
    std::string graph_size(raw_data_path);
    graph_size.append(".size");
    DataLoader::load_data_size(g_vcnt,g_ecnt,graph_size);
    
    int fd = open(raw_data_path.c_str(),O_RDWR);
    mem = static_cast<int*>(mmap(NULL, (g_vcnt+g_ecnt)*sizeof(int), PROT_READ, MAP_SHARED ,fd,0));
    mmp_vertex = (unsigned int *)mem;
    mmp_edge = mem+g_vcnt;
    extra_v_cnt = g_vcnt;
    if (bitmap!=NULL) delete[] bitmap;
    bitmap = new bool[g_vcnt];
    memset(bitmap,0,g_vcnt);
    g_fd = fd;
    return fd;
}

int Graph::memory_lmap() {
    
    std::string graph_size(raw_data_path);
    graph_size.append(".size");
    printf("size path: %s\n",graph_size.c_str());
    DataLoader::load_data_size(g_vcnt,lg_ecnt,graph_size);
    printf("large graph size %d %lld\n",g_vcnt,lg_ecnt);
    
    int fd = open(raw_data_path.c_str(),O_RDWR);
    mem = (int *)(mmap(NULL, g_vcnt*sizeof(size_t)+lg_ecnt*sizeof(int), PROT_READ, MAP_SHARED ,fd,0));
    mmp_l_vertex = (size_t *)mem;
    mmp_edge = (int *)(mmp_l_vertex+g_vcnt);
    extra_v_cnt = g_vcnt;
    if (bitmap!=NULL) delete[] bitmap;
    // bitmap = new bool[g_vcnt];
    // memset(bitmap,0,g_vcnt);
    g_fd = fd;
    return fd;
}

void Graph::free_map(int &fd) {
    release_patch_ptr();
    munmap(mem, (g_vcnt+g_ecnt)*sizeof(int));
    close(fd);
    if (bitmap!=NULL)
        delete[] bitmap;
    bitmap=NULL;
    return;
}

void Graph::free_lmap(int &fd) {
    munmap(mem, g_vcnt*sizeof(size_t)+lg_ecnt*sizeof(int));
    close(fd);
    if (bitmap!=NULL)
        delete[] bitmap;
    bitmap=NULL;
    return;
}

// control by client, invoke when switch to new snapshot
void Graph::init_patch_ptr() {
    patch_read_open = true;
    std::string pi(patch_insert_path);
    std::string pd(patch_delete_path);
    
    pi_fd = open(pi.c_str(), O_RDWR);
    int *insert_info = static_cast<int*>(mmap(NULL, (insert_cnt+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,pi_fd,0));
    mmp_patch_vertex = (unsigned int*)insert_info;
    mmp_patch_edge = insert_info+extra_v_cnt;
    
    pd_fd = open(pd.c_str(), O_RDWR);
    int *delete_info = static_cast<int*>(mmap(NULL, (delete_cnt+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED ,pd_fd,0));
    mmp_patch_delete_vertex = (unsigned int*)delete_info;
    mmp_patch_delete_edge = delete_info+extra_v_cnt;
    return;
}

// control by client, invoke when switch to new snapshot
void Graph::release_patch_ptr() {
    if (!patch_read_open) return;
    munmap(mmp_patch_vertex, (insert_cnt+extra_v_cnt)*sizeof(int));
    close(pi_fd);
    mmp_patch_vertex=NULL;
    mmp_patch_edge=NULL;
    
    munmap(mmp_patch_delete_vertex, (delete_cnt+extra_v_cnt)*sizeof(int));
    mmp_patch_delete_edge=NULL;
    mmp_patch_delete_vertex=NULL;
    close(pd_fd);
    patch_read_open = false;
    return;
}


void Graph::insert_edge(int v1, int v2) {
    if (patch_vi.size() <= v2) {
        patch_vi.resize(v2+1);
        extra_v_cnt = v2+1;
    }
    if (patch_vi.size() <= v1) {
        patch_vi.resize(v1+1);
        extra_v_cnt = v1+1;
    }
    patch_vi[v1].push_back(v2);
    patch_vi[v2].push_back(v1);
}

void Graph::delete_edge(int v1, int v2) {
    patch_vd[v1].push_back(v2);
    patch_vd[v2].push_back(v1);
}

// gather to edge list
void Graph::update(std::vector<int> &u, std::vector<int> &d) {
    // This function is aimed to executed in single thread, overlapped by matching
    patch_recycle();
    double t1 = get_wall_time2();
    for (int i=0; i<u.size();i+=2) {
        insert_edge(u[i],u[i+1]);
    }
    for (int i=0; i<d.size();i+=2) {
        delete_edge(d[i],d[i+1]);
    }
    double t2 = get_wall_time2();
    
    // printf("\nupdating time: %.6lf\n", t2 - t1);
}

// invoke by background generator
void Graph::dump_v(int offset) {
    int *ins_idx = new int[extra_v_cnt];
    int *del_idx = new int[extra_v_cnt];
    int ins=0,del=0;
    
    for (int i=0; i<extra_v_cnt; i++) {
        ins_idx[i]=ins;
        ins+=patch_vi[i].size();
        del_idx[i]=del;
        del+=patch_vd[i].size();
    }
    insert_cnt = ins;
    delete_cnt = del;

    double t1 = get_wall_time2();
    std::string pi(patch_insert_path);
    // pi.append(std::to_string(offset));
    std::string pd(patch_delete_path);
    // pd.append(std::to_string(offset));
    std::string psize(raw_data_path);
    psize+="_patch.size";
    // printf("dump size %s\n",psize.c_str());
    // printf("pi %s\n",pi.c_str());
    // printf("pd %s\n",pd.c_str());
    DataLoader::gen_data_size(ins,del,psize);

    
    int insFile = open(pi.c_str(), O_RDWR| O_CREAT,0777);
    lseek (insFile, (insert_cnt+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (insFile, "", 1);
    int *insBuffer = (int *)mmap(NULL, (insert_cnt+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, insFile, 0);
    
    int delFile = open(pd.c_str(), O_RDWR| O_CREAT,0777);
    lseek (delFile, (delete_cnt+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (delFile, "", 1);
    int *delBuffer = (int *)mmap(NULL, (delete_cnt+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, delFile, 0);


#pragma omp parallel num_threads(available_threads)
{
#pragma omp for nowait
    for (int i=0; i<extra_v_cnt; i++) {
        std::sort(patch_vi[i].begin(),patch_vi[i].end());
        std::sort(patch_vd[i].begin(),patch_vd[i].end());
        // // Do AIO(Normal I/O first) Dump to patch file
        memcpy(insBuffer+i,ins_idx+i,sizeof(int));
        memcpy(insBuffer+ins_idx[i]+extra_v_cnt,patch_vi[i].data(),patch_vi[i].size()*sizeof(int));
        memcpy(delBuffer+i,del_idx+i,sizeof(int));
        memcpy(delBuffer+del_idx[i]+extra_v_cnt,patch_vd[i].data(),patch_vd[i].size()*sizeof(int));
    }
}
    double t2 = get_wall_time2();
    delete[] ins_idx;
    delete[] del_idx;
    // Update Done
    
    munmap(insBuffer, (insert_cnt+extra_v_cnt)*sizeof(int));
    close(insFile);
    
    munmap(delBuffer, (delete_cnt+extra_v_cnt)*sizeof(int));
    close(delFile);
    // init_patch_ptr();
    return;
}

void Graph::dump_coarse(int offset) {
    int *ins_idx = new int[extra_v_cnt];
    int *del_idx = new int[extra_v_cnt];
    int ins=0,del=0;
    
    for (int i=0; i<extra_v_cnt; i++) {
        ins_idx[i]=ins;
        ins+=patch_vi[i].size();
        del_idx[i]=del;
        del+=patch_vd[i].size();
    }
    insert_cnt = ins;
    delete_cnt = del;

    std::string pi(patch_insert_path);
    pi.append(std::to_string(offset));
    std::string pd(patch_delete_path);
    pd.append(std::to_string(offset));

    
    int insFile = open(pi.c_str(), O_RDWR| O_CREAT,0777);
    lseek (insFile, (insert_cnt+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (insFile, "", 1);
    int *insBuffer = (int *)mmap(NULL, (insert_cnt+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, insFile, 0);
    
    int delFile = open(pd.c_str(), O_RDWR| O_CREAT,0777);
    lseek (delFile, (delete_cnt+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (delFile, "", 1);
    int *delBuffer = (int *)mmap(NULL, (delete_cnt+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, delFile, 0);


#pragma omp parallel num_threads(available_threads)
{
#pragma omp for nowait
    for (int i=0; i<extra_v_cnt; i++) {
        std::sort(patch_vi[i].begin(),patch_vi[i].end());
        std::sort(patch_vd[i].begin(),patch_vd[i].end());
        // // Do AIO(Normal I/O first) Dump to patch file
        memcpy(insBuffer+i,ins_idx+i,sizeof(int));
        memcpy(insBuffer+ins_idx[i]+extra_v_cnt,patch_vi[i].data(),patch_vi[i].size()*sizeof(int));
        memcpy(delBuffer+i,del_idx+i,sizeof(int));
        memcpy(delBuffer+del_idx[i]+extra_v_cnt,patch_vd[i].data(),patch_vd[i].size()*sizeof(int));
    }
}
    delete[] ins_idx;
    delete[] del_idx;
    // Update Done
    
    munmap(insBuffer, (insert_cnt+extra_v_cnt)*sizeof(int));
    close(insFile);
    
    munmap(delBuffer, (delete_cnt+extra_v_cnt)*sizeof(int));
    close(delFile);
    init_patch_ptr();
    patch_recycle();
    return;
}

void Graph::patch_recycle() {
    patch_vi.clear();
    patch_vd.clear();
    patch_vi.resize(extra_v_cnt);
    patch_vd.resize(extra_v_cnt);
}

int Graph::init_backup_ptr(int v, unsigned int e) {
    std::string p_back(backup_data_path);

    std::string graph_size(p_back);
    graph_size.append(".size");
    DataLoader::gen_data_size(v,e,graph_size);

    int gfd = open(p_back.c_str(), O_RDWR| O_CREAT,0777);
    long gsize = v+e; 
    lseek(gfd, gsize*sizeof(int)-1, SEEK_SET);
    write (gfd, "", 1);
    g_back_up_V = (unsigned int *)mmap(NULL, gsize*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, gfd, 0);
    g_back_up_E = (int *) g_back_up_V + v;
    return gfd;
}

void Graph::close_backup_ptr(int fd, int v, unsigned int e) {
    munmap(g_back_up_V,(size_t)(v+e)*sizeof(int));
    close(fd);
}

void Graph::update_patch(const std::string& piu, const std::string& pdu, int update_ins, int update_del) {
    // case that not using patch
    if (!patch_read_open) {
        std::string pi(patch_insert_path);
        std::string pd(patch_delete_path);
        // printf("%s\n",piu.c_str());
        // printf("%s\n",pi.c_str());
        // printf("%d\n",piu==pi);
        if (piu==pi && pdu==pd) {
            
            int piuFile = open(piu.c_str(), O_RDWR);
            int *piuBuffer = (int *)mmap(NULL, (update_ins+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED, piuFile, 0);
            
            int pduFile = open(pdu.c_str(), O_RDWR);
            int *pduBuffer = (int *)mmap(NULL, (update_del+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED, pduFile, 0);

            insert_cnt = update_ins;
            delete_cnt = update_del;
            backup_ecnt = g_ecnt + insert_cnt - delete_cnt;
            back_fd=init_backup_ptr(extra_v_cnt, backup_ecnt);
            int ioff=0, doff=0;
            for (int i=0;i<extra_v_cnt;i++) {
                
                unsigned int l,r;
                get_mmp_edge_index(i,l,r);
                g_back_up_V[i] = l+ioff-doff;
                if (i==extra_v_cnt-1) {
                    ioff += update_ins-piuBuffer[i];
                    doff += update_del-pduBuffer[i];
                    continue;
                }
                ioff += piuBuffer[i+1]-piuBuffer[i];
                doff += pduBuffer[i+1]-pduBuffer[i];
            }
            // printf("inner update %d %d\n",piuBuffer[extra_v_cnt-1],doff);
            // assert(ioff==update_ins);
            // assert(doff==update_del);
            munmap(piuBuffer, (update_ins+extra_v_cnt)*sizeof(int));
            close(piuFile);
            munmap(pduBuffer, (update_del+extra_v_cnt)*sizeof(int));
            close(pduFile);
            init_patch_ptr();
            return;
        }

        int piuFile = open(piu.c_str(), O_RDWR);
        int *piuBuffer = (int *)mmap(NULL, (update_ins+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED, piuFile, 0);
        
        int pduFile = open(pdu.c_str(), O_RDWR);
        int *pduBuffer = (int *)mmap(NULL, (update_del+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED, pduFile, 0);

        int insFile = open(pi.c_str(), O_RDWR| O_CREAT,0777);
        lseek (insFile, (update_ins+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
        write (insFile, "", 1);
        int *insBuffer = (int *)mmap(NULL, (update_ins+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, insFile, 0);
        memcpy(insBuffer,piuBuffer,(update_ins+extra_v_cnt)*sizeof(int));
        
        int delFile = open(pd.c_str(), O_RDWR| O_CREAT,0777);
        lseek (delFile, (update_del+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
        write (delFile, "", 1);
        int *delBuffer = (int *)mmap(NULL, (update_del+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, delFile, 0);
        memcpy(delBuffer,pduBuffer,(update_del+extra_v_cnt)*sizeof(int));

        insert_cnt = update_ins;
        delete_cnt = update_del;

        munmap(piuBuffer, (update_ins+extra_v_cnt)*sizeof(int));
        close(piuFile);
        munmap(pduBuffer, (update_del+extra_v_cnt)*sizeof(int));
        close(pduFile);
        
        munmap(insBuffer, (update_ins+extra_v_cnt)*sizeof(int));
        close(insFile);
        
        munmap(delBuffer, (update_del+extra_v_cnt)*sizeof(int));
        close(delFile);
        init_patch_ptr();
        backup_ecnt = g_ecnt + insert_cnt - delete_cnt;
        back_fd=init_backup_ptr(extra_v_cnt, backup_ecnt);
        int ioff=0, doff=0;
        for (int i=0;i<extra_v_cnt;i++) {
            
            unsigned int l,r;
            get_mmp_edge_index(i,l,r);
            g_back_up_V[i] = l+ioff-doff;
            if (i==extra_v_cnt-1) {
                continue;
            }
            ioff += insBuffer[i+1]-insBuffer[i];
            doff += delBuffer[i+1]-delBuffer[i];
        }
        return;
    }
    release_patch_ptr();
    double t1 = get_wall_time2();
    std::string pi(patch_insert_path_backup);
    std::string pd(patch_delete_path_backup);
    
    // Do Open mmap() patch file as write
    int piuFile = open(piu.c_str(), O_RDWR);
    int *piuBuffer = (int *)mmap(NULL, (update_ins+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED, piuFile, 0);
    int *piuEdge = piuBuffer+extra_v_cnt;
    
    int pduFile = open(pdu.c_str(), O_RDWR);
    int *pduBuffer = (int *)mmap(NULL, (update_del+extra_v_cnt)*sizeof(int), PROT_READ, MAP_SHARED, pduFile, 0);
    int *pduEdge = pduBuffer+extra_v_cnt;

    int insFile = open(pi.c_str(), O_RDWR| O_CREAT,0777);
    lseek (insFile, (insert_cnt+update_ins+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (insFile, "", 1);
    int *insBuffer = (int *)mmap(NULL, (insert_cnt+update_ins+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, insFile, 0);
    
    int delFile = open(pd.c_str(), O_RDWR| O_CREAT,0777);
    lseek (delFile, (delete_cnt+update_del+extra_v_cnt)*sizeof(int)-1, SEEK_SET);
    write (delFile, "", 1);
    int *delBuffer = (int *)mmap(NULL, (delete_cnt+update_del+extra_v_cnt)*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, delFile, 0);

    std::vector<int> insert_offset, delete_offset;
    int ioff=0, doff=0;
    backup_ecnt = g_ecnt + insert_cnt - delete_cnt;
    back_fd=init_backup_ptr(extra_v_cnt, backup_ecnt);
    for (int i=0; i< extra_v_cnt; i++) {
        
        unsigned int il,ir,dl,dr,uil,uir,udl,udr;
        get_mmp_patch_edge_index(i,il,ir);
        get_mmp_patch_delete_edge_index(i,dl,dr);
        if (i == extra_v_cnt-1) {
            uir = update_ins;
            udr = update_del;
        }
        else {
            uir = (unsigned int)piuBuffer[i+1];
            udr = (unsigned int)pduBuffer[i+1];
        }
        uil = (unsigned int)piuBuffer[i];
        udl = (unsigned int)pduBuffer[i];

        insert_offset.push_back(ioff);
        delete_offset.push_back(doff);
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        g_back_up_V[i] = l+ioff-doff;
        doff += udr-udl+dr-dl;
        ioff += uir-uil+ir-il;
    }

#pragma omp parallel num_threads(available_threads)
{
#pragma omp for nowait
    for (int i=0; i<extra_v_cnt; i++) {
        unsigned int il,ir,dl,dr,uil,uir,udl,udr;
        get_mmp_patch_edge_index(i,il,ir);
        get_mmp_patch_delete_edge_index(i,dl,dr);
        if (i == extra_v_cnt-1) {
            uir = update_ins;
            udr = update_del;
        }
        else {
            uir = (unsigned int)piuBuffer[i+1];
            udr = (unsigned int)pduBuffer[i+1];
        }
        uil = (unsigned int)piuBuffer[i];
        udl = (unsigned int)pduBuffer[i];
        int ic=il,dc=dl,iuc=uil,duc=udl;
        int *iptr = new int[ir-il+uir-uil];
        int *dptr = new int[dr-dl+udr-udl];
        int icount=0, dcount=0;
        while (ic < ir || iuc < uir)
        {
            if (iuc < uir && ic < ir){
                if (mmp_patch_edge[ic]>piuEdge[iuc]) {
                    iptr[icount] = piuEdge[iuc];
                    iuc++;
                }
                else {
                    iptr[icount] = mmp_patch_edge[ic];
                    ic++;
                }
            }
            else if (ic < ir){
                iptr[icount] = piuEdge[iuc];
                iuc++;
            }
            else {
                iptr[icount] = mmp_patch_edge[ic];
                ic++;
            }
            icount++;
        }
        
        while (dc < dr || duc < udr)
        {
            if (duc < udr && dc < dr){
                if (mmp_patch_delete_edge[dc]>pduEdge[duc]) {
                    dptr[dcount] = pduEdge[duc];
                    duc++;
                }
                else {
                    dptr[dcount] = mmp_patch_edge[dc];
                    dc++;
                }
            }
            else if (dc < dr){
                dptr[dcount] = pduEdge[duc];
                duc++;
            }
            else {
                dptr[dcount] = mmp_patch_edge[dc];
                dc++;
            }
            dcount++;
        }
        assert(icount == ir-il+uir-uil);
        assert(dcount == dr-dl+udr-udl);
        

        memcpy(insBuffer+i,insert_offset.data()+i,sizeof(int));
        memcpy(insBuffer+insert_offset[i]+extra_v_cnt, iptr, icount*sizeof(int));
        memcpy(delBuffer+i,delete_offset.data()+i,sizeof(int));
        memcpy(delBuffer+insert_offset[i]+extra_v_cnt, dptr, dcount*sizeof(int));

        delete[] iptr;
        delete[] dptr;
    }
}
    double t2 = get_wall_time2();

    munmap(piuBuffer, (update_ins+extra_v_cnt)*sizeof(int));
    close(piuFile);
    munmap(pduBuffer, (update_del+extra_v_cnt)*sizeof(int));
    close(pduFile);
    
    munmap(insBuffer, (insert_cnt+update_ins+extra_v_cnt)*sizeof(int));
    close(insFile);
    
    munmap(delBuffer, (delete_cnt+update_del+extra_v_cnt)*sizeof(int));
    close(delFile);
    patch_delete_path.swap(patch_delete_path_backup);
    patch_insert_path.swap(patch_insert_path_backup);
    insert_cnt += update_ins;
    delete_cnt += update_del;
    init_patch_ptr();
    return;
}

void Graph::needle(int *&p, int v, int& size) {
    unsigned int insert_l,insert_r,delete_l,delete_r;
    get_mmp_patch_edge_index(v,insert_l,insert_r);
    get_mmp_patch_delete_edge_index(v,delete_l,delete_r);
    unsigned int l,r;
    get_mmp_edge_index(v,l,r);
    size = r-l+insert_r-insert_l-delete_r+delete_l;
    p = g_back_up_E+g_back_up_V[v];
    int cur = 0;
    // printf("\nneedle %d\n",v);
    // printf("\ninsert: %d\n",insert_r-insert_l);
    // for (int a=insert_l;a<insert_r;a++) {
    //     printf(" %d",mmp_patch_edge[a]);
    // }
    // printf("\ndelete: %d\n",delete_r-delete_l);
    // for (int a=delete_l;a<delete_r;a++) {
    //     printf(" %d",mmp_patch_delete_edge[a]);
    // }
    // printf("\nraw: %d\n",r-l);
    // for (int a=l;a<r;a++) {
    //     printf(" %d",mmp_edge[a]);
    // }
    while (cur < size)
    {
        if (delete_l < delete_r ) {
            if (mmp_patch_delete_edge[delete_l]==mmp_edge[l]) {
                delete_l++;
                l++;
                continue;
            }
        }
        if (insert_l < insert_r) {
            if (l==r || mmp_patch_edge[insert_l] < mmp_edge[l])
            {
                p[cur] = mmp_patch_edge[insert_l];
                insert_l++;
                cur++;
                continue;
            }
        }
        p[cur] = mmp_edge[l];
        l++;
        cur++;
    }
    // assert(l==r);
    // assert(delete_l==delete_r);
    // assert(insert_l==insert_r);
    // assert(cur==degree);
    // printf("\nafter:");
    // for (int i=0;i<size;i++) {
    //     assert(p[i]<extra_v_cnt);
    //     printf(" %d",p[i]);
    // }
    // printf("\n");
}

void Graph::bandaid(int *p, int v, int& size) {
    unsigned int insert_l,insert_r,delete_l,delete_r;
    get_mmp_patch_delete_edge_index(v,delete_l,delete_r);
    unsigned int l,r;
    get_mmp_edge_index(v,l,r);
    int degree = r-l+insert_r-insert_l-delete_r+delete_l;
    p = new int[degree];
    size = degree;
    int cur = 0;
    while (cur < degree || delete_l < delete_r)
    {
        if (delete_l < delete_r ) {
            if (mmp_patch_delete_edge[delete_l]==mmp_edge[l]) {
                delete_l++;
                l++;
                continue;
            }
        }
        if (insert_l < insert_r) {
            if (l==r || mmp_patch_edge[insert_l] < mmp_edge[l])
            {
                p[cur] = mmp_patch_edge[insert_l];
                insert_l++;
                cur++;
                continue;
            }
        }
        p[cur] = mmp_edge[l];
        l++;
        cur++;
    }
    assert(l==r);
    assert(delete_l==delete_r);
    assert(insert_l==insert_r);
    assert(cur==degree);
}

// fsync dump on computing
void Graph::refine_graph(bool preprocess) {
    // patch_vi.clear();
    // patch_vd.clear();
    // set backup graph path

    // Open mmap() backup graph file as write
    unsigned int total_cnt = g_ecnt + insert_cnt - delete_cnt;
    // TODO: path problem
    std::string pi(patch_insert_path);
    std::string pd(patch_delete_path);
    std::string p_back(backup_data_path);
    if (preprocess) p_back=raw_data_path+"_preprocess";

    // printf("used path %s\n",pi.c_str());
    // printf("used path %s\n",pd.c_str());

    std::string graph_size(p_back);
    graph_size.append(".size");
    DataLoader::gen_data_size(extra_v_cnt,total_cnt,graph_size);

    int gfd = open(p_back.c_str(), O_RDWR| O_CREAT,0777);
    long gsize = extra_v_cnt+total_cnt; 
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
    unsigned int *adj_idx = new unsigned int[extra_v_cnt];
    int poi = 0;
    for (int i=0;i<extra_v_cnt; i++) {
        adj_idx[i] = poi;
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        assert(l<=r);
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
#pragma omp parallel num_threads(1)
{
#pragma omp for
    for (int i=0; i<extra_v_cnt; i++) {
        // TODO: asyn I/O
        unsigned int fill_cur = adj_idx[i];
        int fill_end;
        if (i<extra_v_cnt-1) {
            fill_end = adj_idx[i+1];
        }
        else {
            fill_end = total_cnt;
        }
        unsigned int l,r;
        get_mmp_edge_index(i,l,r);
        assert(l<=r);
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
        // assert(rd-ld==patch_vd[i].size());
        // assert(ri-li==patch_vi[i].size());
        while (fill_cur < fill_end || ld < rd)
        {
            if (ld < rd) {
                assert(l<r);
                if (delete_content[ld]<mmp_edge[l]) {
                    printf("%d %d %d %d\n",delete_content[ld],mmp_edge[l],i,l);
                }
                assert(delete_content[ld]>=mmp_edge[l]);
                
                if (delete_content[ld]==mmp_edge[l])
                {
                    ld++;
                    l++;
                    continue;
                }
            }
            if (li < ri)
            {
                if (l==r || insert_content[li] < mmp_edge[l]) {
                    fill_edge_start[fill_cur] = insert_content[li];
                    li++;
                    fill_cur++;
                    continue;
                }
            }
            fill_edge_start[fill_cur] = mmp_edge[l];
            l++;
            fill_cur++;
        }
        assert(ld==rd);
        // printf("%d %d %d\n",i,l,r);
        // printf("%d %d %d\n",i,fill_cur,fill_end);
        assert(li==ri);
        assert(fill_cur == fill_end);
        assert(l==r);
    }
}
    delete []adj_idx;
    munmap(idx_iinfo,(insert_cnt+extra_v_cnt)*sizeof(int));
    close(idx_ifd);
    munmap(idx_dinfo,(delete_cnt+extra_v_cnt)*sizeof(int));
    close(idx_dfd);
    munmap(g_ptr,gsize*sizeof(int));
    close(gfd);
    if (!preprocess) {
        raw_data_path.swap(backup_data_path);
        patch_delete_path.swap(patch_delete_path_backup);
        patch_insert_path.swap(patch_insert_path_backup);
    }
    patch_vi.clear();
    patch_vd.clear();
}

