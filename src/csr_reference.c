/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 1: CRS construction
// Simple two-pass CRS construction using Active Messages


#include "common.h"
#include "csr_reference.h"
#include "aml.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <search.h>
#include <parmetisbin.h>
#define NCON    1

int64_t nverts_known = 0;
int *degrees;
int64_t *column;
float *weights;
int* columndomain;


extern FILE* subgraph;
extern FILE* subgraphFO;
extern FILE* subgraphFI;
extern oned_csr_graph g; //from bfs_reference for isisolated function
//this function is needed for roots generation
int isisolated(int64_t v) {
    if(my_pe()==VERTEX_OWNER(v)) return (g.rowstarts[VERTEX_LOCAL(v)]==g.rowstarts[VERTEX_LOCAL(v)+1]);
    return 0; //locally no evidence, allreduce required
}

struct queryInfo{
    int vcolumnID;
    int queryPE;
};
typedef struct ListNode {
    struct queryInfo val; 
    struct ListNode* next;
}ListNode;
//global
ListNode** queryinfolist;
ListNode* createNode(int vcolumnID, int queryPE) {
    ListNode* temp; 
    temp = (ListNode*)malloc(sizeof(struct ListNode));
    temp -> val.vcolumnID = vcolumnID;
    temp -> val.queryPE = queryPE;
    temp -> next = NULL;
    return temp;
}
ListNode* addNode(ListNode** head, int vcolumnID, int queryPE) {
    ListNode* node = createNode(vcolumnID, queryPE);
    if(*head == NULL) {
        *head = node; 
    }else {
        node -> next = *head;
        *head = node;
    }
}
void freelist(ListNode* head) {
    ListNode* tmp;
    while(head != NULL){
        tmp = head; 
        head = head -> next;
        free(tmp);
    }
}

void halfedgehndl(int from,void* data,int sz)
{  degrees[*(int*)data]++; }

void fulledgehndl(int frompe,void* data,int sz) {
    int vloc = *(int*)data;
    int64_t gtgt = *((int64_t*)(data+4));
    int64_t gsrc = VERTEX_TO_GLOBAL(my_pe(), vloc);
    assert(my_pe() == VERTEX_OWNER(gsrc));
    char edgetuple [100];
    SETCOLUMN(degrees[vloc]++,gtgt);
/*
#ifdef SSSP
    float w = ((float*)data)[3];
    weights[degrees[vloc]-1]=w;
    sprintf(edgetuple, "%lld %lld %f", gsrc, gtgt, w);
    fprintf(rangepartition, "%s\n", edgetuple);
#else
    sprintf(edgetuple, "%lld %lld", gsrc, gtgt);
    fprintf(rangepartition, "%s\n", edgetuple);
#endif
*/
}
void queryrequesthndl(int from, void* data, int sz){
    int64_t v = *(int64_t*)data;
    int vloc = VERTEX_LOCAL(v);
    int columnindex = *((int*)(data+8));
    addNode(&queryinfolist[vloc], columnindex, from); 
}
void queryreplyhndl(int from, void* data, int sz){
    int columnindex = *(int*)data;
    int replySubdomain = *((int*)(data+4));
    columndomain[columnindex] = replySubdomain; 
}
void dumpedgehndl(int from, void* data, int sz) {
    int src = *(int64_t*)data;
    int tgt = *((int64_t*)(data+8));
    int tgtdomain = *((int*)(data+16));
    int srcdomain = my_pe();
    if(srcdomain == tgtdomain){
        char edgetuple[100];
        sprintf(edgetuple, "%lld %lld %d", src, tgt, tgtdomain);
        fprintf(subgraph, "%s\n", edgetuple);
    }else {
        char edgetuple[100];
        sprintf(edgetuple, "%lld %lld %d", src, tgt, tgtdomain);
        fprintf(subgraphFO, "%s\n", edgetuple);
        char edgetuple1[100];
        sprintf(edgetuple1, "%lld %lld %d", tgt, src, tgtdomain);
        fprintf(subgraphFI, "%s\n", edgetuple1);
    }

}

void send_dumpedge(int64_t src, int64_t tgt, int tgtdomain, int pe) {
    int dumpedge[5];
    memcpy(dumpedge, &src, 8);
    memcpy(dumpedge+2, &tgt, 8);
    memcpy(dumpedge+4, &tgtdomain, 4);
    aml_send(&dumpedge, 1, 20, pe);
}
void send_queryrequest(int64_t vertex, int columnindex, int pe) {
    int queryinfo[3];
    memcpy(queryinfo, &vertex, 8);
    memcpy(queryinfo+2, &columnindex, 4);
    aml_send(&queryinfo,1,12,pe);
}
void send_queryreply(int columnindex, int replySubdomain, int pe) {
    int replyinfo[2];
    memcpy(replyinfo, &columnindex, 4);
    memcpy(replyinfo+1, &replySubdomain, 4);
    aml_send(&replyinfo,1,8,pe);
}
void send_half_edge (int64_t src,int64_t tgt) {
    int pe=VERTEX_OWNER(src);
    int vloc=VERTEX_LOCAL(src);
    aml_send(&vloc,1,4,pe);
    if(tgt>=nverts_known) nverts_known=tgt+1;
}
#ifdef SSSP
void send_full_edge (int64_t src,int64_t tgt,float w) {
    int pe=VERTEX_OWNER(src);
    int vloc[4];
    vloc[0]=VERTEX_LOCAL(src);
    memcpy(vloc+1,&tgt,8);
    memcpy(vloc+3,&w,4);
    aml_send(vloc,1,16,pe);
}
#else
void send_full_edge (int64_t src,int64_t tgt) {
    int pe=VERTEX_OWNER(src);
    int vloc[3];
    vloc[0]=VERTEX_LOCAL(src);
    memcpy(vloc+1,&tgt,8);
    aml_send(vloc,1,12,pe);
}
#endif

void ParMetis_GPart(MPI_Comm comm, idx_t* xadj, idx_t xadjlen, idx_t* adjncy, idx_t adjlen, idx_t* vtxdist, idx_t vtxlen, idx_t* part) {
    idx_t ncon, nparts, npes, mype, edgecut;
    /* graph_t graph, mgraph; */
    idx_t numflag = 0, wgtflag = 0, options[10];
    real_t *tpwgts = NULL, ubvec[MAXNCON];
    int i,j;
    //for casting int64_t into int
    //idx_t* adjcny_new = (idx_t *)malloc(adjlen * sizeof(idx_t));
    //for(i =0; i < adjlen; i++) {
    //    adjcny_new[i] = (idx_t)adjncy[i];
    //}
    //xadj = (idx_t*)xadj;
    printf("%d\n", sizeof(idx_t));
    //vtxdist = (idx_t*)vtxdist;
    gkMPI_Comm_size(comm, &npes);
    gkMPI_Comm_rank(comm, &mype);

    int myid, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* ParallelReadGraph(&graph, filename, comm);*/

    gkMPI_Barrier(comm);
    nparts = num_pes();
    ncon = NCON;

//for(nparts = 2; nparts <= num_pes(); nparts *= 2) {
    tpwgts = rmalloc(ncon * nparts * 2, "TestParMetis_V3: tpwgts");
    rset(MAXNCON, 1.05, ubvec);
    /*======================================================================
     *     / ParMETIS_V3_PartKway
     *         /=======================================================================*/
    options[0] = 1;
    options[1] = 3;
    options[2] = 1;
    wgtflag = 0;
    numflag = 0;

    if (mype == 0)
        printf("\nParMETIS_V3_PartKway with ncon: %"PRIDX", nparts: %"PRIDX"\n", ncon, nparts);

    rset(nparts * ncon, 1.0 / (real_t) nparts, tpwgts);
    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, NULL,
            NULL, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    //free(tpwgts);
//}
    if(mype == 0){
        printf("edgecut:%d\n", edgecut);
    }
    
    j = 0;
    for (i = vtxdist[myid]; i < vtxdist[myid + 1] ; i++) {
        //printf("%d %d\n", i, part[j]);
        j++;
    }
    //free(adjcny_new);
}


void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g) {
    g->tg = tg;

    size_t i,j,k;

    int64_t nvert=tg->nglobaledges/2;
    nvert/=num_pes();
    nvert+=1;
    degrees=xcalloc(nvert,sizeof(int));

    aml_register_handler(halfedgehndl,1);
    int numiters=ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg);
    // First pass : calculate degrees of each vertex
    ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
        ptrdiff_t j;
        for (j = 0; j < bufsize; ++j) {
            int64_t v0 = get_v0_from_edge(&buf[j]);
            int64_t v1 = get_v1_from_edge(&buf[j]);
            if(v0==v1) continue;
            send_half_edge(v0, v1);
            send_half_edge(v1, v0);
        }
        aml_barrier();
    } ITERATE_TUPLE_GRAPH_END;

    int64_t nglobalverts = 0;
    aml_long_allmax(&nverts_known);
    nglobalverts=nverts_known+1;
    g->nglobalverts = nglobalverts;
    printf("all max nglobalverts %lld\n", nglobalverts);
    idx_t nlocalverts = (idx_t)DIV_SIZE(nglobalverts + num_pes() - 1 - my_pe());
    g->nlocalverts = nlocalverts;

    //graph stats printing
#ifdef DEBUGSTATS
    long maxdeg=0,isolated=0,totaledges=0,originaledges;
    long maxlocaledges,minlocaledges;
    for(i=0;i<g->nlocalverts;i++) {
        long deg = degrees[i];
        totaledges+=deg;
        if(maxdeg<deg) maxdeg=deg;
        if(!deg) isolated++;
    }
    originaledges=totaledges;
    maxlocaledges=totaledges;
    minlocaledges=totaledges;
    aml_long_allmax(&maxdeg);
    aml_long_allsum(&isolated);
    aml_long_allsum(&totaledges);
    aml_long_allmin(&minlocaledges);
    aml_long_allmax(&maxlocaledges);
    long averageedges = totaledges/num_pes();
    double disbalance = (double)(maxlocaledges-minlocaledges)/(double)averageedges * 100.0;
    if(!my_pe()) printf("\n maxdeg %lld verts %lld, isolated %lld edges %lld\n\t A max %ld min %ld ave %ld delta %ld percent %3.2f\n ",
            maxdeg,g->nglobalverts,isolated,totaledges,maxlocaledges,minlocaledges,averageedges,maxlocaledges-minlocaledges,disbalance);

    // finished stats printing

    g->notisolated=g->nglobalverts-isolated;
#endif
    printf("size of rowstarts %d\n", nlocalverts);  
    idx_t *rowstarts = xmalloc((nlocalverts + 1) * sizeof(idx_t));
    g->rowstarts = rowstarts;

    rowstarts[0] = 0;
    for (i = 0; i < nlocalverts; ++i) {
        rowstarts[i + 1] = rowstarts[i] + (i >= nlocalverts ? 0 : degrees[i]);
        degrees[i] = rowstarts[i];
    }

    idx_t nlocaledges = rowstarts[nlocalverts];
    g->nlocaledges = nlocaledges;

    idx_t colalloc = BYTES_PER_VERTEX*nlocaledges;
    colalloc += (4095);
    colalloc /= 4096;
    colalloc *= 4096;
    column = xmalloc(colalloc);
    aml_barrier();
#ifdef SSSP
    weights = xmalloc(4*nlocaledges);
    g->weights = weights;
    aml_barrier();
#endif
    //long allocatededges=colalloc;
    g->column = column;

    aml_register_handler(fulledgehndl,1);
    //Next pass , actual data transfer: placing edges to its places in column and hcolumn
    ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
        ptrdiff_t j;
        for (j = 0; j < bufsize; ++j) {
            int64_t v0 = get_v0_from_edge(&buf[j]);
            int64_t v1 = get_v1_from_edge(&buf[j]);
            if(v0==v1) continue;
#ifdef SSSP
            send_full_edge(v0, v1,wbuf[j]);
            send_full_edge(v1, v0,wbuf[j]);
#else               
            send_full_edge(v0, v1);
            send_full_edge(v1, v0);
#endif
        }
        aml_barrier();
    } ITERATE_TUPLE_GRAPH_END;

    free(degrees);
    //use parmetis and csr array partition graph
    idx_t mype, npes;
    MPI_Comm comm_temp;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_temp);
    gkMPI_Comm_size(comm_temp, &npes);
    gkMPI_Comm_rank(comm_temp, &mype);
    idx_t* vtxdist = malloc((num_pes() + 1)*sizeof(idx_t));
    idx_t* part = (idx_t *) malloc((nlocalverts + 3) * sizeof(idx_t));
    for(j = 0; j <= num_pes(); ++j) {
        vtxdist[j] = j*SIZE_PER_RANK;    
    }
    ParMetis_GPart(comm_temp, rowstarts, nlocalverts + 1,(idx_t*)column, nlocaledges, vtxdist, num_pes() + 1, part);
    printf("end partition\n");
    //two rounds one-side message 
    queryinfolist = (ListNode**)malloc(nlocalverts * sizeof(ListNode*));
    for(j = 0; j < nlocalverts; ++j) {
        queryinfolist[j] = NULL;
    }
    aml_register_handler(queryrequesthndl,1);
    //go through adjcny and construct query package then send to target
    for(j = 0; j < nlocaledges; ++j) {
        int target_pe = VERTEX_OWNER(column[j]);
	//printf("send_queryrequest pe:%d\n", target_pe);
        send_queryrequest(column[j], j, target_pe);
    }
    aml_barrier(); 
    columndomain = (int*)malloc(nlocaledges * sizeof(int));
    aml_register_handler(queryreplyhndl,1);
    //go through linklist array and send the reply back to query process
    for(j = 0; j < nlocalverts; ++j) {
        ListNode* temp = queryinfolist[j];
        if(temp == NULL) {
            //printf("this vertex no process query: %d, my pe is %d\n", j, my_pe());
        }else {
            while(temp != NULL){
                int columnindex = temp->val.vcolumnID;
                int querype = temp->val.queryPE;
                int replySubdomain = part[j];
	//	        printf("send_queryreply pe:%d\n", querype);
                send_queryreply(columnindex, replySubdomain, querype);
                temp = temp -> next;
            }
        }
    }
    aml_barrier();
    //free all the linked list
    for(j = 0; j < nlocalverts; ++j) {
        freelist(queryinfolist[j]);
    }
    free(queryinfolist);
    aml_register_handler(dumpedgehndl, 1);
    for(j = 0; j < nlocalverts; ++j) {
        for(k = rowstarts[j]; k < rowstarts[j+1]; ++k) {
            int64_t src = VERTEX_TO_GLOBAL(my_pe(), j);
            int64_t tgt = column[k];
            int srcdomain = part[j];
            int tgtdomain = columndomain[k];
	 //       printf("send_dumpedge pe:%d\n", srcdomain);
            send_dumpedge(src, tgt, tgtdomain, srcdomain);
        }
    }
}

void free_oned_csr_graph(oned_csr_graph* const g) {
    if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
    if (g->column != NULL) {free(g->column); g->column = NULL;}
#ifdef SSSP
    if (g->weights != NULL) {free(g->weights); g->weights = NULL;}
#endif
}
