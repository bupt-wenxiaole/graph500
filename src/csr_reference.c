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
extern oned_csr_graph g; //from bfs_reference for isisolated function
extern FILE* rangepartition;
extern FILE* csrformat;

//this function is needed for roots generation
int isisolated(int64_t v) {
	if(my_pe()==VERTEX_OWNER(v)) return (g.rowstarts[VERTEX_LOCAL(v)]==g.rowstarts[VERTEX_LOCAL(v)+1]);
	return 0; //locally no evidence, allreduce required
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
#ifdef SSSP
	float w = ((float*)data)[3];
	weights[degrees[vloc]-1]=w;
	sprintf(edgetuple, "%lld %lld %f", gsrc, gtgt, w);
	fprintf(rangepartition, "%s\n", edgetuple);
#else
	sprintf(edgetuple, "%lld %lld", gsrc, gtgt);
	fprintf(rangepartition, "%s\n", edgetuple);
#endif
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

void ParMetis_GPart(MPI_Comm comm, int* xadj, int xadjlen, int64_t* adjcny, int adjlen, int* vtxdist, int vtxlen, idx_t* part) {
    idx_t ncon, nparts, npes, mype, edgecut;
    /* graph_t graph, mgraph; */
    idx_t numflag = 0, wgtflag = 0, options[10];
    real_t *tpwgts = NULL, ubvec[MAXNCON];
    int i,j;
    //for casting int64_t into int
    idx_t* adjcny_new = (idx_t *)malloc(adjlen * sizeof(idx_t));
    for(i =0; i < adjlen; i++) {
        adjcny_new[i] = (idx_t)adjcny[i];
    }
    //xadj = (idx_t*)xadj;
    printf("%d\n", sizeof(idx_t));
    //vtxdist = (idx_t*)vtxdist;
    gkMPI_Comm_size(comm, &npes);
    gkMPI_Comm_rank(comm, &mype);

    int myid, numprocs, namelen;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    /* ParallelReadGraph(&graph, filename, comm);*/

    gkMPI_Barrier(comm);
    nparts = num_pes();
    ncon = NCON;

    /* tpwgts = (real_t*)malloc((ncon * nparts + 3) * sizeof(real_t));*/
    tpwgts = rmalloc(ncon * nparts * 2, "TestParMetis_V3: tpwgts");
    rset(MAXNCON, 1.05, ubvec);

    /*graph.vwgt = ismalloc(graph.nvtxs * 5, 1, "TestParMetis_GPart: vwgt");*/


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
    ParMETIS_V3_PartKway((idx_t*)vtxdist, (idx_t*)xadj, (idx_t*)adjcny_new, NULL,
            NULL, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    if(mype == 0){
    	printf("edgecut:%d\n", edgecut);
    }
    
    j = 0;
    for (i = vtxdist[myid]; i < vtxdist[myid + 1] ; i++) {
        printf("%d %d\n", i, part[j]);
        j++;
    }
    free(adjcny_new);
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
	size_t nlocalverts = (size_t)DIV_SIZE(nglobalverts + num_pes() - 1 - my_pe());
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
	unsigned int *rowstarts = xmalloc((nlocalverts + 1) * sizeof(int));
	g->rowstarts = rowstarts;

	rowstarts[0] = 0;
	for (i = 0; i < nlocalverts; ++i) {
		rowstarts[i + 1] = rowstarts[i] + (i >= nlocalverts ? 0 : degrees[i]);
		degrees[i] = rowstarts[i];
	}

	size_t nlocaledges = rowstarts[nlocalverts];
	g->nlocaledges = nlocaledges;

	int64_t colalloc = BYTES_PER_VERTEX*nlocaledges;
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
    int* vtxdist = malloc((num_pes() + 1)*sizeof(int));
    idx_t* part = (idx_t *) malloc((nlocalverts + 3) * sizeof(idx_t));
    for(j = 0; j <= num_pes(); ++j) {
		vtxdist[j] = j*SIZE_PER_RANK;
	}
    ParMetis_GPart(comm_temp, rowstarts, nlocalverts + 1, column, nlocaledges, vtxdist, num_pes() + 1, part);
    printf("end partition\n");
    int subdomainNum = num_pes();
    int64_t subdomainVertexNum[subdomainNum];
    memset(subdomainVertexNum, 0, sizeof(subdomainVertexNum));
    for(j = 0; j < nlocalverts; ++j) {
        subdomainVertexNum[part[j]]++; 
    }
    for(j = 0; j < subdomainNum; ++j) {
        aml_long_allsum(subdomainVertexNum[j]);
    }
/*
	fprintf(csrformat, "xadj\n");
	fprintf(csrformat, "%d\n", nlocalverts + 1);
	for (j = 0; j < nlocalverts + 1; ++j) {
		fprintf(csrformat, "%d\n", rowstarts[j]);
	}
	fprintf(csrformat, "adjncy\n");
	fprintf(csrformat, "%d\n", nlocaledges);
	for(j = 0; j < nlocaledges; ++j) {
		fprintf(csrformat, "%lld\n ", column[j]);
	}
	fprintf(csrformat, "vtxdist\n");
	fprintf(csrformat, "%d\n", num_pes() + 1);
    for(j = 0; j <= num_pes(); ++j) {
        fprintf(csrformat, "%lld\n", j*SIZE_PER_RANK);
    }
*/
    //aml_register_handler(queryrequesthndl,1);
    //third pass, send query src and tgt partition's request to other rank
    //first, go through all the CSR adjcny's elements and caculate the owner of each elements by vtxdist
    //第一轮，每个rank遍历自己CSR中adjcny数组中的每个元素，查的就是这个数组中每个元素的归属
    //这轮pass只能发送查请求，那么对面只是接受下这个查请求，这时候需要保存（哪个rank 查哪个点这样的pair）
    //注意这里的内存组织结构，这里是主要的内存消耗区域  

    //开始下一轮pass, 下一轮pass根据上一轮保存的那个点的pair查自己的part数组，然后向发起查询请求的rank发送回去发起rank的查询
    //这时候每个rank的接受handler接到这个reply后保存自己要查的点的rank。
    //开始下一轮pass, 下一轮pass已经知道了自己这边所有src和tgt的归属，所以直接把这条边发给需要的归属地，并在发送时标识最后dump在那个文件中
    //这轮pass接收到边之后直接写文件，根据标识往进去写（这个思路同样用在hash分图里）



}

void free_oned_csr_graph(oned_csr_graph* const g) {
	if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
	if (g->column != NULL) {free(g->column); g->column = NULL;}
#ifdef SSSP
	if (g->weights != NULL) {free(g->weights); g->weights = NULL;}
#endif
}
