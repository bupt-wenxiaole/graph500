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

int64_t nverts_known = 0;
//number of edges of cross rank
int64_t nedges_boarder = 0;
int *degrees;
int64_t *column;
float *weights;
extern oned_csr_graph g; //from bfs_reference for isisolated function
extern FILE* subgraph;
extern FILE* subgraphFO;
extern FILE* subgraphFI;

struct crossedge{
	int64_t gsrc;
	int64_t gtgt;
	int srcowner;
};
//struct crossedge* tmpbuf;
float *crossweights;
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
	SETCOLUMN(degrees[vloc]++,gtgt);
	int64_t gsrc = VERTEX_TO_GLOBAL(my_pe(), vloc);
	int tgtowner = VERTEX_OWNER(gtgt);
	#ifdef SSSP
	float w = ((float*)data)[3];
	weights[degrees[vloc]-1]=w;
	#endif
	assert(my_pe() == VERTEX_OWNER(gsrc));
	int srcowner = my_pe();
	if(tgtowner == my_pe()){
		char edgetuple [100];
		#ifdef SSSP
		sprintf(edgetuple, "%lld %lld %d %f", gsrc, gtgt, tgtowner, w);
		fprintf(subgraph, "%s\n", edgetuple);
		#else
		sprintf(edgetuple, "%lld %lld %d", gsrc, gtgt, tgtowner);
		fprintf(subgraph, "%s\n", edgetuple);
		#endif
	}
	else {
		//dump this edge tuple to Frank.O
		char edgetuple[100];
		#ifdef SSSP
		sprintf(edgetuple, "%lld %lld %d %f", gsrc, gtgt, tgtowner, w);
		fprintf(subgraphFO, "%s\n", edgetuple);
		/*struct crossedge edge;
		edge.gsrc = gsrc;
		edge.gtgt = gtgt;
		edge.srcowner = srcowner;
		tmpbuf[nedges_boarder] = edge;
		crossweights[nedges_boarder++] = w;*/
		#else
		sprintf(edgetuple, "%lld %lld %d", gsrc, gtgt, tgtowner);
		fprintf(subgraphFO, "%s\n", edgetuple);
		/*struct crossedge edge;
		edge.gsrc = gsrc;
		edge.gtgt = gtgt;
		edge.srcowner = srcowner;
		tmpbuf[nedges_boarder++] = edge;*/
		#endif
	}


}
//edgedumphndl的内容merge到fulledgehndl中去
//this function is handler for dump edge which belongs to F.I
void dumphndl(int frompe,void* data,int sz) {
	int64_t gsrc = *(int64_t*)data;
	int64_t gtgt =  *((int64_t*)(data+8));
	int srcowner = *((int*)(data+16));
	char edgetuple[100];
	#ifdef SSSP
	float w = *((float*)(data+20));
	sprintf(edgetuple, "%lld %lld %d %f", gsrc, gtgt, srcowner, w);
	fprintf(subgraphFI, "%s\n", edgetuple);
	#else
	sprintf(edgetuple, "%lld %lld %d", gsrc, gtgt, srcowner);
	fprintf(subgraphFI, "%s\n", edgetuple);
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
	size_t nlocalverts = VERTEX_LOCAL(nglobalverts + num_pes() - 1 - my_pe());
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
	//tmpbuf=xcalloc(nlocaledges+1,sizeof(struct crossedge));
	aml_barrier();
#ifdef SSSP
	//crossweights=xcalloc(nlocaledges+1,sizeof(float));
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
			//if(VERTEX_OWNER(v0) == VERTEX_OWNER(v1)){
			//	printf("two vertex belongs to one partition!\n");
			//}
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
	aml_barrier();
	//Third pass, transfer the across rank edge to it's owner rank according with target index
	aml_register_handler(dumphndl,1);
	int srcowner = my_pe();
	//assert(nedges_boarder <= nlocaledges);
	for(i = 0; i < nlocalverts; i++) {
		for(j = g->rowstarts[i]; j < g->rowstarts[i + 1]; j++) {
			int64_t gsrc, gdst;
			gsrc = VERTEX_TO_GLOBAL(srcowner, i);
			gdst = COLUMN(j);

			printf("src:%lld  dst:%lld\n", gsrc, gdst);
			int tgtowner = VERTEX_OWNER(gdst);
#ifdef SSSP
//printf("send edge gsrc:%lld gtgt %lld srcowner %lld weight %f, the %lldth send,nedges_boarder %lld, my rank is %d\n",
// tmpbuf[k].gsrc, tmpbuf[k].gtgt, tmpbuf[k].srcowner, crossweights[k], k, nedges_boarder, my_pe());
		    int vgolable[6];
		    memcpy(vgolable,&gsrc,8);
		    memcpy(vgolable+2,&gdst,8);
		    memcpy(vgolable+4,&srcowner,4);
		    memcpy(vgolable+5,&g->weights[j],4);
		    aml_send(vgolable,1,24,tgtowner);
#else
			int vgolable[5];
			memcpy(vgolable,&gsrc,8);
			memcpy(vgolable+2,&gdst,8);
			memcpy(vgolable+4,&srcowner,4);
			aml_send(vgolable,1,20,tgtowner);
#endif
		}
	}
	aml_barrier();
}

void free_oned_csr_graph(oned_csr_graph* const g) {
	if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
	if (g->column != NULL) {free(g->column); g->column = NULL;}
#ifdef SSSP
	if (g->weights != NULL) {free(g->weights); g->weights = NULL;}
#endif
}
