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
int *degrees;
int64_t *column;
float *weights;

int *tgt_degrees; //save the queryed tgt degree
int *edge_part;  //recorde the edge's owner partition
int* vertex_replica_num; //recorde the 

extern oned_csr_graph g; //from bfs_reference for isisolated function

extern FILE* subgraphedgesfile;
extern FILE* isolatedfile;
extern FILE* masterfile;
extern FILE* mirrorfile;
//this function is needed for roots generation
int isisolated(int64_t v) {
	if(my_pe()==VERTEX_OWNER(v)) return (g.rowstarts[VERTEX_LOCAL(v)]==g.rowstarts[VERTEX_LOCAL(v)+1]);
	return 0; //locally no evidence, allreduce required
}


typedef struct vertex_spawn_info_node{
	int location;
	int isMaster;
	struct vertex_spawn_info_node* next;
}vertex_spawn_info_node;

vertex_spawn_info_node** vertex_spawn_info;

vertex_spawn_info_node* create_vsi_node(int location){
	vertex_spawn_info_node* temp;
	temp = (vertex_spawn_info_node*)malloc(sizeof(struct vertex_spawn_info_node));
	temp -> location = location;
	temp -> isMaster = 0;// mirror by default
	temp -> next = NULL;
	return temp;
} 
vertex_spawn_info_node* add_vsi_node(vertex_spawn_info_node**head, int location){
	vertex_spawn_info_node* node = create_vsi_node(location);
	if(*head == NULL){
        *head = node; 
    }else {
    	vertex_spawn_info_node* temp = *head;
    	while(temp -> next != NULL){
    		if(temp -> location == node -> location)
    			break;
    		else
    			temp = temp -> next;
    	}
    	if(temp -> location == node -> location)
    		return *head;
    	temp -> next = node;
    }
    return *head;
}

int list_length(vertex_spawn_info_node* head){
	vertex_spawn_info_node* temp = head;
	int len = 0;
	while(temp != NULL){
		len++;
		temp = temp -> next;
	}
	return len; 
}

void free_vsi_list(vertex_spawn_info_node* head){
	vertex_spawn_info_node* tmp;
    while(head != NULL){
        tmp = head; 
        head = head -> next;
        free(tmp);
    }
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
	//size_t len = g->nlocaledges + 1000;
	//printf("set column vloc %d\n", degrees[vloc]);
	SETCOLUMN(degrees[vloc]++,gtgt);
#ifdef SSSP
	float w = ((float*)data)[3];
	weights[degrees[vloc]-1]=w;
#endif
}
void queryrequesthndl(int from, void* data, int sz){
    //printf("receive message from %d\n", from);
    int64_t v = *(int64_t*)data;
    int vloc = VERTEX_LOCAL(v);
    int columnindex = *((int*)(data+8));
    addNode(&queryinfolist[vloc], columnindex, from); 
}
void queryreplyhndl(int from, void* data, int sz){
    int columnindex = *(int*)data;
    int reply_degree = *((int*)(data+4));
    tgt_degrees[columnindex] = reply_degree; 
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

void dump_edge_hndl(int frompe,void* data,int sz) {
	int64_t gsrc = *(int64_t*)data;
	int64_t gtgt =  *((int64_t*)(data+8));
	char edgetuple[100];
	#ifdef SSSP
	float w = *((float*)(data+16));
	sprintf(edgetuple, "%lld %lld %f", gsrc, gtgt, w);
	fprintf(subgraphedgesfile, "%s\n", edgetuple);
	#else
	sprintf(edgetuple, "%lld %lld", gsrc, gtgt);
	fprintf(subgraphedgesfile, "%s\n", edgetuple);
	#endif
}

void dump_vertex_hndl(int from, void* data, int sz) {
	int isMaster = *(int*)data;
	int64_t vertexid = *((int64_t*)(data+4));
	if(isMaster){
		int other_location_num  = (sz - 12)/4;
		char* vertextuple = malloc((sz - 12) * 2 + 20);
		char* offset = vertextuple;
		offset += sprintf(offset, "%lld ", vertexid);
		int i = 0;
		for(; i < other_location_num; ++i){
			int location = *((int*)(data + 12 +4*i));
			if(i == other_location_num - 1)
				offset += sprintf(offset, "%d", location);
			else
				offset += sprintf(offset, "%d ", location);
		}
		fprintf(masterfile, "%s\n", vertextuple);
		free(vertextuple);
	}else{
		int master_location = *((int*)(data+12));
		char* vertextuple = malloc(15);
		char* offset = vertextuple; 
		offset += sprintf(offset, "%lld ", vertexid);
		offset += sprintf(offset, "%d", master_location);
		fprintf(mirrorfile, "%s\n", vertextuple);
		free(vertextuple);
	}
	
}
void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g) {
	g->tg = tg;

	size_t i,j,k;
	int64_t nvert= 23947347;
	nvert/=num_pes();
	nvert+=1;
	degrees=xcalloc(nvert,sizeof(int));

	aml_register_handler(halfedgehndl,1);
	int numiters=ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg);
	// First pass : calculate degrees of each vertex

	// This part for reading from graph500 kronecker generator
	/*ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
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
	
    */
	//This part for reading from real-world graphs which are edge-oriented
	
	{
		char rank_partion_path[100];
		ssize_t read = 0; 
		size_t len = 0;
		char* line = NULL;
		sprintf(rank_partion_path, "/mnt/nfs/xwen/partition.%d", my_pe());
		FILE* partitionfile = fopen(rank_partion_path, "r");
		if(partitionfile == NULL)
			exit(EXIT_FAILURE);
		printf("start send\n");
		while((read = getline(&line, &len, partitionfile)) != -1){
			//preprocessing input file to each line contains two elements
			char* token1 = strtok(line, " ");
			char* token2 = strtok(NULL, " ");
			int64_t v0 = atoll(token1);
			int64_t v1 = atoll(token2);
                        v0 -= 1;
                        v1 -= 1;
			if(v0==v1) continue;
			send_half_edge(v0, v1);
			
    	}
    	if(line) free(line);
    	fclose(partitionfile);
    	printf("end send\n");
	}
	
    aml_barrier();

    //This part for reading from real-world graph which are vertex-oriented
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

	int64_t colalloc = BYTES_PER_VERTEX*(nlocaledges);
	colalloc += (4095);
	colalloc /= 4096;
	colalloc *= 4096;
	//colalloc += 100;
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
	//This for reading from real-world graphs which are edge-oriented
/*	ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
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
*/
	printf("start 2th send\n");
	{
		char rank_partion_path[100];
		ssize_t read = 0; 
		size_t len = 0;
		char* line = NULL;
		sprintf(rank_partion_path, "/mnt/nfs/xwen/partition.%d", my_pe());
		FILE *partitionfile = fopen(rank_partion_path, "r");
		if(partitionfile == NULL)
			exit(EXIT_FAILURE);
		while((read = getline(&line, &len, partitionfile)) != -1){
			//preprocessing input file to each line contains two elements
			char* token1 = strtok(line, " ");
			char* token2 = strtok(NULL, " ");
			char* token3 = strtok(NULL, " ");

			int64_t v0 = atoll(token1);
			int64_t v1 = atoll(token2);
			//for making vertexID start from 0
                        v0 -= 1;
                        v1 -= 1;
			float weight = atof(token3);
			if(v0==v1) continue;
#ifdef SSSP

			send_full_edge(v0, v1,weight);
#else				
			send_full_edge(v0, v1);
#endif
    	}
    	if(line) free(line);
    	fclose(partitionfile);
    	printf("end 2th send\n");

	}
	
    aml_barrier();
	//two rounds one-side message 
    queryinfolist = (ListNode**)malloc(nlocalverts * sizeof(ListNode*));
    for(j = 0; j < nlocalverts; ++j) {
        queryinfolist[j] = NULL;
    }
    aml_register_handler(queryrequesthndl, 1);
    //go through adjcny and construct query package then send to target
    for(j = 0; j < nlocaledges; ++j) {
        int target_pe = VERTEX_OWNER(COLUMN(j));
		//printf("send_queryrequest pe:%d\n", target_pe);
        send_queryrequest(COLUMN(j), j, target_pe);
    }
    aml_barrier();
    tgt_degrees = xmalloc(4*nlocaledges);
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
                int reply_degree = rowstarts[j+1] - rowstarts[j];
	//	        printf("send_queryreply pe:%d\n", querype);
                send_queryreply(columnindex, reply_degree, querype);
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
    //DBH for assigning edges
    //start request
	edge_part = xmalloc(4*nlocaledges);
	vertex_spawn_info = (vertex_spawn_info_node**)malloc(nlocalverts * sizeof(vertex_spawn_info_node*));
	for(j = 0; j < nlocalverts; ++j) {
        vertex_spawn_info[j] = NULL;
    }
	for(j = 0; j < nlocalverts; ++j){
		int k = 0;
		int src_degree = rowstarts[j+1] - rowstarts[j];
		if(rowstarts[j] == rowstarts[j+1]){ //this means isolated vertex
			//add_vsi_node(&vertex_spawn_info[j], my_pe());
			//vertex_spawn_info[j] -> isMaster = 1;
			//printf("dump the isolated vertices\n");
			int64_t gsrc = VERTEX_TO_GLOBAL(my_pe(), j);
			char isolated_vertex[8];
			sprintf(isolated_vertex, "%lld", gsrc);
			fprintf(isolatedfile, "%s\n", isolated_vertex);
		}
		for(k = rowstarts[j]; k < rowstarts[j+1]; ++k){
			int tgt_degree = tgt_degrees[k];
			if(src_degree > tgt_degree){
				edge_part[k] = VERTEX_OWNER(COLUMN(k));
			}else if(src_degree = tgt_degree){
				edge_part[k] = VERTEX_TO_GLOBAL(my_pe(), j) > COLUMN(k) ? my_pe() : VERTEX_OWNER(COLUMN(k));
			}else {
				edge_part[k] = my_pe();
			}
			add_vsi_node(&vertex_spawn_info[j], edge_part[k]); 
		}
	}
	free(tgt_degrees);


	vertex_replica_num = xmalloc((nlocalverts) * sizeof(int));


	for(j = 0; j < nlocalverts; ++j){
		int len = list_length(vertex_spawn_info[j]);
		vertex_replica_num[j] = len;
		if(len <= 1)
			continue;
		int random_jump = random() % len;
		vertex_spawn_info_node* temp = vertex_spawn_info[j];
		while(random_jump--){
			temp = temp -> next;
		}
		temp -> isMaster = 1;
	}
	free(degrees);
	aml_barrier();

	//next round, send edges to it's partition
	aml_register_handler(dump_edge_hndl,1);

	for(j = 0; j < nlocalverts; ++j){
		int k = 0;
		for(k = rowstarts[j]; k < rowstarts[j+1]; ++k){
			int64_t gsrc, gdst;
			gsrc = VERTEX_TO_GLOBAL(my_pe(), j);
			gdst = COLUMN(k);
#ifdef SSSP 
			int vglobal[5];
			memcpy(vglobal,&gsrc,8);
			memcpy(vglobal+2,&gdst,8);
			memcpy(vglobal+4,&weights[k],4);
			aml_send(vglobal, 1, 20, edge_part[k]);

#else
			int vglobal[4];
			memcpy(vglobal,&gsrc,8);
			memcpy(vglobal+2,&gdst,8);
			aml_send(vglobal, 1, 16, edge_part[k]);

#endif

	   }

   }
   aml_barrier();
   //next pass, dump the vertices and compute statistics for vertex-cut communication
   float  cutted_vertex_num = 0;
   aml_register_handler(dump_vertex_hndl, 1);
   for(j = 0; j < nlocalverts; ++j){
   		if(vertex_replica_num[j] <= 1)
   			continue;  // no-cut vertex
   		else{
   			//these vertices envolve in communication
   			int master_location;
   			int* all_other_location = malloc(4*(vertex_replica_num[j]-1));
			cutted_vertex_num++;
   			vertex_spawn_info_node* temp = vertex_spawn_info[j];
   			int i = 0;
   			while(temp != NULL){
   				if(temp -> isMaster)
   					master_location = temp -> location;
   				else
   					all_other_location[i++] = temp -> location;
   				temp = temp -> next;

   			}
   			assert(i = vertex_replica_num[j] - 1);
   			//construct master message
   			int* completeMessage = malloc(4*(vertex_replica_num[j] + 3));
   			int isMaster = 1;
   			int64_t vertexID = VERTEX_TO_GLOBAL(my_pe(), j);
   			memcpy(completeMessage, &isMaster, 4);
   			memcpy(completeMessage + 1, &vertexID, 8);
   			memcpy(completeMessage + 3, all_other_location, 4*(vertex_replica_num[j]-1));
   			//int isMastertest = *(int*)completeMessage;
			//int64_t vertexidtest = *((int64_t*)(completeMessage+1));
			//printf("isMastertest %d\n", isMastertest);
			//printf("vertexidtest %lld\n", vertexidtest);
   			aml_send(completeMessage, 1, 4*(vertex_replica_num[j] + 2), master_location);
   			//construct mirror message
   			isMaster = 0;
   			for(i = 0; i < vertex_replica_num[j] - 1; ++i){
   				int message[4];
   				memcpy(message, &isMaster, 4);
   				memcpy(message+1, &vertexID, 8);
   				memcpy(message+3, &master_location, 4);
   				aml_send(message, 1, 16, all_other_location[i]);
   			}
   			free(all_other_location);
   			free(completeMessage);
		}
   }
   aml_barrier();
   float global_cutted_vertex_num;
   MPI_Reduce(&cutted_vertex_num, &global_cutted_vertex_num, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
   if(my_pe() == 0){
	printf("cutted vertex percent = %f\n", cutted_vertex_num/global_cutted_vertex_num);
   }
}

void free_oned_csr_graph(oned_csr_graph* const g) {
	if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
	if (g->column != NULL) {free(g->column); g->column = NULL;}
}
