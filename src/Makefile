CFLAGS = -Drestrict=__restrict__ -O3 -DGRAPH_GENERATOR_MPI -DREUSE_CSR_FOR_VALIDATION -I../aml
LDFLAGS = -lpthread
MPICC = mpicc
PARMETISCFLAGS = -I/home/xwen/parmetis-4.0.3/include -I/home/xwen/parmetis-4.0.3/metis/GKlib -I/home/xwen/parmetis-4.0.3/metis/include -I/home/xwen/parmetis-4.0.3/programs/. \
-DLINUX -D_FILE_OFFSET_BITS=64 -fPIC -DNDEB -DHAVE_EXECINFO_H -DHAVE_GETLINE 


all: graph500_reference_bfs_sssp graph500_reference_bfs 
#graph500_custom_bfs graph500_custom_bfs_sssp

GENERATOR_SOURCES = ../generator/graph_generator.c ../generator/make_graph.c ../generator/splittable_mrg.c ../generator/utils.c
SOURCES = main.c utils.c validate.c ../aml/aml.c
HEADERS = common.h csr_reference.h bitmap_reference.h

graph500_reference_bfs_sssp: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)  csr_reference.c sssp_reference.c
	$(MPICC) $(CFLAGS) $(PARMETISCFLAGS) -DSSSP $(LDFLAGS) -o graph500_reference_bfs_sssp bfs_reference.c sssp_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -rdynamic /home/xwen/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a /home/xwen/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a   -lm

graph500_reference_bfs: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS)  $(PARMETISCFLAGS) $(LDFLAGS) -o graph500_reference_bfs bfs_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -rdynamic /home/xwen/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a  /home/xwen/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a  -lm

graph500_custom_bfs: bfs_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS)) $(LDFLAGS) -o graph500_custom_bfs bfs_custom.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_custom_bfs_sssp: bfs_custom.c sssp_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o graph500_custom_bfs_sssp bfs_custom.c sssp_custom.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm

clean:
	-rm -f graph500_*
