
CC=gcc
CXX=g++
CFLAGS=-Wall -Wno-unused-function -std=c99 -O3 -msse4.1 -fopenmp
CXXFLAGS=-Wall -Wno-unused-function -std=gnu++11 -O3 -msse4.1 -fopenmp

BENCH_SRCS=main.cc blast.cc simdblast.cc adaptive.cc scalar.cc vertical.cc diagonal.cc striped.cc
BENCH_MODULES=wave/DB.o wave/QV.o wave/align.o ssw.o parasail/cpuid.o parasail/io.o parasail/matrix_lookup.o parasail/memory.o parasail/memory_sse.o parasail/time.o sg_striped_sse41_128_16.o full.o

all: bench

$(BENCH_MODULES):
	$(CC) $(CFLAGS) -c -o wave/DB.o wave/DB.c
	$(CC) $(CFLAGS) -c -o wave/QV.o wave/QV.c
	$(CC) $(CFLAGS) -c -o wave/align.o wave/align.c
	$(CC) $(CFLAGS) -c -o ssw.o ssw.c
	$(CC) $(CFLAGS) -c -o parasail/cpuid.o -I. parasail/cpuid.c
	$(CC) $(CFLAGS) -c -o parasail/io.o -I. parasail/io.c
	$(CC) $(CFLAGS) -c -o parasail/matrix_lookup.o -I. parasail/matrix_lookup.c
	$(CC) $(CFLAGS) -c -o parasail/memory.o -I. parasail/memory.c
	$(CC) $(CFLAGS) -c -o parasail/memory_sse.o -I. parasail/memory_sse.c
	$(CC) $(CFLAGS) -c -o parasail/time.o -I. parasail/time.c
	$(CC) $(CFLAGS) -c -o sg_striped_sse41_128_16.o -I. sg_striped_sse41_128_16.c
	$(CC) $(CFLAGS) -c -o full.o -I. full.c

bench: $(BENCH_MODULES)
	$(CXX) $(CXXFLAGS) -o bin/bench -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

clean:
	rm -rf *.o bin/*

