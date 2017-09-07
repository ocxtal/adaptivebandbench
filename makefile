
CC=gcc
CXX=g++
CFLAGS=-Wall -std=c99 -O3 -msse4.1 -fopenmp
CXXFLAGS=-Wall -std=gnu++11 -O3 -msse4.1 -fopenmp

BENCH_SRCS=main.cc blast.cc simdblast.cc
# BENCH_SRCS=main.cc alinear.cc aaffine.cc blast.cc simdblast.cc
BENCH_MODULES=wave/DB.o wave/QV.o wave/align.o ssw.o parasail/cpuid.o parasail/io.o parasail/matrix_lookup.o parasail/memory.o parasail/memory_sse.o parasail/time.o sg_striped_sse41_128_16.o full.o
ABAND_MODULES=$(shell seq -f'aband.%g.o ' 32 8 256)

all: recall bench

recall:
	mkdir -p bin
	for f in aband blast simdblast; do for bw in 8 16 24 32 40 48 56 64; do $(CXX) $(CXXFLAGS) -DMAIN -DBW=$$bw $$f.cc -o bin/$$f-$$bw; done; done
	for bw in 8 16 24 32 40 48 56 64; do $(CC) $(CFLAGS) -DMAIN -DBW=$$bw full.c -o bin/full-$$bw; done

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
	$(CC) $(CFLAGS) -c -o full.o full.c

$(ABAND_MODULES):
	$(CXX) $(CXXFLAGS) -c -o $@ -DBW=`echo $@ | cut -d'.' -f2` aband.cc

bench: $(BENCH_MODULES) $(ABAND_MODULES)
	$(CXX) $(CXXFLAGS) -o bin/bench -DBENCH $(BENCH_SRCS) $(BENCH_MODULES) $(ABAND_MODULES)

clean:
	rm *.o bin/*
