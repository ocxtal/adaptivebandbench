
CC=gcc
CXX=g++
CFLAGS=-Wall -std=c99 -O3 -msse4.1 -fopenmp
CXXFLAGS=-Wall -std=gnu++11 -O3 -msse4.1 -fopenmp

BENCH_SRCS=main.cc aband.cc blast.cc simdblast.cc
# BENCH_SRCS=main.cc alinear.cc aaffine.cc blast.cc simdblast.cc
BENCH_MODULES=wave/DB.o wave/QV.o wave/align.o ssw.o parasail/cpuid.o parasail/io.o parasail/matrix_lookup.o parasail/memory.o parasail/memory_sse.o parasail/time.o sg_striped_sse41_128_16.o full.o

all: recall bench

recall:
	mkdir -p bin
	for f in aband blast simdblast; do for bw in 8 16 24 32 40 48 56 64; do $(CXX) $(CXXFLAGS) -DMAIN -DBW=$$bw $$f.cc -o bin/$$f-$$bw; done; done
	for bw in 8 16 24 32 40 48 56 64; do $(CC) $(CFLAGS) -DMAIN -DBW=$$bw full.c -o bin/full-$$bw; done



bench.modules:
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

bench: bench.32 bench.64 bench.aband.32 bench.aband.64 bench.blast.32 bench.blast.64

bench.32: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.32 -DBW=32 -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.64: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.64 -DBW=64 -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.96: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.96 -DBW=96 -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.128: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.128 -DBW=128 -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.aband.32: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.aband.32 -DBW=32 -DDEBUG_PATH -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.aband.64: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.blast.64 -DBW=64 -DDEBUG_PATH -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.blast.32: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.blast.32 -DBW=32 -DDEBUG_PATH -DDEBUG_BLAST -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

bench.blast.64: bench.modules
	$(CXX) $(CXXFLAGS) -o bin/bench.blast.64 -DBW=64 -DDEBUG_PATH -DDEBUG_BLAST -DBENCH $(BENCH_SRCS) $(BENCH_MODULES)

