
CC=gcc
CXX=g++
CFLAGS=-Wall -O3 -march=native -fopenmp


all: recall bench

recall:
	mkdir -p bin
	for f in aband blast simdblast; do for bw in 8 16 24 32 40 48 56 64; do $(CXX) $(CFLAGS) -std=c++11 -DMAIN -DBW=$$bw $$f.cc -o bin/$$f-$$bw; done; done
	for bw in 8 16 24 32 40 48 56 64; do $(CC) $(CFLAGS) -std=c99 -DMAIN -DBW=$$bw full.c -o bin/full-$$bw; done

bench:
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
	$(CXX) $(CFLAGS) -std=c++11 -o bin/bench -DBW=64 main.cc aband.cc blast.cc simdblast.cc wave/DB.o wave/QV.o wave/align.o ssw.o parasail/cpuid.o parasail/io.o parasail/matrix_lookup.o parasail/memory.o parasail/memory_sse.o parasail/time.o sg_striped_sse41_128_16.o full.o -DBENCH


