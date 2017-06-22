
CC=gcc
CXX=g++
CFLAGS=-Wall -O3 -march=native -fopenmp


all: recall bench

recall:
	mkdir -p bin
	for f in aband blast simdblast; do for bw in 8 16 24 32 40 48 56 64; do $(CXX) $(CFLAGS) -DMAIN -DBW=$$bw $$f.cc -o bin/$$f-$$bw; done; done
	for bw in 8 16 24 32 40 48 56 64; do $(CC) $(CFLAGS) -std=c99 -DMAIN -DBW=$$bw full.c -o bin/full-$$bw; done

bench:
	$(CC) $(CFLAGS) -c -o wave/DB.o wave/DB.c
	$(CC) $(CFLAGS) -c -o wave/QV.o wave/QV.c
	$(CC) $(CFLAGS) -c -o wave/align.o wave/align.c
	$(CC) $(CFLAGS) -c -o ssw.o ssw.c
	$(CC) $(CFLAGS) -c -o full.o full.c
	$(CXX) $(CFLAGS) -o bin/bench main.cc aband.cc blast.cc simdblast.cc wave/DB.o wave/QV.o wave/align.o ssw.o full.o -DBENCH -DBW=32


