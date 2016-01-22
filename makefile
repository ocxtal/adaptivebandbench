
CC=gcc
CXX=g++
CFLAGS=""

all:
	$(CC) -Wall -O3 -msse4.1 -std=c99 -c full.c
	$(CC) -Wall -O3 -msse4.1 -std=c99 -c ssw.c
	$(CXX) -o bench -Wall -O3 -msse4.1 main.cc rognes.cc ddiag.cc diag.cc blast.cc simdblast.cc full.o ssw.o -DBENCH -DALL -DSSW
