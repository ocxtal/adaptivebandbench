
CC=gcc
CXX=g++
CFLAGS=

all:
	$(CC) $(CFLAGS) -Wall -O3 -msse4.1 -std=c99 -c full.c
	$(CC) $(CFLAGS) -Wall -O3 -msse4.1 -std=c99 -c ssw.c
	$(CXX) $(CFLAGS) -o bench -Wall -O3 -msse4.1 main.cc rognes.cc ddiag.cc diag.cc blast.cc simdblast.cc full.o ssw.o -DBENCH -DALL -DSSW
