
CC  = gcc
CFLAGS = -Wall -Wno-unused-function -std=gnu99 -DBENCH
OFLAGS = -O3 -march=native

# CXX = g++
# CXXFLAGS = -Wall -Wno-unused-function -std=gnu++11 -O3 -msse4.1 -fopenmp

BENCH_MODULES = full staticband ssw parasail wavefront # adaptiveband
LIBS = $(join $(BENCH_MODULES), $(BENCH_MODULES:%=/%.a))
SRCS = main.c
TGTS = bench

all: $(TGTS)

$(LIBS):
	$(MAKE) -C `echo $@ | sed 's/\/.*//g'` CFLAGS="$(CFLAGS)" OFLAGS="$(OFLAGS)" all

$(TGTS): $(SRCS) $(LIBS)
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ $(SRCS) $(LIBS)

clean:
	rm -rf $(TGTS) $(join $(BENCH_MODULES), $(BENCH_MODULES:%=/%.a)) $(BENCH_MODULES:%=%/*.o)

