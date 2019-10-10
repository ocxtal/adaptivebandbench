
CC  = gcc
CFLAGS = -Wall -Wno-unused-function -std=c99 -O3 -msse4.1 -fopenmp

# CXX = g++
# CXXFLAGS = -Wall -Wno-unused-function -std=gnu++11 -O3 -msse4.1 -fopenmp

BENCH_MODULES = full staticband ssw parasail wavefront # adaptiveband
TGTS = bench

all: $(TGTS)

$(BENCH_MODULES):
	$(MAKE) -C $@ all

$(TGTS): $(BENCH_MODULES)
	$(CC) $(CFLAGS) -o $@ $(join $^, $(^:%=/%.a))

clean:
	rm -rf $(TGTS) $(join $(BENCH_MODULES), $(BENCH_MODULES:%=/%.a)) $(BENCH_MODULES:%=%/*.o)

