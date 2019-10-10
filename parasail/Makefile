
# CC  = gcc
# CFLAGS = -Wall -Wno-unused-function -std=c99 -O3 -march=native

SRCS = cpuid.c io.c matrix_lookup.c memory.c memory_sse.c sg_striped_sse41_128_16.c time.c
OBJS = $(SRCS:.c=.o)
TGTS = parasail.a

all: $(TGTS)
clean:
	$(RM) -rf $(TGTS) $(OBJS)

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(OFLAGS) -I. -c -o $@ $(@:.o=.c)

$(TGTS): $(OBJS)
	$(AR) -crs $(TGTS) $^


