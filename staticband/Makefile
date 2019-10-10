
# CC  = gcc
# CFLAGS = -Wall -Wno-unused-function -std=c99 -O3 -march=native

SCALAR_SRCS = scalar.c
SIMD_SRCS   = vertical.c diagonal.c striped.c
OBJS = $(SCALAR_SRCS:.c=.o) $(SIMD_SRCS:.c=.sse.o) $(SIMD_SRCS:.c=.avx.o)
TGTS = staticband.a

all: $(TGTS)
clean:
	$(RM) -rf $(TGTS) $(OBJS)

$(OBJS): $(SCALAR_SRCS) $(SIMD_SRCS)
	$(CC) $(CFLAGS) $(OFLAGS) `if echo $@ | grep avx >/dev/null; then echo "-DBENCH_ARCH_AVX"; else echo "-DBENCH_ARCH_SSE"; fi` -c -o $@ `echo $@ | sed 's/\..*o$$/\.c/g'`

$(TGTS): $(OBJS)
	$(AR) -crs $(TGTS) $^


