
# CC  = gcc
# CFLAGS = -Wall -Wno-unused-function -std=c99 -O3 -march=native

SRCS = align.c DB.c QV.c
OBJS = $(SRCS:.c=.o)
TGTS = wavefront.a

all: $(TGTS)
clean:
	$(RM) -rf $(TGTS) $(OBJS)

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(OFLAGS) -I. -c -o $@ $(@:.o=.c)

$(TGTS): $(OBJS)
	$(AR) -crs $(TGTS) $^


