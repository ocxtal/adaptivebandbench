
# CC  = gcc
# CFLAGS = -Wall -Wno-unused-function -std=c99 -O3 -march=native

SRCS = full.c
OBJS = $(SRCS:.c=.o)
TGTS = full.a

all: $(TGTS)
clean:
	$(RM) -rf $(TGTS) $(OBJS)

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $(@:.o=.c)

$(TGTS): $(OBJS)
	$(AR) -crs $(TGTS) $^


