
# CC  = gcc
# CFLAGS = -Wall -Wno-unused-function -std=c99 -O3 -march=native

SRCS = ssw.c
OBJS = $(SRCS:.c=.o)
TGTS = ssw.a

all: $(TGTS)
clean:
	$(RM) -rf $(TGTS) $(OBJS)

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $(@:.o=.c)

$(TGTS): $(OBJS)
	$(AR) -crs $(TGTS) $^


