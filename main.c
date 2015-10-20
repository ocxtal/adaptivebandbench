
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "util.h"

int rognes(
	uint8_t const *a,
	uint64_t alen,
	uint8_t const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge);
int diag(
	uint8_t const *a,
	uint64_t alen,
	uint8_t const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge);



/**
 * random sequence generator, modifier.
 * rseq generates random nucleotide sequence in ascii,
 * mseq takes ascii sequence, modifies the sequence in given rate.
 */
static char rbase(void)
{
	switch(rand() % 4) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'A';
	}
}

char *rseq(int len)
{
	int i;
	char *seq;
	seq = (char *)malloc(sizeof(char) * (len + 1));

	for(i = 0; i < len; i++) {
		seq[i] = rbase();
	}
	seq[len] = '\0';
	return(seq);
}

char *mseq(char const *seq, int x, int ins, int del)
{
	int i;
	int len = strlen(seq);
	char *mod, *ptr;
	mod = (char *)malloc(sizeof(char) * 2 * (len + 1));

	ptr = mod;
	for(i = 0; i < len; i++) {
		if(rand() % x == 0) { *ptr++ = rbase(); }
		else if(rand() % ins == 0) { *ptr++ = rbase(); i--; }
		else if(rand() % del == 0) { /* skip a base */ }
		else { *ptr++ = seq[i]; }
	}
	*ptr = '\0';
	return(mod);
}

int main(int argc, char *argv[])
{
	int i;
	char *a, *b, *at, *bt;
	int len = (argc > 1) ? atoi(argv[1]) : 1000;
	int cnt = (argc > 2) ? atoi(argv[2]) : 1000;
	bench_t r, d;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	unsigned long s = (argc > 3) ? atoi(argv[3]) : tv.tv_usec;
	srand(s);
	printf("%lu\n", s);

	a = rseq(len * 9 / 10);
	b = mseq(a, 10, 40, 40);
	at = rseq(len / 10);
	bt = rseq(len / 10);
	a = realloc(a, 2*len); strcat(a, at); free(at);
	b = realloc(b, 2*len); strcat(b, bt); free(bt);

//	printf("%s\n%s\n", a, b);

	bench_init(r);
	bench_init(d);

	bench_start(r);
	for(i = 0; i < cnt; i++) {
		rognes((uint8_t const *)a, strlen(a), (uint8_t const *)b, strlen(b), 2, -3, -5, -1);
	}
	bench_end(r);
	printf("rognes: %lld\n", bench_get(r) / 1000);

	bench_start(d);
	for(i = 0; i < cnt; i++) {
		diag((uint8_t const *)a, strlen(a), (uint8_t const *)b, strlen(b), 2, -3, -5, -1);
	}
	bench_end(d);
	printf("diag:   %lld\n", bench_get(d) / 1000);

	free(a);
	free(b);
	return 0;
}
