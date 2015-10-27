
/**
 * @file main.cc
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "util.h"

int rognes_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge);
int rognes_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge);
int diag_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge);
int diag_affine(
	char const *a,
	uint64_t alen,
	char const *b,
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
	bench_t rl, ra, dl, da;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	unsigned long s = (argc > 3) ? atoi(argv[3]) : tv.tv_usec;
	srand(s);
	printf("%lu\n", s);

	a = rseq(len * 9 / 10);
	b = mseq(a, 10, 40, 40);
	at = rseq(len / 10);
	bt = rseq(len / 10);
	a = (char *)realloc(a, 2*len); strcat(a, at); free(at);
	b = (char *)realloc(b, 2*len); strcat(b, bt); free(bt);

//	printf("%s\n%s\n", a, b);

	printf("len:\t%d\ncnt:\t%d\n", len, cnt);

	/* rognes */
	bench_init(rl);
	bench_init(ra);
	bench_start(rl);
	for(i = 0; i < cnt; i++) {
		rognes_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	}
	bench_end(rl);
	bench_start(ra);
	for(i = 0; i < cnt; i++) {
		rognes_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	}
	bench_end(ra);
	printf("rognes:\t%lld\t%lld\n",
		bench_get(rl) / 1000,
		bench_get(ra) / 1000);

	/* diag */
	bench_init(dl);
	bench_init(da);
	bench_start(dl);
	for(i = 0; i < cnt; i++) {
		diag_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	}
	bench_end(dl);
	bench_start(da);
	for(i = 0; i < cnt; i++) {
		diag_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	}
	bench_end(da);
	printf("diag:\t%lld\t%lld\n",
		bench_get(dl) / 1000,
		bench_get(da) / 1000);

	free(a);
	free(b);
	return 0;
}

/**
 * end of main.cc
 */
