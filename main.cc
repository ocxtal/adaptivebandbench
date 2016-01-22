
/**
 * @file main.cc
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <gperftools/profiler.h>
#include <sys/time.h>
#include "util.h"
#include "full.h"
#include "ssw.h"

// #define ALL

int rognes_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int rognes_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);

int blast_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int blast_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);

int simdblast_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int simdblast_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);

int diag_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int diag_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);

int ddiag_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int ddiag_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);


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
	int const m = 1, x = -1, gi = -1, ge = -1;
	int const xt = 30;
	char *a, *b, *at, *bt;
	int len = (argc > 1) ? atoi(argv[1]) : 1000;
	int cnt = (argc > 2) ? atoi(argv[2]) : 1000;
	bench_t fl, fa, rl, ra, bl, ba, sl, sa, dl, da, ddl, dda, wa, wat;
	volatile int64_t sfl = 0, sfa = 0, srl = 0, sra = 0;
	volatile int64_t sbl = 0, sba = 0, ssl = 0, ssa = 0;
	volatile int64_t sdl = 0, sda = 0, sddl = 0, sdda = 0;
	volatile int64_t swa = 0, swat = 0;
	struct timeval tv;

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, m, x);

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
	printf("%p, %p\n", a, b);

	// printf("%s\n%s\n", a, b);

	printf("len:\t%d\ncnt:\t%d\n", len, cnt);
	printf("m: %d\tx: %d\tgi: %d\tge: %d\n", m, x, gi, ge);
	printf("alg\tlinear\taffine\tsc(l)\tsc(a)\n");

#ifdef FULL
	/* full */
	bench_init(fl);
	bench_init(fa);
	bench_start(fl);
	for(i = 0; i < cnt; i++) {
		sw_result_t r = sw_linear(a, strlen(a), b, strlen(b), score_matrix, ge);
		sfl += r.score;
		free(r.path);
	}
	bench_end(fl);
	bench_start(fa);
	for(i = 0; i < cnt; i++) {
		sw_result_t r = sw_affine(a, strlen(a), b, strlen(b), score_matrix, gi, ge);
		sfa += r.score;
		free(r.path);
	}
	bench_end(fa);
	printf("full:\t%lld\t%lld\t%lld\t%lld\n",
		bench_get(fl) / 1000,
		bench_get(fa) / 1000,
		sfl, sfa);
#endif

#ifdef ALL
	/* rognes */
	bench_init(rl);
	bench_init(ra);
	bench_start(rl);
	for(i = 0; i < cnt; i++) {
		srl += rognes_linear(a, strlen(a), b, strlen(b), score_matrix, ge, xt);
	}
	bench_end(rl);
	bench_start(ra);
	for(i = 0; i < cnt; i++) {
		sra += rognes_affine(a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
	}
	bench_end(ra);
	printf("rognes:\t%lld\t%lld\t%lld\t%lld\n",
		bench_get(rl) / 1000,
		bench_get(ra) / 1000,
		srl, sra);
#endif

	/* blast */
	bench_init(bl);
	bench_init(ba);
	bench_start(bl);
	for(i = 0; i < cnt; i++) {
		sbl += blast_linear(a, strlen(a), b, strlen(b), score_matrix, ge, xt);
	}
	bench_end(bl);
	bench_start(ba);
	for(i = 0; i < cnt; i++) {
		sba += blast_affine(a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
	}
	bench_end(ba);
	printf("blast:\t%lld\t%lld\t%lld\t%lld\n",
		bench_get(bl) / 1000,
		bench_get(ba) / 1000,
		sbl, sba);

#ifdef ALL
	/* simdblast */
	bench_init(sl);
	bench_init(sa);
	bench_start(sl);
	for(i = 0; i < cnt; i++) {
		ssl += simdblast_linear(a, strlen(a), b, strlen(b), score_matrix, ge, xt);
	}
	bench_end(sl);
	bench_start(sa);
	for(i = 0; i < cnt; i++) {
		ssa += simdblast_affine(a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
	}
	bench_end(sa);
	printf("simd:\t%lld\t%lld\t%lld\t%lld\n",
		bench_get(sl) / 1000,
		bench_get(sa) / 1000,
		ssl, ssa);

	/* diag */
	bench_init(dl);
	bench_init(da);
	bench_start(dl);
	for(i = 0; i < cnt; i++) {
		sdl += diag_linear(a, strlen(a), b, strlen(b), score_matrix, ge, xt);
	}
	bench_end(dl);
	bench_start(da);
	for(i = 0; i < cnt; i++) {
		sda += diag_affine(a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
	}
	bench_end(da);
	printf("diag:\t%lld\t%lld\t%lld\t%lld\n",
		bench_get(dl) / 1000,
		bench_get(da) / 1000,
		sdl, sda);
#endif

	/* dynamic diag */
	bench_init(ddl);
	bench_init(dda);
	bench_start(ddl);
	for(i = 0; i < cnt; i++) {
		sddl += ddiag_linear(a, strlen(a), b, strlen(b), score_matrix, ge, xt);
	}
	bench_end(ddl);
	bench_start(dda);
	for(i = 0; i < cnt; i++) {
		sdda += ddiag_affine(a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
	}
	bench_end(dda);
	printf("ddiag:\t%lld\t%lld\t%lld\t%lld\n",
		bench_get(ddl) / 1000,
		bench_get(dda) / 1000,
		sddl, sdda);

#ifdef SSW
	/* SSW library */
	/* convert sequence to number string */
	int8_t *na = (int8_t *)malloc(strlen(a));
	int8_t *nb = (int8_t *)malloc(strlen(b));
	for(i = 0; i < (int)strlen(a); i++) {
		na[i] = encode(a[i]);
	}
	for(i = 0; i < (int)strlen(b); i++) {
		nb[i] = encode(b[i]);
	}

	bench_init(wa);
	bench_init(wat);
	s_profile *sp = ssw_init(na, strlen(a), score_matrix, 4, 1);
	bench_start(wa);
	for(i = 0; i < cnt; i++) {
		s_align *r = ssw_align(sp, nb, strlen(b), -gi-ge, -ge, 0, 0, 0, 30);
		swa += r->score1;
		align_destroy(r);
	}
	bench_end(wa);
	init_destroy(sp);
	bench_start(wat);
	for(i = 0; i < cnt; i++) {
		s_profile *sp = ssw_init(na, strlen(a), score_matrix, 4, 1);
		s_align *r = ssw_align(sp, nb, strlen(b), -gi-ge, -ge, 0, 0, 0, 30);
		swat += r->score1;
		align_destroy(r);
		init_destroy(sp);
	}
	bench_end(wat);
	printf("ssw (w/o prep):\t%lld\t\t%lld\n",
		bench_get(wa) / 1000,
		swa);
	printf("ssw (w/ prep):\t%lld\t\t%lld\n",
		bench_get(wat) / 1000,
		swat);
	free(na);
	free(nb);
#endif

	free(a);
	free(b);
	return 0;
}

/**
 * end of main.cc
 */
