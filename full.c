
/**
 * @file full.c
 *
 * @brief full Smith-Waterman
 */
#include "full.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define _a(p, q, plen)		( (q) * ((plen) + 1) + (p) )
#define MAX2(p, q)			( ((p) < (q)) ? (q) : (p) )
#define MAX3(p, q, r)		( MAX2(p, MAX2(q, r)) )
#define MAX4(p, q, r, s)	( MAX2(MAX2(p, q), MAX2(r, s)) )

typedef struct _maxpos {
	int16_t score;
	uint64_t apos;
	uint64_t bpos;
} maxpos_t;

/**
 * @fn sw_linear
 */
sw_result_t sw_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge)
{
	/* util macros */
	#define a(p, q)		_a(p, q, alen)
	#define s(p, q)		( (a[(p) - 1] == b[(q) - 1]) ? m : x )

	int16_t *mat = (int16_t *)malloc(
		(alen + 1) * (blen + 1) * sizeof(int16_t));

	/* init */
	maxpos_t max = { 0, 0, 0 };
	mat[a(0, 0)] = 0;
	for(uint64_t i = 1; i < alen+1; i++) { mat[a(i, 0)] = i * ge; }
	for(uint64_t j = 1; j < blen+1; j++) { mat[a(0, j)] = j * ge; }

	for(uint64_t j = 1; j < blen+1; j++) {
		for(uint64_t i = 1; i < alen+1; i++) {
			int16_t score = mat[a(i, j)] = MAX3(
				mat[a(i - 1, j - 1)] + s(i, j),
				mat[a(i, j - 1)] + ge,
				mat[a(i - 1, j)] + ge);
			if(score >= max.score) { max = (maxpos_t){ score, i, j }; }
		}
	}
	if(max.score == 0) { max = (maxpos_t){ 0, 0, 0 }; }

	sw_result_t result;
	result = (sw_result_t){
		.score = max.score,
		.apos = max.apos,
		.bpos = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + 1)
	};
	uint32_t path_index = max.apos + max.bpos + 1;
	while(max.apos != 0 && max.bpos != 0) {
		if(mat[a(max.apos, max.bpos)] == mat[a(max.apos, max.bpos - 1)] + ge) {
			max.bpos--;
			result.path[--path_index] = 'I';
		} else if(mat[a(max.apos, max.bpos)] == mat[a(max.apos - 1, max.bpos)] + ge) {
			max.apos--;
			result.path[--path_index] = 'D';
		} else {
			if(mat[a(max.apos, max.bpos)] == mat[a(max.apos - 1, max.bpos - 1)] + x) {
				result.path[--path_index] = 'X';
			} else {
				result.path[--path_index] = 'M';
			}
			max.apos--;
			max.bpos--;
		}
	}

	result.path_length -= path_index;
	for(uint64_t i = 0; i < result.path_length; i++) {
		result.path[i] = result.path[path_index++];
	}
	result.path[result.path_length] = '\0';

	free(mat);

	#undef a
	#undef s
	return(result);
}

/**
 * @fn sw_affine
 */
sw_result_t sw_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge)
{
	/* utils */
	#define a(p, q)		_a(p, 3*(q), alen)
	#define f(p, q)		_a(p, 3*(q)+1, alen)
	#define e(p, q)		_a(p, 3*(q)+2, alen)
	#define s(p, q)		( (a[(p) - 1] == b[(q) - 1]) ? m : x )

	int16_t *mat = (int16_t *)malloc(
		3 * (alen + 1) * (blen + 1) * sizeof(int16_t));

	/* init */
	maxpos_t max = { 0, 0, 0 };
	mat[a(0, 0)] = mat[f(0, 0)] = mat[e(0, 0)] = 0;
	for(uint64_t i = 1; i < alen+1; i++) {
		mat[a(i, 0)] = mat[f(i, 0)] = gi + (i - 1) * ge;
		mat[e(i, 0)] = gi + (i - 1) * ge + gi - ge;
	}
	for(uint64_t j = 1; j < blen+1; j++) {
		mat[a(0, j)] = mat[e(0, j)] = gi + (j - 1) * ge;
		mat[f(0, j)] = gi + (j - 1) * ge + gi - ge;
	}

	for(uint64_t j = 1; j < blen+1; j++) {
		for(uint64_t i = 1; i < alen+1; i++) {
			int16_t score_f = mat[f(i, j)] = MAX2(
				mat[a(i - 1, j)] + gi,
				mat[f(i - 1, j)] + ge);
			int16_t score_e = mat[e(i, j)] = MAX2(
				mat[a(i, j - 1)] + gi,
				mat[e(i, j - 1)] + ge);
			int16_t score = mat[a(i, j)] = MAX3(
				mat[a(i - 1, j - 1)] + s(i, j),
				score_f, score_e);
			if(score >= max.score) { max = (maxpos_t){ score, i, j }; }
		}
	}
	if(max.score == 0) { max = (maxpos_t){ 0, 0, 0 }; }

	sw_result_t result;
	result = (sw_result_t){
		.score = max.score,
		.apos = max.apos,
		.bpos = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + 1)
	};
	uint32_t path_index = max.apos + max.bpos + 1;
	while(max.apos != 0 && max.bpos != 0) {
		if(mat[a(max.apos, max.bpos)] == mat[e(max.apos, max.bpos)]) {
			while(mat[e(max.apos, max.bpos)] == mat[e(max.apos, max.bpos - 1)] + ge) {
				max.bpos--;
				result.path[--path_index] = 'I';
			}
			max.bpos--;
			result.path[--path_index] = 'I';
		} else if(mat[a(max.apos, max.bpos)] == mat[f(max.apos, max.bpos)]) {
			while(mat[f(max.apos, max.bpos)] == mat[f(max.apos - 1, max.bpos)] + ge) {
				max.apos--;
				result.path[--path_index] = 'D';
			}
			max.apos--;
			result.path[--path_index] = 'D';
		} else {
			if(mat[a(max.apos, max.bpos)] == mat[a(max.apos - 1, max.bpos - 1)] + x) {
				result.path[--path_index] = 'X';
			} else {
				result.path[--path_index] = 'M';
			}
			max.apos--;
			max.bpos--;
		}
	}

	result.path_length -= path_index;
	for(uint64_t i = 0; i < result.path_length; i++) {
		result.path[i] = result.path[path_index++];
	}
	result.path[result.path_length] = '\0';

	free(mat);

	#undef a
	#undef f
	#undef e
	#undef s
	return(result);
}

#ifdef TEST
#include <stdio.h>
#include <assert.h>
#define _f(f, a, b, s)		f(a, strlen(a), b, strlen(b), s[0], s[1], s[2], s[3])
#define _linear(a, b, s)	_f(sw_linear, a, b, s)
#define _affine(a, b, s)	_f(sw_affine, a, b, s)
static int8_t const t[][5] = {
	{1, -1, -1, -1},
	{1, -1, -2, -1},
	{1, -3, -5, -1},
	{2, -3, -5, -1}
};

void test_linear_1_1_1(void)
{
	#define l(s, ap, bp, aln, p, q) { \
		sw_result_t r = _linear(p, q, t[0]); \
		assert(r.score == (s)); \
		assert(r.apos == (ap)); \
		assert(r.bpos == (bp)); \
		assert(r.path_length == strlen(aln)); \
		assert(strcmp(r.path, aln) == 0); \
		if(r.path != NULL) { free(r.path); } \
	}
	l( 0,  0,  0, "", "", "");
	l( 0,  0,  0, "", "A", "");
	l( 1,  1,  1, "M", "A", "A");
	l( 3,  3,  3, "MMM", "AAA", "AAA");
	l( 0,  0,  0, "", "AAA", "TTT");
	l( 3,  6,  9, "MMM", "GGGAAAGGG", "TTTTTTAAATTTTTT");

	l( 5, 10, 13, "MMMDMMM", "GGGAAACAAAGGG", "TTTTTTTAAAAAATTTTTTT");
	l( 4, 11, 13, "MMMDDMMM", "GGGAAACCAAAGGG", "TTTTTTTAAAAAATTTTTTT");
	#undef l
}

void test_affine_1_1_1(void)
{
	#define a(s, ap, bp, aln, p, q) { \
		sw_result_t r = _affine(p, q, t[0]); \
		assert(r.score == (s)); \
		assert(r.apos == (ap)); \
		assert(r.bpos == (bp)); \
		assert(r.path_length == strlen(aln)); \
		assert(strcmp(r.path, aln) == 0); \
		if(r.path != NULL) { free(r.path); } \
	}
	a( 0,  0,  0, "", "", "");
	a( 0,  0,  0, "", "A", "");
	a( 1,  1,  1, "M", "A", "A");
	a( 3,  3,  3, "MMM", "AAA", "AAA");
	a( 0,  0,  0, "", "AAA", "TTT");
	a( 3,  6,  9, "MMM", "GGGAAAGGG", "TTTTTTAAATTTTTT");

	a( 4, 10, 13, "MMMDMMM", "GGGAAACAAAGGG", "TTTTTTTAAAAAATTTTTTT");
	a( 3, 11, 13, "MMMDDMMM", "GGGAAACCAAAGGG", "TTTTTTTAAAAAATTTTTTT");
	#undef a
}

int main(void)
{
	test_linear_1_1_1();
	test_affine_1_1_1();
	return(0);
}
#endif


#ifdef MAIN
#include <stdio.h>
int main_ext(int argc, char *argv[])
{
	if(strcmp(argv[1], "linear") == 0) {
		sw_result_t result = sw_linear(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]));
		printf("%d\t%llu\t%llu\t%s\n",
			result.score,
			result.apos,
			result.bpos,
			result.path);
		free(result.path);
	} else if(strcmp(argv[1], "affine") == 0) {
		sw_result_t result = sw_affine(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]));
		printf("%d\t%llu\t%llu\t%s\n",
			result.score,
			result.apos,
			result.bpos,
			result.path);
		free(result.path);
	} else {
		printf("./a.out linear AAA AAA 2 -3 -5 -1 30\n");
	}
	return(0);
}

int main(int argc, char *argv[])
{
	char const *a = "aabbcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	if(argc > 1) { return(main_ext(argc, argv)); }

	sw_result_t sl = sw_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	printf("%d\t%llu\t%llu\t%s\n",
		sl.score,
		sl.apos,
		sl.bpos,
		sl.path);
	free(sl.path);

	sw_result_t sa = sw_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	printf("%d\t%llu\t%llu\t%s\n",
		sa.score,
		sa.apos,
		sa.bpos,
		sa.path);
	free(sa.path);

	return(0);
}
#endif

/**
 * end of full.c
 */
