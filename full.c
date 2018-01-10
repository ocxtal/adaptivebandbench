
/**
 * @file full.c
 *
 * @brief full Smith-Waterman
 */
#include "full.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define _a(p, q, plen)		( (q) * ((plen) + 1) + (p) )
typedef struct { int64_t score; uint64_t apos, bpos; } sw_maxpos_t;

/**
 * @fn sw_affine
 */
sw_result_t sw_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge)
{
	/* utils */
	#define a(p, q)		_a(p, 3*(q), alen)
	#define f(p, q)		_a(p, 3*(q)+1, alen)
	#define e(p, q)		_a(p, 3*(q)+2, alen)
	// #define s(p, q)		( (a[(p) - 1] == b[(q) - 1]) ? m : x )
	#define s(p, q)		( score_matrix[encode_a(a[(p) - 1]) | encode_b(b[(q) - 1])] )

	int16_t const min = INT16_MIN - extract_min_score(score_matrix) - gi;

	int16_t *mat = (int16_t *)malloc(
		3 * (alen + 1) * (blen + 1) * sizeof(int16_t));

	/* fix gi */
	gi += ge;

	/* init */
	sw_maxpos_t max = { 0, 0, 0 };
	mat[a(0, 0)] = mat[f(0, 0)] = mat[e(0, 0)] = 0;
	for(uint64_t i = 1; i < alen+1; i++) {
		mat[a(i, 0)] = mat[f(i, 0)] = MAX2(min, gi + (i - 1) * ge);
		mat[e(i, 0)] = MAX2(min, gi + (i - 1) * ge + gi - ge - 1);
	}
	for(uint64_t j = 1; j < blen+1; j++) {
		mat[a(0, j)] = mat[e(0, j)] = MAX2(min, gi + (j - 1) * ge);
		mat[f(0, j)] = MAX2(min, gi + (j - 1) * ge + gi - ge - 1);
	}

	for(uint64_t j = 1; j < blen+1; j++) {
		for(uint64_t i = 1; i < alen+1; i++) {
			int16_t score_f = mat[f(i, j)] = MAX2(
				mat[a(i - 1, j)] + gi,
				mat[f(i - 1, j)] + ge);
			int16_t score_e = mat[e(i, j)] = MAX2(
				mat[a(i, j - 1)] + gi,
				mat[e(i, j - 1)] + ge);
			int16_t score = mat[a(i, j)] = MAX4(min,
				mat[a(i - 1, j - 1)] + s(i, j),
				score_f, score_e);
			if(score > max.score) { max = (sw_maxpos_t){ score, i, j }; }
		}
	}
	if(max.score == 0) { max = (sw_maxpos_t){ 0, 0, 0 }; }

	sw_result_t result;
	result = (sw_result_t){
		.score = max.score,
		.apos = max.apos,
		.bpos = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + 1)
	};
	uint32_t path_index = max.apos + max.bpos + 1;
	while(max.apos != 0 || max.bpos != 0) {
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
			// if(mat[a(max.apos, max.bpos)] == mat[a(max.apos - 1, max.bpos - 1)] + x) {
			if(a[max.apos - 1] != b[max.bpos - 1]) {
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
#define _linear(a, b, s)	sw_linear(a, strlen(a), b, strlen(b), score_matrix, s[3])
#define _affine(a, b, s)	sw_affine(a, strlen(a), b, strlen(b), score_matrix, s[2], s[3])
static int8_t const t[][5] = {
	{1, -1, -1, -1},
	{1, -1, -2, -1},
	{1, -3, -5, -1},
	{2, -3, -5, -1}
};

void test_linear_1_1_1(void)
{
	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, t[0][0], t[0][1]);

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
	l( 3,  3,  3, "MMM", "AAAGGG", "AAATTTTTT");

	l( 5,  7,  6, "MMMDMMM", "AAACAAAGGG", "AAAAAATTTTTTT");
	l( 4,  8,  6, "MMMDDMMM", "AAACCAAAGGG", "AAAAAATTTTTTT");
	#undef l
}

void test_affine_1_1_1(void)
{
	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, t[0][0], t[0][1]);

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
	a( 3,  3,  3, "MMM", "AAAGGG", "AAATTTTTT");

	a( 5,  8,  7, "MMMDMMMM", "AAACACGTGGG", "AAAACGTTTTTTTT");
	a( 4,  9,  7, "MMMDDMMMM", "AAACCACGTGGG", "AAAACGTTTTTTTT");
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
	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, atoi(argv[4]), atoi(argv[5]));

	if(strcmp(argv[1], "linear") == 0) {
		sw_result_t result = sw_linear(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			score_matrix,
			atoi(argv[5]));
		printf("%d\t%lu\t%lu\t%s\n",
			result.score,
			result.apos,
			result.bpos,
			result.path);
		free(result.path);
	} else if(strcmp(argv[1], "affine") == 0) {
		sw_result_t result = sw_affine(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			score_matrix,
			atoi(argv[5]),
			atoi(argv[6]));
		printf("%d\t%lu\t%lu\t%s\n",
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
	char const *a = "aattcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	if(argc > 1) { return(main_ext(argc, argv)); }

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, 1, -1);

	sw_result_t sl = sw_linear(a, strlen(a), b, strlen(b), score_matrix, -1);
	printf("%d\t%lu\t%lu\t%s\n",
		sl.score,
		sl.apos,
		sl.bpos,
		sl.path);
	free(sl.path);

	sw_result_t sa = sw_affine(a, strlen(a), b, strlen(b), score_matrix, -1, -1);
	printf("%d\t%lu\t%lu\t%s\n",
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
