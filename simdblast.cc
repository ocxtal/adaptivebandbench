#define DEBUG
/**
 * @file simdblast.cc
 *
 * @brief simd parallel blast_SemiGappedAlign algorithm
 */
#include <string.h>
#include "sse.h"
/*
#define DEBUG
#include "log.h"
#undef DEBUG
*/
#include "util.h"

// #define BW		( 128 )
#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn simdblast_affine
 */
int
simdblast_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt,
	uint32_t bw)		/* unused */
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	uint64_t i, a_size, first_a_index, last_a_index, a_index, b_index;
	vec const xtv(xt), zv, ofsv(OFS);
	vec const giv(-gi), gev(-ge), gev2(-2*ge), gev4(-4*ge), gev8(-8*ge);
	int16_t const acc_ge[vec::LEN] __attribute__(( aligned(16) )) = {
		(int16_t)(0),
		(int16_t)(-ge),
		(int16_t)(-2*ge),
		(int16_t)(-3*ge),
		(int16_t)(-4*ge),
		(int16_t)(-5*ge),
		(int16_t)(-6*ge),
		(int16_t)(-7*ge)
	};
	vec acc_gev; acc_gev.load(acc_ge);

	int8_t sc_min = extract_min_score(score_matrix);

	struct _dp {
		int16_t best[vec::LEN];
		int16_t best_gap[vec::LEN];
	};
	struct _dp *ptr = (struct _dp *)work, *prev = (struct _dp *)work + roundup(alen, vec::LEN) / vec::LEN;

	debug("%llu, %llu, %llu, %llu, %lu", alen, roundup((alen + 1), vec::LEN), roundup((alen + 1), vec::LEN) / vec::LEN + 1, blen, sizeof(struct _dp));

	/* init top row */
	vec init_v = ofsv - (giv<<1) - acc_gev;
	init_v.store(&ptr[0].best);
	(init_v - giv).store(&ptr[0].best_gap);
	init_v = ofsv - giv - gev8 - acc_gev;
	for(i = 1; i < roundup(alen + 1, vec::LEN) / vec::LEN; i++) {
		if((ofsv - init_v > xtv) == 0xffff) { break; }
		init_v.store(&ptr[i].best);
		(init_v - giv).store(&ptr[i].best_gap);
		init_v -= gev8;
	}

	debug("%llu", i);

	a_size = last_a_index = i;
	first_a_index = 0;
	vec best_score_v = ofsv;
	vec next_score_v;
	vec smv; smv.load(score_matrix);

	for(b_index = 1; b_index <= blen; b_index++) {
		char_vec bv(encode_b(b[b_index-1]));

		debug("%llu", last_a_index + 1 - first_a_index);
		prev = ptr; ptr += (last_a_index + 1 - first_a_index);
		last_a_index = first_a_index;

		vec score_v = zv;
		vec score_gap_row_v = zv;
		vec prev_pv = giv + gev4;
		char_vec prev_av((uint64_t)0);
		vec best_temp_v = best_score_v;

		debug("loop b(%llu), first_a_index(%llu), a_size(%llu)", b_index, first_a_index, a_size);
		for(a_index = first_a_index; a_index < a_size; a_index++) {

			/* load prev vector */
			vec pv; pv.load(&prev[a_index].best);
			vec pf; pf.load(&prev[a_index].best_gap);
			pv.print("pv");

			/* calc f */
			vec score_gap_col_v = vec::max(pv - giv, pf) - gev;

			/* calc e */
			score_gap_row_v = (vec::max(score_v - giv, score_gap_row_v) - gev)>>7;

			/* calc d */
			score_v = pv.dsl(prev_pv);
			char_vec av; av.load_encode_a(&a[a_index * vec::LEN], a_size - a_index * vec::LEN);
			score_v += smv.shuffle(av.dsl(prev_av) | bv);
			score_v = vec::max(score_v, score_gap_col_v);

			/* chain e */
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = (vec::max(score_v - giv, score_gap_row_v) - gev)<<1;
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = (vec::max(score_v - giv, score_gap_row_v) - gev2)<<2;
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = (vec::max(score_v - giv, score_gap_row_v) - gev4)<<4;
			score_v = vec::max(score_v, score_gap_row_v);

			/* xdrop test */
			if((best_score_v - score_v > xtv) == 0xffff) {
				vec const xxv(-sc_min);
				xxv.store(&ptr[a_index].best);
				xxv.store(&ptr[a_index].best_gap);
				if(a_index == first_a_index) {
					next_score_v = score_v;
					first_a_index++;
				}
			} else {
				last_a_index = a_index;
				score_v.store(&ptr[a_index].best);
				score_gap_col_v.store(&ptr[a_index].best_gap);
				best_temp_v = vec::max(best_temp_v, score_v);
			}

			prev_pv = pv;
			prev_av = av;
		}

		debug("xdrop, first_a_index(%llu), last_a_index(%llu), a_size(%llu)", first_a_index, last_a_index, a_size);
		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - 1) {
			a_size = last_a_index + 1;
		} else {
			score_v.print("pv (raw)"); score_gap_row_v.print("gap_v (raw)");
			while((best_score_v - score_v <= xtv) && a_size <= roundup(alen + 1, vec::LEN) / vec::LEN) {
				score_v = (vec::max(score_v - giv, score_gap_row_v) - gev)>>7;
				score_v.set(score_v[0]);
				score_v -= acc_gev;
				score_gap_row_v = score_v - giv;
				score_v.store(&ptr[a_size].best);
				score_gap_row_v.store(&ptr[a_size].best_gap);
				a_size++;
			}
			// printf("a(%d), b(%d)\n", a_size, b_index);
			last_a_index = a_size;
		}
/*
		if(a_size <= alen) {
			zv.store(ptr[a_size].best);
			zv.store(ptr[a_size].best_gap);
			a_size += vec::LEN;
		}
*/
		best_score_v.set(best_temp_v.hmax());
	}
	// free(mat);

	return(best_score_v.hmax() - OFS);
}

#ifdef MAIN
#include <assert.h>
#include <stdlib.h>
int main_ext(int argc, char *argv[])
{
	uint64_t alen = strlen(argv[2]);
	uint64_t blen = strlen(argv[3]);
	char *a = (char *)malloc(alen + vec::LEN + 1);
	char *b = (char *)malloc(blen + vec::LEN + 1);

	memcpy(a, argv[2], alen);
	memset(a + alen, 0, vec::LEN + 1);
	memcpy(b, argv[3], blen);
	memset(b + blen, 0x80, vec::LEN + 1);

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, atoi(argv[4]), atoi(argv[5]));


	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	if(0) {
		printf("./a.out AAA AAA 2 -3 -5 -1 30\n");
	}

	int score = simdblast_affine(
		work,
		a, alen, b, blen,
		score_matrix,
		atoi(argv[6]),
		atoi(argv[7]),
		atoi(argv[8]));
	printf("%d\n", score);
	free(a); free(b); free(work);
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


	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	#define a(s, p, q) { \
		assert(simdblast_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10) == (s)); \
	}
	a( 0, "", "");
	a( 0, "A", "");
	a( 1, "A", "A");
	a( 3, "AAA", "AAA");
	a( 0, "AAA", "TTT");
	a( 3, "AAAGGG", "AAATTTTTT");
	a( 3, "TTTGGGGGAAAA", "TTTCCCCCCCCAAAA");
	a( 4, "AAACAAAGGG", "AAAAAATTTTTTT");
	a( 3, "AAACCAAAGGG", "AAAAAATTTTTTT");

	int sa = simdblast_affine(work, a, strlen(a), b, strlen(b), score_matrix, -1, -1, 30);
	printf("%d\n", sa);

	free(work);
	return(0);
}
#endif

/**
 * end of simdblast.cc
 */
