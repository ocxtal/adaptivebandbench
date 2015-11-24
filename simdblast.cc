
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

#define BW		( 128 )
#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn simdblast_linear
 */
int
simdblast_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	uint64_t i, a_size, first_a_index, last_a_index, a_index, b_index;
	vec const mv(m), xv(x), gv(-gi), gv2(-2*gi), gv4(-4*gi), xtv(xt), zv(0), ofsv(OFS);
	int16_t const acc_g[vec::LEN] __attribute__(( aligned(16) )) = {
		0,
		-gi,
		-2*gi,
		-3*gi,
		-4*gi,
		-5*gi,
		-6*gi,
		-7*gi,
	};

	int16_t *mat = (int16_t *)aligned_malloc(
		(roundup((alen + 1), vec::LEN) + 2*vec::LEN) * (blen + 1) * sizeof(int16_t),
		vec::SIZE);
	debug("%llu, %llu, %llu, %llu", alen, roundup((alen + 2), vec::LEN), (blen + 1), roundup((alen + 2), vec::LEN) * (blen + 1) * sizeof(int16_t));
	int16_t *ptr = mat, *prev;

	/* init top row */
	int16_t const init[vec::LEN] __attribute__(( aligned(16) )) = {
		OFS,
		OFS + gi,
		OFS + 2*gi,
		OFS + 3*gi,
		OFS + 4*gi,
		OFS + 5*gi,
		OFS + 6*gi,
		OFS + 7*gi,
	};
	vec init_v; init_v.load(init);
	for(i = 0; i < roundup(alen, vec::LEN); i += vec::LEN) {
		if((ofsv - init_v > xtv) == 0xffff) { break; }
		init_v.store(&ptr[i]);
		init_v -= (gv4 + gv4);
	}

	a_size = last_a_index = i;
	first_a_index = 0;
	vec best_score_v = ofsv;
	vec next_score_v(gv4);

	for(b_index = 1; b_index <= blen; b_index++) {
		vec bv(b[b_index-1]);

		prev = ptr; ptr += (last_a_index + vec::LEN - first_a_index);
		last_a_index = first_a_index;

		vec score_v = zv;
		vec prev_pv = gv4;
		vec prev_av = zv;
		vec best_temp_v = best_score_v;
		if(first_a_index != 0) { prev_pv = next_score_v; }

		debug("loop b(%llu), a_size(%llu)", b_index, a_size);
		for(a_index = first_a_index; a_index < a_size; a_index += vec::LEN) {

			/* load prev vector */
			vec pv; pv.load(&prev[a_index]);
			pv.print();

			/* calc f */
			vec score_gap_col_v = pv - gv;

			/* calc e */
			vec score_gap_row_v = (score_v - gv)>>7;		/* vec::LEN == 8 */

			/* calc d */
			score_v = prev_pv>>7 | pv<<1;
			vec av; av.load_expand(&a[a_index]);
			score_v += vec::comp(prev_av>>7 | av<<1, bv).select(mv, xv);
			score_v = vec::max(score_v, score_gap_col_v);

			/* chain e */
			score_v = vec::max(score_v, score_gap_row_v);
			score_v = vec::max(score_v, (score_v - gv)<<1);
			score_v = vec::max(score_v, (score_v - gv2)<<2);
			score_v = vec::max(score_v, (score_v - gv4)<<4);

			score_v.print();

			/* xdrop test */
			if((best_score_v - score_v > xtv) == 0xffff) {
				vec const xxv(-x);
				xxv.store(&ptr[a_index]);
				if(a_index == first_a_index) {
					next_score_v = score_v;
					first_a_index += vec::LEN;
				}
			} else {
				last_a_index = a_index;
				score_v.store(&ptr[a_index]);
				best_temp_v = vec::max(best_temp_v, score_v);
			}

			prev_pv = pv;
			prev_av = av;
		}

		debug("xdrop, first_a_index(%llu), last_a_index(%llu), a_size(%llu)", first_a_index, last_a_index, a_size);
		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - vec::LEN) {
			a_size = last_a_index + vec::LEN;
		} else {
			best_score_v.print();
			score_v.print();
			vec acc_gv; acc_gv.load(acc_g);
// #if 0
			while((best_score_v - score_v <= xtv) && a_size <= roundup(alen + 1, vec::LEN)) {
				score_v = (score_v>>7) - gv;		/* vec::LEN == 8 */
				score_v.set(score_v[0]);
				score_v -= acc_gv;
/*
				score_v = vec::max(score_v, (score_v - gv)<<1);
				score_v = vec::max(score_v, (score_v - gv2)<<2);
				score_v = vec::max(score_v, (score_v - gv4)<<4);*/
				score_v.store(&ptr[a_size]);
				a_size += vec::LEN;
				score_v.print();
			}
// #endif
			// printf("a(%d), b(%d)\n", a_size, b_index);
			last_a_index = a_size;
		}
/*
		if(a_size <= alen) {
			zv.store(&ptr[a_size]); a_size += vec::LEN;
		}
*/
		best_score_v.set(best_temp_v.hmax());
	}
	free(mat);

	return(best_score_v.hmax() - OFS);
}

/**
 * @fn simdblast_affine
 */
int
simdblast_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	uint64_t i, a_size, first_a_index, last_a_index, a_index, b_index;
	/* gvはaffine gap costにあわせて書き換える */
	vec const mv(m), xv(x), xtv(xt), zv(0), ofsv(OFS);
	vec const giv(-gi), gev(-ge), giev(ge-gi), gev3(-3*ge), gev8(-8*ge);
	int16_t const acc_ge[vec::LEN] __attribute__(( aligned(16) )) = {
		0,
		-ge,
		-2*ge,
		-3*ge,
		-4*ge,
		-5*ge,
		-6*ge,
		-7*ge
	};

	struct _dp {
		int16_t best[vec::LEN];
		int16_t best_gap[vec::LEN];
	};
	struct _dp *mat = (struct _dp *)aligned_malloc(
		(roundup(alen + 1, vec::LEN) / vec::LEN + 2) * (blen + 1) * sizeof(struct _dp),
		vec::SIZE);
	struct _dp *ptr = mat, *prev = mat + roundup(alen, vec::LEN) / vec::LEN;
	// memset(mat, 0, roundup(alen, vec::LEN) / vec::LEN * roundup(blen, vec::LEN) * sizeof(struct _dp));

	debug("%llu, %llu, %llu, %llu, %lu", alen, roundup((alen + 1), vec::LEN), roundup((alen + 1), vec::LEN) / vec::LEN + 1, blen, sizeof(struct _dp));

	/* init top row */
	int16_t const init[vec::LEN] __attribute__(( aligned(16) )) = {
		OFS,
		OFS + gi,
		OFS + gi + ge,
		OFS + gi + 2*ge,
		OFS + gi + 3*ge,
		OFS + gi + 4*ge,
		OFS + gi + 5*ge,
		OFS + gi + 6*ge
	};
	int16_t const init_g[vec::LEN] __attribute__(( aligned(16) )) = {
		OFS + gi - ge,
		OFS + gi,
		OFS + gi + ge,
		OFS + gi + 2*ge,
		OFS + gi + 3*ge,
		OFS + gi + 4*ge,
		OFS + gi + 5*ge,
		OFS + gi + 6*ge,
	};
	vec init_v; init_v.load(init);
	vec init_gv; init_gv.load(init_g);
	for(i = 0; i < roundup(alen + 1, vec::LEN) / vec::LEN; i++) {
		if((ofsv - init_v > xtv) == 0xffff) { break; }
		init_v.store(&ptr[i].best);
		(init_v - giev).store(&ptr[i].best_gap);
		init_gv -= gev8; init_v = init_gv;
	}

	debug("%llu", i);

	a_size = last_a_index = i;
	first_a_index = 0;
	vec best_score_v = ofsv;
	vec next_score_v;

	for(b_index = 1; b_index <= blen; b_index++) {
		vec bv(b[b_index-1]);

		debug("%llu", last_a_index + 1 - first_a_index);
		prev = ptr; ptr += (last_a_index + 1 - first_a_index);
		last_a_index = first_a_index;

		vec score_v = zv;
		vec score_gap_row_v = zv;
		vec prev_pv = giv + gev3;
		vec prev_av = zv;
		vec best_temp_v = best_score_v;

		debug("loop b(%llu), first_a_index(%llu), a_size(%llu)", b_index, first_a_index, a_size);
		for(a_index = first_a_index; a_index < a_size; a_index++) {

			/* load prev vector */
			vec pv; pv.load(&prev[a_index].best);
			vec pf; pf.load(&prev[a_index].best_gap);
			pv.print();

			/* calc f */
			vec score_gap_col_v = vec::max(pv - giv, pf - gev);

			/* calc e */
			score_gap_row_v = vec::max(score_v - giv, score_gap_row_v - gev)>>7;

			/* calc d */
			score_v = prev_pv>>7 | pv<<1;
			// vec av('A');
			vec av; av.load_expand(&a[a_index * vec::LEN]);
			score_v += vec::comp(prev_av>>7 | av<<1, bv).select(mv, xv);
			score_v = vec::max(score_v, score_gap_col_v);

			/* chain e */
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = vec::max(score_v - giv, score_gap_row_v - gev)<<1;
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = vec::max(score_v - giv - gev, score_gap_row_v - gev - gev)<<2;
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = vec::max(score_v - giv - gev3, score_gap_row_v - gev - gev3)<<4;
			score_v = vec::max(score_v, score_gap_row_v); score_v.print();

			/* xdrop test */
			if((best_score_v - score_v > xtv) == 0xffff) {
				vec const xxv(-x);
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
			vec acc_gev; acc_gev.load(acc_ge);
			score_v.print(); score_gap_row_v.print();
			while((best_score_v - score_v <= xtv) && a_size <= roundup(alen + 1, vec::LEN) / vec::LEN) {
				score_v = vec::max(score_v - giv, score_gap_row_v - gev)>>7; score_v.print();
				score_v.set(score_v[0]); score_v.print();
				score_v -= acc_gev; score_v.print();
				score_gap_row_v = score_v - giev;
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
	free(mat);

	return(best_score_v.hmax() - OFS);
}

#ifdef MAIN
#include <stdlib.h>
int main_ext(int argc, char *argv[])
{
	if(strcmp(argv[1], "linear") == 0) {
		int score = simdblast_linear(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]),
			atoi(argv[8]));
		printf("%d\n", score);
	} else if(strcmp(argv[1], "affine") == 0) {
		int score = simdblast_affine(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]),
			atoi(argv[8]));
		printf("%d\n", score);
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

	int sl = simdblast_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1, 30);
	printf("%d\n", sl);

	int sa = simdblast_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1, 20);
	printf("%d\n", sa);

	return(0);
}
#endif

/**
 * end of simdblast.cc
 */
