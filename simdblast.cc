
/**
 * @file simdblast.cc
 *
 * @brief simd parallel blast_SemiGappedAlign algorithm
 */
#include <string.h>
#include "sse.h"
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
	int i, a_size, first_a_index, last_a_index, a_index, b_index;
	vec const mv(m), xv(x), gv(-gi), gv2(-2*gi), gv4(-4*gi), xtv(xt), zv(0), ofsv(OFS);
	int16_t acc_g[vec::LEN] __attribute__(( aligned(16) )) = {
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
		roundup(alen, vec::LEN) * roundup(blen, vec::LEN) * sizeof(int16_t),
		vec::SIZE);
	int16_t *ptr = mat, *prev;

	/* init top row */
	int16_t init[vec::LEN] __attribute__(( aligned(16) )) = {
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
	vec next_score_v;

	for(b_index = 1; b_index <= blen; b_index++) {
		vec bv(b[b_index-1]);

		prev = ptr; ptr += (last_a_index + vec::LEN - first_a_index);
		last_a_index = first_a_index;

		vec score_v = zv;
		vec prev_pv = gv4;
		vec prev_av = zv;
		vec best_temp_v = best_score_v;
		if(first_a_index != 0) { prev_pv = next_score_v; }

		debug("loop b(%d), a_size(%d)", b_index, a_size);
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
				zv.store(&ptr[a_index]);
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

		debug("xdrop, first_a_index(%d), last_a_index(%d), a_size(%d)", first_a_index, last_a_index, a_size);
		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - vec::LEN) {
			a_size = last_a_index + vec::LEN;
		} else {
			best_score_v.print();
			score_v.print();
			vec acc_gv; acc_gv.load(acc_g);
			while((best_score_v - score_v <= xtv) && a_size <= roundup(alen, vec::LEN)) {
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
	int i, a_size, first_a_index, last_a_index, a_index, b_index;
	/* gvはaffine gap costにあわせて書き換える */
	vec const mv(m), xv(x), xtv(xt), zv(0), ofsv(OFS);
	vec const giv(-gi), gev(-ge), giev(ge-gi), gev3(-3*ge), gev8(-8*ge);
	int16_t acc_ge[vec::LEN] __attribute__(( aligned(16) )) = {
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
		roundup(alen, vec::LEN) / vec::LEN * roundup(blen, vec::LEN) * sizeof(struct _dp),
		vec::SIZE);
	struct _dp *ptr = mat, *prev;

	/* init top row */
	int16_t init[vec::LEN] __attribute__(( aligned(16) )) = {
		OFS,
		OFS + gi,
		OFS + gi + ge,
		OFS + gi + 2*ge,
		OFS + gi + 3*ge,
		OFS + gi + 4*ge,
		OFS + gi + 5*ge,
		OFS + gi + 6*ge
	};
	int16_t init_g[vec::LEN] __attribute__(( aligned(16) )) = {
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
	for(i = 0; i < roundup(alen, vec::LEN); i += vec::LEN) {
		if((ofsv - init_v > xtv) == 0xffff) { break; }
		init_v.store(ptr[i].best);
		(init_v - giev).store(ptr[i].best_gap);
		init_gv -= gev8; init_v = init_gv;
	}

	a_size = last_a_index = i;
	first_a_index = 0;
	vec best_score_v = ofsv;
	vec next_score_v;

	for(b_index = 1; b_index <= blen; b_index++) {
		vec bv(b[b_index-1]);

		prev = ptr; ptr += (last_a_index + 1 - first_a_index);
		last_a_index = first_a_index;

		vec score_v = zv;
		vec score_gap_row_v = zv;
		vec prev_pv = giv + gev3;
		vec prev_av = zv;
		vec best_temp_v = best_score_v;

		debug("loop b(%d), a_size(%d)", b_index, a_size);
		for(a_index = first_a_index; a_index < a_size; a_index += vec::LEN) {

			/* load prev vector */
			vec pv; pv.load(prev[a_index].best);
			vec pf; pf.load(prev[a_index].best_gap);
			pv.print();

			/* calc f */
			vec score_gap_col_v = vec::max(pv - giv, pf - gev);

			/* calc e */
			score_gap_row_v = vec::max(score_v - giv, score_gap_row_v - gev)>>7;

			/* calc d */
			score_v = prev_pv>>7 | pv<<1;
			vec av; av.load_expand(&a[a_index]);
			score_v += vec::comp(prev_av>>7 | av<<1, bv).select(mv, xv);
			score_v = vec::max(score_v, score_gap_col_v);

			/* chain e */
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = vec::max(score_v - giv, score_gap_row_v - gev)<<1;
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = vec::max(score_v - giv - gev, score_gap_row_v - gev - gev)<<2;
			score_v = vec::max(score_v, score_gap_row_v);
			score_gap_row_v = vec::max(score_v - giv - gev3, score_gap_row_v - gev - gev3)<<4;
			score_v = vec::max(score_v, score_gap_row_v);

			/* xdrop test */
			if((best_score_v - score_v > xtv) == 0xffff) {
				zv.store(ptr[a_index].best);
				zv.store(ptr[a_index].best_gap);
				if(a_index == first_a_index) {
					next_score_v = score_v;
					first_a_index += vec::LEN;
				}
			} else {
				last_a_index = a_index;
				score_v.store(ptr[a_index].best);
				score_gap_col_v.store(ptr[a_index].best_gap);
				best_temp_v = vec::max(best_temp_v, score_v);
			}

			prev_pv = pv;
			prev_av = av;
		}

		debug("xdrop, first_a_index(%d), last_a_index(%d), a_size(%d)", first_a_index, last_a_index, a_size);
		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - vec::LEN) {
			a_size = last_a_index + vec::LEN;
		} else {
			vec acc_gev; acc_gev.load(acc_ge);
			score_v.print(); score_gap_row_v.print();
			while((best_score_v - score_v <= xtv) && a_size <= roundup(alen, vec::LEN)) {
				score_v = vec::max(score_v - giv, score_gap_row_v - gev)>>7; score_v.print();
				score_v.set(score_v[0]); score_v.print();
				score_v -= acc_gev; score_v.print();
				score_gap_row_v = score_v - giev;
				score_v.store(ptr[a_size].best);
				score_gap_row_v.store(ptr[a_size].best_gap);
				a_size += vec::LEN;
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
int main(void)
{
	char const *a = "aabbcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

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
