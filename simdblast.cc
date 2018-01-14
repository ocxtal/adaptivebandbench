
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

	vec const xtv(xt), zv, ofsv(OFS), smv((uint16_t const *)score_matrix);
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
	vec acc_gev((uint16_t const *)acc_ge);
	uint64_t vblen = roundup(blen, vec::LEN) / vec::LEN;
	uint64_t first_b_index = 0, last_b_index = vblen;		/* [first_b_index, last_b_index) */

	struct _dp {
		uint16_t s[vec::LEN], e[vec::LEN], f[vec::LEN];
	};
	struct _dp *ptr = (struct _dp *)work + sizeof(maxpos_t), *prev;
	debug("%llu, %llu, %llu, %llu, %lu", alen, roundup((alen + 1), vec::LEN), roundup((alen + 1), vec::LEN) / vec::LEN + 1, blen, sizeof(struct _dp));

	/* init the top row */
	vec init_pv = ofsv - giv - gev - acc_gev, init_ev = zv;
	init_pv.print("init_pv");
	init_pv.store(ptr[0].s);
	init_ev.store(ptr[0].e);
	init_pv.store(ptr[0].f);
	for(uint64_t i = 0; i < vblen - 1; i++) {
		init_pv -= gev8;
		if((init_pv < ofsv - xtv) == 0xffff) { last_b_index = i + 1; }
		init_pv.print("init_pv");
		init_pv.store(ptr[i + 1].s);
		init_ev.store(ptr[i + 1].e);
		init_pv.store(ptr[i + 1].f);
	}
	ptr[last_b_index].s[0] = OFS;
	ptr[last_b_index].s[1] = last_b_index - first_b_index;
	ptr[last_b_index].s[2] = first_b_index;

	#define _gap(_i)		( ((_i) > 0 ? gi : 0) + (_i) * ge )
	vec max(OFS);
	uint64_t amax = 0;
	struct _dp *bmax = ptr + last_b_index;
	for(uint64_t a_index = 0; a_index < alen; a_index++) {
		debug("a_index(%llu), ch(%c), b_range(%llu, %llu)", a_index, a[a_index], first_b_index, last_b_index);

		prev = ptr; ptr += last_b_index + 1 - first_b_index;
		char_vec av(encode_a(a[a_index]));

		#define _update_vector(_i) { \
			char_vec bv; bv.load_encode_b(&b[(_i) * vec::LEN], blen - (_i) * vec::LEN); \
			vec th(prev[_i].s), te(prev[_i].e); \
			/* calc tentative s and e */ \
			pe = vec::max(te, th - giv) - gev; \
			pv = vec::max(th.dsl(ch) + smv.shuffle(av | bv), pe); ch = th; \
			av.print("av"); bv.print("bv"); smv.shuffle(av | bv).print("score"); \
			/* calc f from the tentative s */ \
			pf = vec::max(pv - giv, (pf>>7) - gev); \
			pf = vec::max(pf, (pf<<1) - gev); \
			pf = vec::max(pf, (pf<<2) - gev2); \
			pf = vec::max(pf, (pf<<4) - gev4); \
			/* fixup s */ \
			pv = vec::max(pv, pf); \
			/* update max */ \
			mv = vec::max(mv, pv); \
			pv.print("pv"); pe.print("pe"); pf.print("pf"); \
		}

		vec ch(OFS + _gap(a_index)), pv, pe, pf, mv(max[0]);
		while(1) {
			_update_vector(first_b_index);
			if((pv < max - xt) != 0xffff) { break; }
			first_b_index++; ptr--;
		}

		pv.store(ptr[first_b_index].s);
		pe.store(ptr[first_b_index].e);
		pf.store(ptr[first_b_index].f);

		debug("ptr(%p), b_range(%llu, %llu)", ptr, first_b_index, last_b_index);
		uint64_t next_last_b_index = last_b_index;
		for(uint64_t b_index = first_b_index + 1; b_index < last_b_index; b_index++) {
			_update_vector(b_index);
			if((pv < max - xt) != 0xffff) { next_last_b_index = b_index + 1; }

			pv.store(ptr[b_index].s);
			pe.store(ptr[b_index].e);
			pf.store(ptr[b_index].f);
		}
		last_b_index = next_last_b_index;
		debug("updated b_range(%llu, %llu)", first_b_index, last_b_index);

		if(last_b_index < vblen) {
			char_vec bv; bv.load_encode_b(&b[last_b_index * vec::LEN], blen - last_b_index * vec::LEN);
			vec d = gev8.dsl(ch) + smv.shuffle(av | bv);
			mv = vec::max(mv, d);

			pf = vec::max(d - giv, (pf>>7) - gev);
			pf.set(pf[0]); pf -= acc_gev;
			pv = vec::max(d, pf);
			while(last_b_index < vblen) {
				if((pv < max - xt) == 0xffff) { break; }

				pv.print("pv (tail)");
				pf.print("pf (tail)");

				pv.store(ptr[last_b_index].s);
				zv.store(ptr[last_b_index].e);
				pf.store(ptr[last_b_index].f);
				last_b_index++;

				pv -= gev8;
				pf -= gev8;
			}
		}

		int32_t m = mv.hmax();
		if(m > max[0]) { amax = a_index + 1; bmax = ptr + last_b_index; }
		max.set(m);

		ptr[last_b_index].s[0] = m;
		ptr[last_b_index].s[1] = last_b_index - first_b_index;
		ptr[last_b_index].s[2] = first_b_index;
	}

	debug("bmax(%p, %d, %d, %d)", bmax, bmax->s[0] - OFS, bmax->s[1], bmax->s[2]);
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	r->apos = 0;
	r->bpos = 0;
	struct _dp *bbase = bmax - bmax->s[1];
	for(uint64_t b_index = 0; b_index < bmax->s[1]; b_index++) {
		vec pv(bbase[b_index].s);
		pv.print("pv (search)");
		if(pv == max) {
			r->apos = amax;
			r->bpos = vec::LEN * (bmax->s[2] + b_index) + (tzcnt(pv == max)>>1) + 1;
			break;
		}
	}
	return(max[0] - OFS);
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
