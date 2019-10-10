
/**
 * @file vertical.cc
 *
 * @brief vertical parallelization of the standard banded matrix
 */
#include <string.h>
#include "util.h"
#include "log.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn vertical_affine
 */
int
cat(vertical_affine, SUFFIX)(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, uint32_t bw)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("start, %s, %s", a, b);

	/* s: score vector, e: horizontal gap, f: verticalical gap */
	#define _s(_p, _i)		( (_p)[         (_i)] )
	#define _e(_p, _i)		( (_p)[2 * bw + (_i)] )
	#define _f(_p, _i)		( (_p)[4 * bw + (_i)] )
	#define _vlen()			( 6 * bw )
	uint8_t c[2 * bw + VLEN];

	/* init the leftmost vector (verticalically placed) */
	uint16_t *base = (uint16_t *)((uint8_t *)work + sizeof(maxpos_t)), *curr = base, *prev = base;
	#define _gap(_i)		( ((_i) > 0 ? gi : 0) + (_i) * ge )
	for(uint64_t i = 0; i < 2 * bw; i++) {
		if(i < bw) {
			_s(curr, i) = _e(curr, i) = _f(curr, i) = 0;
			c[i + 1] = 0;
		} else {
			_s(curr, i) = _f(curr, i) = OFS + _gap(i - bw);
			_e(curr, i) = 0;
			c[i + 1] = (i - bw) < blen ? encode_b(b[i - bw]) : encode_n();
		}
	}
	_e(curr, bw) = OFS + gi;					/* fix gap cells at (0, 0) */
	_f(curr, bw) = OFS + gi;
	uint64_t smax = OFS, amax = 0;				/* max score and its position */

	vdp_t const z = zero_vdp();
	vdp_t const giv  = seta_vdp(-1 * gi);
	vdp_t const gev  = seta_vdp(-1 * ge);
	vdp_t const gev2 = seta_vdp(-2 * ge);
	vdp_t const gev4 = seta_vdp(-4 * ge);
	vmat_t const smv = loadu_vmat(score_matrix);
	vdp_t max = seta_vdp(OFS);

	for(uint64_t apos = 0; apos < alen; apos++) {
		debug("apos(%llu)", apos);
		vchar_t const va = seta_vchar(encode_a(a[apos]));
		prev = curr; curr += _vlen();

		/* fetch the next base */
		c[2 * bw] = (apos + bw - 1) < blen ? encode_b(b[apos + bw - 1]) : encode_n();

		/* init f */
		vdp_t pf = zero_vdp();
		vdp_t ce = load_vdp(&_e(prev, 0));

		/* bpos = apos + bofs - bw/2 */
		for(uint64_t bofs = 0; bofs < 2 * bw; bofs += VLEN) {
			debug("bofs(%llu)", bofs);
			/* load prev vectors */
			vdp_t pv = load_vdp(&_s(prev, bofs));
			vdp_t te = load_vdp(&_e(prev, bofs + VLEN));
			if(bofs + VLEN >= 2 * bw) { te = z; }
			vdp_t pe = bsrd_vdp(te, ce);
			ce = te;

			/* calc score */
			vchar_t const vb = loadu_vchar(&c[bofs + 1]);
			store_vchar(&c[bofs], vb);			/* shift by one */
			vchar_t const xt = or_vchar(va, vb);
			vmat_t const yt  = shuffle_vmat(smv, cvt_vchar_vmat(xt));

			print_vdp(pe);
			print_vdp(add_vdp(cvt_vmat_vdp(yt), seta_vdp(32768)));

			pv = max_vdp(pv, giv);				/* ensure not rounded by adding the score profile */
			pv = add_vdp(pv, cvt_vmat_vdp(yt));

			/* calc e */
			pe = sub_vdp(pe, gev);
			pv = max_vdp(pv, pe);

			/* calc f: verticalical gap propagation */
			pf = sub_vdp(bsr_vdp(pf, 7), gev);
			pf = max_vdp(pf, sub_vdp(pv, giv));
			pf = max_vdp(pf, sub_vdp(bsl_vdp(pf, 1), gev));
			pf = max_vdp(pf, sub_vdp(bsl_vdp(pf, 2), gev2));
			pf = max_vdp(pf, sub_vdp(bsl_vdp(pf, 4), gev4));

			/* fixup s */
			pv = max_vdp(pv, pf);
			pe = max_vdp(pe, sub_vdp(pv, giv));

			store_vdp(&_s(curr, bofs), pv);
			store_vdp(&_e(curr, bofs), pe);
			store_vdp(&_f(curr, bofs), pf);

			print_vdp(pv); print_vdp(pe); print_vdp(pf);

			/* update max */
			max = max_vdp(max, pv);
		}

		/* update maxpos */
		uint16_t const m = hmax_vdp(max);
		debug("m(%u), pm(%u)", m, smax);
		if(m > smax) { smax = m; amax = apos + 1; }
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	#ifdef DEBUG_CNT
		r->ccnt = alen * 2 * bw;
		r->fcnt = 0;
	#endif

	base += _vlen() * amax;
	uint16_t const m = hmax_vdp(max);
	debug("m(%u), amax(%llu)", m, amax);
	for(uint64_t bofs = 0; bofs < bw; bofs += VLEN) {
		vdp_t const s = load_vdp(&_s(base, bofs));
		vdp_t const t = seta_vdp(m);
		if(eq_vdp(s, t)) {
			r->apos = amax;
			r->bpos = tzcnt(eq_vdp(s, t)) / 2 + bofs - bw + amax;
			debug("amax(%llu), bmax(%llu), mask(%hx), bofs(%llu)", r->apos, r->bpos, eq_vdp(s, t), bofs);
			break;
		}
	}
	debug("score(%d)", m - OFS);
	return(m - OFS);
}

#ifdef MAIN
#include <assert.h>
#include <stdlib.h>
int main_ext(int argc, char *argv[])
{
	uint64_t alen = strlen(argv[2]);
	uint64_t blen = strlen(argv[3]);
	char *a = (char *)malloc(alen + VLEN + 1);
	char *b = (char *)malloc(blen + VLEN + 1);

	memcpy(a, argv[2], alen);
	memset(a + alen, 0, VLEN + 1);
	memcpy(b, argv[3], blen);
	memset(b + blen, 0x80, VLEN + 1);

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, atoi(argv[4]), atoi(argv[5]));


	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	if(0) {
		printf("./a.out AAA AAA 2 -3 -5 -1 30\n");
	}

	int score = vertical_affine(
		work,
		a, alen, b, blen,
		score_matrix,
		atoi(argv[6]),
		atoi(argv[7]),
		atoi(argv[8]),
		32);
	printf("%d\n", score);
	free(a); free(b); free(work);
	return(0);
}

int main(int argc, char *argv[])
{
	if(argc > 1) { return(main_ext(argc, argv)); }

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, 1, -1);
	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	#define a(s, p, q) { \
		assert(vertical_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10, 32) == (s)); \
	}
	a( 0, "", "");
	a( 0, "A", "");
	a( 1, "A", "A");
	a( 3, "AAA", "AAA");
	a( 0, "AAA", "TTT");
	a( 3, "AAAGGG", "AAATTTTTT");
	a( 3, "TTTGGGGGAAAA", "TTTCCCCCCCCAAAA");
	a( 4, "TTTAAAA", "TTTCCAAAA");
	a( 4, "GGGCCCC", "GGGAACCCC");
	a( 6, "AAACAAAAAGGG", "AAAAAAAATTTTTTT");
	a( 6, "AAAAAAAATTTTTTT", "AAACAAAAAGGG");
	a( 5, "AAACCAAAAAGGG", "AAAAAAAATTTTTTT");
	a( 5, "AAAAAAAATTTTTTT", "AAACCAAAAAGGG");

	free(work);
	return(0);
}
#endif

/**
 * end of vertical.cc
 */
