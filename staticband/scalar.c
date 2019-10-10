
/**
 * @file scalar.cc
 *
 * @brief scalar parallelization of the standard banded matrix
 */
#include <string.h>
#include "util.h"
#include "log.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn scalar_affine
 */
int
scalar_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, uint32_t bw)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	/* s: score vector, e: horizontal gap, f: scalarical gap */
	#define _s(_p, _i)		( (_p)[         (_i)] )
	#define _e(_p, _i)		( (_p)[2 * bw + (_i)] )
	#define _f(_p, _i)		( (_p)[4 * bw + (_i)] )
	#define _vlen()			( 6 * bw )
	uint8_t c[2 * bw + VLEN];

	/* init the leftmost vector */
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
	int32_t max = OFS;
	uint64_t amax = 0, bmax = 0;					/* max score and its position */
	for(uint64_t apos = 0; apos < alen; apos++) {
		debug("apos(%llu)", apos);
		int8_t ach = encode_a(a[apos]);
		prev = curr; curr += _vlen();

		/* fetch the next base */
		c[2 * bw] = (apos + bw - 1) < blen ? encode_b(b[apos + bw - 1]) : encode_n();

		/* init f */
		int32_t pf = 0;

		/* bpos = apos + bofs - bw/2 */
		for(size_t bofs = 0; bofs < 2 * bw; bofs++) {
			/* load prev vectors */
			int32_t pv = _s(prev, bofs);
			int32_t pe = bofs < (2 * bw - 1) ? _e(prev, bofs + 1) : 0;
			int32_t score = score_matrix[ach | c[bofs + 1]];
			debug("(%c, %c), pd(%d), ph(%d), pe(%d), score(%d)", "ACGT"[ach], "ACGT"[c[bofs + 1]>>2], pd - OFS, ph - OFS, pe - OFS, score - OFS);

			c[bofs] = c[bofs + 1];				/* shift by one */

			pe = MAX2(0, pe + ge);
			pf = MAX2(0, pf + ge);
			pv = MAX4(-gi, pv + score, pe, pf);

			pe = MAX2(pv + gi, pe);
			pf = MAX2(pv + gi, pf);
			_s(curr, bofs) = pv;
			_e(curr, bofs) = pe;
			_f(curr, bofs) = pf;
			if(pv > max) { max = pv; amax = apos + 1; bmax = bofs - bw + apos + 1; }

			debug("apos(%llu), bofs(%llu), e(%d), f(%d), s(%d), max(%d)", apos, bofs, pe - OFS, pf - OFS, pv - OFS, max - OFS);
		}

		/* update maxpos */
		debug("max(%d)", max - OFS);
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	r->apos = amax;
	r->bpos = bmax;
	#ifdef debug
		r->ccnt = alen * 2 * bw;
		r->fcnt = 0;
	#endif
	return(max - OFS);
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

	int score = scalar_affine(
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
		assert(scalar_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10, 32) == (s)); \
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
 * end of scalar.cc
 */
