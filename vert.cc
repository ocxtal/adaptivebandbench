
/**
 * @file vert.cc
 *
 * @brief simd parallel blast_SemiGappedAlign algorithm
 */
#include <string.h>
#include "sse.h"
#include "util.h"
#include "log.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn vert_affine
 */
int
vert_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, int bw)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	/* s: score vector, e: horizontal gap, f: vertical gap */
	uint16_t *base = (uint16_t *)((uint8_t *)work + sizeof(maxpos_t)), *curr = base, *prev = base;
	#define _s(_p, _i)		( (_p)[         (_i)] )
	#define _e(_p, _i)		( (_p)[    bw + (_i)] )
	#define _f(_p, _i)		( (_p)[2 * bw + (_i)] )
	#define _m(_p)			( (_p)[3 * bw       ] )
	#define _vlen()			( 3 * bw + vec::LEN )
	uint8_t c[bw + 1];

	/* init the leftmost vector (vertically placed) */
	int8_t margin = -2 * extract_min_score(score_matrix);
	#define _gap(_i)		( ((_i) > 0 ? gi : 0) + (_i) * ge )
	for(uint64_t i = 0; i < bw; i++) {
		if(i < bw / 2) {
			_s(curr, i) = _e(curr, i) = _f(curr, i) = margin;
			c[i + 1] = 0;
		} else {
			_s(curr, i) = _f(curr, i) = OFS + _gap(i - bw/2);
			_e(curr, i) = margin;
			c[i + 1] = (i - bw / 2) < blen ? encode_b(b[i - bw / 2]) : 0;
		}
	}
	_m(curr) = OFS;								/* clear the max vector */
	uint64_t amax = 0;							/* max score position */

	vec const giv(-gi), gev(-ge), gev2(-2*ge), gev4(-4*ge), smv((uint16_t const *)score_matrix);
	vec max(OFS);
	for(uint64_t apos = 0; apos < alen; apos++) {
		debug("apos(%llu)", apos);
		char_vec av(encode_a(a[apos]));
		prev = curr; curr += _vlen();

		/* fetch the next base */
		c[bw] = (apos + bw / 2) < blen ? encode_b(b[apos + bw / 2]) : 0;

		/* init f */
		vec pf;

		/* bpos = apos + bofs - bw/2 */
		for(uint64_t bofs = 0; bofs < bw; bofs += vec::LEN) {
			/* load prev vectors */
			vec pv(&_s(prev, bofs));
			vec pe(&_e(prev, bofs + 1));		/* unaligned */

			/* calc score */
			char_vec bv((int8_t const *)&c[bofs + 1]);
			bv.store(&c[bofs]);					/* shift by one */
			pv += smv.shuffle(av | bv);
			av.print("av"); bv.print("bv"); smv.shuffle(av | bv).print("score");

			/* calc e */
			pe = pe + gev;
			pe.store(&_e(curr, bofs));

			/* calc tentative s */
			pv = vec::max(pv, pe); pv.print("pv (tentative)");

			/* calc f */
			pf = vec::max(pv - giv, (pf>>7) - gev); pf.print("pf1");
			pf = vec::max(pf, (pf<<1) - gev); pf.print("pf2");
			pf = vec::max(pf, (pf<<2) - gev2); pf.print("pf3");
			pf = vec::max(pf, (pf<<4) - gev4); pf.print("pf4");

			/* fixup s */
			pv = vec::max(pv, pf);

			pv.print("pv"); pe.print("pe"); pf.print("pf");
			pv.store(&_s(curr, bofs));
			pf.store(&_f(curr, bofs));

			/* update max */
			max = vec::max(max, pv);
		}

		/* update maxpos */
		uint16_t m = max.hmax();
		debug("m(%u), pm(%u)", m, _m(prev));
		if(m > _m(prev)) { amax = apos + 1; }
		_m(curr) = m;
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	base += _vlen() * amax;
	uint16_t m = max.hmax();
	debug("m(%u), amax(%llu)", m, amax);
	for(uint64_t bofs = 0; bofs < bw; bofs += vec::LEN) {
		vec p(&_s(base, bofs)), q(m);
		p.print("p");
		if(p == q) {
			r->apos = amax;
			r->bpos = tzcnt(p == q) + bofs - bw / 2 + amax;
			debug("amax(%llu), bmax(%llu)", r->apos, r->bpos);
			break;
		}
	}
	printf("score(%d)\n", m - OFS);
	return(m - OFS);
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

	int score = vert_affine(
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
		assert(vert_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10, 32) == (s)); \
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

	free(work);
	return(0);
}
#endif

/**
 * end of vert.cc
 */
