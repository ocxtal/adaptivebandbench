
/**
 * @file striped.cc
 *
 * @brief striped parallelization of the standard banded matrix
 */
#include <string.h>
#include "sse.h"
#include "util.h"
#include "log.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn striped_affine
 */
int
striped_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, uint32_t bw)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	/* s: score vector, e: horizontal gap, f: stripedical gap */
	#define _s(_p, _s, _t)	( (_p)[         (_s) * vec::LEN + (_t)] )
	#define _e(_p, _s, _t)	( (_p)[2 * bw + (_s) * vec::LEN + (_t)] )
	#define _f(_p, _s, _t)	( (_p)[4 * bw + (_s) * vec::LEN + (_t)] )
	#define _vlen()			( 6 * bw )
	#define _blen()			( 2 * bw / vec::LEN )	/* block height */

	/* construct score profile vector */
	uint64_t tlen = (blen + _blen() - 1) / _blen();
	if(tlen < vec::LEN) { tlen = vec::LEN; }

	uint16_t *scv = (uint16_t *)((uint8_t *)work + sizeof(maxpos_t));
	uint16_t *base = scv + 4 * _blen() * tlen;

	#define _scv(_s, _t, _k)		( scv[(((_s) * 4 + (_k)) * tlen) + (_t)] )
	for(uint64_t s = 0; s < _blen(); s++) {
		for(uint64_t t = 0; t < tlen; t++) {
			uint64_t i = t * _blen() + s;
			int8_t ch = (i - bw) < blen ? encode_b(b[i - bw]) : -1;
			debug("s(%llu), t(%llu), idx(%llu), pos(%llu), ch(%c)", s, t, t + s * tlen, i, (i - bw) < blen ? b[i - bw] : '-');
			_scv(s, t, 0) = ch != -1 ? score_matrix[encode_a('A') | ch] : -1;
			_scv(s, t, 1) = ch != -1 ? score_matrix[encode_a('C') | ch] : -1;
			_scv(s, t, 2) = ch != -1 ? score_matrix[encode_a('G') | ch] : -1;
			_scv(s, t, 3) = ch != -1 ? score_matrix[encode_a('T') | ch] : -1;
		}
	}

	/* init the leftmost vector */
	uint16_t *curr = base, *prev = base;
	#define _gap(_i)		( ((_i) > 0 ? gi : 0) + (_i) * ge )
	for(uint64_t i = 0; i < 2 * bw; i++) {
		uint64_t s = i % _blen(), t = i / _blen();
		if(i < bw) {
			_s(curr, s, t) = _e(curr, s, t) = _f(curr, s, t) = 0;
		} else {
			_s(curr, s, t) = _f(curr, s, t) = OFS + _gap(i - bw);
			_e(curr, s, t) = 0;
		}
		debug("i(%llu), s(%llu), t(%llu), ofs(%llu), s(%d, %d, %d)", i, s, t, s * vec::LEN + t, _s(curr, s, t) - 32768, _e(curr, s, t) - 32768, _f(curr, s, t) - 32768);
	}
	/* fix gap cells at (0, 0) */
	_e(curr, bw % _blen(), bw / _blen()) = OFS + gi;
	_f(curr, bw % _blen(), bw / _blen()) = OFS + gi;
	uint64_t smax = OFS, amax = 0;				/* max score and its position */

	vec const giv(-gi), gev(-ge);
	vec max(OFS);
	#ifdef DEBUG
		uint64_t fcnt = 0;						/* lazy-f chain length */
	#endif
	for(uint64_t apos = 0; apos < alen; apos++) {
		debug("apos(%llu)", apos);
		char_vec av(encode_a(a[apos]));
		prev = curr; curr += _vlen();

		/* fetch the next base */
		int8_t ch = encode_a(a[apos]);

		/* init f */
		vec pf, pv, cv(&_s(prev, 0, 0));

		/* speculative score calculation */
		for(uint64_t bofs = 0; bofs < 2 * bw; bofs += vec::LEN) {
			debug("bofs(%llu)", bofs);
			/* load prev vectors */
			vec ph(&_s(prev, 1, bofs)), pe(&_e(prev, 1, bofs));
			if(bofs + vec::LEN >= 2 * bw) {
				ph.load(&_s(prev, 0, 0)); ph >>= 1;
				pe.load(&_e(prev, 0, 0)); pe >>= 1;
			}

			/* load score */
			uint64_t s = (apos + bofs / vec::LEN) % _blen(), t = (apos + bofs / vec::LEN) / _blen();

			vec sc; sc.loadu(&_scv(s, t, ch)); sc.print("sc");
			cv = vec::max(cv, giv);				/* ensure not rounded by adding the score profile */
			sc += cv; cv = ph; sc.print("pv (tentative)q");

			/* calc e and f */
			pe = vec::max(ph - giv, pe) - gev;
			pf = vec::max(pv - giv, pf) - gev;
			pe.store(&_e(curr, 0, bofs)); pe.print("pe");
			pf.store(&_f(curr, 0, bofs)); pf.print("pf");

			/* calc s */
			pv = vec::max(sc, vec::max(pe, pf));
			pv.store(&_s(curr, 0, bofs)); pv.print("pv");

			/* update max */
			max = vec::max(max, pv);
		}
		pf <<= 1;								/* shift by one to move to the next block */

		/* lazy-f loop */
		debug("lazy-f");
		for(uint64_t bofs = 0; bofs < 2 * bw; bofs += vec::LEN) {
			vec pv(&_s(curr, 0, bofs));
			pv.print("pv"); pf.print("pf");
			if(!(pv < (pf -= gev))) { break; }
			#ifdef DEBUG
				fcnt++;
			#endif
			pv = vec::max(pv, pf);				/* max score cannot be updated here because max is always  */
			pv.store(&_s(curr, 0, bofs));
			pf.store(&_f(curr, 0, bofs));
		}

		/* update maxpos */
		uint16_t m = max.hmax();
		debug("m(%u), pm(%u)", m, smax);
		if(m > smax) { smax = m; amax = apos + 1; }
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	#ifdef DEBUG
		r->ccnt = alen * 2 * bw;
		r->fcnt = fcnt;
	#endif

	base += _vlen() * amax;
	uint16_t m = max.hmax();
	debug("m(%u), amax(%llu)", m, amax);
	for(uint64_t bofs = 0; bofs < bw; bofs += vec::LEN) {
		vec s(&_s(base, 0, bofs)), t(m);
		s.print("s");
		if(s == t) {
			r->apos = amax;
			r->bpos = (tzcnt(s == t) / 2) * _blen() + bofs / vec::LEN - bw + amax;
			debug("amax(%llu), bmax(%llu), mask(%hx), bofs(%llu)", r->apos, r->bpos, s == t, bofs);
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

	int score = striped_affine(
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
		assert(striped_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10, 32) == (s)); \
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
 * end of striped.cc
 */
