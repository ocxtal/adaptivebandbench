
/**
 * @file striped.cc
 *
 * @brief striped parallelization of the standard banded matrix
 */
#include <string.h>
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
	#define _s(_p, _s, _t)	( (_p)[         (_s) * VLEN + (_t)] )
	#define _e(_p, _s, _t)	( (_p)[2 * bw + (_s) * VLEN + (_t)] )
	#define _f(_p, _s, _t)	( (_p)[4 * bw + (_s) * VLEN + (_t)] )
	#define _vlen()			( 6 * bw )
	#define _blen()			( 2 * bw / VLEN )	/* block height */

	/* construct score profile vector */
	uint64_t tlen = (blen + 2 * bw - 1) / _blen();
	if(tlen < VLEN) { tlen = VLEN; }

	uint16_t *scv = (uint16_t *)((uint8_t *)work + sizeof(maxpos_t));
	uint16_t *base = scv + 4 * _blen() * tlen;

	#define _scv(_s, _t, _k)		( scv[(((_s) * 4 + (_k)) * tlen) + (_t)] )
	for(uint64_t s = 0; s < _blen(); s++) {
		for(uint64_t t = 0; t < tlen; t++) {
			uint64_t i = t * _blen() + s;
			int8_t ch = (i - bw) < blen ? encode_b(b[i - bw]) : encode_n();
			debug("s(%llu), t(%llu), idx(%llu), pos(%llu), ch(%c)", s, t, t + s * tlen, i, (i - bw) < blen ? b[i - bw] : '-');
			_scv(s, t, 0) = score_matrix[encode_a('A') | ch];
			_scv(s, t, 1) = score_matrix[encode_a('C') | ch];
			_scv(s, t, 2) = score_matrix[encode_a('G') | ch];
			_scv(s, t, 3) = score_matrix[encode_a('T') | ch];
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
		debug("i(%llu), s(%llu), t(%llu), ofs(%llu), s(%d, %d, %d)", i, s, t, s * VLEN + t, _s(curr, s, t) - 32768, _e(curr, s, t) - 32768, _f(curr, s, t) - 32768);
	}
	/* fix gap cells at (0, 0) */
	_e(curr, bw % _blen(), bw / _blen()) = OFS + gi;
	_f(curr, bw % _blen(), bw / _blen()) = OFS + gi;
	uint64_t smax = OFS, amax = 0;				/* max score and its position */

	vdp_t const giv  = seta_vdp(-1 * gi);
	vdp_t const gev  = seta_vdp(-1 * ge);
	vdp_t const gebv = seta_vdp(-1 * ge * _blen());
	vdp_t max = seta_vdp(OFS);
	#ifdef DEBUG_CNT
		uint64_t fcnt = 0;						/* lazy-f chain length */
		bench_t body, lazyf;
		bench_init(body);
		bench_init(lazyf);
	#endif
	for(uint64_t apos = 0; apos < alen; apos++) {
		#ifdef DEBUG_CNT
			bench_start(body);
		#endif

		debug("apos(%llu)", apos);
		prev = curr; curr += _vlen();

		/* fetch the next base */
		int8_t ch = encode_a(a[apos]);

		/* init f */
		vdp_t pf = zero_vdp();

		/* speculative score calculation */
		for(size_t bofs = 0; bofs < 2 * bw; bofs += VLEN) {
			debug("bofs(%llu)", bofs);
			/* load prev vectors */
			vdp_t pv = load_vdp(&_s(prev, 0, bofs));
			vdp_t pe = load_vdp(&_e(prev, 1, bofs));
			if(bofs + VLEN >= 2 * bw) {
				pe = load_vdp(&_e(prev, 0, 0));
				pe = bsr_vdp(pe, 1);
			}

			print_vdp(pv);
			print_vdp(pe);

			/* load score */
			size_t const s = (apos + bofs / VLEN) % _blen();
			size_t const t = (apos + bofs / VLEN) / _blen();

			vdp_t const sc = loadu_vdp(&_scv(s, t, ch));
			pv = max_vdp(pv, giv);				/* ensure not rounded by adding the score profile */
			pv = add_vdp(pv, sc);

			pe = sub_vdp(pe, gev);
			pf = sub_vdp(pf, gev);
			pv = max_vdp(pv, max_vdp(pe, pf));
			pe = max_vdp(pe, sub_vdp(pv, giv));
			pf = max_vdp(pf, sub_vdp(pv, giv));

			print_vdp(add_vdp(sc, seta_vdp(32768)));
			print_vdp(pv);
			print_vdp(pe);
			print_vdp(pf);

			store_vdp(&_s(curr, 0, bofs), pv);
			store_vdp(&_e(curr, 0, bofs), pe);
			store_vdp(&_f(curr, 0, bofs), pf);

			/* update max */
			max = max_vdp(max, pv);
		}

		/* propagate gap to the end */
		for(size_t i = 0; i < VLEN - 1; i++) {
			vdp_t const tf = sub_vdp(bsl_vdp(pf, 1), gebv);
			pf = max_vdp(pf, tf);
		}

		#ifdef DEBUG_CNT
			bench_end(body);
			bench_start(lazyf);
		#endif

		/* lazy-f loop */
		debug("lazy-f");
		while(1) {
			pf = bsl_vdp(pf, 1);				/* shift by one to move to the next block */
			for(size_t bofs = 0; bofs < 2 * bw; bofs += VLEN) {
				vdp_t pv = load_vdp(&_s(curr, 0, bofs));
				pf = sub_vdp(pf, gev);
				if(gt_vdp(pf, sub_vdp(pv, giv)) == 0) { goto _tail; }
				// if((pv - giv < pf) == 0) { goto _tail; }

				#ifdef DEBUG_CNT
					fcnt++;
				#endif
				pv = max_vdp(pv, pf);			/* max score cannot be updated here because max is always  */
				store_vdp(&_s(curr, 0, bofs), pv);
				store_vdp(&_f(curr, 0, bofs), pf);
			}
		}
	_tail:;

		/* update maxpos */
		uint16_t m = hmax_vdp(max);
		debug("m(%d), pm(%d)", m - OFS, smax - OFS);
		if(m > smax) { smax = m; amax = apos + 1; }

		#ifdef DEBUG_CNT
			bench_end(lazyf);
		#endif
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	#ifdef DEBUG_CNT
		r->ccnt = alen * 2 * bw;
		r->fcnt = fcnt * VLEN;
		r->time[0] = bench_get(body) / 1000;
		r->time[1] = bench_get(lazyf) / 1000;
	#endif

	base += _vlen() * amax;
	uint16_t const m = hmax_vdp(max);
	debug("m(%u), amax(%llu)", m, amax);
	for(uint64_t bofs = 0; bofs < bw; bofs += VLEN) {
		vdp_t const s = load_vdp(&_s(base, 0, bofs));
		vdp_t const t = seta_vdp(m);
		if(eq_vdp(s, t)) {
			r->apos = amax;
			r->bpos = (tzcnt(eq_vdp(s, t)) / 2) * _blen() + bofs / VLEN - bw + amax;
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
