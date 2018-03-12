
/**
 * @file adaptive.cc
 */
#include <string.h>
#include "sse.h"
#include "util.h"
#include "log.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )
#define roundup(a, bound)		( (((a) + (bound) - 1) / (bound)) * (bound) )

/**
 * @fn adaptive_affine
 */
int
adaptive_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, uint32_t bw)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	/* s: score vector, e: horizontal gap, f: vertical gap */
	#define _s(_p, _i)		( (_p)[         (_i)] )
	#define _e(_p, _i)		( (_p)[    bw + (_i)] )
	#define _f(_p, _i)		( (_p)[2 * bw + (_i)] )
	#define _vlen()			( 3 * bw )
	uint8_t abuf[bw + vec::LEN], bbuf[bw + vec::LEN];

	/* init the leftmost vector (vertically placed) */
	uint16_t *base = (uint16_t *)((uint8_t *)work + sizeof(maxpos_t));
	uint16_t *curr = base + _vlen();
	int8_t max_match = extract_max_score(score_matrix);
	#define _gap(_i)		( ((_i) > 0 ? gi : 0) + (_i) * ge )
	for(uint64_t i = 0; i < bw / 2; i++) {
		uint16_t *prev = curr - _vlen();

		/* vector at p = 0 */
		_s(prev, bw/2 - 1 - i) = OFS - (i + 1) * max_match + _gap((i + 1)*2);
		_s(prev, bw/2     + i) = OFS -  i      * max_match + _gap( i     *2);
		_e(prev, bw/2 - 1 - i) = 0;
		_e(prev, bw/2     + i) = 0;
		_f(prev, bw/2 - 1 - i) = 0;
		_f(prev, bw/2     + i) = 0;

		/* vectors at p = 1 */
		_s(curr, bw/2 - 1 - i) = OFS - i * max_match + _gap(i*2 + 1);
		_s(curr, bw/2     + i) = OFS - i * max_match + _gap(i*2 + 1);
		_e(curr, bw/2 - 1 - i) = i == 0 ? OFS + _gap(1) : 0;
		_e(curr, bw/2     + i) = 0;
		_f(curr, bw/2 - 1 - i) = 0;
		_f(curr, bw/2     + i) = i == 0 ? OFS + _gap(1) : 0;

		/* char buffers */
		abuf[bw/2 - 1 - i] = i < alen ? encode_a(a[i]) : encode_n();
		abuf[bw/2     + i] = 0;
		bbuf[bw/2 - 1 - i] = 0;
		bbuf[bw/2     + i] = i < blen ? encode_b(b[i]) : encode_n();
	}

	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = _s(base, bw - 1) <= _s(base, 0) ? RR : RD;		/* always RD */
	uint64_t apos = bw / 2;
	uint64_t bpos = bw / 2;

	int32_t max = OFS;
	uint64_t amax = 0, bmax = 0;					/* max score and its position */
	for(uint64_t p = 0; p < (uint64_t)(alen+blen-1); p++) {
		curr += _vlen();
		uint16_t *pprv = curr - 2 * _vlen(), *prev = curr - _vlen();

		debug("%lld, %d, (%d, %d), max(%d), (%lu, %lu)", dir, _s(prev, bw - 1) > _s(prev, 0), _s(prev, bw - 1), _s(prev, 0), max, amax, bmax);
		dir = dir_trans[_s(prev, bw - 1) > _s(prev, 0)][dir];
		switch(dir) {
		case RR: {
			int8_t ach = apos < alen ? encode_a(a[apos]) : encode_n(); apos++;
			for(uint64_t i = 0; i < bw; i++) {
				int32_t ph = _s(prev, i), pv = _s(prev, i - 1), pe = _e(prev, i), pf = _f(prev, i - 1), pd = _s(pprv, i - 1);
				if(i == 0) { pv = 0; pf = 0; pd = 0; }
				int32_t score = score_matrix[ach | bbuf[i]];
				int8_t at = abuf[i]; abuf[i] = ach; ach = at;		/* shift by one */

				pe = MAX2(0, MAX2(ph + gi, pe) + ge);
				pf = MAX2(0, MAX2(pv + gi, pf) + ge);
				pv = MAX4(0, pd + score, pe, pf);
				_s(curr, i) = pv; _e(curr, i) = pe; _f(curr, i) = pf;

				if(pv > max) { max = pv; amax = apos - i; bmax = bpos + i; }
			}
		} break;
		case RD: {
			bbuf[bw] = bpos < blen ? encode_b(b[bpos]) : encode_n(); bpos++;
			for(uint64_t i = 0; i < bw; i++) {
				int32_t ph = _s(prev, i + 1), pv = _s(prev, i), pe = _e(prev, i + 1), pf = _f(prev, i), pd = _s(pprv, i);
				if(i == bw - 1) { ph = 0; pe = 0; }
				int32_t score = score_matrix[abuf[i] | bbuf[i + 1]];
				bbuf[i] = bbuf[i + 1];								/* shift by one */

				pe = MAX2(0, MAX2(ph + gi, pe) + ge);
				pf = MAX2(0, MAX2(pv + gi, pf) + ge);
				pv = MAX4(0, pd + score, pe, pf);
				_s(curr, i) = pv; _e(curr, i) = pe; _f(curr, i) = pf;

				if(pv > max) { max = pv; amax = apos - i; bmax = bpos + i; }
			}
		} break;
		case DR: {
			int8_t ach = apos < alen ? encode_a(a[apos]) : encode_n(); apos++;
			for(uint64_t i = 0; i < bw; i++) {
				int32_t ph = _s(prev, i), pv = _s(prev, i - 1), pe = _e(prev, i), pf = _f(prev, i - 1), pd = _s(pprv, i);
				if(i == 0) { pv = 0; pf = 0; }
				int32_t score = score_matrix[ach | bbuf[i]];
				int8_t at = abuf[i]; abuf[i] = ach; ach = at;		/* shift by one */

				pe = MAX2(0, MAX2(ph + gi, pe) + ge);
				pf = MAX2(0, MAX2(pv + gi, pf) + ge);
				pv = MAX4(0, pd + score, pe, pf);
				_s(curr, i) = pv; _e(curr, i) = pe; _f(curr, i) = pf;

				if(pv > max) { max = pv; amax = apos - i; bmax = bpos + i; }
			}
		} break;
		case DD: {
			bbuf[bw] = bpos < blen ? encode_b(b[bpos]) : encode_n(); bpos++;
			for(uint64_t i = 0; i < bw; i++) {
				int32_t ph = _s(prev, i + 1), pv = _s(prev, i), pe = _e(prev, i + 1), pf = _f(prev, i), pd = _s(pprv, i + 1);
				if(i == bw - 1) { ph = 0; pe = 0; pd = 0; }
				int32_t score = score_matrix[abuf[i] | bbuf[i + 1]];
				bbuf[i] = bbuf[i + 1];								/* shift by one */

				pe = MAX2(0, MAX2(ph + gi, pe) + ge);
				pf = MAX2(0, MAX2(pv + gi, pf) + ge);
				pv = MAX4(0, pd + score, pe, pf);
				_s(curr, i) = pv; _e(curr, i) = pe; _f(curr, i) = pf;

				if(pv > max) { max = pv; amax = apos - i; bmax = bpos + i; }
			}
		} break;
		}
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	r->apos = amax;
	r->bpos = bmax - bw + 1;
	#ifdef DEBUG_CNT
		r->ccnt = (apos + bpos - bw) * bw;
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

	int score = adaptive_affine(
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
		assert(adaptive_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10, 32) == (s)); \
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
 * end of adaptive.cc
 */
