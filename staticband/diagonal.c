
/**
 * @file diagonal.cc
 *
 * @brief diagonal parallelization of the standard banded matrix
 */
#include <string.h>
#include "util.h"
#include "log.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )

/**
 * @fn diagonal_affine
 */
int
cat(diagonal_affine, SUFFIX)(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t score_matrix[16], int8_t gi, int8_t ge, int16_t xt, uint32_t bw)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	/* s: score vector, e: horizontal gap, f: vertical gap */
	#define _s(_p, _i)		( (_p)[         (_i)] )
	#define _e(_p, _i)		( (_p)[    bw + (_i)] )
	#define _f(_p, _i)		( (_p)[2 * bw + (_i)] )
	#define _vlen()			( 3 * bw )
	uint8_t abuf[bw + VLEN], bbuf[bw + VLEN];

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
	uint64_t smax = OFS, pmax = 0;				/* max score and its position */

	vdp_t const z = zero_vdp();
	vdp_t const giv = seta_vdp(-1 * gi);
	vdp_t const gev = seta_vdp(-1 * ge);
	vmat_t const smv = loadu_vmat(score_matrix);
	vdp_t max = seta_vdp(OFS);

	for(uint64_t p = 2; p < alen + blen + 1; p++) {
		curr += _vlen();
		uint16_t *pprv = curr - 2 * _vlen(), *prev = curr - _vlen();

		if(p & 0x01) {
			debug("D");
			bbuf[bw] = (p / 2 + bw / 2 - 1) < blen ? encode_b(b[p / 2 + bw / 2 - 1]) : encode_n();

			vchar_t cb = loadu_vchar(&bbuf[0]);
			vdp_t ce = load_vdp(&_e(prev, 0));
			for(size_t i = 0; i < bw; i += VLEN) {
				debug("loop: %llu", i);

				/* load the previous buffers */
				vchar_t const va = loadu_vchar(&abuf[i]);
				vchar_t const tb = loadu_vchar(&bbuf[i + VLEN]);

				vdp_t pf = load_vdp(&_f(prev, i));
				vdp_t te = load_vdp(&_e(prev, i + VLEN));

				/* clear if the tail */
				if(i + VLEN >= bw) { te = z; }
				vdp_t pe = bsrd_vdp(te, ce);
				ce = te;

				/* rotate the vectors to align the vectors */
				vchar_t const vb = bsrd_vchar(tb, cb);
				storeu_vchar(&bbuf[i], vb);				/* store shifted */
				cb = tb;

				/* load the second previous */
				vdp_t pv = load_vdp(&_s(pprv, i));
				vchar_t const xt = or_vchar(va, vb);
				vmat_t const yt  = shuffle_vmat(smv, cvt_vchar_vmat(xt));

				print_vdp(pv);
				print_vdp(add_vdp(cvt_vmat_vdp(yt), seta_vdp(32768)));

				pv = add_vdp(pv, cvt_vmat_vdp(yt));

				pe = sub_vdp(pe, gev);
				pf = sub_vdp(pf, gev);
				pv = max_vdp(pv, max_vdp(pe, pf));
				pe = max_vdp(pe, sub_vdp(pv, giv));
				pf = max_vdp(pf, sub_vdp(pv, giv));

				store_vdp(&_s(curr, i), pv);
				store_vdp(&_e(curr, i), pe);
				store_vdp(&_f(curr, i), pf);

				print_vdp(pv); print_vdp(pe); print_vdp(pf);

				max = max_vdp(max, pv);
			}
		} else {
			debug("R");
			vchar_t ca = seta_vchar((int8_t)((p / 2 + bw / 2 - 1) < alen ? encode_a(a[p / 2 + bw / 2 - 1]) : encode_n()));
			vdp_t cf = zero_vdp();
			for(size_t i = 0; i < bw; i += VLEN) {
				debug("loop: %llu", i);

				/* load the previous buffers */
				vchar_t const ta = loadu_vchar(&abuf[i]);
				vchar_t const vb = loadu_vchar(&bbuf[i]);

				vdp_t pe = load_vdp(&_e(prev, i));
				vdp_t tf = load_vdp(&_f(prev, i));
				vdp_t pf = bsld_vdp(tf, cf);
				cf = tf;

				/* rotate the vectors */
				vchar_t const va = bsld_vchar(ta, ca);
				storeu_vchar(&abuf[i], va);				/* store shifted */
				ca = ta;

				/* load the second previous */
				vdp_t pv = load_vdp(&_s(pprv, i));
				vchar_t const xt = or_vchar(va, vb);
				vmat_t const yt  = shuffle_vmat(smv, cvt_vchar_vmat(xt));

				print_vdp(pv);
				print_vdp(add_vdp(cvt_vmat_vdp(yt), seta_vdp(32768)));

				pv = add_vdp(pv, cvt_vmat_vdp(yt));

				pe = sub_vdp(pe, gev);
				pf = sub_vdp(pf, gev);
				pv = max_vdp(pv, max_vdp(pe, pf));
				pe = max_vdp(pe, sub_vdp(pv, giv));
				pf = max_vdp(pf, sub_vdp(pv, giv));

				store_vdp(&_s(curr, i), pv);
				store_vdp(&_e(curr, i), pe);
				store_vdp(&_f(curr, i), pf);

				print_vdp(pv); print_vdp(pe); print_vdp(pf);

				max = max_vdp(max, pv);
			}
		}

		/* update maxpos */
		uint16_t const m = hmax_vdp(max);
		debug("m(%u), pm(%llu)", m, smax);
		if(m > smax) { smax = m; pmax = p; }
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	#ifdef DEBUG_CNT
		r->ccnt = bw * (alen + blen - 1);
		r->fcnt = 0;
	#endif

	base += _vlen() * pmax;
	uint16_t const m = hmax_vdp(max);
	debug("m(%u), amax(%llu)", m, smax);
	for(uint64_t i = 0; i < bw; i += VLEN) {
		vdp_t const s = load_vdp(&_s(base, i));
		vdp_t const t = seta_vdp(m);
		if(eq_vdp(s, t)) {
			uint64_t q = i + tzcnt(eq_vdp(s, t)) / 2 - bw / 2;
			r->apos =  pmax      / 2 - q;
			r->bpos = (pmax + 1) / 2 + q;
			debug("amax(%llu), bmax(%llu)", r->apos, r->bpos);
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
	int score = diagonal_affine(
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
	char const *a = "aattcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	if(argc > 1) { return(main_ext(argc, argv)); }

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, 1, -1);

	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	#define a(s, p, q) { \
		assert(diagonal_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10, 32) == (s)); \
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
 * end of diagonal.cc
 */
