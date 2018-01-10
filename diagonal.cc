
/**
 * @file diagonal.cc
 *
 * @brief diagonal parallelization of the standard banded matrix
 */
#include <string.h>
#include "sse.h"
#include "util.h"

#define MIN 	( 0 )
#define OFS 	( 32768 )

/**
 * @fn diagonal_affine
 */
int
diagonal_affine(
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

		/* vectors at p = 1 */
		_s(curr, bw/2 - 1 - i) = OFS - i * max_match + _gap(i*2 + 1);
		_s(curr, bw/2     + i) = OFS - i * max_match + _gap(i*2 + 1);

		_e(curr, bw/2 - 1 - i) = OFS - i * max_match + _gap(i*2 + 1);
		_e(curr, bw/2     + i) = OFS - i * max_match + _gap(i*2 + 1);
		_f(curr, bw/2 - 1 - i) = OFS - i * max_match + _gap(i*2 + 1);
		_f(curr, bw/2     + i) = OFS - i * max_match + _gap(i*2 + 1);

		/* char buffers */
		abuf[bw/2 - 1 - i] = i < alen ? encode_a(a[i]) : 0;
		abuf[bw/2     + i] = 0;
		bbuf[bw/2 - 1 - i] = 0;
		bbuf[bw/2     + i] = i < blen ? encode_b(b[i]) : 0;
	}
	uint64_t smax = OFS, pmax = 0;				/* max score and its position */

	char_vec const cz;							/* zero */
	vec const z, giv(-gi), gev(-ge), smv((uint16_t const *)score_matrix);
	vec max(OFS);
	for(uint64_t p = 2; p < alen + blen + 1; p++) {
		curr += _vlen();
		uint16_t *pprv = curr - 2 * _vlen(), *prev = curr - _vlen();

		if(p & 0x01) {
			debug("D");
			bbuf[bw] = (p / 2 + bw / 2) < blen ? encode_b(b[p / 2 + bw / 2]) : 0;

			char_vec cb((int8_t const *)&bbuf[0]);
			vec ch(&_s(curr, 0)), ce(&_e(curr, 0)); ch -= giv;
			for(uint64_t i = 0; i < bw; i += vec::LEN) {
				debug("loop: %llu", i);

				/* load the previous buffers */
				char_vec va((int8_t const *)&abuf[i]), tb((int8_t const *)&bbuf[i + vec::LEN]);
				vec th(&_s(prev, i + vec::LEN)), te(&_e(prev, i + vec::LEN)), pf(&_f(prev, i));
				th -= giv;						/* add gap-open, common for the two directions */

				/* clear if the tail */
				if(i + vec::LEN >= bw) { tb = cz; th = z; te = z; }
				th.print("th"); te.print("te");

				/* rotate the vectors to align the vectors */
				char_vec vb = tb.dsr(cb);
				vec pv = ch, ph = th.dsr(ch), pe = te.dsr(ce);
				cb = tb; ch = th; ce = te;
				vb.store(&bbuf[i]);				/* store shifted */
				pv.print("pv"); pe.print("pe"); pf.print("pf");
				va.print("a"); vb.print("b"); smv.shuffle(va | vb).print("score");

				/* load the second previous */
				vec ppv(&_s(pprv, i));
				ppv += smv.shuffle(va | vb);

				/* calc next */
				vec ne = vec::max(ph, pe) - gev;
				vec nf = vec::max(pv, pf) - gev;
				ne.store(&_e(curr, i)); ne.print("ne");
				nf.store(&_f(curr, i)); nf.print("nf");

				/* update s */
				vec nv = vec::max(vec::max(ne, nf), ppv);
				nv.store(&_s(curr, i)); nv.print("nv");

				max = vec::max(max, nv);
			}
		} else {
			debug("R");
			char_vec ca((int8_t)((p / 2 + bw / 2) < alen ? encode_a(a[p / 2 + bw / 2]) : 0));
			vec cv, cf;
			for(uint64_t i = 0; i < bw; i += vec::LEN) {
				debug("loop: %llu", i);

				/* load the previous buffers */
				char_vec ta((int8_t const *)&abuf[i]), vb((int8_t const *)&bbuf[i]);
				vec tv(&_s(prev, i)), pe(&_e(prev, i)), tf(&_s(prev, i));
				tv -= giv;
				tv.print("tv"); tf.print("tf");

				/* rotate the vectors */
				char_vec va = ta.dsl(ca);
				vec ph = tv, pv = tv.dsl(cv), pf = tf.dsl(cf);
				ca = ta; cv = tv; cf = pf;
				va.store(&abuf[i]);				/* store shifted */
				pv.print("pv"); pe.print("pe"); pf.print("pf");
				va.print("a"); vb.print("b"); smv.shuffle(va | vb).print("score");

				/* load the second previous */
				vec ppv(&_s(pprv, i));
				ppv += smv.shuffle(va | vb);

				/* calc next */
				vec ne = vec::max(ph, pe) - gev;
				vec nf = vec::max(pv, pf) - gev;
				ne.store(&_e(curr, i)); ne.print("ne");
				nf.store(&_f(curr, i)); nf.print("nf");

				/* update s */
				vec nv = vec::max(vec::max(ne, nf), ppv);
				nv.store(&_s(curr, i)); nv.print("nv");

				max = vec::max(max, nv);
			}
		}

		/* update maxpos */
		uint16_t m = max.hmax();
		debug("m(%u), pm(%llu)", m, smax);
		if(m > smax) { smax = m; pmax = p; }
	}

	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	#ifdef debug
		r->ccnt = bw * (alen + blen - 1);
		r->fcnt = 0;
	#endif

	base += _vlen() * pmax;
	uint16_t m = max.hmax();
	debug("m(%u), amax(%llu)", m, smax);
	for(uint64_t i = 0; i < bw; i += vec::LEN) {
		vec s(&_s(base, i)), t(m);
		s.print("s");
		if(s == t) {
			uint64_t q = i + tzcnt(s == t) / 2 - bw / 2;
			r->apos = pmax / 2 - q;
			r->bpos = pmax / 2 + q;
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
