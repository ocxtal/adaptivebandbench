
/**
 * @file ddiag.cc
 *
 * @brief SIMD dynamic banded
 */
#include <string.h>
#include <smmintrin.h>
#include "sse.h"
#include "util.h"

#ifndef BW
#define BW		( 32 )
#endif

#define MIN 	( 0 )
#define OFS 	( 32768 )

/**
 * @fn ddiag_linear
 */
int
ddiag_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt)
//	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	uint16_t *mat = (uint16_t *)aligned_malloc(
		(alen+blen-1) * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	/* extract max and min */
	int8_t sc_max = extract_max_score(score_matrix);
	int8_t sc_min = extract_min_score(score_matrix);

	struct _w {
		int8_t b[BW];
		int8_t pad1[vec::SIZE];
		int8_t a[BW];
		uint16_t pv[BW];
		uint16_t pad2[vec::LEN];
		uint16_t cv[BW];
		uint16_t pad3[vec::LEN];
		uint16_t max[BW];
	} w __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 + i] = 0x80;
		w.b[i] = 0xff;
	}
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 - i - 1] = encode_a(a[i]);
		w.b[BW/2 + i] = encode_b(b[i]);
	}

	/* init vec */
	#define _Q(x)		( (x) - BW/2 )
	for(int i = 0; i < (int)BW; i++) {
		w.pv[i] =      (_Q(i) < 0 ? -_Q(i)   : _Q(i)) * (2*ge - sc_max) + OFS;
		w.cv[i] = ge + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)) * (2*ge - sc_max) + OFS;
		debug("pv(%d), cv(%d)", w.pv[i], w.cv[i]);
	}
	#undef _Q

	/* init pad */
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = -sc_min; w.pad3[i] = -ge;
	}

	/* init maxv */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) {
		w.max[i] = w.pv[i];
	}

	/* direction determiner */
	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = w.pv[BW-1] <= w.pv[0] ? RR : RD;
	debug("%u, %u", w.pv[BW-1], w.pv[0]);

	uint64_t apos = BW/2;
	uint64_t bpos = BW/2;
	uint64_t const L = vec::LEN;
	// vec mv(m), xv(x), giv(-gi);
	vec smv, gev(-ge); smv.load(score_matrix);
	for(uint64_t p = 0; p < (uint64_t)(alen+blen-1); p++) {
		debug("%lld, %d, %d", dir, w.cv[BW-1], w.cv[0]);
		dir = dir_trans[w.cv[BW-1] > w.cv[0]][dir];

		// dump(w.pv, sizeof(uint16_t) * BW);
		// dump(w.cv, sizeof(uint16_t) * BW);

		switch(dir & 0x03) {
			case DD: {
				debug("DD");
				w.pad1[0] = encode_b(b[bpos++]);

				char_vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				vec cd; cd.load(w.pv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va; va.load(&w.a[L*i]);
					char_vec tb; tb.load(&w.b[L*(i+1)]);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(&w.b[L*i]);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					/* load diag */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*(i+1)]);
					vec vd = td.dsr(cd);
					cd = td;

					/* load horizontal and vertical */
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = th.dsr(ch);
					ch = th;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
			case RD: {
				debug("RD");
				w.pad1[0] = encode_b(b[bpos++]);

				char_vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va; va.load(&w.a[L*i]);
					char_vec tb; tb.load(&w.b[L*(i+1)]);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(&w.b[L*i]);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					vec vd; vd.load(&w.pv[L*i]);
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = th.dsr(ch);
					ch = th;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
			case DR: {
				debug("DR");
				w.pad1[vec::SIZE-1] = encode_a(a[apos++]);

				char_vec ca; ca.load(&w.pad1[vec::SIZE/2]);
				vec cv(-ge);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta; ta.load(&w.a[L*i]);
					char_vec va = ta.dsl(ca);
					char_vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					vec vd; vd.load(&w.pv[L*i]);
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = tv.dsl(cv);
					cv = tv;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
			case RR: {
				debug("RR");
				w.pad1[vec::SIZE-1] = encode_a(a[apos++]);

				char_vec ca; ca.load(&w.pad1[vec::SIZE/2]);
				vec cv(-ge);
				vec cd(-sc_min);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta; ta.load(&w.a[L*i]);
					char_vec va = ta.dsl(ca);
					char_vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					/* load diag */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*i]);
					vec vd = td.dsl(cd);
					cd = td;

					/* load horizontal and vertical */
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = tv.dsl(cv);
					cv = tv;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
		}
		ptr += BW;

		if(w.cv[BW/2] < w.max[BW/2] - xt) { break; }
	}
	free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
		vec t; t.load(&w.max[L*i]);
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
	return(max - OFS);
}

/**
 * @fn ddiag_affine
 */
int
ddiag_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t score_matrix[16], int8_t gi, int8_t ge, int16_t xt)
//	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	uint16_t *mat = (uint16_t *)aligned_malloc(
		(alen+blen-1) * 3 * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	/* extract max and min */
	int8_t sc_max = extract_max_score(score_matrix);
	int8_t sc_min = extract_min_score(score_matrix);
	/* fix gap open penalty */
	gi += ge;

	struct _w {
		int8_t b[BW];
		int8_t pad1[vec::SIZE];
		int8_t a[BW];
		uint16_t pv[BW];
		uint16_t pad2[vec::LEN];
		uint16_t cv[BW];
		uint16_t pad3[vec::LEN];
		uint16_t ce[BW];
		uint16_t pad4[vec::LEN];
		uint16_t cf[BW];
		uint16_t max[BW];
	} w __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 + i] = 0x80;
		w.b[i] = 0xff;
	}
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 - i - 1] = encode_a(a[i]);
		w.b[BW/2 + i] = encode_b(b[i]);
	}

	/* init vec */
	#define _Q(x)		( (x) - BW/2 )
	for(int i = 0; i < (int)BW; i++) {
		w.pv[i] =      (_Q(i) < 0 ? -_Q(i)   : _Q(i)) * (2*gi - sc_max) + OFS;
		w.cv[i] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)) * (2*gi - sc_max) + OFS;
		w.ce[i] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)+1) * (2*gi - sc_max) + OFS;
		w.cf[i] = gi + (_Q(i) < 0 ? -_Q(i) : _Q(i)) * (2*gi - sc_max) + OFS;
		debug("pv(%d), cv(%d)", w.pv[i], w.cv[i]);
	}
	#undef _Q

	/* init pad */
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = -sc_min;
		w.pad3[i] = -gi; w.pad4[i] = -ge;
	}

	/* init maxv */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) {
		w.max[i] = w.pv[i];
	}

	/* direction determiner */
	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = w.pv[BW-1] <= w.pv[0] ? RR : RD;

	uint64_t apos = BW/2;
	uint64_t bpos = BW/2;
	uint64_t const L = vec::LEN;
	// vec mv(m), xv(x), giv(-gi), gev(-ge);
	vec smv, giv(-gi), gev(-ge); smv.load(score_matrix);
	for(uint64_t p = 0; p < (uint64_t)(alen+blen-1); p++) {
		debug("%lld, %d, %d", dir, w.cv[BW-1], w.cv[0]);
		dir = dir_trans[w.cv[BW-1] > w.cv[0]][dir];

		// dump(w.pv, sizeof(uint16_t) * BW);
		// dump(w.cv, sizeof(uint16_t) * BW);
		// dump(w.ce, sizeof(uint16_t) * BW);
		// dump(w.cf, sizeof(uint16_t) * BW);

		switch(dir & 0x03) {
			case DD: {
				debug("DD");
				w.pad1[0] = encode_b(b[bpos++]);

				char_vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				vec ce; ce.load(w.ce);
				vec cd; cd.load(w.pv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va; va.load(&w.a[L*i]);
					char_vec tb; tb.load(&w.b[L*(i+1)]);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(&w.b[L*i]);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*(i+1)]);
					vec vd = td.dsr(cd);
					cd = td;

					/* load v and h */
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = th.dsr(ch);
					ch = th;

					/* load f and e */
					vec te; te.load(&w.ce[L*(i+1)]);
					vec vf; vf.load(&w.cf[L*i]);
					vec ve = te.dsr(ce);
					ce = te;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]); ne.print();
					nf.store(&w.cf[L*i]); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]); nv.print();
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
			case RD: {
				debug("RD");
				w.pad1[0] = encode_b(b[bpos++]);

				char_vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				vec ce; ce.load(w.ce);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va; va.load(&w.a[L*i]);
					char_vec tb; tb.load(&w.b[L*(i+1)]);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(&w.b[L*i]);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					vec vd; vd.load(&w.pv[L*i]);

					/* load v and h */
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = th.dsr(ch);
					ch = th;

					/* load f and e */
					vec te; te.load(&w.ce[L*(i+1)]);
					vec vf; vf.load(&w.cf[L*i]);
					vec ve = te.dsr(ce);
					ce = te;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]); ne.print();
					nf.store(&w.cf[L*i]); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]); nv.print();
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
			case DR: {
				debug("DR");
				w.pad1[vec::SIZE-1] = encode_a(a[apos++]);

				char_vec ca; ca.load(&w.pad1[vec::SIZE/2]);
				vec cv(-gi);
				vec cf(-ge);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta; ta.load(&w.a[L*i]);
					char_vec va = ta.dsl(ca);
					char_vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					vec vd; vd.load(&w.pv[L*i]);

					/* load v and h */
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = tv.dsl(cv);
					cv = tv;

					/* load f and e */
					vec ve; ve.load(&w.ce[L*i]);
					vec tf; tf.load(&w.cf[L*i]);
					vec vf = tf.dsl(cf);
					cf = tf;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]); ne.print();
					nf.store(&w.cf[L*i]); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]); nv.print();
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
			case RR: {
				debug("RR");
				w.pad1[vec::SIZE-1] = encode_a(a[apos++]);

				char_vec ca; ca.load(&w.pad1[vec::SIZE/2]);
				vec cv(-gi);
				vec cf(-ge);
				vec cd(-sc_min);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta; ta.load(&w.a[L*i]);
					char_vec va = ta.dsl(ca);
					char_vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*i]);
					vec vd = td.dsl(cd);
					cd = td;

					/* load v and h */
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = tv.dsl(cv);
					cv = tv;

					/* load f and e */
					vec ve; ve.load(&w.ce[L*i]);
					vec tf; tf.load(&w.cf[L*i]);
					vec vf = tf.dsl(cf);
					cf = tf;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]); ne.print();
					nf.store(&w.cf[L*i]); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]); nv.print();
					nv.store(&ptr[L*i]);

					vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
					t.store(&w.max[L*i]);
				}
			} break;
		}
		ptr += BW;

		if(w.cv[BW/2] < w.max[BW/2] - xt) { break; }
	}
	free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
		vec t; t.load(&w.max[L*i]);
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
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

	if(strcmp(argv[1], "linear") == 0) {
		int score = ddiag_linear(
			a, alen, b, blen,
			score_matrix,
			atoi(argv[6]),
			atoi(argv[7]));
		printf("%d\n", score);
	} else if(strcmp(argv[1], "affine") == 0) {
		int score = ddiag_affine(
			a, alen, b, blen,
			score_matrix,
			atoi(argv[6]),
			atoi(argv[7]),
			atoi(argv[8]));
		printf("%d\n", score);
	} else {
		printf("./a.out linear AAA AAA 2 -3 -5 -1 30\n");
	}

	free(a); free(b);
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

	#define l(s, p, q) { \
		assert(ddiag_linear(p, strlen(p), q, strlen(q), score_matrix, -1, 10) == (s)); \
	}
	l( 0, "", "");
	l( 0, "A", "");
	l( 1, "A", "A");
	l( 3, "AAA", "AAA");
	l( 0, "AAA", "TTT");
	l( 3, "AAAGGG", "AAATTTTTT");
	l( 3, "TTTGGGGGAAAA", "TTTCCCCCCCCAAAA");
	l( 5, "AAACAAAGGG", "AAAAAATTTTTTT");
	l( 4, "AAACCAAAGGG", "AAAAAATTTTTTT");

	#define a(s, p, q) { \
		assert(ddiag_affine(p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10) == (s)); \
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

	int sl = ddiag_linear(a, strlen(a), b, strlen(b), score_matrix, -1, 30);
	printf("%d\n", sl);

	int sa = ddiag_affine(a, strlen(a), b, strlen(b), score_matrix, -1, -1, 30);
	printf("%d\n", sa);

	return(0);
}
#endif

/**
 * end of ddiag.cc
 */
