
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
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt)
//	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	// uint16_t *mat = (uint16_t *)aligned_malloc(
		// (alen+blen-1) * BW * sizeof(uint16_t),
		// sizeof(__m128i));
	uint16_t *ptr = (uint16_t *)work;

	/* extract max and min */
	int8_t sc_max = extract_max_score(score_matrix);
	int8_t sc_min = extract_min_score(score_matrix);

	uint64_t const L = vec::LEN;
	struct _w {
		int8_t b[vec::LEN];
		int8_t a[vec::LEN];
		uint16_t pv[vec::LEN];
		uint16_t cv[vec::LEN];
		uint16_t max[vec::LEN];
	} w[BW/L+1] __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w[(BW/2 + i)/L].a[(BW/2 + i)%L] = 0x80;
		w[i/L].b[i%L] = 0xff;
	}
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w[(BW/2 - i - 1)/L].a[(BW/2 - i - 1)%L] = encode_a(a[i]);
		w[(BW/2 + i)/L].b[(BW/2 + i)%L] = encode_b(b[i]);
	}

	/* init vec */
	#define _Q(x)		( (x) - BW/2 )
	for(int i = 0; i < (int)BW; i++) {
		w[i/L].pv[i%L] =      (_Q(i) < 0 ? -_Q(i)   : _Q(i)) * (2*ge - sc_max) + OFS;
		w[i/L].cv[i%L] = ge + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)) * (2*ge - sc_max) + OFS;
		debug("pv(%d), cv(%d)", w[i/L].pv[i%L], w[i/L].cv[i%L]);
	}
	#undef _Q

	/* init pad */
	/*
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = -sc_min; w.pad3[i] = -ge;
	}
	*/
	for(uint64_t i = 0; i < 8; i++) {
		w[BW/L].b[i] = 0;
		w[BW/L].a[i] = 0;
		w[BW/L].pv[i] = -sc_min;
		w[BW/L].cv[i] = -ge;
	}

	/* init maxv */
	for(uint64_t i = 0; i < (uint64_t)BW/L; i++) {
		vec t(w[i].pv);
		t.store(w[i].max);
	}

	/* direction determiner */
	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = w[BW/L-1].pv[L-1] <= w[0].pv[0] ? RR : RD;
	debug("%u, %u", w[BW/L].pv[L-1], w[0].pv[0]);

	uint64_t apos = BW/2;
	uint64_t bpos = BW/2;
	// vec mv(m), xv(x), giv(-gi);
	vec smv, gev(-ge); smv.load(score_matrix);
	for(uint64_t p = 0; p < (uint64_t)(alen+blen-1); p++) {
		debug("%lld, %d, %d", dir, w[BW/L-1].cv[L-1], w[0].cv[0]);
		dir = dir_trans[w[BW/L-1].cv[L-1] > w[0].cv[0]][dir];

		// dump(w.pv, sizeof(uint16_t) * BW);
		// dump(w.cv, sizeof(uint16_t) * BW);

		switch(dir & 0x03) {
			case DD: {
				debug("DD");
				w[BW/L].b[0] = encode_b(b[bpos++]);

				char_vec cb(w[0].b);
				vec ch(w[0].cv);
				vec cd(w[0].pv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va(w[i].a);
					char_vec tb(w[i+1].b);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(w[i].b);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					/* load diag */
					// vec vd(w[i].pv);
					vec td(w[i+1].pv);
					vec vd = td.dsr(cd);
					cd = td;

					/* load horizontal and vertical */
					vec th(w[i+1].cv);
					vec vv = ch; ch.store(w[i].pv);
					vec vh = th.dsr(ch);
					ch = th;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(w[i].cv);
					nv.store(ptr);

					vec t(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
			case RD: {
				debug("RD");
				w[BW/L].b[0] = encode_b(b[bpos++]);

				char_vec cb(w[0].b);
				vec ch(w[0].cv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va(w[i].a);
					char_vec tb(w[i+1].b);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(w[i].b);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					vec vd(w[i].pv);
					vec th(w[i+1].cv);
					vec vv = ch; ch.store(w[i].pv);
					vec vh = th.dsr(ch);
					ch = th;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(w[i].cv);
					nv.store(ptr);

					vec t(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
			case DR: {
				debug("DR");

				char_vec ca(encode_a(a[apos++]));
				vec cv(-ge);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta(w[i].a);
					char_vec va = ta.dsl(ca);
					char_vec vb(w[i].b);
					ca = ta; va.store(w[i].a);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					vec vd(w[i].pv);
					vec tv(w[i].cv);
					vec vh = tv; tv.store(w[i].pv);
					vec vv = tv.dsl(cv);
					cv = tv;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(w[i].cv);
					nv.store(ptr);

					vec t(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
			case RR: {
				debug("RR");

				char_vec ca(encode_a(a[apos++]));
				vec cv(-ge);
				vec cd(-sc_min);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta(w[i].a);
					char_vec va = ta.dsl(ca);
					char_vec vb(w[i].b);
					ca = ta; va.store(w[i].a);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					/* load diag */
					// vec vd(w[i].pv);
					vec td(w[i].pv);
					vec vd = td.dsl(cd);
					cd = td;

					/* load horizontal and vertical */
					vec tv(w[i].cv);
					vec vh = tv; tv.store(w[i].pv);
					vec vv = tv.dsl(cv);
					cv = tv;

					vec nv = vec::max(vec::max(vh, vv) - gev, vd + scv);
					nv.store(w[i].cv);
					nv.store(ptr);

					vec t(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
		}
		ptr += BW;

		if(w[BW/2/L].cv[0] < w[BW/2/L].max[0] - xt) { break; }
	}
	// free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
		vec t(w[i].max);
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
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t score_matrix[16], int8_t gi, int8_t ge, int16_t xt)
//	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	// uint16_t *mat = (uint16_t *)aligned_malloc(
		// (alen+blen-1) * 3 * BW * sizeof(uint16_t),
		// sizeof(__m128i));
	uint16_t *ptr = (uint16_t *)work;

	/* extract max and min */
	int8_t sc_max = extract_max_score(score_matrix);
	int8_t sc_min = extract_min_score(score_matrix);
	/* fix gap open penalty */
	gi += ge;

	/*
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
	*/
	uint64_t const L = vec::LEN;
	struct _w {
		int8_t b[vec::LEN];
		int8_t a[vec::LEN];
		uint16_t pv[vec::LEN];
		uint16_t cv[vec::LEN];
		uint16_t ce[vec::LEN];
		uint16_t cf[vec::LEN];
		uint16_t max[vec::LEN];
	} w[BW/L+1] __attribute__(( aligned(16) ));


	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w[(BW/2 + i)/L].a[(BW/2 + i)%L] = 0x80;
		w[i/L].b[i%L] = 0xff;
	}
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w[(BW/2 - i - 1)/L].a[(BW/2 - i - 1)%L] = encode_a(a[i]);
		w[(BW/2 + i)/L].b[(BW/2 + i)%L] = encode_b(b[i]);
	}

	/* init vec */
	#define _Q(x)		( (x) - BW/2 )
	for(int i = 0; i < (int)BW; i++) {
		w[i/L].pv[i%L] =      (_Q(i) < 0 ? -_Q(i)   : _Q(i)) * (2*gi - sc_max) + OFS;
		w[i/L].cv[i%L] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)) * (2*gi - sc_max) + OFS;
		w[i/L].ce[i%L] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)+1) * (2*gi - sc_max) + OFS;
		w[i/L].cf[i%L] = gi + (_Q(i) < 0 ? -_Q(i) : _Q(i)) * (2*gi - sc_max) + OFS;
		debug("pv(%d), cv(%d)", w[i/L].pv[i%L], w[i/L].cv[i%L]);
	}
	#undef _Q

	/* init pad */
	/*
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = -sc_min;
		w.pad3[i] = -gi; w.pad4[i] = -ge;
	}
	*/
	for(uint64_t i = 0; i < L; i++) {
		w[BW/L].b[i] = 0;
		w[BW/L].a[i] = 0;
		w[BW/L].pv[i] = -sc_min;
		w[BW/L].cv[i] = -gi;
		w[BW/L].ce[i] = -ge;
		w[BW/L].cf[i] = -ge;
		w[BW/L].max[i] = 0;
	}

	/* init maxv */
	for(uint64_t i = 0; i < (uint64_t)BW/L; i++) {
		vec t(w[i].pv);
		t.store(w[i].max);
	}

	/* direction determiner */
	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = w[BW/L-1].pv[L-1] <= w[0].pv[0] ? RR : RD;

	uint64_t apos = BW/2;
	uint64_t bpos = BW/2;
	// vec mv(m), xv(x), giv(-gi), gev(-ge);
	vec smv, giv(-gi), gev(-ge); smv.load(score_matrix);
	for(uint64_t p = 0; p < (uint64_t)(alen+blen-1); p++) {
		debug("%lld, %d, %d", dir, w[BW/L-1].cv[L-1], w[0].cv[0]);
		dir = dir_trans[w[BW/L-1].cv[L-1] > w[0].cv[0]][dir];

		// dump(w.pv, sizeof(uint16_t) * BW);
		// dump(w.cv, sizeof(uint16_t) * BW);
		// dump(w.ce, sizeof(uint16_t) * BW);
		// dump(w.cf, sizeof(uint16_t) * BW);

		switch(dir & 0x03) {
			case DD: {
				debug("DD");
				w[BW/L].b[0] = encode_b(b[bpos++]);

				char_vec cb(w[0].b);
				vec ch(w[0].cv);
				vec ce(w[0].ce);
				vec cd(w[0].pv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va(w[i].a);
					char_vec tb(w[i+1].b);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(w[i].b);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td(w[i+1].pv);
					vec vd = td.dsr(cd);
					cd = td;

					/* load v and h */
					vec th(w[i+1].cv);
					vec vv = ch; ch.store(w[i].pv);
					vec vh = th.dsr(ch);
					ch = th;

					/* load f and e */
					vec te(w[i+1].ce);
					vec vf(w[i].cf);
					vec ve = te.dsr(ce);
					ce = te;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(w[i].ce); ne.print();
					nf.store(w[i].cf); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(w[i].cv); nv.print();
					nv.store(&ptr[L*i]);

					vec t; t.load(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
			case RD: {
				debug("RD");
				w[BW/L].b[0] = encode_b(b[bpos++]);

				char_vec cb(w[0].b);
				vec ch(w[0].cv);
				vec ce(w[0].ce);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec va(w[i].a);
					char_vec tb(w[i+1].b);
					char_vec vb = tb.dsr(cb);
					cb = tb; vb.store(w[i].b);

					// vec scv = vec::comp(va, vb).select(mv, xv);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					vec vd(w[i].pv);

					/* load v and h */
					vec th(w[i+1].cv);
					vec vv = ch; ch.store(w[i].pv);
					vec vh = th.dsr(ch);
					ch = th;

					/* load f and e */
					vec te(w[i+1].ce);
					vec vf(w[i].cf);
					vec ve = te.dsr(ce);
					ce = te;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(w[i].ce); ne.print();
					nf.store(w[i].cf); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(w[i].cv); nv.print();
					nv.store(&ptr[L*i]);

					vec t; t.load(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
			case DR: {
				debug("DR");
				char_vec ca(encode_a(a[apos++]));
				vec cv(-gi);
				vec cf(-ge);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta(w[i].a);
					char_vec va = ta.dsl(ca);
					char_vec vb(w[i].b);
					ca = ta; va.store(w[i].a);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					vec vd(w[i].pv);

					/* load v and h */
					vec tv(w[i].cv);
					vec vh = tv; tv.store(w[i].pv);
					vec vv = tv.dsl(cv);
					cv = tv;

					/* load f and e */
					vec ve(w[i].ce);
					vec tf(w[i].cf);
					vec vf = tf.dsl(cf);
					cf = tf;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(w[i].ce); ne.print();
					nf.store(w[i].cf); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(w[i].cv); nv.print();
					nv.store(&ptr[L*i]);

					vec t(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
			case RR: {
				debug("RR");

				char_vec ca(encode_a(a[apos++]));
				vec cv(-gi);
				vec cf(-ge);
				vec cd(-sc_min);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %llu", i);
					char_vec ta(w[i].a);
					char_vec va = ta.dsl(ca);
					char_vec vb(w[i].b);
					ca = ta; va.store(w[i].a);

					// vec scv = vec::comp(va, vb).select(m, x);
					vec scv = smv.shuffle(va | vb);

					/* load pv */
					// vec vd(w.pv[L*i]);
					vec td(w[i].pv);
					vec vd = td.dsl(cd);
					cd = td;

					/* load v and h */
					vec tv(w[i].cv);
					vec vh = tv; tv.store(w[i].pv);
					vec vv = tv.dsl(cv);
					cv = tv;

					/* load f and e */
					vec ve(w[i].ce);
					vec tf(w[i].cf);
					vec vf = tf.dsl(cf);
					cf = tf;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(w[i].ce); ne.print();
					nf.store(w[i].cf); nf.print();

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(w[i].cv); nv.print();
					nv.store(&ptr[L*i]);

					vec t(w[i].max); t = vec::max(t, nv);
					t.store(w[i].max);
				}
			} break;
		}
		ptr += BW;

		if(w[BW/2/L].cv[0] < w[BW/2/L].max[0] - xt) {
			debug("xdrop");
			break;
		}
	}
	// free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
		vec t(w[i].max);
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


	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	if(strcmp(argv[1], "linear") == 0) {
		int score = ddiag_linear(
			work,
			a, alen, b, blen,
			score_matrix,
			atoi(argv[6]),
			atoi(argv[7]));
		printf("%d\n", score);
	} else if(strcmp(argv[1], "affine") == 0) {
		int score = ddiag_affine(
			work,
			a, alen, b, blen,
			score_matrix,
			atoi(argv[6]),
			atoi(argv[7]),
			atoi(argv[8]));
		printf("%d\n", score);
	} else {
		printf("./a.out linear AAA AAA 2 -3 -5 -1 30\n");
	}

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

	#define l(s, p, q) { \
		assert(ddiag_linear(work, p, strlen(p), q, strlen(q), score_matrix, -1, 10) == (s)); \
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
		assert(ddiag_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10) == (s)); \
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

	int sl = ddiag_linear(work, a, strlen(a), b, strlen(b), score_matrix, -1, 30);
	printf("%d\n", sl);

	int sa = ddiag_affine(work, a, strlen(a), b, strlen(b), score_matrix, -1, -1, 30);
	printf("%d\n", sa);


	free(work);
	return(0);
}
#endif

/**
 * end of ddiag.cc
 */
