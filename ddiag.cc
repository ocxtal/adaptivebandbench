
/**
 * @file ddiag.cc
 *
 * @brief SIMD dynamic banded
 */
#include <string.h>
#include <smmintrin.h>
#include "sse.h"
#include "util.h"

#define BW		( 32 )
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
	int8_t m, int8_t x, int8_t gi, int8_t ge)
{
	uint16_t *mat = (uint16_t *)aligned_malloc(
		(MIN3(2*alen, alen+blen, 2*blen)) * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	debug("%lld, %lld",
		MIN3(2*alen, alen+blen, 2*blen),
		MIN3(2*alen, alen+blen, 2*blen) * BW * sizeof(uint16_t));

	struct _w {
		uint16_t b[BW];
		uint16_t pad1[8];
		uint16_t a[BW];
		uint16_t pv[BW];
		uint16_t pad2[8];
		uint16_t cv[BW];
		uint16_t pad3[8];
	} w __attribute__(( aligned(16) ));
	uint16_t maxv[BW] __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 + i] = 0x80;
		w.b[i] = 0xff;
	}
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 - i - 1] = a[i];
		w.b[BW/2 + i] = b[i];
	}

	/* init vec */
	#define _Q(x)		( (x) - BW/2 )
	for(int i = 0; i < (int)BW; i++) {
		w.pv[i] =      (_Q(i) < 0 ? -_Q(i)   : _Q(i)) * (2*gi - m) + OFS;
		w.cv[i] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)) * (2*gi - m) + OFS;
		debug("pv(%d), cv(%d)", w.pv[i], w.cv[i]);
	}
	#undef _Q

	/* init pad */
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = -x; w.pad3[i] = 0;
	}

	/* init maxv */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) {
		maxv[i] = 0;
	}

	/* direction determiner */
	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = w.pv[BW-1] <= w.pv[0] ? RR : RD;
	debug("%u, %u", w.pv[BW-1], w.pv[0]);

	uint64_t apos = BW/2;
	uint64_t bpos = BW/2;
	uint64_t const L = vec::LEN;
	vec mv(m), xv(x), giv(-gi);
	for(uint64_t p = 0; p < (uint64_t)MIN3(2*alen, alen+blen, 2*blen)-1; p++) {
		debug("%lld, %d, %d", dir, w.cv[BW-1], w.cv[0]);
		dir = dir_trans[w.cv[BW-1] > w.cv[0]][dir];
		switch(dir & 0x03) {
			case DD: {
				debug("DD");
				w.pad1[0] = b[bpos++];

				vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				vec cd; cd.load(w.pv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec va; va.load(&w.a[L*i]);
					vec tb; tb.load(&w.b[L*(i+1)]);
					vec vb = (tb<<7) | (cb>>1);
					cb = tb; vb.store(&w.b[L*i]);

					vec scv = vec::comp(va, vb).select(mv, xv);

					/* load diag */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*(i+1)]);
					vec vd = (td<<7) | (cd>>1);
					cd = td;

					/* load horizontal and vertical */
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = (th<<7) | (ch>>1);
					ch = th;

					vec nv = vec::max(vec::max(vh, vv) - giv, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
			case RD: {
				debug("RD");
				w.pad1[0] = b[bpos++];

				vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec va; va.load(&w.a[L*i]);
					vec tb; tb.load(&w.b[L*(i+1)]);
					vec vb = (tb<<7) | (cb>>1);
					cb = tb; vb.store(&w.b[L*i]);

					vec scv = vec::comp(va, vb).select(mv, xv);

					vec vd; vd.load(&w.pv[L*i]);
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = (th<<7) | (ch>>1);
					ch = th;

					vec nv = vec::max(vec::max(vh, vv) - giv, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
			case DR: {
				debug("DR");
				w.pad1[7] = a[apos++];

				vec ca; ca.load(w.pad1);
				vec cv; cv.zero();
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec ta; ta.load(&w.a[L*i]);
					vec va = (ta<<1) | (ca>>7);
					vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					vec scv = vec::comp(va, vb).select(m, x);

					vec vd; vd.load(&w.pv[L*i]);
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = (tv<<1) | (cv>>7);
					cv = tv;

					vec nv = vec::max(vec::max(vh, vv) - giv, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
			case RR: {
				debug("RR");
				w.pad1[7] = a[apos++];

				vec ca; ca.load(w.pad1);
				vec cv; cv.zero();
				vec cd(-x);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec ta; ta.load(&w.a[L*i]);
					vec va = (ta<<1) | (ca>>7);
					vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					vec scv = vec::comp(va, vb).select(m, x);

					/* load diag */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*i]);
					vec vd = (td<<1) | (cd>>7);
					cd = td;

					/* load horizontal and vertical */
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = (tv<<1) | (cv>>7);
					cv = tv;

					vec nv = vec::max(vec::max(vh, vv) - giv, vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
		}
		ptr += BW;
	}
	free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
		vec t; t.load(&maxv[L*i]);
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
	int8_t m, int8_t x, int8_t gi, int8_t ge)
{
	uint16_t *mat = (uint16_t *)aligned_malloc(
		(MIN3(2*alen, alen+blen, 2*blen)) * 3 * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	struct _w {
		uint16_t b[BW];
		uint16_t pad1[8];
		uint16_t a[BW];
		uint16_t pv[BW];
		uint16_t pad2[8];
		uint16_t cv[BW];
		uint16_t pad3[8];
		uint16_t ce[BW];
		uint16_t pad4[8];
		uint16_t cf[BW];
	} w __attribute__(( aligned(16) ));
	uint16_t maxv[BW] __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 + i] = 0x80;
		w.b[i] = 0xff;
	}
	for(uint64_t i = 0; i < (uint64_t)BW/2; i++) {
		w.a[BW/2 - i - 1] = a[i];
		w.b[BW/2 + i] = b[i];
	}

	/* init vec */
	#define _Q(x)		( (x) - BW/2 )
	for(int i = 0; i < (int)BW; i++) {
		w.pv[i] =      (_Q(i) < 0 ? -_Q(i)   : _Q(i)) * (2*gi - m) + OFS;
		w.cv[i] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)) * (2*gi - m) + OFS;
		w.ce[i] = gi + (_Q(i) < 0 ? -_Q(i)-1 : _Q(i)+1) * (2*gi - m) + OFS;
		w.cf[i] = gi + (_Q(i) < 0 ? -_Q(i) : _Q(i)) * (2*gi - m) + OFS;
		debug("pv(%d), cv(%d)", w.pv[i], w.cv[i]);
	}
	#undef _Q

	/* init pad */
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = 0;
		w.pad3[i] = 0; w.pad4[i] = 0;
	}

	/* init maxv */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) {
		maxv[i] = 0;
	}

	/* direction determiner */
	uint64_t const RR = 0, RD = 1, DR = 2, DD = 3;
	uint64_t const dir_trans[2][4] = {{RR, DR, RR, DR}, {RD, DD, RD, DD}};
	uint64_t dir = w.pv[BW-1] <= w.pv[0] ? RR : RD;

	uint64_t apos = BW/2;
	uint64_t bpos = BW/2;
	uint64_t const L = vec::LEN;
	vec mv(m), xv(x), giv(-gi), gev(-ge);
	for(uint64_t p = 0; p < (uint64_t)MIN3(2*alen, alen+blen, 2*blen)-1; p++) {
		debug("%lld, %d, %d", dir, w.cv[BW-1], w.cv[0]);
		dir = dir_trans[w.cv[BW-1] > w.cv[0]][dir];
		switch(dir & 0x03) {
			case DD: {
				debug("DD");
				w.pad1[0] = b[bpos++];

				vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				vec ce; ce.load(w.ce);
				vec cd; cd.load(w.pv);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec va; va.load(&w.a[L*i]);
					vec tb; tb.load(&w.b[L*(i+1)]);
					vec vb = (tb<<7) | (cb>>1);
					cb = tb; vb.store(&w.b[L*i]);

					vec scv = vec::comp(va, vb).select(mv, xv);

					/* load pv */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*(i+1)]);
					vec vd = (td<<7) | (cd>>1);
					cd = td;

					/* load v and h */
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = (th<<7) | (ch>>1);
					ch = th;

					/* load f and e */
					vec te; te.load(&w.ce[L*(i+1)]);
					vec vf; vf.load(&w.cf[L*i]);
					vec ve = (te<<7) | (ce>>1);
					ce = te;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]);
					nf.store(&w.cf[L*i]);

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
			case RD: {
				debug("RD");
				w.pad1[0] = b[bpos++];

				vec cb; cb.load(w.b);
				vec ch; ch.load(w.cv);
				vec ce; ce.load(w.ce);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec va; va.load(&w.a[L*i]);
					vec tb; tb.load(&w.b[L*(i+1)]);
					vec vb = (tb<<7) | (cb>>1);
					cb = tb; vb.store(&w.b[L*i]);

					vec scv = vec::comp(va, vb).select(mv, xv);

					/* load pv */
					vec vd; vd.load(&w.pv[L*i]);

					/* load v and h */
					vec th; th.load(&w.cv[L*(i+1)]);
					vec vv = ch; ch.store(&w.pv[L*i]);
					vec vh = (th<<7) | (ch>>1);
					ch = th;

					/* load f and e */
					vec te; te.load(&w.ce[L*(i+1)]);
					vec vf; vf.load(&w.cf[L*i]);
					vec ve = (te<<7) | (ce>>1);
					ce = te;

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]);
					nf.store(&w.cf[L*i]);

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
			case DR: {
				debug("DR");
				w.pad1[7] = a[apos++];

				vec ca; ca.load(w.pad1);
				vec cv; cv.zero();
				vec cf; cf.zero();
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec ta; ta.load(&w.a[L*i]);
					vec va = (ta<<1) | (ca>>7);
					vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					vec scv = vec::comp(va, vb).select(m, x);

					/* load pv */
					vec vd; vd.load(&w.pv[L*i]);

					/* load v and h */
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = (tv<<1) | (cv>>7);
					cv = tv;

					/* load f and e */
					vec ve; ve.load(&w.ce[L*i]);
					vec tf; tf.load(&w.cf[L*i]);
					vec vf = (tf<<1) | (cf>>7);

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]);
					nf.store(&w.cf[L*i]);

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
			case RR: {
				debug("RR");
				w.pad1[7] = a[apos++];

				vec ca; ca.load(w.pad1);
				vec cv; cv.zero();
				vec cf; cf.zero();
				vec cd(-x);
				for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
					debug("loop: %d", i);
					vec ta; ta.load(&w.a[L*i]);
					vec va = (ta<<1) | (ca>>7);
					vec vb; vb.load(&w.b[L*i]);
					ca = ta; va.store(&w.a[L*i]);

					vec scv = vec::comp(va, vb).select(m, x);

					/* load pv */
					// vec vd; vd.load(&w.pv[L*i]);
					vec td; td.load(&w.pv[L*i]);
					vec vd = (td<<1) | (cv>>7);
					cd = td;

					/* load v and h */
					vec tv; tv.load(&w.cv[L*i]);
					vec vh = tv; tv.store(&w.pv[L*i]);
					vec vv = (tv<<1) | (cv>>7);
					cv = tv;

					/* load f and e */
					vec ve; ve.load(&w.ce[L*i]);
					vec tf; tf.load(&w.cf[L*i]);
					vec vf = (tf<<1) | (cf>>7);

					/* update e and f */
					vec ne = vec::max(vh - giv, ve - gev);
					vec nf = vec::max(vv - giv, vf - gev);
					ne.store(&w.ce[L*i]);
					nf.store(&w.cf[L*i]);

					/* update s */
					vec nv = vec::max(vec::max(ne, nf), vd + scv);
					nv.store(&w.cv[L*i]);
					nv.store(&ptr[L*i]);

					vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
					t.store(&maxv[L*i]);
				}
			} break;
		}
	}
	free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(BW / L); i++) {
		vec t; t.load(&maxv[L*i]);
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
	return(max - OFS);
}

#ifdef MAIN
#include <stdlib.h>
int main_ext(int argc, char *argv[])
{
	if(strcmp(argv[1], "linear") == 0) {
		int score = ddiag_linear(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]));
		printf("%d\n", score);
	} else if(strcmp(argv[1], "affine") == 0) {
		int score = ddiag_affine(
			argv[2], strlen(argv[2]),
			argv[3], strlen(argv[3]),
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]));
		printf("%d\n", score);
	} else {
		printf("./a.out linear AAA AAA 2 -3 -5 -1 30\n");
	}
	return(0);
}

int main(int argc, char *argv[])
{
	char const *a = "aabbcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	if(argc > 1) { return(main_ext(argc, argv)); }

	int sl = ddiag_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	printf("%d\n", sl);

	int sa = ddiag_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	printf("%d\n", sa);

	return(0);
}
#endif

/**
 * end of ddiag.cc
 */
