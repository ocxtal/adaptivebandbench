
/**
 * @file rognes.cc
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
 * @fn rognes_linear
 */
int
rognes_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	uint16_t *mat = (uint16_t *)aligned_malloc(
		alen * 2 * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	struct _w {
		uint16_t b[2*BW];
		uint16_t pad1[8];
		uint16_t pv[2*BW];
		uint16_t pad2[8];
		uint16_t max[2*BW];
	} w __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { w.b[i] = 0; }
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { w.b[i+BW] = b[i]; }

	/* init vec */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { w.pv[i] = -BW*gi; }
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { w.pv[i+BW] = i * gi + OFS; }

	/* init pad */
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = 0;
	}

	/* init maxv */
	for(uint64_t i = 0; i < 2*BW; i++) {
		w.max[i] = w.pv[i];
	}

	vec mv(m), xv(x), giv(-gi);

	int16_t gsum[BW*3] __attribute__(( aligned (sizeof(__m128i)) ));
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { gsum[i] = -gi; }
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { gsum[i+BW] = -2*gi; }
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { gsum[i+2*BW] = -4*gi; }

	/* fill-in loop */
	uint64_t const L = vec::LEN;
	for(uint64_t p = 0; p < (uint64_t)alen; p++) {
		vec va; va.set(a[p]);
		w.pad1[0] = b[p+BW];

		debug("p: %llu", p);

		vec cb; cb.load(w.b);
		vec ch; ch.load(w.pv);
		vec cv; cv.zero();
		for(uint64_t i = 0; i < (uint64_t)(2*BW / L); i++) {
			debug("i: %llu", i);
			vec scv = vec::comp(va, cb).select(mv, xv); scv.print();

			/* load char vec */
			vec tb; tb.load(&w.b[L*(i+1)]);
			vec vb = (tb<<7) | (cb>>1);
			cb = tb; vb.store(&w.b[L*i]);

			/* load pv */
			vec th; th.load(&w.pv[L*(i+1)]);
			vec vd = ch;
			vec vh = (th<<7) | (ch>>1);
			ch = th;

			vec nv = vec::max(vec::max(vh - giv, cv>>7), vd + scv); nv.print();
			/* chain vertical */
			vec g1; g1.load(&gsum[0]);
			nv = vec::max((nv - g1)<<1, nv);
			vec g2; g2.load(&gsum[BW]);
			nv = vec::max((nv - g2)<<2, nv);
			vec g4; g4.load(&gsum[2*BW]);
			nv = vec::max((nv - g4)<<4, nv);
			cv = nv - giv;

			nv.store(&w.pv[L*i]); nv.print();
			nv.store(&ptr[L*i]);

			vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
			t.store(&w.max[L*i]);
		}
		ptr += 2*BW;

		if(w.pv[BW] < w.max[BW] - xt) { break; }
	}
	free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(2*BW / L); i++) {
		vec t; t.load(&w.max[L*i]);
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
	return(max - OFS);
}

/**
 * @fn rognes_affine
 */
int
rognes_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	/* init vectors */
	uint16_t *mat = (uint16_t *)aligned_malloc(
		alen * 3 * 2 * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	struct _w {
		uint16_t b[2*BW];
		uint16_t pad1[8];
		uint16_t pv[2*BW];
		uint16_t pad2[8];
		uint16_t pe[2*BW];
		uint16_t pad3[8];
		uint16_t max[2*BW];
	} w __attribute__(( aligned(16) ));

	/* init char vec */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { w.b[i] = 0; }
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { w.b[i+BW] = b[i]; }

	/* init vec */
	for(uint64_t i = 0; i < (uint64_t)BW; i++) {
		w.pv[i] = (BW-i) * (2*gi - m) + OFS;
		w.pe[i] = (BW-i) * (2*gi - m) + OFS;
	}
	for(uint64_t i = 0; i < (uint64_t)BW; i++) {
		w.pv[i+BW] = i * gi + OFS;
		w.pe[i+BW] = (i+1) * gi + OFS;
	}

	/* init pad */
	for(uint64_t i = 0; i < 8; i++) {
		w.pad1[i] = 0; w.pad2[i] = 0; w.pad3[i] = 0;
	}

	/* init maxv */
	for(uint64_t i = 0; i < 2*BW; i++) {
		w.max[i] = w.pv[i];
	}

	vec mv(m), xv(x), giv(-gi), gev(-ge);

	int16_t gsum[BW*3] __attribute__(( aligned (sizeof(__m128i)) ));
	for(uint64_t i = 0; i < (uint64_t)BW; i++) { gsum[i] = -3*ge; }

	/* fill-in loop */
	uint64_t const L = vec::LEN;
	for(uint64_t p = 0; p < (uint64_t)alen; p++) {
		vec va; va.set(a[p]);
		w.pad1[0] = b[p+BW];

		debug("p: %llu", p);

		vec cb; cb.load(w.b);
		vec ch; ch.load(w.pv);
		vec ce; ce.load(w.pe);
		vec cv; cv.zero();
		vec cf; cf.zero();
		for(uint64_t i = 0; i < (uint64_t)(2*BW / L); i++) {
			debug("i: %llu", i);
			vec scv = vec::comp(va, cb).select(mv, xv); scv.print();

			/* load char vec */
			vec tb; tb.load(&w.b[L*(i+1)]);
			vec vb = (tb<<7) | (cb>>1);
			cb = tb; vb.store(&w.b[L*i]);

			/* load pv */
			vec th; th.load(&w.pv[L*(i+1)]);
			vec vd = ch;
			vec vh = (th<<7) | (ch>>1);
			ch = th;

			/* load pe */
			vec te; te.load(&w.pe[L*(i+1)]);
			vec ve = (te<<7) | (ce>>1);
			ce = te;
			vec ne = vec::max(vh - giv, ve - gev);
			ne.store(&w.pe[L*i]);
			ne.store(&ptr[L*i + 2*BW]);

			// vec nv = vec::max(vec::max(vh - giv, cv>>7), vd + scv);
			vec nv = vec::max(ne, vd + scv); nv.print();
			/* chain vertical */
			vec nf = vec::max(cf, cv)>>7;
			nv = vec::max(nf, nv);
			nf = vec::max(nv - giv, nf - gev)<<1;
			nv = vec::max(nf, nv);
			nf = vec::max(nv - giv - gev, nf - gev - gev)<<2;
			nv = vec::max(nf, nv);
			vec g4; g4.load(&gsum[0]);
			nf = vec::max(nv - giv - g4, nf - gev - g4)<<4;
			nv = vec::max(nf, nv); nv.print();

			cv = nv - giv;
			cf = nf - gev;

			nv.store(&w.pv[L*i]);
			nv.store(&ptr[L*i]);
			nf.store(&ptr[L*i + 2*2*BW]);

			vec t; t.load(&w.max[L*i]); t = vec::max(t, nv);
			t.store(&w.max[L*i]);
		}
		ptr += 3*2*BW;
		debug("%p", ptr);

		if(w.pv[BW] < w.max[BW] - xt) { break; }
	}
	free(mat);

	int32_t max = 0;
	for(uint64_t i = 0; i < (uint64_t)(2*BW / L); i++) {
		vec t; t.load(&w.max[L*i]);
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
	return(max - OFS);
}

#ifdef MAIN
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

	if(strcmp(argv[1], "linear") == 0) {
		int score = rognes_linear(
			a, alen, b, blen,
			atoi(argv[4]),
			atoi(argv[5]),
			atoi(argv[6]),
			atoi(argv[7]),
			atoi(argv[8]));
		printf("%d\n", score);
	} else if(strcmp(argv[1], "affine") == 0) {
		int score = rognes_affine(
			a, alen, b, blen,
			atoi(argv[4]),
			atoi(argv[5]),
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
	char const *a = "aabbcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	if(argc > 1) { return(main_ext(argc, argv)); }

	int sl = rognes_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1, 30);
	printf("%d\n", sl);

	int sa = rognes_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1, 30);
	printf("%d\n", sa);

	return(0);
}
#endif

/**
 * end of rognes.cc
 */
