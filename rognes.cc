
/**
 * @file rognes.cc
 */
#include <string.h>
#include <smmintrin.h>
#include "sse.h"
#include "util.h"

#define BW		( 8 )
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
	int8_t m, int8_t x, int8_t gi, int8_t ge)
{
	/* init vectors */
	uint16_t *mat = (uint16_t *)aligned_malloc(
		MIN2(alen, blen) * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	struct _w {
		uint16_t b[2*BW];
		uint16_t pad1[8];
		uint16_t pv[2*BW];
		uint16_t pad2[8];
	} w __attribute__(( aligned(16) ));
	uint16_t maxv[2*BW] __attribute__(( aligned(16) ));

	/* init char vec */
	for(int i = 0; i < BW; i++) { w.b[i] = 0; }
	for(int i = 0; i < BW; i++) { w.b[i+BW] = b[i]; }

	/* init vec */
	for(int i = 0; i < BW; i++) { w.pv[i] = (BW-i) * (2*gi - m) + OFS; }
	for(int i = 0; i < BW; i++) { w.pv[i+BW] = i * gi + OFS; }

	/* init maxv */
	for(int i = 0; i < 2*BW; i++) {
		maxv[i] = 0;
	}

	vec mv(m), xv(x), giv(gi);

	int16_t gsum[BW*3] __attribute__(( aligned (sizeof(__m128i)) ));
	for(int i = 0; i < BW; i++) { gsum[i] = gi; }
	for(int i = 0; i < BW; i++) { gsum[i+BW] = 2*gi; }
	for(int i = 0; i < BW; i++) { gsum[i+2*BW] = 4*gi; }

	/* fill-in loop */
	uint64_t const L = vec::LEN;
	for(int p = 0; p < MIN2(alen, blen); p++) {
		vec va; va.set(a[p]);
		w.pad1[0] = b[p+BW];
		w.pad2[0] = MIN - 8*gi;

		printf("p: %d\n", p);

		vec cb; cb.load(w.b);
		vec ch; ch.load(w.pv);
		vec cv; cv.zero();
		for(int i = 0; i < 2*BW / L; i++) {
			printf("i: %d\n", i);
			vec scv = vec::comp(va, cb).select(mv, xv);
			scv.print();

			/* load char vec */
			vec tb; tb.load(&w.b[L*(i+1)]);
			vec vb = (tb<<7) | (cb>>1);
			cb = tb; vb.store(&w.b[L*i]);

			/* load pv */
			vec th; th.load(&w.pv[L*(i+1)]);
			vec vd = ch;
			vec vh = (th<<7) | (ch>>1);
			ch = th;
			vd.print(); vh.print();

			vec nv = vec::max(vec::max(vh + giv, cv>>7), vd + scv); nv.print();
			/* chain vertical */
			vec g1; g1.load(&gsum[0]);
			nv = vec::max((nv + g1)<<1, nv); nv.print();
			vec g2; g2.load(&gsum[BW]);
			nv = vec::max((nv + g2)<<2, nv); nv.print();
			vec g4; g4.load(&gsum[2*BW]);
			nv = vec::max((nv + g4)<<4, nv); nv.print();
			cv = nv + g1;

			nv.store(&w.pv[L*i]);
			nv.store(&ptr[L*i]);

			vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
			t.store(&maxv[L*i]); t.print();
		}
		ptr += L;
	}
	free(mat);

	int32_t max = 0;
	for(int i = 0; i < 2*BW / L; i++) {
		vec t; t.load(&maxv[L*i]);
		t.print();
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
	return(max - OFS);
}

/**
 * @fn rognes_affine
 */
int
rognes_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge)
{
	/* init vectors */
	uint16_t *mat = (uint16_t *)aligned_malloc(
		MIN2(alen, blen) * BW * sizeof(uint16_t),
		sizeof(__m128i));
	uint16_t *ptr = mat;

	struct _w {
		uint16_t b[2*BW];
		uint16_t pad1[8];
		uint16_t pv[2*BW];
		uint16_t pad2[8];
	} w __attribute__(( aligned(16) ));
	uint16_t maxv[2*BW] __attribute__(( aligned(16) ));

	/* init char vec */
	for(int i = 0; i < BW; i++) { w.b[i] = 0; }
	for(int i = 0; i < BW; i++) { w.b[i+BW] = b[i]; }

	/* init vec */
	for(int i = 0; i < BW; i++) { w.pv[i] = (BW-i) * (2*gi - m) + OFS; }
	for(int i = 0; i < BW; i++) { w.pv[i+BW] = i * gi + OFS; }

	/* init maxv */
	for(int i = 0; i < 2*BW; i++) {
		maxv[i] = 0;
	}

	vec mv(m), xv(x), giv(gi);

	int16_t gsum[BW*3] __attribute__(( aligned (sizeof(__m128i)) ));
	for(int i = 0; i < BW; i++) { gsum[i] = gi; }
	for(int i = 0; i < BW; i++) { gsum[i+BW] = 2*gi; }
	for(int i = 0; i < BW; i++) { gsum[i+2*BW] = 4*gi; }

	/* fill-in loop */
	uint64_t const L = vec::LEN;
	for(int p = 0; p < MIN2(alen, blen); p++) {
		vec va; va.set(a[p]);
		w.pad1[0] = b[p+BW];
		w.pad2[0] = MIN - 8*gi;

		printf("p: %d\n", p);

		vec cb; cb.load(w.b);
		vec ch; ch.load(w.pv);
		vec cv; cv.zero();
		for(int i = 0; i < 2*BW / L; i++) {
			printf("i: %d\n", i);
			vec scv = vec::comp(va, cb).select(mv, xv);
			scv.print();

			/* load char vec */
			vec tb; tb.load(&w.b[L*(i+1)]);
			vec vb = (tb<<7) | (cb>>1);
			cb = tb; vb.store(&w.b[L*i]);

			/* load pv */
			vec th; th.load(&w.pv[L*(i+1)]);
			vec vd = ch;
			vec vh = (th<<7) | (ch>>1);
			ch = th;
			vd.print(); vh.print();

			vec nv = vec::max(vec::max(vh + giv, cv>>7), vd + scv); nv.print();
			/* chain vertical */
			vec g1; g1.load(&gsum[0]);
			nv = vec::max((nv + g1)<<1, nv); nv.print();
			vec g2; g2.load(&gsum[BW]);
			nv = vec::max((nv + g2)<<2, nv); nv.print();
			vec g4; g4.load(&gsum[2*BW]);
			nv = vec::max((nv + g4)<<4, nv); nv.print();
			cv = nv + g1;

			nv.store(&w.pv[L*i]);
			nv.store(&ptr[L*i]);

			vec t; t.load(&maxv[L*i]); t = vec::max(t, nv);
			t.store(&maxv[L*i]); t.print();
		}
		ptr += L;
	}
	free(mat);

	int32_t max = 0;
	for(int i = 0; i < 2*BW / L; i++) {
		vec t; t.load(&maxv[L*i]);
		t.print();
		debug("%d", t.hmax());
		if(t.hmax() > max) { max = t.hmax(); }
	}
	return(max - OFS);
}

int main(void)
{
	char const *a = "abcdefpppijkl";
	char const *b = "abefpppghijkl";

	int sl = rognes_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1);
	printf("%d\n", sl);

	return(0);
}

/**
 * end of rognes.cc
 */
