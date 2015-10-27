
/**
 * @file sse.h
 *
 * @brief class implementation
 */
#ifndef _SSE_H_INCLUDED
#define _SSE_H_INCLUDED

#include <smmintrin.h>
#include <stdint.h>
#include <stdio.h>

/**
 * @class vec
 *
 * @brief SSE4.1 16bit 8cell
 */
class vec {

private:
	__m128i v;

public:
	/* consts */
	static uint16_t const MAX = 65535;
	static uint16_t const MIN = 0;
	static uint64_t const SIZE = sizeof(__m128i);
	static uint64_t const LEN = sizeof(__m128i) / sizeof(uint16_t);

	/* constructors */
	vec(void) {}
	vec(uint16_t k) {
		set(k);
	}
	vec(__m128i i) {
		v = i;
	}

	/* setter */
	void zero(void) {
		v = _mm_setzero_si128();
	}
	void set(int16_t k) {
		v = _mm_set1_epi16(k);
	}

	/* getter */
	__m128i const &get(void) const { return(v); }

	/* assign */
	vec operator=(vec const &b) {
		return(vec(v = b.get()));
	}
	/* cast */

	/* add */
	vec operator+(vec const &b) const {
		return(vec(_mm_add_epi16(v, b.get())));
	}
	/* sub */
	vec operator-(vec const &b) const {
		return(vec(_mm_subs_epu16(v, b.get())));
	}
	/* and */
	vec operator&(vec const &b) const {
		return(vec(_mm_and_si128(v, b.get())));
	}
	/* or */
	vec operator|(vec const &b) const {
		return(vec(_mm_or_si128(v, b.get())));
	}
	/* compare */
	uint16_t operator<(vec const &b) const {
		return(_mm_movemask_epi8(_mm_cmplt_epi16(v, b.get())));
	}
	uint16_t operator>(vec const &b) const {
		return(_mm_movemask_epi8(_mm_cmpgt_epi16(v, b.get())));
	}
	uint16_t operator<=(vec const &b) const { return(~operator>(b)); }
	uint16_t operator>=(vec const &b) const { return(~operator<(b)); }
	uint16_t operator==(vec const &b) const {
		return(_mm_movemask_epi8(_mm_cmpeq_epi16(v, b.get())));
	}
	uint16_t operator!=(vec const &b) const {
		return(~_mm_movemask_epi8(_mm_cmpeq_epi16(v, b.get())));
	}
	/* shift left */
	vec operator<<(int s) const {
		switch(s) {
		#define l(n) c(n+1) c(n+2) c(n+3) c(n+4) c(n+5) c(n+6) c(n+7)
			case 0: return(vec(v));
		#define c(n) case n: return(vec(_mm_slli_si128(v, 2*(n))));
			l(0);
		#undef c
			default: return(vec(_mm_setzero_si128()));
		#undef l
		}
	}
	vec operator>>(int s) const {
		switch(s) {
		#define l(n) c(n+1) c(n+2) c(n+3) c(n+4) c(n+5) c(n+6) c(n+7)
			case 0: return(vec(v));
		#define c(n) case n: return(vec(_mm_srli_si128(v, 2*(n))));
			l(0);
		#undef c
			default: return(vec(_mm_setzero_si128()));
		#undef l
		}
	}
	/* binary assign */
	vec operator+=(vec const &b) { return(operator=(operator+(b))); }
	vec operator-=(vec const &b) { return(operator=(operator-(b))); }
	vec operator&=(vec const &b) { return(operator=(operator&(b))); }
	vec operator|=(vec const &b) { return(operator=(operator|(b))); }
	vec operator<<=(int s) { return(operator=(operator<<(s))); }
	vec operator>>=(int s) { return(operator=(operator>>(s))); }
	/* array */
	uint16_t operator[](uint64_t i) const {
		if(i < 8) {
			return((uint16_t)_mm_extract_epi16(v, i));
		} else {
			return(0);
		}
	}
	uint16_t lsb(void) const { return(operator[](0)); }
	uint16_t center(void) const { return(operator[](4)); }
	uint16_t msb(void) const { return(operator[](7)); }
	uint16_t ins(uint16_t k, uint64_t i) {
		if(i < 8) {
			v = _mm_insert_epi16(v, k, i);
		} else {
			return(0);
		}
		return(k);
	}

	/* compare and select */
	vec static comp(vec const &a, vec const &b) {
		return(vec(_mm_cmpeq_epi16(a.get(), b.get())));
	}
	vec select(uint16_t m, uint16_t x) const {
		__m128i mv = _mm_set1_epi16(m);
		__m128i xv = _mm_set1_epi16(x);
		return(vec(_mm_blendv_epi8(xv, mv, v)));
	}
	vec select(vec const &m, vec const &x) const {
		return(vec(_mm_blendv_epi8(x.get(), m.get(), v)));
	}
	/* make mask */
	int64_t mask(void) const {
		return(_mm_movemask_epi8(v));
	}
	/* max */
	vec static max(vec const &a, vec const &b) {
		return(vec(_mm_max_epu16(a.get(), b.get())));
	}
	/* horizontal max */
	uint16_t hmax(void) const {
		__m128i t = _mm_max_epu16(v, _mm_srli_si128(v, 2));
		t = _mm_max_epu16(t, _mm_srli_si128(t, 4));
		t = _mm_max_epu16(t, _mm_srli_si128(t, 8));
		return((uint16_t)_mm_extract_epi16(t, 0));
	}
	/* load and store */
	void load(void const *ptr) {
		v = _mm_load_si128((__m128i *)ptr);
	}
	void loadu(void const *ptr) {
		v = _mm_loadu_si128((__m128i *)ptr);
	}
	void load_expand(void const *ptr) {
		uint64_t a = *((uint64_t *)ptr);
		v = _mm_cvtepu8_epi16(_mm_cvtsi64_si128(a));
	}
	void store(void *ptr) const {
		_mm_store_si128((__m128i *)ptr, v);
	}
	void storeu(void *ptr) const {
		_mm_storeu_si128((__m128i *)ptr, v);
	}
	/* print */
	#ifdef DEBUG
	void print(void) const {
		print(stderr);
	}
	void print(FILE *fp) const {
		uint16_t b[8] __attribute__(( aligned (16) ));
		store(b);
		fprintf(fp,
			"[%04x %04x %04x %04x %04x %04x %04x %04x]\n",
			b[7], b[6], b[5], b[4], b[3], b[2], b[1], b[0]);
	}
	#else
	void print(void) const {}
	void print(FILE *fp) const {}
	#endif
};
/**
 * end of SSE4.1 16bit 8cell
 */

#endif
/**
 * end of sse.h
 */
