
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

#define DEBUG2

/**
 * @class char_vec
 */
typedef uint64_t char_vec_reg;

class char_vec {

private:
	uint64_t v;

public:
	/* consts */
	static int8_t const MAX = 127;
	static int8_t const MIN = -128;
	static uint64_t const SIZE = sizeof(uint64_t);
	static uint64_t const LEN = sizeof(uint64_t);

	/* constructors */
	char_vec(void) {
		set(0);
	}
	char_vec(int8_t k) {
		set(k);
	}
	char_vec(uint64_t i) {
		v = i;
	}
	char_vec(uint64_t const *p) {
		v = *p;
	}
	char_vec(int8_t const *p) {
		v = *((uint64_t *)p);
	}

	/* setter */
	inline void zero(void) {
		v = 0;
	}
	inline void set(int8_t k) {
		v = ((uint64_t)k<<8) | (uint64_t)k;
		v = (v<<16) | v;
		v = (v<<32) | v;
	}

	/* getter */
	inline uint64_t const &get(void) const { return(v); }

	/* assign */
	inline char_vec operator=(char_vec const &b) {
		return(char_vec(v = b.get()));
	}

	/* and */
	inline char_vec operator&(char_vec const &b) const {
		return(char_vec(v & b.get()));
	}
	/* or */
	inline char_vec operator|(char_vec const &b) const {
		return(char_vec(v | b.get()));
	}
	/* shift left */
	inline char_vec operator<<(int s) const {
		switch(s) {
		#define l(n) c(n+1) c(n+2) c(n+3) c(n+4) c(n+5) c(n+6) c(n+7)
			case 0: return(char_vec(v));
		#define c(n) case n: return(char_vec(v<<(8*(n))));
			l(0);
		#undef c
			default: return(char_vec((uint64_t)0));
		#undef l
		}
	}
	inline char_vec operator>>(int s) const {
		switch(s) {
		#define l(n) c(n+1) c(n+2) c(n+3) c(n+4) c(n+5) c(n+6) c(n+7)
			case 0: return(char_vec(v));
		#define c(n) case n: return(char_vec(v>>(8*(n))));
			l(0);
		#undef c
			default: return(char_vec((uint64_t)0));
		#undef l
		}
	}
	/* double shift: (a<<7) | (b>>1) */
	inline char_vec dsr(char_vec const &b) const {
		return(char_vec((v<<56) | (b.get()>>8)));
	}
	/* double shift: (a<<1) | (b>>7) */
	inline char_vec dsl(char_vec const &b) const {
		return(char_vec((v<<8) | (b.get()>>56)));
	}
	/* binary assign */
	inline char_vec operator&=(char_vec const &b) { return(operator=(operator&(b))); }
	inline char_vec operator|=(char_vec const &b) { return(operator=(operator|(b))); }
	inline char_vec operator<<=(int s) { return(operator=(operator<<(s))); }
	inline char_vec operator>>=(int s) { return(operator=(operator>>(s))); }
	/* array */
	inline int8_t operator[](uint64_t i) const {
		if(i < 8) {
			return(operator>>(i).get() & 0xff);
		} else {
			return(0);
		}
	}
	inline int8_t lsb(void) const { return(operator[](0)); }
	inline int8_t center(void) const { return(operator[](4)); }
	inline int8_t msb(void) const { return(operator[](7)); }
	inline int8_t ins(int8_t k, uint64_t i) {
		if(i < 8) {
			v = (v & ~(0xff<<(8*i))) | (k<<(8*i));
		} else {
			return(0);
		}
		return(k);
	}

	/* load and store */
	inline void load(void const *ptr) {
		v = *((uint64_t *)ptr);
	}
	inline void loadu(void const *ptr) {
		load(ptr);
	}
	inline void load_encode_a(void const *ptr, uint64_t len) {
		load(ptr);
		len = len > 8 ? 8 : len; len--;
		v = (0x0303030303030303>>8*(7 - len)) & ((v>>1) ^ (v>>2));
		v |= 0x0909090909090900<<8*len;
	}
	inline void load_encode_b(void const *ptr, uint64_t len) {
		load(ptr);
		len = len > 8 ? 8 : len; len--;
		v = (0x0c0c0c0c0c0c0c0c>>8*(7 - len)) & ((v<<1) ^ v);
		v |= 0x0909090909090900<<8*len;
	}
	inline void store(void *ptr) const {
		*((uint64_t *)ptr) = v;
	}
	inline void storeu(void *ptr) const {
		store(ptr);
	}
	/* print */
	#ifdef DEBUG
	void print(void) const {
		print(stderr, NULL);
	}
	void print(char const *msg) const {
		print(stderr, msg);
	}
	void print(FILE *fp, char const *msg) const {
		fprintf(fp, "%s%s[%016llx]\n", msg == NULL ? "" : msg, msg == NULL ? "" : " ", v);
	}
	#else
	void print(void) const {}
	void print(char const *msg) const {}
	void print(FILE *fp, char const *msg) const {}
	#endif
};

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
	vec(void) {
		set(0);
	}
	vec(uint16_t k) {
		set(k);
	}
	vec(__m128i i) {
		v = i;
	}
	vec(__m128i const *p) {
		v = _mm_load_si128(p);
	}
	vec(uint16_t const *p) {
		v = _mm_load_si128((__m128i const *)p);
	}

	/* setter */
	inline void zero(void) {
		v = _mm_setzero_si128();
	}
	inline void set(int16_t k) {
		v = _mm_set1_epi16(k);
	}

	/* getter */
	inline __m128i const &get(void) const { return(v); }

	/* assign */
	inline vec operator=(vec const &b) {
		return(vec(v = b.get()));
	}

	/* add */
	inline vec operator+(vec const &b) const {
		return(vec(_mm_add_epi16(v, b.get())));
	}
	/* sub */
	inline vec operator-(vec const &b) const {
		return(vec(_mm_subs_epu16(v, b.get())));
	}
	/* and */
	inline vec operator&(vec const &b) const {
		return(vec(_mm_and_si128(v, b.get())));
	}
	/* or */
	inline vec operator|(vec const &b) const {
		return(vec(_mm_or_si128(v, b.get())));
	}
	/* compare */
	inline uint16_t operator<(vec const &b) const {
		__m128i _ofs = _mm_set1_epi16(32768);
		__m128i _a = _mm_sub_epi16(v, _ofs);
		__m128i _b = _mm_sub_epi16(b.get(), _ofs);
		return(_mm_movemask_epi8(_mm_cmplt_epi16(_a, _b)));
	}
	inline uint16_t operator>(vec const &b) const {
		__m128i _ofs = _mm_set1_epi16(32768);
		__m128i _a = _mm_sub_epi16(v, _ofs);
		__m128i _b = _mm_sub_epi16(b.get(), _ofs);
		return(_mm_movemask_epi8(_mm_cmpgt_epi16(_a, _b)));
	}
	inline uint16_t operator<=(vec const &b) const { return(~operator>(b)); }
	inline uint16_t operator>=(vec const &b) const { return(~operator<(b)); }
	inline uint16_t operator==(vec const &b) const {
		return(_mm_movemask_epi8(_mm_cmpeq_epi16(v, b.get())));
	}
	inline uint16_t operator!=(vec const &b) const {
		return(~_mm_movemask_epi8(_mm_cmpeq_epi16(v, b.get())));
	}
	/* shift left */
	inline vec operator<<(int s) const {
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
	inline vec operator>>(int s) const {
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
	/* double shift: (a<<7) | (b>>1) */
	inline vec dsr(vec const &b) const {
		return(vec(_mm_alignr_epi8(v, b.get(), 2)));
	}
	/* double shift: (a<<1) | (b>>7) */
	inline vec dsl(vec const &b) const {
		return(vec(_mm_alignr_epi8(v, b.get(), 14)));
	}
	/* binary assign */
	inline vec operator+=(vec const &b) { return(operator=(operator+(b))); }
	inline vec operator-=(vec const &b) { return(operator=(operator-(b))); }
	inline vec operator&=(vec const &b) { return(operator=(operator&(b))); }
	inline vec operator|=(vec const &b) { return(operator=(operator|(b))); }
	inline vec operator<<=(int s) { return(operator=(operator<<(s))); }
	inline vec operator>>=(int s) { return(operator=(operator>>(s))); }
	/* array */
	inline uint16_t operator[](uint64_t const i) const {
		#ifdef DEBUG2
			switch(i) {
				case 0: return((uint16_t)_mm_extract_epi16(v, 0));
				case 1: return((uint16_t)_mm_extract_epi16(v, 1));
				case 2: return((uint16_t)_mm_extract_epi16(v, 2));
				case 3: return((uint16_t)_mm_extract_epi16(v, 3));
				case 4: return((uint16_t)_mm_extract_epi16(v, 4));
				case 5: return((uint16_t)_mm_extract_epi16(v, 5));
				case 6: return((uint16_t)_mm_extract_epi16(v, 6));
				case 7: return((uint16_t)_mm_extract_epi16(v, 7));
				default: return(0);
			}
		#else
			if(i < 8) {
				return((uint16_t)_mm_extract_epi16(v, i));
			} else {
				return(0);
			}
		#endif
	}
	inline uint16_t lsb(void) const { return(operator[](0)); }
	inline uint16_t center(void) const { return(operator[](4)); }
	inline uint16_t msb(void) const { return(operator[](7)); }
	inline uint16_t ins(uint16_t k, uint64_t const i) {
		#ifdef DEBUG2
			switch(i) {
				case 0: v = _mm_insert_epi16(v, k, 0); return(k);
				case 1: v = _mm_insert_epi16(v, k, 1); return(k);
				case 2: v = _mm_insert_epi16(v, k, 2); return(k);
				case 3: v = _mm_insert_epi16(v, k, 3); return(k);
				case 4: v = _mm_insert_epi16(v, k, 4); return(k);
				case 5: v = _mm_insert_epi16(v, k, 5); return(k);
				case 6: v = _mm_insert_epi16(v, k, 6); return(k);
				case 7: v = _mm_insert_epi16(v, k, 7); return(k);
				default: return(0);
			}
		#else
			if(i < 8) {
				v = _mm_insert_epi16(v, k, i);
			} else {
				return(0);
			}
		#endif
		return(k);
	}

	/* compare and select */
	inline vec static comp(vec const &a, vec const &b) {
		return(vec(_mm_cmpeq_epi16(a.get(), b.get())));
	}
	inline vec select(uint16_t m, uint16_t x) const {
		__m128i mv = _mm_set1_epi16(m);
		__m128i xv = _mm_set1_epi16(x);
		return(vec(_mm_blendv_epi8(xv, mv, v)));
	}
	inline vec select(vec const &m, vec const &x) const {
		return(vec(_mm_blendv_epi8(x.get(), m.get(), v)));
	}
	inline vec shuffle(char_vec const &a) const {
		__m128i index = _mm_cvtsi64_si128(a.get());
		return(vec(_mm_cvtepi8_epi16(_mm_shuffle_epi8(v, index))));
	}
	/* make mask */
	inline int64_t mask(void) const {
		return(_mm_movemask_epi8(v));
	}
	/* max */
	inline vec static max(vec const &a, vec const &b) {
		return(vec(_mm_max_epu16(a.get(), b.get())));
	}
	/* horizontal max */
	inline uint16_t hmax(void) const {
		__m128i t = _mm_max_epu16(v, _mm_srli_si128(v, 2));
		t = _mm_max_epu16(t, _mm_srli_si128(t, 4));
		t = _mm_max_epu16(t, _mm_srli_si128(t, 8));
		return((uint16_t)_mm_extract_epi16(t, 0));
	}
	/* load and store */
	inline void load(void const *ptr) {
		v = _mm_load_si128((__m128i *)ptr);
	}
	inline void loadu(void const *ptr) {
		v = _mm_loadu_si128((__m128i *)ptr);
	}
	inline void load_expand(void const *ptr) {
		uint64_t a = *((uint64_t *)ptr);
		v = _mm_cvtepu8_epi16(_mm_cvtsi64_si128(a));
	}
	inline void store(void *ptr) const {
		_mm_store_si128((__m128i *)ptr, v);
	}
	inline void storeu(void *ptr) const {
		_mm_storeu_si128((__m128i *)ptr, v);
	}
	/* print */
	#ifdef DEBUG
	void print(void) const {
		print(stderr, NULL);
	}
	void print(char const *msg) const {
		print(stderr, msg);
	}
	void print(FILE *fp, char const *msg) const {
		uint16_t b[8] __attribute__(( aligned (16) ));
		store(b);
		fprintf(fp,
			"%s%s[%d, %d, %d, %d, %d, %d, %d, %d]\n",
			msg == NULL ? "" : msg, msg == NULL ? "" : " ",
			b[7] - 32768, b[6] - 32768, b[5] - 32768, b[4] - 32768, b[3] - 32768, b[2] - 32768, b[1] - 32768, b[0] - 32768);
	}
	#else
	void print(void) const {}
	void print(char const *msg) const {}
	void print(FILE *fp, char const *msg) const {}
	#endif
};
/**
 * end of SSE4.1 16bit 8cell
 */

#endif
/**
 * end of sse.h
 */
