
/*
 * @file util.h
 *
 * @brief a collection of utilities
 */

#ifndef _UTIL_H_INCLUDED
#define _UTIL_H_INCLUDED

#include <stdio.h>
#include <stdint.h>				/** uint32_t, uint64_t, ... */
#include <stdlib.h>
#include <stddef.h>				/** offsetof */
#include <string.h>				/** memset, memcpy */
#include <x86intrin.h>


/* result container */
typedef struct maxpos_s {
	uint64_t apos, bpos;
	uint64_t alen, blen;
	uint64_t ccnt;				/** #cells calculated */
	uint64_t fcnt;				/** lazy-f count for debugging */
	uint64_t time[2];
	char *path;
	uint64_t path_length;
} maxpos_t;

/**
 * aligned malloc
 */
static inline
void *aligned_malloc(
	size_t size,
	size_t align)
{
	void *ptr = NULL;
	if(posix_memalign(&ptr, align, size) != 0) {
		return(NULL);
	}
	return(ptr);
}

/**
 * misc
 */
static inline
int8_t encode(char a)
{
	return(0x03 & ((a>>1) ^ (a>>2)));
}

static inline
int8_t encode_a(char a)
{
	return(0x03 & ((a>>1) ^ (a>>2)));
}

static inline
int8_t encode_b(char b)
{
	return(0x0c & ((b<<1) ^ b));
}
#define encode_n()		( 0x09 )

static inline
void build_score_matrix(int8_t *matrix, int8_t m, int8_t x)
{
	int i = 0;
	matrix[i++] = m;	// (A, A)
	matrix[i++] = x;	// (C, A)
	matrix[i++] = x;	// (G, A)
	matrix[i++] = x;	// (T, A)
	matrix[i++] = x;	// (A, C)
	matrix[i++] = m;	// (C, C)
	matrix[i++] = x;	// (G, C)
	matrix[i++] = x;	// (T, C)
	matrix[i++] = x;	// (A, G)
	matrix[i++] = x;	// (C, G)
	matrix[i++] = m;	// (G, G)
	matrix[i++] = x;	// (T, G)
	matrix[i++] = x;	// (A, T)
	matrix[i++] = x;	// (C, T)
	matrix[i++] = x;	// (G, T)
	matrix[i++] = m;	// (T, T)
	return;
}

static inline
int8_t extract_max_score(int8_t *matrix)
{
	__m128i a = _mm_load_si128((__m128i *)matrix);
	a = _mm_max_epi8(a, _mm_srli_si128(a, 1));
	a = _mm_max_epi8(a, _mm_srli_si128(a, 2));
	a = _mm_max_epi8(a, _mm_srli_si128(a, 4));
	a = _mm_max_epi8(a, _mm_srli_si128(a, 8));
	return(_mm_extract_epi8(a, 0));
}

static inline
int8_t extract_min_score(int8_t *matrix)
{
	__m128i a = _mm_load_si128((__m128i *)matrix);
	a = _mm_min_epi8(a, _mm_srli_si128(a, 1));
	a = _mm_min_epi8(a, _mm_srli_si128(a, 2));
	a = _mm_min_epi8(a, _mm_srli_si128(a, 4));
	a = _mm_min_epi8(a, _mm_srli_si128(a, 8));
	return(_mm_extract_epi8(a, 0));
}


/**
 * @macro popcnt
 */
#ifdef __POPCNT__
	#define popcnt(x)		_mm_popcnt_u64(x)
#else
	static inline
	int popcnt(uint64_t n)
	{
		uint64_t c = 0;
		c = (n & 0x5555555555555555) + ((n>>1) & 0x5555555555555555);
		c = (c & 0x3333333333333333) + ((c>>2) & 0x3333333333333333);
		c = (c & 0x0f0f0f0f0f0f0f0f) + ((c>>4) & 0x0f0f0f0f0f0f0f0f);
		c = (c & 0x00ff00ff00ff00ff) + ((c>>8) & 0x00ff00ff00ff00ff);
		c = (c & 0x0000ffff0000ffff) + ((c>>16) & 0x0000ffff0000ffff);
		c = (c & 0x00000000ffffffff) + ((c>>32) & 0x00000000ffffffff);
		return(c);
	}
#endif

/**
 * @macro tzcnt
 * @brief trailing zero count (count #continuous zeros from LSb)
 */
#ifdef __BMI__
	/** immintrin.h is already included */
	#define tzcnt(x)		_tzcnt_u64(x)
#else
	static inline
	int tzcnt(uint64_t n)
	{
		n |= n<<1;
		n |= n<<2;
		n |= n<<4;
		n |= n<<8;
		n |= n<<16;
		n |= n<<32;
		return(64-popcnt(n));
	}
#endif

/**
 * @macro lzcnt
 * @brief leading zero count (count #continuous zeros from MSb)
 */
#ifdef __LZCNT__
	#define lzcnt(x)		_lzcnt_u64(x)
#else
	static inline
	int lzcnt(uint64_t n)
	{
		n |= n>>1;
		n |= n>>2;
		n |= n>>4;
		n |= n>>8;
		n |= n>>16;
		n |= n>>32;
		return(64-popcnt(n));
	}
#endif

/**
 * coordinate conversion macros
 */
#define cox(p, q, band)				( ((p)>>1) - (q) )
#define coy(p, q, band)				( (((p)+1)>>1) + (q) )
#define cop(x, y, band)				( (x) + (y) )
#define coq(x, y, band) 			( ((y)-(x))>>1 )

/**
 * max and min
 */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z) 	( MAX2(x, MAX2(y, z)) )
#define MAX4(w,x,y,z) 	( MAX2(MAX2(w, x), MAX2(y, z)) )

#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z) 	( MIN2(x, MIN2(y, z)) )
#define MIN4(w,x,y,z) 	( MIN2(MIN2(w, x), MIN2(y, z)) )


/* split_foreach */
#define mm_split_foreach(_ptr, _delims, _body) { \
	char const *_q = (_ptr); \
	int64_t i = 0; \
	__m128i _dv = _mm_loadu_si128((__m128i const *)_delims); \
	uint16_t _m, _mask = 0x02<<tzcnt(		/* reserve space for '\0' */ \
		_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_set1_epi8('\0'), _dv)) \
	); \
	_dv = _mm_slli_si128(_dv, 1); _mask--;						/* push '\0' at the head of the vector */ \
	do { \
		char const *_p = _q; \
		/* test char one by one until dilimiter found */ \
		while(((_m = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_set1_epi8(*_q), _dv))) & _mask) == 0) { _q++; } \
		/* delimiter found, pass to _body */ \
		char const *p = _p; \
		uint64_t l = _q++ - _p; \
		if(l > 0) { _body; i++; } \
	} while((_m & 0x01) == 0); \
}

static
char *mm_strdup(char const *p)
{
	if(!p) { return(NULL); }
	uint64_t l = 0;
	p--; while(*++p) { l++; }
	char *a = (char *)malloc(l + 1);
	memcpy(a, p - l, l);
	a[l] = '\0';
	return(a);
}

#endif /* #ifndef _UTIL_H_INCLUDED */

/*
 * end of util.h
 */
