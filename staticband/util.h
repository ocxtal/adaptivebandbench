
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


#define cat_intl(x, y)		x##_##y
#define cat(x, y)			cat_intl(x, y)


#if defined(BENCH_ARCH_SSE)

#define SUFFIX				sse
#define VLEN				( 8 )

typedef struct { __m128i v; } vdp_t;
typedef struct { __m128i v; } vchar_t;		/* lower half */
typedef struct { __m128i v; } vmat_t;

#define zero_vdp()			((vdp_t){ .v = _mm_setzero_si128() })
#define seta_vdp(x)			((vdp_t){ .v = _mm_set1_epi16((uint16_t)(x)) })
#define load_vdp(p)			((vdp_t){ .v = _mm_load_si128((__m128i const *)(p)) })
#define loadu_vdp(p)		((vdp_t){ .v = _mm_loadu_si128((__m128i const *)(p)) })
#define store_vdp(p, x)		{ _mm_store_si128((__m128i *)(p), (x).v); }
#define storeu_vdp(p, x)	{ _mm_storeu_si128((__m128i *)(p), (x).v); }
#define bsl_vdp(x, y)		((vdp_t){ .v = _mm_slli_si128((x).v, 2 * (y)) })
#define bsr_vdp(x, y)		((vdp_t){ .v = _mm_srli_si128((x).v, 2 * (y)) })
#define bsld_vdp(x, y)		((vdp_t){ .v = _mm_alignr_epi8((x).v, (y).v, 14) })
#define bsrd_vdp(x, y)		((vdp_t){ .v = _mm_alignr_epi8((x).v, (y).v, 2) })
#define add_vdp(x, y)		((vdp_t){ .v = _mm_add_epi16((x).v, (y).v) })
#define sub_vdp(x, y)		((vdp_t){ .v = _mm_subs_epu16((x).v, (y).v) })
#define max_vdp(x, y)		((vdp_t){ .v = _mm_max_epu16((x).v, (y).v) })
#define eq_vdp(x, y)		( (uint64_t)_mm_movemask_epi8(_mm_cmpeq_epi16((x).v, (y).v)) )
#define gt_vdp(x, y)		( (uint64_t)_mm_movemask_epi8(_mm_cmpgt_epi16(_mm_sub_epi16((x).v, _mm_set1_epi16(32768)), _mm_sub_epi16((y).v, _mm_set1_epi16(32768)))) )

static inline
uint16_t hmax_vdp(vdp_t v)
{
	__m128i t = v.v;
	t = _mm_max_epu16(t, _mm_srli_si128(t, 8));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 4));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 2));
	return((uint16_t)_mm_extract_epi16(t, 0));
}

#if defined(DEBUG)
#define print_vdp(x) { \
	uint16_t buf[8]; \
	storeu_vdp(buf, x); \
	fprintf(stderr, \
		"%s[%d, %d, %d, %d, %d, %d, %d, %d]\n", #x, \
		buf[7] - 32768, \
		buf[6] - 32768, \
		buf[5] - 32768, \
		buf[4] - 32768, \
		buf[3] - 32768, \
		buf[2] - 32768, \
		buf[1] - 32768, \
		buf[0] - 32768 \
	); \
}
#else
#define print_vdp(x)		;
#endif


#define zero_vchar()		((vchar_t){ .v = _mm_setzero_si128() })
#define seta_vchar(x)		((vchar_t){ .v = _mm_set1_epi16(x) })
#define load_vchar(p)		((vchar_t){ .v = _mm_loadl_epi64((__m128i const *)(p)) })
#define loadu_vchar(p)		((vchar_t){ .v = _mm_loadl_epi64((__m128i const *)(p)) })
#define store_vchar(p, x)	{ _mm_storel_epi64((__m128i *)(p), (x).v); }
#define storeu_vchar(p, x)	{ _mm_storel_epi64((__m128i *)(p), (x).v); }
// #define bsl_vchar(x, y)		((vchar_t){ .v = _mm_slli_si128((x).v, (y)) })
// #define bsr_vchar(x, y)		((vchar_t){ .v = _mm_srli_si128((x).v, (y)) })
#define bsld_vchar(x, y)	((vchar_t){ .v = _mm_srli_si128(_mm_unpack_epi64((y).v, (x).v), 7) })
#define bsrd_vchar(x, y)	((vchar_t){ .v = _mm_srli_si128(_mm_unpack_epi64((y).v, (x).v), 1) })


#define loadu_vmat(p)		((vmat_t){ .v = _mm_loadu_si128((__m128i const *)(p)) })
#define shuffle_vmat(x, y)	((vmat_t){ .v = _mm_shuffle_epi8((x).v, (y).v); })


#define cvt_vchar_vmat(x)	((vmat_t){ .v = (x).v })
#define cvt_vmat_vdp(x)		((vdp_t){ .v = _mm_cvtepi8_epi16((x).v) })




#elif defined(BENCH_ARCH_AVX)

#define SUFFIX				avx
#define VLEN				( 16 )

typedef struct { __m256i v; } vdp_t;
typedef struct { __m128i v; } vchar_t;		/* lower half */
typedef struct { __m128i v; } vmat_t;

#define zero_vdp()			((vdp_t){ .v = _mm256_setzero_si256() })
#define seta_vdp(x)			((vdp_t){ .v = _mm256_set1_epi16((uint16_t)(x)) })
#define load_vdp(p)			((vdp_t){ .v = _mm256_load_si256((__m256i const *)(p)) })
#define loadu_vdp(p)		((vdp_t){ .v = _mm256_loadu_si256((__m256i const *)(p)) })
#define store_vdp(p, x)		{ _mm256_store_si256((__m256i *)(p), (x).v); }
#define storeu_vdp(p, x)	{ _mm256_storeu_si256((__m256i *)(p), (x).v); }
/*
#define bsl_vdp(x, y)		((vdp_t){ .v = _mm256_slli_si256((x).v, 2 * (y)) })
#define bsr_vdp(x, y)		((vdp_t){ .v = _mm256_srli_si256((x).v, 2 * (y)) })
#define bsld_vdp(x, y)		((vdp_t){ .v = _mm256_alignr_epi8((x).v, (y).v, 14) })
#define bsrd_vdp(x, y)		((vdp_t){ .v = _mm256_alignr_epi8((x).v, (y).v, 2) })
*/
#define bsl_vdp(x, imm) ( \
	(imm) >= 8 ? (vdp_t){ .v = _mm256_slli_si256(_mm256_inserti128_si256(_mm256_setzero_si256(), _mm256_castsi256_si128((x).v), 1), 2 * (imm) - 16) } \
	           : (vdp_t){ .v = _mm256_alignr_epi8((x).v, _mm256_permute2x128_si256((x).v, (x).v, 0x08), 16 - 2 * (imm)) } \
)
#define bsr_vdp(x, imm) ( \
	(imm) >= 8 ? (vdp_t){ .v = _mm256_srli_si256(_mm256_castsi128_si256(_mm256_extracti128_si256((x).v, 1)), 2 * (imm) - 16) } \
	           : (vdp_t){ .v = _mm256_alignr_epi8(_mm256_castsi128_si256(_mm256_extracti128_si256((x).v, 1)), (x).v, 2 * (imm)) } \
)

#define bsld_vdp(x, y)		((vdp_t){ .v = _mm256_alignr_epi8((x).v, _mm256_permute2x128_si256((x).v, (y).v, 0x03), 14) })
#define bsrd_vdp(x, y)		((vdp_t){ .v = _mm256_alignr_epi8(_mm256_permute2x128_si256((x).v, (y).v, 0x03), (y).v, 2) })

#define add_vdp(x, y)		((vdp_t){ .v = _mm256_add_epi16((x).v, (y).v) })
#define sub_vdp(x, y)		((vdp_t){ .v = _mm256_subs_epu16((x).v, (y).v) })
#define max_vdp(x, y)		((vdp_t){ .v = _mm256_max_epu16((x).v, (y).v) })
#define eq_vdp(x, y)		( (uint64_t)_mm256_movemask_epi8(_mm256_cmpeq_epi16((x).v, (y).v)) )
#define gt_vdp(x, y)		( (uint64_t)_mm256_movemask_epi8(_mm256_cmpgt_epi16(_mm256_sub_epi16((x).v, _mm256_set1_epi16(32768)), _mm256_sub_epi16((y).v, _mm256_set1_epi16(32768)))) )

static inline
uint16_t hmax_vdp(vdp_t v)
{
	__m128i t = _mm_max_epu16(
		_mm256_castsi256_si128(v.v),
		_mm256_extracti128_si256(v.v, 1)
	);
	t = _mm_max_epu16(t, _mm_srli_si128(t, 8));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 4));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 2));
	return((uint16_t)_mm_extract_epi16(t, 0));
}

#if defined(DEBUG)
#define print_vdp(x) { \
	uint16_t buf[16]; \
	storeu_vdp(buf, x); \
	fprintf(stderr, \
		"%s[%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", #x, \
		buf[8 + 7] - 32768, \
		buf[8 + 6] - 32768, \
		buf[8 + 5] - 32768, \
		buf[8 + 4] - 32768, \
		buf[8 + 3] - 32768, \
		buf[8 + 2] - 32768, \
		buf[8 + 1] - 32768, \
		buf[8 + 0] - 32768, \
		buf[7] - 32768, \
		buf[6] - 32768, \
		buf[5] - 32768, \
		buf[4] - 32768, \
		buf[3] - 32768, \
		buf[2] - 32768, \
		buf[1] - 32768, \
		buf[0] - 32768 \
	); \
}
#else
#define print_vdp(x)		;
#endif


#define zero_vchar()		((vchar_t){ .v = _mm_setzero_si128() })
#define seta_vchar(x)		((vchar_t){ .v = _mm_set1_epi16(x) })
#define load_vchar(p)		((vchar_t){ .v = _mm_load_si128((__m128i const *)(p)) })
#define loadu_vchar(p)		((vchar_t){ .v = _mm_loadu_si128((__m128i const *)(p)) })
#define store_vchar(p, x)	{ _mm_store_si128((__m128i *)(p), (x).v); }
#define storeu_vchar(p, x)	{ _mm_storeu_si128((__m128i *)(p), (x).v); }
// #define bsl_vchar(x, y)		((vchar_t){ .v = _mm_slli_si128((x).v, (y)) })
// #define bsr_vchar(x, y)		((vchar_t){ .v = _mm_srli_si128((x).v, (y)) })
#define bsld_vchar(x, y)	((vchar_t){ .v = _mm_alignr_epi8((x).v, (y).v, 15) })
#define bsrd_vchar(x, y)	((vchar_t){ .v = _mm_alignr_epi8((x).v, (y).v, 1) })


#define loadu_vmat(p)		((vmat_t){ .v = _mm_loadu_si128((__m128i const *)(p)) })
#define shuffle_vmat(x, y)	((vmat_t){ .v = _mm_shuffle_epi8((x).v, (y).v); })


#define cvt_vchar_vmat(x)	((vmat_t){ .v = (x).v })
#define cvt_vmat_vdp(x)		((vdp_t){ .v = _mm256_cvtepi8_epi16((x).v) })




#endif


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
