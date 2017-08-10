
/**
 * @file ddiag.h
 */
#ifndef _DDIAG_H_INCLUDED
#define _DDIAG_H_INCLUDED
#include <limits.h>
#include "sea.h"

#ifndef _SIMD_H_INCLUDED
#define _SIMD_H_INCLUDED

/**
 * Constants representing algorithms
 *
 * Notice: This constants must be consistent with the sea_flags_alg in sea.h.
 */
#define SW 								( 1 )
#define SEA 							( 2 )
#define XSEA 							( 3 )
#define NW 								( 6 )

/**
 * max and min
 */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z) 	( MAX2(x, MAX2(y, z)) )
#define MAX4(w,x,y,z) 	( MAX2(MAX2(w, x), MAX2(y, z)) )

#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z) 	( MIN2(x, MIN2(y, z)) )
#define MIN4(w,x,y,z) 	( MIN2(MIN2(w, x), MIN2(y, z)) )

/**
 * architecture flag definitions
 */
#define SSE 		( 1 )
#define AVX 		( 2 )

/**
 * @struct pos
 * @brief a struct containing a position
 */
struct pos {
	sea_int_t i, j, p, q;
};

/**
 * @struct mpos
 * @biref contains multiple position, max pos and end pos
 */
struct mpos {
	struct pos m;			/** max score pos */
	struct pos e;			/** end pos */
};

/**
 * bitwidth selection
 */
#if BIT_WIDTH == 8
/*
	#define CELL_TYPE					char
	#define BYTES_PER_CELL				( sizeof(CELL_TYPE) )
	#define CELL_MAX					( SCHAR_MAX )
	#define CELL_MIN					( SCHAR_MIN )
*/
	#define CELL_TYPE					short
	#define BYTES_PER_CELL				( sizeof(CELL_TYPE) )
	#define CELL_MAX					( SHRT_MAX )
	#define CELL_MIN					( SHRT_MIN + 1000 )

#elif BIT_WIDTH == 16
/*
	#define CELL_TYPE					short
	#define BYTES_PER_CELL				( sizeof(CELL_TYPE) )
	#define CELL_MAX					( SHRT_MAX )
	#define CELL_MIN					( SHRT_MIN )
*/
	#define CELL_TYPE					int
	#define BYTES_PER_CELL				( sizeof(CELL_TYPE) )
	#define CELL_MAX					( INT_MAX )
	#define CELL_MIN					( INT_MIN + 1000 )

#else
 	#error "the BIT_WIDTH must be 8 or 16 in diag algorithms."
#endif

/**
 * In the diag algorithms, the bit width of SIMD packed variables and 
 * the bit width of a cell in the memory are always the same.
 */
//#define SIMD_BIT_WIDTH 				BIT_WIDTH
#define SIMD_BIT_WIDTH 				BIT_WIDTH * 2
#define SIMD_BAND_WIDTH 			BAND_WIDTH
#define BYTES_PER_LINE 				( sizeof(CELL_TYPE)*BAND_WIDTH )
#define SCORE(c)					( *((CELL_TYPE *)(c)) )

/**
 * wrapper for the affine-gap cost algorithms
 */
#define ASCOREV(c)					SCORE(c)
#define ASCOREF(c)					SCORE((c) +   BYTES_PER_LINE)
#define ASCOREE(c)					SCORE((c) + 2*BYTES_PER_LINE)

/**
 * score saturation macro
 */
#define SAT(a)		( ((a) > CELL_MAX) ? CELL_MAX \
									   : (((a) < CELL_MIN) \
									   ? CELL_MIN \
									   : (a)) )

/**
 * char vector shift operations
 */
#define PUSHQ(x, y)					{ VEC_CHAR_SHIFT_L(y); VEC_CHAR_INSERT_LSB(y, x); }
#define PUSHT(x, y)					{ VEC_CHAR_SHIFT_R(y); VEC_CHAR_INSERT_MSB(y, x); }

/**
 * vector initialization macro (replace of VEC_INIT_PV)
 */
#define VEC_INIT_PVN(v, m, g, i) { \
	for(i = -bw/2; i < bw/2; i++) { \
		VEC_SHIFT_R(v); \
		VEC_INSERT_MSB(v, \
			i == 0 ? 0 \
				   : i < 0 \
				   ? SAT((g)-i*(2*(g)-(m))) \
				   : SAT((g)+i*(2*(g)-(m)))); \
	} \
}

/**
 * coordinate conversion macros (common for all algorithms)
 */
#define COX(p, q)				( ((p)>>1) - (q) )
#define COY(p, q)				( (((p)+1)>>1) + (q) )
#define COP(x, y)				( (x) + (y) )
#define COQ(x, y) 				( ((y)-(x))>>1 )

#define DIR_V 					( 0x01 )
#define DIR_H 					( 0 )
#define DIR_VV					( 0x03 )
#define DIR_HH					( 0 )

/**
 * address calculation macros for the linear-gap cost algorithms
 */
#define	ADDR(p, q)				( ((BAND_WIDTH)*(p)+(q)+(BAND_WIDTH)/2) * BYTES_PER_CELL )
#define TOPQ(p, q) 				( - !((p)&0x01) * BYTES_PER_CELL )
#define LEFTQ(p, q) 			( ((p)&0x01) * BYTES_PER_CELL )
#define TOP(p, q)				( -(BYTES_PER_LINE) + TOPQ(p, q) )
#define LEFT(p, q)				( -(BYTES_PER_LINE) + LEFTQ(p, q) )
#define TOPLEFT(p, q) 			( -2*(BYTES_PER_LINE) )

#define DTOPQ(dir) 				( - !((dir)&0x01) * BYTES_PER_CELL )
#define DLEFTQ(dir) 			( ((dir)&0x01) * BYTES_PER_CELL )
#define DTOP(dir)				( -(BYTES_PER_LINE) + DTOPQ(dir) )
#define DLEFT(dir)				( -(BYTES_PER_LINE) + DLEFTQ(dir) )
#define DTOPLEFT(dir) 			( DTOP(dir) + DLEFT((dir)>>1) )

/**
 * address calculation macros for the affine-gap cost algorithms
 */
#define	AADDR(p, q)				( (3 * (BAND_WIDTH)*(p)+(q)+(BAND_WIDTH)/2) * BYTES_PER_CELL )
#define ATOPQ(p, q) 			( - !((p)&0x01) * sizeof(CELL_TYPE) )
#define ALEFTQ(p, q) 			( ((p)&0x01) * sizeof(CELL_TYPE) )
#define ATOP(p, q)				( -(3 * BYTES_PER_LINE) + ATOPQ(p, q) )
#define ALEFT(p, q)				( -(3 * BYTES_PER_LINE) + ALEFTQ(p, q) )
#define ATOPLEFT(p, q) 			( -2*(3 * BYTES_PER_LINE) )

#define DATOPQ(dir) 			( - !((dir)&0x01) * sizeof(CELL_TYPE) )
#define DALEFTQ(dir) 			( ((dir)&0x01) * sizeof(CELL_TYPE) )
#define DATOP(dir)				( -(3 * BYTES_PER_LINE) + DATOPQ(dir) )
#define DALEFT(dir)				( -(3 * BYTES_PER_LINE) + DALEFTQ(dir) )
#define DATOPLEFT(dir) 			( DATOP(dir) + DALEFT((dir)>>1) )


/**
 * include SIMD intrinsic macros, depending on the value of BIT_WIDTH.
 */
#if defined(__AVX2__)
	#include "x86_64/avx.h"
#elif defined(__SSE4_1__)
	#include "x86_64/sse.h"
#else
 	#error "unsupported architecture. check definition of the 'ARCH' constant."
#endif

#endif /* _SIMD_H_INCLUDED */

struct mpos
diag_linear_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat);

struct mpos
diag_linear_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diag_linear_dynamic_banded_trace(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diag_linear_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);


struct mpos
diag_affine_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat);

struct mpos
diag_affine_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diag_affine_dynamic_banded_trace(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
diag_affine_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);


#endif
/**
 * end of ddiag.h
 */
