
/**
 * @file sse.h
 *
 * @brief a top level header for macros for x86_64 SSE4 instructions.
 *
 * @detail
 * Requirements:
 * x86_64 processors, SSE4 instructions.
 * a compiler must have 'smmintrin.h'.
 *
 * Description:
 * The SSE instructions on x86_64 processors handles the 128-bit wide registers,
 * which can be divided into sixteen 8-bit variables or eight 16-bit variables.
 * Therefore, the wider SIMD operations, which is requered in the implementation
 * of the wide-band diag and diff algorithms, must be implemented with a concatenation
 * of two or more SSE4(xmm) registers. This header provides the automatic selection of 
 * headers which implement concatenaton macros, based on SIMD_BAND_WIDTH and 
 * SIMD_BIT_WIDTH constant given to the file.
 *
 * The selection strategy is shown in the table below. The top row of the cell
 * (e.g. 2 regs parallel) indicates the composition of a band, the middle of the
 * cell is the name of the header, and the bottom is an algorithmic capability of
 * the header.
 * +------------------------------------------------+
 * |            |            band width             |
 * |            |        16       |        32       |
 * +------------+-----------------+-----------------+
 * |            | single register | 2 regs parallel |
 * |          8 |   sse_b8_r1.h   |   sse_b8_r2.h   |
 * | bit        |   diag / diff   |   diag / diff   |
 * | width      +-----------------+-----------------+
 * |            | 2 regs parallel | 4 regs parallel |
 * |         16 |  sse_b16_r2.h   |  sse_b16_r4.h   |
 * |            |   diag / diff   |   diag / diff   |
 * +------------+-----------------+-----------------+
 *
 * Complier flags:
 * To use macros in this header, the SIMD_BAND_WIDTH and the SIMD_BIT_WIDTH
 * constants must be defined. These constants are normally defined in the
 * diag_variant.h or diff_variant.h, based on the BAND_WIDTH and BIT_WIDTH 
 * value given to the files. The details of the BAND_WIDTH and BIT_WIDTH are
 * described in these headers.
 */
#ifndef _SSE_H_INCLUDED
#define _SSE_H_INCLUDED


#if SIMD_BAND_WIDTH == 16
	#if SIMD_BIT_WIDTH == 8
		#include "sse_b8_r1.h"
	#elif SIMD_BIT_WIDTH == 16
		#include "sse_b16_r2.h"
 	#elif SIMD_BIT_WIDTH == 32
 		#include "sse_b32_r4.h"
	#else
 		#error "The SIMD_BIT_WIDTH must be 8 or 16 in the 16-cell wide SSE4 variants."
 	#endif
#elif SIMD_BAND_WIDTH == 32
	#if SIMD_BIT_WIDTH == 8
		#include "sse_b8_r2.h"
	#elif SIMD_BIT_WIDTH == 16
		#include "sse_b16_r4.h"
 	#elif SIMD_BIT_WIDTH == 32
 		#include "sse_b32_r8.h"
	#else
 		#error "The SIMD_BIT_WIDTH must be 8 or 16 in the 32-cell wide SSE4 variants."
 	#endif
#elif SIMD_BAND_WIDTH == 64
	#if SIMD_BIT_WIDTH == 8
		#include "sse_b8_r4.h"
	#elif SIMD_BIT_WIDTH == 16
		#include "sse_b16_r8.h"
	#else
 		#error "The SIMD_BIT_WIDTH must be 8 or 16 in the 64-cell wide SSE4 variants."
 	#endif
#else
 	#error "The SIMD_BAND_WIDTH must be 16, 32, or 64 in the SSE4 variants."
#endif


#endif /* #ifndef _SSE_H_INCLUDED */
/**
 * end of sse.h
 */
