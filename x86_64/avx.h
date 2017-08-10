
/**
 * @file avx.h
 *
 * @brief a top level header for macros for x86_64 AVX2 instructions.
 *
 * @detail
 * Requirements:
 * x86_64 processors, AVX2 instructions.
 * a compiler must have 'immintrin.h'.
 *
 * Description:
 * The AVX2 instructions on x86_64 processors handles the 256-bit wide registers,
 * which can be divided into thirty two 8-bit variables or sixteen 16-bit variables.
 * Therefore, the wider SIMD operations, which is requered in the implementation
 * of the wide-band diag and diff algorithms, must be implemented with a concatenation
 * of two or more AVX2(ymm) registers. This header provides the automatic selection of 
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
 * |            |  unsupported    | single register |
 * |          8 |    use sse4     |   avx_b8_r1.h   |
 * | bit        |    instead      |   diag / diff   |
 * | width      +-----------------+-----------------+
 * |            | single register | 2 regs parallel |
 * |         16 |  avx_b16_r1.h   |  avx_b16_r2.h   |
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
#ifndef _AVX_H_INCLUDED
#define _AVX_H_INCLUDED


#if SIMD_BAND_WIDTH == 16
 	#if SIMD_BIT_WIDTH == 8
 		/** if register width is 128-bit, use SSE4.1 instructions instead */
 		#include "sse_b8_r1.h"
	#elif SIMD_BIT_WIDTH == 16
		#include "avx_b16_r1.h"
	#elif SIMD_BIT_WIDTH == 32
		#include "avx_b32_r2.h"
	#else
 		#error "The SIMD_BIT_WIDTH must be 8 or 16 in the 16-cell wide AVX2 variants."
 	#endif
#elif SIMD_BAND_WIDTH == 32
	#if SIMD_BIT_WIDTH == 8
		#include "avx_b8_r1.h"
	#elif SIMD_BIT_WIDTH == 16
		#include "avx_b16_r2.h"
	#elif SIMD_BIT_WIDTH == 32
		#include "avx_b32_r4.h"
	#else
 		#error "The SIMD_BIT_WIDTH must be 8 or 16 in the 16-cell wide AVX2 variants."
 	#endif
#else 
 	#error "The SIMD_BAND_WIDTH must be 16 or 32 in the AVX2 variants."
#endif


#endif /* #ifndef _AVX_H_INCLUDED */
/**
 * end of avx.h
 */