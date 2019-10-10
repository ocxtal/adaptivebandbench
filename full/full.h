
/**
 * @file full.h
 */
#ifndef _FULL_H_INCLUDED
#define _FULL_H_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _sw_result {
	int score;
	uint32_t path_length;
	uint64_t apos, bpos;
	char *path;
} sw_result_t;

sw_result_t sw_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge);
	// int8_t m, int8_t x, int8_t gi, int8_t ge);

sw_result_t sw_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge);
	// int8_t m, int8_t x, int8_t gi, int8_t ge);

#ifdef __cplusplus
}
#endif

#endif
/**
 * end of full.h
 */
