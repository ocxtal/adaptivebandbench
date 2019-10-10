
/*
 * staticband.h
 */

#define _base_signature			void *work, char const *a, uint64_t alen, char const *b, uint64_t blen, int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, uint32_t bw
int scalar_affine(_base_signature);
int vertical_affine_sse(_base_signature);
int vertical_affine_avx(_base_signature);
int diagonal_affine_sse(_base_signature);
int diagonal_affine_avx(_base_signature);
int striped_affine_sse(_base_signature);
int striped_affine_avx(_base_signature);

/*
 * end of staticband.h
 */
