
/**
 * @file blast.cc
 *
 * @brief reimplementation of the blast_SemiGappedAlign algorithm
 * (original implementation found in <ncbi-blast>/c++/src/algo/blast/core/blast_gapalign.c)
 */
#include <string.h>
#include "util.h"

#define MIN			( -32768 + 10 )

/**
 * @fn blast_linear
 */
int
blast_linear(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	int i, a_size, first_a_index, last_a_index, a_index, b_index;
	int16_t best_score, score, next_score, score_gap_row, score_gap_col;

	int16_t *mat = (int16_t *)malloc(alen * blen * sizeof(int16_t));
	int16_t *ptr = mat, *prev;

	/* initialize top row */
	ptr[0] = (score = 0);
	for(i = 1; i < alen; i++) {
		if(0 - score > xt) { break; }
		ptr[i] = (score += gi);
	}

	a_size = last_a_index = i;
	first_a_index = 0;
	next_score = MIN;
	best_score = 0;

	for(b_index = 1; b_index <= blen; b_index++) {
		prev = ptr; ptr += (last_a_index + 1 - first_a_index);
		last_a_index = first_a_index;

		score = MIN;
		if(first_a_index != 0) {
			score = next_score + (a[0] == b[b_index-1] ? m : x);
		}
		debug("next_score(%d)", next_score);
		next_score = MIN;

		for(a_index = first_a_index; a_index < a_size; a_index++) {

			/* calc f and d */
			int16_t prev_best = prev[a_index];
			score_gap_col = prev_best + gi;
			if(score < score_gap_col) { score = score_gap_col; }

			debug("score(%d)", score);

			/* xdrop test */
			if(best_score - score > xt) {
				ptr[a_index] = MIN;
				if(a_index == first_a_index) {
					next_score = score;
					first_a_index++;
				}
			} else {
				last_a_index = a_index;
				ptr[a_index] = score;
				/* update best_score */
				if(score > best_score) { best_score = score; }
			}

			/* calc next score */
			score_gap_row = score + gi;
			score = prev_best + (a[a_index] == b[b_index-1] ? m : x);
			if(score < score_gap_row) { score = score_gap_row; }

			debug("(%d, %d), row(%d), col(%d), diag(%d)", a_index, b_index, score_gap_row, score_gap_col, score);
		}

		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - 1) {
			a_size = last_a_index + 1;
		} else {
			while((best_score - score <= xt) && a_size <= alen) {
				ptr[a_size] = (score += gi); a_size++;
			}
			last_a_index = a_size;
		}

		if(a_size <= alen) {
			ptr[a_size] = MIN; a_size++;
		}
	}
	free(mat);

	return(best_score);
}

/**
 * @fn blast_affine
 */
int
blast_affine(
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	int i, a_size, first_a_index, last_a_index, a_index, b_index;
	int16_t best_score, score, next_score, score_gap_row, score_gap_col;

	struct _dp {
		int16_t best;
		int16_t best_gap;
	};
	struct _dp *mat = (struct _dp *)malloc(alen * blen * sizeof(struct _dp));
	struct _dp *ptr = mat, *prev;

	/* initialize top row */
	ptr[0].best = (score = 0); ptr[0].best_gap = (gi - ge);
	for(i = 1; i < alen; i++) {
		if(0 - score > xt) { break; }
		ptr[i].best = score; ptr[i].best_gap = score + (gi - ge);
		score += ge;
	}

	a_size = last_a_index = i;
	first_a_index = 0;
	next_score = MIN;
	best_score = 0;

	for(b_index = 1; b_index <= blen; b_index++) {
		prev = ptr; ptr += (last_a_index + 1 - first_a_index);
		last_a_index = first_a_index;

		score = MIN;
		score_gap_row = MIN;
		if(first_a_index != 0) {
			score = next_score + (a[0] == b[b_index-1] ? m : x);
		}
		debug("first_next_score(%d)", next_score);
		next_score = MIN;

		for(a_index = first_a_index; a_index < a_size; a_index++) {

			/* calc f */
			int16_t prev_best = prev[a_index].best;
			score_gap_col = prev[a_index].best_gap + ge;
			if(score_gap_col < prev_best + gi) { score_gap_col = prev_best + gi; }
			ptr[a_index].best_gap = score_gap_col;

			/* update d */
			if(score < score_gap_col) { score = score_gap_col; }

			debug("(%d, %d), row(%d), col(%d), diag(%d)", a_index, b_index, score_gap_row, score_gap_col, score);

			debug("updated score(%d)", score);

			/* xdrop test */
			if(best_score - score > xt) {
				ptr[a_index].best = MIN;
				if(a_index == first_a_index) {
					next_score = score;
					first_a_index++;
				}
			} else {
				last_a_index = a_index;
				ptr[a_index].best = score;
				/* update best score */
				if(score > best_score) { best_score = score; }
			}

			/* calc e */
			score_gap_row += ge;
			if(score_gap_row < score + gi) { score_gap_row = score + gi; }

			/* calc next score */
			score = prev_best + (a[a_index] == b[b_index-1] ? m : x);
			if(score < score_gap_row) { score = score_gap_row; }
		}

		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - 1) {
			a_size = last_a_index + 1;
		} else {
			while((best_score - score <= xt) && a_size <= alen) {
				ptr[a_size].best = (score += gi);
				ptr[a_size].best_gap = (score_gap_row += ge);
				if(score < score_gap_row) { score = score_gap_row; }
				a_size++;
			}
			last_a_index = a_size;
		}

		if(a_size <= alen) {
			ptr[a_size].best = MIN;
			ptr[a_size].best_gap = MIN;
			a_size++;
		}
	}
	free(mat);

	return(best_score);
}

#ifdef MAIN
int main(void)
{
	char const *a = "aabbcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	int sl = blast_linear(a, strlen(a), b, strlen(b), 2, -3, -5, -1, 30);
	printf("%d\n", sl);

	int sa = blast_affine(a, strlen(a), b, strlen(b), 2, -3, -5, -1, 20);
	printf("%d\n", sa);

	return(0);
}
#endif

/**
 * end of blast.cc
 */
