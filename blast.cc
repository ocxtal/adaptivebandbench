
/**
 * @file blast.cc
 *
 * @brief reimplementation of the blast_SemiGappedAlign algorithm
 * (original implementation found in <ncbi-blast>/c++/src/algo/blast/core/blast_gapalign.c)
 */
#include <string.h>
#include "util.h"

#define MIN			( -32768 + 30 )

/**
 * @fn blast_affine
 */
int
blast_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt)
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt)
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	uint64_t i, a_size, first_a_index, last_a_index, a_index, b_index;
	int16_t best_score, score, next_score, score_gap_row, score_gap_col;

	struct _dp {
		int16_t best;
		int16_t best_gap;
	};
	// struct _dp *mat = (struct _dp *)malloc((alen + 2) * (blen + 1) * sizeof(struct _dp));
	struct _dp *ptr = (struct _dp *)work, *prev;

	/* initialize top row */
	ptr[0].best = (score = 0); ptr[0].best_gap = gi;
	score += gi + ge;
	for(i = 1; i <= alen; i++) {
		if(0 - score > xt) { break; }
		ptr[i].best = score; ptr[i].best_gap = score + gi;
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
			// score = next_score + (a[0] == b[b_index-1] ? m : x);
			score = next_score + score_matrix[encode_a(a[0]) | encode_b(b[b_index-1])];
		}
		debug("first_next_score(%d)", next_score);
		next_score = MIN;

		for(a_index = first_a_index; a_index < a_size; a_index++) {

			/* calc f */
			int16_t prev_best = prev[a_index].best;
			score_gap_col = prev[a_index].best_gap + ge;
			if(score_gap_col < prev_best + gi + ge) { score_gap_col = prev_best + gi + ge; }
			ptr[a_index].best_gap = score_gap_col;

			/* update d */
			if(score < score_gap_col) { score = score_gap_col; }

			debug("(%llu, %llu), row(%d), col(%d), diag(%d)", a_index, b_index, score_gap_row, score_gap_col, score);

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
			if(score_gap_row < score + gi + ge) { score_gap_row = score + gi + ge; }

			/* calc next score */
			// score = prev_best + (a[a_index] == b[b_index-1] ? m : x);
			score = prev_best + score_matrix[encode_a(a[a_index]) | encode_b(b[b_index-1])];
			if(score < score_gap_row) { score = score_gap_row; }
		}

		if(first_a_index == a_size) { break; }
		if(last_a_index < a_size - 1) {
			a_size = last_a_index + 1;
		} else {
			while((best_score - score <= xt) && a_size <= alen) {
				ptr[a_size].best = (score += gi + ge);
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
	// free(mat);

	return(best_score);
}

#ifdef MAIN
#include <assert.h>
#include <stdlib.h>
int main_ext(int argc, char *argv[])
{
	uint64_t alen = strlen(argv[2]);
	uint64_t blen = strlen(argv[3]);
	char *a = argv[2];
	char *b = argv[3];

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, atoi(argv[4]), atoi(argv[5]));

	void *work = aligned_malloc(128 * 1024 * 1024, 16);
	
	if(0) {
		printf("./a.out linear AAA AAA 2 -3 -5 -1 30\n");
	}
	int score = blast_affine(
		work,
		a, alen, b, blen,
		score_matrix,
		atoi(argv[6]),
		atoi(argv[7]),
		atoi(argv[8]));
	printf("%d\n", score);
	free(work);
	return(0);
}

int main(int argc, char *argv[])
{
	char const *a = "aattcccccc";
	char const *b = "aacccccc";
	// char const *a = "abefpppbbqqqqghijkltttt";
	// char const *b = "abcdefpppqqqqijklggtttt";

	if(argc > 1) { return(main_ext(argc, argv)); }

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, 1, -1);


	void *work = aligned_malloc(128 * 1024 * 1024, 16);

	#define a(s, p, q) { \
		assert(blast_affine(work, p, strlen(p), q, strlen(q), score_matrix, -1, -1, 10) == (s)); \
	}
	a( 0, "", "");
	a( 0, "A", "");
	a( 1, "A", "A");
	a( 3, "AAA", "AAA");
	a( 0, "AAA", "TTT");
	a( 3, "AAAGGG", "AAATTTTTT");
	a( 3, "TTTGGGGGAAAA", "TTTCCCCCCCCAAAA");
	a( 4, "AAACAAAGGG", "AAAAAATTTTTTT");
	a( 3, "AAACCAAAGGG", "AAAAAATTTTTTT");

	int sa = blast_affine(work, a, strlen(a), b, strlen(b), score_matrix, -1, -1, 30);
	printf("%d\n", sa);

	free(work);
	return(0);
}
#endif

/**
 * end of blast.cc
 */
