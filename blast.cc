
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
	int8_t *score_matrix, int8_t _gi, int8_t _ge, int16_t _xt,
	uint32_t bw)		/* unused */
{
	if(alen == 0 || blen == 0) { return(0); }
	debug("%s, %s", a, b);

	int32_t gi = _gi, ge = _ge, xt = _xt, best_score = 0;
	uint64_t first_a_index = 0, last_a_index = alen + 1;		/* [first_a_index, last_a_index) */
	uint64_t amax = 0, bmax = 0;

	struct _dp { int16_t s, e, f; };
	struct _dp *ptr = (struct _dp *)work + sizeof(maxpos_t), *prev;

	/* initialize the top row */
	ptr[0].s = 0; ptr[0].e = gi; ptr[0].f = gi;
	for(uint64_t i = 0; i < alen; i++) {
		if(ptr[i].s + ge < -xt) { last_a_index = i + 1; break; }
		ptr[i + 1].s = ptr[i].e + ge;
		ptr[i + 1].e = ptr[i].e + ge;
		ptr[i + 1].f = MIN;
	}

	#define _sc(_i, _bch)		( score_matrix[encode_a(a[_i]) | _bch] )

	ptr[last_a_index].s = best_score;
	ptr[last_a_index].e = last_a_index - first_a_index;
	ptr[last_a_index].f = first_a_index;
	prev = ptr;
	for(uint64_t b_index = 0; b_index < blen; b_index++) {
		int8_t bch = encode_b(b[b_index]);
		debug("b_index(%llu), ch(%d), a_range(%llu, %llu)", b_index, bch, first_a_index, last_a_index);
		prev = ptr;

		int32_t e = MIN;
		int32_t f = MAX2(MIN, MAX2(prev[first_a_index].f, prev[first_a_index].s + gi) + ge);
		int32_t s = f;
		debug("initial, a_index(%llu), score(%d, %d, %d)", first_a_index, s, e, f);
		while(s < best_score - xt) {
			e = MAX2(e, s + gi) + ge;
			f = MAX2(prev[first_a_index + 1].f, prev[first_a_index + 1].s + gi) + ge;
			s = MAX3(prev[first_a_index].s + _sc(first_a_index, bch), e, f);
			first_a_index++;
			debug("forward head, a_index(%llu), ch(%x), score(%d, %d, %d)", first_a_index, encode_a(a[first_a_index]) | bch, s, e, f);
		}

		ptr += last_a_index - first_a_index;
		debug("ptr(%p, %p), a_range(%llu, %llu)", prev, ptr, first_a_index, last_a_index);
		ptr[first_a_index].s = s;
		ptr[first_a_index].e = e;
		ptr[first_a_index].f = f;

		uint64_t next_last_a_index = first_a_index;
		for(uint64_t a_index = first_a_index + 1; a_index < last_a_index; a_index++) {
			e = MAX2(e, s + gi) + ge;
			f = MAX2(prev[a_index].f, prev[a_index].s + gi) + ge;
			s = MAX3(prev[a_index - 1].s + _sc(a_index - 1, bch), e, f);
			if(s > best_score) { best_score = s; amax = a_index; bmax = b_index + 1; }
			if(s >= best_score - xt) { next_last_a_index = a_index + 1; }

			ptr[a_index].s = s;
			ptr[a_index].e = e;
			ptr[a_index].f = f;
			debug("fill, a_index(%llu, %p, %p), ch(%x), score(%d, %d, %d)", a_index, &ptr[a_index].s, &prev[a_index].s, encode_a(a[a_index - 1]) | bch, s, e, f);
		}
		last_a_index = next_last_a_index;

		if(last_a_index <= alen) {
			int32_t d = prev[last_a_index - 1].s + _sc(last_a_index - 1, bch);
			if(d > best_score) { best_score = d; amax = last_a_index; bmax = b_index + 1; }
			while(last_a_index <= alen) {
				e = MAX2(e, s + gi) + ge;
				f = MIN;
				s = MAX2(d, e);
				if(s < best_score - xt) { debug("xdrop failed, s(%d, %d)", s, best_score - xt); break; }

				d = MIN;
				ptr[last_a_index].s = s;
				ptr[last_a_index].e = e;
				ptr[last_a_index].f = f;
				last_a_index++;
				debug("forward tail, a_index(%llu, %p), score(%d, %d, %d)", last_a_index, &ptr[last_a_index - 1].s, s, e, f);
			}
		}
	}
	/* save the maxpos */
	maxpos_t *r = (maxpos_t *)work;
	r->alen = alen;
	r->blen = blen;
	r->apos = amax;
	r->bpos = bmax;

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
