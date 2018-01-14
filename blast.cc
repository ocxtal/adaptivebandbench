
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
	uint64_t first_b_index = 0, last_b_index = blen + 1;		/* [first_b_index, last_b_index) */
	uint64_t amax = 0, bmax = 0;

	struct _dp { int16_t s, e, f; };
	struct _dp *ptr = (struct _dp *)work + sizeof(maxpos_t), *prev;

	/* initialize the top row */
	ptr[0].s = 0; ptr[0].e = gi; ptr[0].f = gi;
	for(uint64_t i = 0; i < blen; i++) {
		if(ptr[i].s + ge < -xt) { last_b_index = i + 1; break; }
		ptr[i + 1].s = ptr[i].f + ge;
		ptr[i + 1].e = MIN;
		ptr[i + 1].f = ptr[i].f + ge;
		debug("init, b_index(%llu), score(%d)", i + 1, ptr[i + 1].s);
	}
	ptr[last_b_index].s = best_score;
	ptr[last_b_index].e = last_b_index - first_b_index;
	ptr[last_b_index].f = first_b_index;

	#define _sc(_ach, _i)		( score_matrix[_ach | encode_b(b[_i])] )
	prev = ptr;
	for(uint64_t a_index = 0; a_index < alen; a_index++) {
		prev = ptr; ptr += last_b_index + 1 - first_b_index;
		int8_t ach = encode_a(a[a_index]);
		debug("a_index(%llu), ch(%d), b_range(%llu, %llu)", a_index, ach, first_b_index, last_b_index);

		int32_t e = MAX2(MIN, MAX2(prev[first_b_index].e, prev[first_b_index].s + gi) + ge);
		int32_t f = MIN;
		int32_t s = e;
		debug("initial, b_index(%llu), score(%d, %d, %d)", first_b_index, s, e, f);
		while(s < best_score - xt) {
			e = MAX2(prev[first_b_index + 1].e, prev[first_b_index + 1].s + gi) + ge;
			f = MAX2(f, s + gi) + ge;
			s = MAX3(prev[first_b_index].s + _sc(ach, first_b_index), e, f);
			first_b_index++; ptr--;
			debug("forward head, b_index(%llu), ch(%x), score(%d, %d, %d)", first_b_index, encode_b(b[first_b_index]) | ach, s, e, f);
		}

		debug("ptr(%p, %p), b_range(%llu, %llu)", prev, ptr, first_b_index, last_b_index);
		ptr[first_b_index].s = s;
		ptr[first_b_index].e = e;
		ptr[first_b_index].f = f;

		uint64_t next_last_b_index = last_b_index;
		for(uint64_t b_index = first_b_index + 1; b_index < last_b_index; b_index++) {
			e = MAX2(prev[b_index].e, prev[b_index].s + gi) + ge;
			f = MAX2(f, s + gi) + ge;
			s = MAX3(prev[b_index - 1].s + _sc(ach, b_index - 1), e, f);
			if(s > best_score) { best_score = s; amax = a_index + 1; bmax = b_index; }
			if(s >= best_score - xt) { next_last_b_index = b_index + 1; }

			ptr[b_index].s = s;
			ptr[b_index].e = e;
			ptr[b_index].f = f;
			debug("fill, b_index(%llu, %p, %p), ch(%x), score(%d, %d, %d)", b_index, &ptr[b_index].s, &prev[b_index].s, encode_b(b[b_index - 1]) | ach, s, e, f);
		}
		last_b_index = next_last_b_index;

		if(last_b_index <= blen) {
			int32_t d = prev[last_b_index - 1].s + _sc(ach, last_b_index - 1);
			debug("tail, b_index(%llu, %p, %p), ch(%x), d(%d, %d)", last_b_index, &ptr[last_b_index].s, &prev[last_b_index].s, encode_b(b[last_b_index - 1]) | ach, prev[last_b_index - 1].s, d);

			if(d > best_score) { best_score = d; amax = a_index + 1; bmax = last_b_index; }
			while(last_b_index <= blen) {
				e = MIN;
				f = MAX2(f, s + gi) + ge;
				s = MAX2(d, f);
				if(s < best_score - xt) { debug("xdrop failed, s(%d, %d)", s, best_score - xt); break; }

				d = MIN;
				ptr[last_b_index].s = s;
				ptr[last_b_index].e = e;
				ptr[last_b_index].f = f;
				last_b_index++;
				debug("forward tail, a_index(%llu, %p), score(%d, %d, %d)", last_b_index, &ptr[last_b_index - 1].s, s, e, f);
			}
		}
		ptr[last_b_index].s = best_score;
		ptr[last_b_index].e = last_b_index - first_b_index;
		ptr[last_b_index].f = first_b_index;
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
