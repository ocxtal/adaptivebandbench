
/**
 * @file diag_linear_dynamic_banded.c
 *
 * @brief a linear-gap cost, banded SIMD implementation.
 *
 * @detail
 * This implementation is a vectorized variant of naive_linear_banded.c.
 * The consistency of the results are tested in the test_cross_diag_linear_dynamic_banded
 * function.
 *
 * This file requires some constant definition to compile. See naive_linear_full.c
 * and table.th for details of required definitions.
 */
#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>						/* for definition of the NULL */

#define BIT_WIDTH 			8
#define BAND_WIDTH 			32

#include "seqreader.h"
#include "ddiag.h"

#define ALG 				SEA

/**
 * function declarations
 */
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

/**
 * @fn diag_linear_dynamic_banded
 *
 * @brief a linear-gap cost banded SIMD implementation.
 *
 * @param[ref] aln : a pointer to a sea_result structure. aln must NOT be NULL. (a structure which contains an alignment result.)
 * @param[in] aln->a : a pointer to a query sequence a.
 * @param[inout] aln->apos : (input) the start position of the search section on sequence a. (0 <= apos < length(sequence a)) (output) the start position of the alignment.
 * @param[inout] aln->alen : (input) the length of the search section on sequence a. (0 < alen) (output) the length of the alignment.
 * @param[in] aln->b : a pointer to a query sequence b.
 * @param[inout] aln->bpos : same as apos.
 * @param[inout] aln->blen : same as alen.
 *
 * @param[out] aln->aln : a pointer to the alignment string.
 * @param[inout] aln->len : (input) the reserved length of the aln->aln. (output) the length of the alignment string.
 * @param[out] aln->score : an alignment score.
 *
 * @param[in] scmat : unused parameter.
 * @param[in] xdrop : unused parameter. 
 * @param[in] bandwidth : unused parameter.
 * @param[in] mat : a pointer to the DP matrix.
 * @param[in] matsize : the size of the DP matrix.
 *
 * @return status. zero means success. see sea.h for the other error code.
 */
sea_int_t
diag_linear_dynamic_banded(
	struct sea_result *aln,
	struct sea_params param,
	void *mat)
{
	sea_int_t retval = SEA_ERROR;
	struct mpos o;

	/**
	 * fill in a matrix
	 */
	o = diag_linear_dynamic_banded_fill(
		aln, param,
		(char *)mat + BYTES_PER_LINE);		/** barriar */
	/**
	 * search maximum score position
	 */
	o = diag_linear_dynamic_banded_search(
		aln, param,
		(char *)mat + BYTES_PER_LINE,
		o);

#ifdef DEBUG
	sea_int_t const bw = BAND_WIDTH;
	fprintf(stderr, "ldmat, %d, %d\n", aln->alen+BAND_WIDTH, aln->blen+BAND_WIDTH);
	for(int p = 0; p < o.e.p; p++) {
		for(int q = -bw/2; q < bw/2; q++) {
			fprintf(stderr, "%d, %d, %d\n", p, q, SCORE((char *)mat + ADDR(p, q)));
		}
	}
	fprintf(stderr, "\n");
#endif
	/**
	 * if aln->aln is not NULL and score did not overflow, do traceback.
	 */
	if(o.m.p < 0) { return o.m.p; }
	if(aln->aln == NULL) { return SEA_SUCCESS; }
	retval = diag_linear_dynamic_banded_trace(
		aln, param,
		(char *)mat + BYTES_PER_LINE,
		o);
	return(retval);
}

/**
 * @fn diag_linear_dynamic_banded_fill
 *
 * @brief a matrix fill-in function.
 *
 * @param[ref] aln : a pointer to a sea_result struct.
 * @param[in] xdrop : currently unused.
 * @param[in] mat : a pointer to a matrix.
 * @param[in] matlen : the length of the p-direction of the matrix.
 *
 * @return an end position of the extension.
 */
struct mpos
diag_linear_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat)
{
	sea_int_t const bw = BAND_WIDTH;
	sea_int_t i, j;
	sea_int_t m = param.m,
			  x = param.x,
			  g = param.gi,
			  xdrop = param.xdrop;
	sea_int_t apos = 0,				/** seq.a read position */
			  bpos = 0;				/** seq.b read position */
	sea_int_t alen = aln->alen,
			  blen = aln->blen;
	sea_int_t alim = alen + bw/2,
			  blim = blen + bw/2;
	sea_int_t dir = -1;				/** direction of the band (downward or leftward, in diag_variant.h) */
	struct mpos o = {{0, 0, 0, 0}, {0, 0, 0, 0}};

	/**
	 * declarations of SIMD registers.
	 */
	DECLARE_VEC_CELL_REG(mv);		/** (m, m, m, ..., m) */
	DECLARE_VEC_CELL_REG(xv);		/** (x, x, x, ..., x) */
	DECLARE_VEC_CELL_REG(gv);		/** (g, g, g, ..., g) */
	DECLARE_VEC_CHAR_REG(wq);		/** a buffer for seq.a */
	DECLARE_VEC_CHAR_REG(wt);		/** a buffer for seq.b */
	DECLARE_VEC_CELL_REG(v);		/** score vector */
	DECLARE_VEC_CELL_REG(pv);		/** previous score vector */
	DECLARE_VEC_CELL_REG(tmp1);		/** temporary */
	DECLARE_VEC_CELL_REG(tmp2);		/** temporary */
	DECLARE_VEC_CELL_REG(maxv);		/** a vector which holds maximum scores */

	/**
	 * seqreaders
	 */
	DECLARE_SEQ(a);
	DECLARE_SEQ(b);

	CLEAR_SEQ(a, aln->a, aln->apos, aln->alen);
	CLEAR_SEQ(b, aln->b, aln->bpos, aln->blen);

	xdrop -= g;						/** compensation */
	VEC_SET(mv, m);					/** initialize mv vector with m */
	VEC_SET(xv, x);					/** xv with x */
	VEC_SET(gv, g);					/** gv with g */
	VEC_SET(maxv, CELL_MIN);		/** max score vector */

	VEC_SET(pv, CELL_MIN);			/** phantom vector at p = -1 */
	VEC_INIT_PVN(v, m, g, i);		/** initialize vector at p = 0 */
	VEC_STORE(mat, v);
	VEC_MAX(maxv, maxv, v);

	/**
	 * pre feed of sequences
	 */
	VEC_CHAR_SETZERO(wq);			/** padding of the seq.a is 0 */
	VEC_CHAR_SETONES(wt);			/** padding of the seq.b is 0xff */
	for(apos = 0; apos < bw/2; apos++) {
		FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0,    wq);
	}
	for(bpos = 0; bpos < bw/2-1; bpos++) {
		FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt);
	}

	/**
	 * @macro DIAG_LINEAR_DBANDED_REGISTER_UPDATE
	 * @brief a set of band advancing operations.
	 */
	#define DIAG_LINEAR_DBANDED_REGISTER_UPDATE() { \
	 	VEC_ADD(tmp1, v, gv); \
		VEC_ADD(tmp2, tmp2, gv); \
		VEC_MAX(tmp1, tmp1, tmp2); \
		VEC_COMPARE(tmp2, wq, wt); \
		VEC_SELECT(tmp2, xv, mv, tmp2); \
		VEC_ADD(tmp2, pv, tmp2); \
		VEC_MAX(tmp1, tmp1, tmp2); \
		VEC_ASSIGN(pv, v); \
		VEC_ASSIGN(v, tmp1); \
		VEC_MAX(maxv, maxv, v); \
	}

	i = 0; j = 0;							/** the center cell of the init vector */
	while(i < alim && j < blim) {
		VEC_ASSIGN(tmp2, v);
		if(VEC_MSB(v) > VEC_LSB(v)) {
			if(dir == DIR_V) { VEC_SHIFT_R(pv); VEC_INSERT_MSB(pv, CELL_MIN); }
			VEC_SHIFT_R(tmp2); VEC_INSERT_MSB(tmp2, CELL_MIN);
			j++; dir = DIR_V;				/** go downward */
			FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt); bpos++;
		} else {
			if(dir == DIR_H) { VEC_SHIFT_L(pv); VEC_INSERT_LSB(pv, CELL_MIN); }
			VEC_SHIFT_L(tmp2); VEC_INSERT_LSB(tmp2, CELL_MIN);
			i++; dir = DIR_H;				/** go rightward */
			FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0, wq); apos++;
		}
		DIAG_LINEAR_DBANDED_REGISTER_UPDATE();
		VEC_STORE(mat, v);
		if(ALG == NW && COP(i, j) == COP(alen, blen)) { break; }
		if(ALG == XSEA && VEC_CENTER(v) + xdrop - VEC_CENTER(maxv) < 0) { break; }
	}
	#undef VEC_SW_MAX
	#undef VEC_SEA_MAX
	#undef DIAG_LINEAR_DBANDED_REGISTER_UPDATE

	aln->len = COP(alim, blim);
	o.e.i = i; o.e.j = j;
	o.e.p = COP(i, j); o.e.q = COQ(i, j);	/** e.q == 0; */
	if(ALG != NW) {
		VEC_STORE(mat, maxv);
		int16_t max;
		VEC_HMAX(max, maxv);
		*((int16_t *)mat) = max;
	}
	return(o);
}

/**
 * @fn diag_linear_dynamic_banded_search
 *
 * @brief search a cell with maximal score.
 *
 * @return a struct mpos which contains maximal score position.
 */
struct mpos
diag_linear_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t const bw = BAND_WIDTH,
					bc = BYTES_PER_CELL,
					bl = BYTES_PER_LINE;

	sea_int_t i, j, p, q;
	char *lmat = mat + ADDR(o.e.p+1, -bw/2),
		 *tmat, *pl, *pu;					
	sea_int_t max = SCORE(lmat + bl);

	o.m.i = o.m.j = o.m.p = o.m.q = 0;
	if(max == CELL_MAX) { o.m.p = SEA_ERROR_OVERFLOW; return(o); }
	if(max == 0) { return(o); } /** wind back to (0, 0) if no meaningful score was found */
	for(q = -bw/2; q < bw/2; q++, lmat += bc) {
		if(SCORE(lmat) == max) {
			i = o.e.i - q; j = o.e.j + q;
			pl = mat + ADDR(o.e.p-1, bw/2-1);
			pu = mat + ADDR(o.e.p-1, -bw/2);
			for(p = o.e.p, tmat = lmat-bl;
				p >= 0 && SCORE(tmat) != max;
				p--, tmat -= bl) {
				if(SCORE(pl) > SCORE(pu)) { j--; } else { i--; }
				pl -= bl; pu -= bl;
			}
			if(p >= o.m.p) { o.m.i = i; o.m.j = j; o.m.p = p; o.m.q = q; }
		}
	}
	aln->alen = o.m.i;
	aln->blen = o.m.j;
	aln->score = SCORE(mat + ADDR(o.m.p, o.m.q));
	return(o);
}

/**
 * @fn diag_linear_dynamic_banded_trace
 *
 * @brief traceback function.
 *
 * @param aln : a pointer to struct sea_result.
 * @param mat : a pointer to dynamic programming matrix.
 * @param matlen : the length of the p-coordinate of the matrix.
 * @param mpos : the start position of the trace.
 */
sea_int_t
diag_linear_dynamic_banded_trace(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t const bw = BAND_WIDTH,
					bl = BYTES_PER_LINE;
	sea_int_t mi = o.m.i,
			  mj = o.m.j,
			  mp = o.m.p,
			  mq = o.m.q;
	sea_int_t g = param.gi;
	char *tmat = (char *)mat + ADDR(mp, mq),
		 *pu = (char *)mat + ADDR(mp-1, -bw/2),
		 *pl = (char *)mat + ADDR(mp-1, bw/2-1),
		 *p = (char *)aln->aln + aln->len - 1;
	sea_int_t score, v, h, dir;
	sea_int_t icnt = 0, ibases = 0, dcnt = 0, dbases = 0;
	char prev = 0;

	#define DET_DIR(dir, pl, pu) { \
		(dir) = (SCORE(pl) > SCORE(pu)) ? DIR_V : DIR_H; \
		(pl) -= bl; (pu) -= bl; \
	}

	score = SCORE(tmat);
	dir = 0;
	while(mi > 0 || mj > 0) {
		DET_DIR(dir, pl, pu);
		if(score == ((v = SCORE(tmat + DTOP(dir))) + g)) {
			tmat += DTOP(dir); score = v;
			mj--; icnt++; ibases += prev != 'I';
			*p-- = 'I'; prev = 'I';
		} else if(score == ((h = SCORE(tmat + DLEFT(dir))) + g)) {
			tmat += DLEFT(dir); score = h;
			mi--; dcnt++; dbases += prev != 'D';
			*p-- = 'D'; prev = 'D';
		} else {
			tmat += DLEFT(dir); DET_DIR(dir, pl, pu);
			tmat += DTOP(dir);
			mi--;
			mj--;
			*p-- = 'M'; prev = 'M';
			score = SCORE(tmat);
		}
	}
	aln->icnt = icnt;
	aln->ibases = ibases;
	aln->dcnt = dcnt;
	aln->dbases = dbases;

	char *b = (char *)aln->aln;
	int64_t l = 0;
	while(p < (char *)aln->aln + aln->len - 1) {
		*b++ = *++p;
		l++;
	}
	*b = '\0';
	aln->len = l;

	#undef DET_DIR

	return SEA_SUCCESS;
}

/**
 * @fn diag_linear_dynamic_banded_matsize
 *
 * @brief returns the size of matrix for diag_linear_dynamic_banded, in bytes.
 *
 * @param[in] alen : the length of the input sequence a (in base pairs).
 * @param[in] blen : the length of the input sequence b (in base pairs).
 * @param[in] bandwidth : unused. give zero for practice.
 *
 * @return the size of a matrix in bytes.
 */
sea_int_t
diag_linear_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth)
{
	/** barriar (1) + matrix (alen+blen+bw) + max vectors (2) */
	return((1 + alen + blen + bandwidth + 4) * BYTES_PER_LINE);
}

/**
 * @fn ddiag_linear
 * @brief wrapper to main.cc
 */
int ddiag_linear(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt)
{

	struct sea_result *aln = (struct sea_result *)work;
	aln->score = 0;
//	aln->ctx = NULL;

	void *base = (void *)(aln + 1);
	memset(base, 0, 32 * 6 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 6 * sizeof(int16_t);

	aln->score = 0;
	aln->a = a; aln->apos = 0; aln->alen = alen;
	aln->b = b; aln->bpos = 0; aln->blen = blen;
	aln->len = 1024 * 1024 * 1024;


	struct sea_params param;
	param.m = score_matrix[0];
	param.x = score_matrix[1];
	param.gi = 0;
	param.ge = ge;
	param.xdrop = xt;

	struct mpos o = diag_linear_dynamic_banded_fill(aln, param, mat);
	diag_linear_dynamic_banded_search(aln, param, mat, o);
	return(aln->score);
}

#if 0
/**
 * @fn main
 *
 * @brief gets two ascii strings from stdin, align strings with naive_affine_full, and print the result.
 */
#ifdef MAIN

int main(int argc, char *argv[])
{
	int matlen, alnlen;
	void *mat;
	struct sea_result aln;
	struct sea_params param;

	param.flags = 0;
	param.m = 2;
	param.x = -3;
	param.gi = -4;
	param.ge = -1;
	param.xdrop = 12;
	param.bandwidth = 32;

	aln.a = argv[1];
	aln.alen = strlen(aln.a);
	aln.apos = 0;
	aln.b = argv[2];
	aln.blen = strlen(aln.b);
	aln.bpos = 0;
	alnlen = aln.alen + aln.blen;
	aln.len = alnlen;

	alnlen = aln.alen + aln.blen;
	matlen = diag_linear_dynamic_banded_matsize(
		aln.alen, aln.blen, param.bandwidth);

	aln.aln = (void *)malloc(alnlen);
	mat = (void *)malloc(matlen);
	diag_linear_dynamic_banded(&aln, param, mat);

	printf("%d, %d, %s\n", aln.score, aln.len, aln.aln);

	free(mat);
	free(aln.aln);
	return 0;
}

#endif /* #ifdef MAIN */

/**
 * unittest functions.
 *
 * give a definition of TEST to a compiler, e.g. -DTEST=1, to make a library for tests.
 */
#ifdef TEST
#if SEQ == ascii && ALN == ascii
#if !(BAND_WIDTH == 32 && BIT_WIDTH == 8 && ALG == XSEA)

extern sea_int_t
naive_linear_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);

extern sea_int_t
naive_linear_banded(
	struct sea_result *aln,
	struct sea_params param,
	void *mat);

/**
 * @fn test_2_diag_linear_dynamic_banded
 *
 * @brief a unittest function of diag_linear_dynamic_banded.
 *
 * @detail
 * This function is an aggregation of simple fixed ascii queries.
 * In this function, a sea_assert_align macro in tester.h is called. This macro
 * calls sea_align function with given context, checks if the score and the alignment
 * string is the same as the given score and string. If both or one of the results
 * are different, the macro prints the failure status (filename, lines, input sequences,
 * and the content of sea_result) and dumps a content of dynamic programming matrix
 * to dumps.log.
 */
void
test_2_diag_linear_dynamic_banded(
	void)
{
	sea_int_t m = 2,
			  x = -3,
			  g = -4;
	struct sea_context *ctx;

	ctx = sea_init_fp(
		SEA_BANDWIDTH_64,
		diag_linear_dynamic_banded,
		diag_linear_dynamic_banded_matsize,
		m, x, g, 0,			/** the default blast scoring scheme */
		12);				/** xdrop threshold */

	/**
	 * when both sequences are empty
	 */
	sea_assert_align(ctx, "", 				"", 			0,			"");
	/**
	 * when one sequence is empty
	 */
	sea_assert_align(ctx, "AAA", 			"", 			0,			"");
	sea_assert_align(ctx, "TTTTTTT", 		"", 			0,			"");
	sea_assert_align(ctx, "", 				"AAA", 			0,			"");
	sea_assert_align(ctx, "", 				"TTTTTTT", 		0,			"");
	sea_assert_align(ctx, "CTAG",			"", 			0,			"");

	/**
	 * a match
	 */
	sea_assert_align(ctx, "A",				"A", 			m,			"M");
	sea_assert_align(ctx, "C", 				"C", 			m,			"M");
	sea_assert_align(ctx, "G", 				"G", 			m,			"M");
	sea_assert_align(ctx, "T", 				"T", 			m,			"M");

	/**
	 * a mismatch
	 */
	if(ALG == NW) {
		/**
		 * the Needleman-Wunsch algorithm
		 */
		sea_assert_align(ctx, "A", 				"C", 			x,			"X");
		sea_assert_align(ctx, "C", 				"A", 			x,			"X");
		sea_assert_align(ctx, "G", 				"T", 			x,			"X");
		sea_assert_align(ctx, "T", 				"A", 			x,			"X");
	} else if(ALG == SEA || ALG == SW || ALG == XSEA) {
		/**
		 * the Smith-Waterman algorithm and the seed-and-extend algorithm
		 */
		sea_assert_align(ctx, "A", 				"C", 			0,			"");
		sea_assert_align(ctx, "C", 				"A", 			0,			"");
		sea_assert_align(ctx, "G", 				"T", 			0,			"");
		sea_assert_align(ctx, "T", 				"A", 			0,			"");
	}

	/**
	 * homopolymers with different lengths.
	 */
	if(ALG == NW) {
		/**
		 * the Needleman-Wunsch algorithm
		 */
		sea_assert_align(ctx, "A", 				"AA", 			g+m,			"IM");
		sea_assert_align(ctx, "A", 				"AAA", 			g+g+m,			"IIM");
		sea_assert_align(ctx, "AAAA", 			"AA", 			g+g+m+m,		"DDMM");
		sea_assert_align(ctx, "TTTT", 			"TTTTTTTT", 	4*g+4*m,		"IIIIMMMM");
	} else if(ALG == SEA || ALG == SW || ALG == XSEA) {
		/**
		 * the Smith-Waterman algorithm and the seed-and-extend algorithm
		 */
		sea_assert_align(ctx, "A", 				"AA", 			m,				"M");
		sea_assert_align(ctx, "A", 				"AAA", 			m,				"M");
		sea_assert_align(ctx, "AAAA", 			"AA", 			m+m,			"MM");
		sea_assert_align(ctx, "TTTT", 			"TTTTTTTT", 	4*m,			"MMMM");
	}

	/**
	 * when mismatches occurs.
	 */
	sea_assert_align(ctx, "AAAAAAAA", 		"AAAATAAA", 	7*m+x,			"MMMMXMMM");
	sea_assert_align(ctx, "TTTTTTTT", 		"TTTCTTTT", 	7*m+x,			"MMMXMMMM");
	sea_assert_align(ctx, "CCCCCCCC", 		"CCTCCCCC", 	7*m+x,			"MMXMMMMM");
	sea_assert_align(ctx, "GGGGGGGG", 		"GGCGGTGG", 	6*m+2*x,		"MMXMMXMM");

	/**
	 * when gaps with a base occurs on seq a (insertion).
	 */
	sea_assert_align(ctx, "AAAAATTTT", 		"AAAAAGTTTT", 	9*m+g,			"MMMMMIMMMM");
	sea_assert_align(ctx, "TTTTCCCCC", 		"TTTTACCCCC", 	9*m+g,			"MMMMIMMMMM");
	sea_assert_align(ctx, "CCCGGGGGG", 		"CCCTGGGGGG", 	9*m+g,			"MMMIMMMMMM");
	sea_assert_align(ctx, "GGGAATTT", 		"GGGCAAGTTT", 	8*m+2*g,		"MMMIMMIMMM");

	/**
	 * when gaps with a base occurs on seq b (deletion).
	 */
	sea_assert_align(ctx, "AAAAAGTTTT", 	"AAAAATTTT", 	9*m+g,			"MMMMMDMMMM");
	sea_assert_align(ctx, "TTTTACCCCC", 	"TTTTCCCCC", 	9*m+g,			"MMMMDMMMMM");
	sea_assert_align(ctx, "CCCTGGGGGG", 	"CCCGGGGGG", 	9*m+g,			"MMMDMMMMMM");
	sea_assert_align(ctx, "GGGCAAGTTT", 	"GGGAATTT", 	8*m+2*g,		"MMMDMMDMMM");

	/**
	 * when a gap longer than two bases occurs on seq a.
	 */
	sea_assert_align(ctx, "AAAAATTTTT", 	"AAAAAGGTTTTT", 10*m+2*g,		"MMMMMIIMMMMM");
	sea_assert_align(ctx, "TTTAAAAGGG", 	"TTTCAAAACGGG", 10*m+2*g,		"MMMIMMMMIMMM");

	/**
	 * when a gap longer than two bases occurs on seq b.
	 */
	sea_assert_align(ctx, "AAAAAGGTTTTT",	"AAAAATTTTT", 	10*m+2*g,		"MMMMMDDMMMMM");
	sea_assert_align(ctx, "TTTCAAAACGGG", 	"TTTAAAAGGG", 	10*m+2*g,	 	"MMMDMMMMDMMM");

	/**
	 * when outer gaps occurs.
	 */
	if(ALG == NW) {
	 	sea_assert_align(ctx, "TTAAAAAAAATT", 	"CCAAAAAAAACC",	8*m+4*x, 		"XXMMMMMMMMXX");
	} else if(ALG == SEA || ALG == XSEA) {
	 	sea_assert_align(ctx, "TTAAAAAAAATT", 	"CCAAAAAAAACC",	8*m+2*x, 		"XXMMMMMMMM");
	} else if(ALG == SW) {
	 	sea_assert_align(ctx, "TTAAAAAAAATT", 	"CCAAAAAAAACC",	8*m,	 		"MMMMMMMM");
	}

	/**
	 * X-drop test (here xdrop threshold is 12)
	 */
	if(ALG == XSEA) {
		sea_assert_align(ctx, "TTTTTTAAAAAAAA",	"GGGGGGAAAAAAAA",	0, 			"");
		sea_assert_align(ctx, "TTTTTAAAAAAAAA",	"GGGGGAAAAAAAAA",	0, 			"");
		sea_assert_align(ctx, "TTTTAAAAAAAA",	"GGGGAAAAAAAA",		8*m+4*x, 	"XXXXMMMMMMMM");
		sea_assert_align(ctx, "TTTAAAAAAAAA",	"GGGAAAAAAAAA",		9*m+3*x, 	"XXXMMMMMMMMM");

		sea_assert_align(ctx, "AAAATTTTTAAAAAAA","AAAAGGGGGAAAAAAA",	4*m, 		"MMMM");
		sea_assert_align(ctx, "AAAATTTAAAAAAAAA","AAAAGGGAAAAAAAAA",	13*m+3*x, 	"MMMMXXXMMMMMMMMM");
	}

	sea_clean(ctx);
	return;
}

/**
 * @fn test_8_cross_diag_linear_dynamic_banded
 *
 * @brief cross test between naive_linear_banded and diag_linear_dynamic_banded
 */
#if HAVE_NAIVE_BANDED

void
test_8_cross_diag_linear_dynamic_banded(
	void)
{
	int i;
	int const cnt = 5;
	sea_int_t m = 2,
			  x = -3,
			  g = -4;
	char *a, *b;
	struct sea_context *full, *band;

	#if BIT_WIDTH == 8
		int const len = 50;
	#elif BIT_WIDTH == 16
		int const len = 1000;
	#else
		#error "bit width must be 8 or 16."
	#endif

	full = sea_init_fp(
		SEA_BANDWIDTH_64,
		naive_linear_banded,
		naive_linear_banded_matsize,
		m, x, g, 0,
		10000);
	band = sea_init_fp(
		SEA_BANDWIDTH_64,
		diag_linear_dynamic_banded,
		diag_linear_dynamic_banded_matsize,
		m, x, g, 0,
		10000);


	for(i = 0; i < cnt; i++) {
		a = rseq(len);
		b = mseq(a, 20, 100, 100);
		sea_assert_cross(
			full, band, a, b);
		free(a); free(b);
	}

	sea_clean(full);
	sea_clean(band);
	return;
}

#endif /* #if HAVE_NAIVE_BANDED */

#endif 	/* !(BAND_WIDTH == 32 && BIT_WIDTH == 8 && ALG == XSEA) */
#endif	/* #if SEQ == ascii && ALN == ascii */
#endif	/* #ifdef TEST */
#endif
/**
 * end of diag_linear_dynamic_banded.c
 */
