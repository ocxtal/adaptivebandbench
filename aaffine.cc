
/**
 * @file diag_affine_dynamic_banded.c
 *
 * @brief a affine-gap cost, banded SIMD implementation.
 *
 * @detail
 * This implementation is a vectorized variant of naive_affine_dynamic_banded.c.
 * The consistency of the results are tested in the test_cross_diag_affine_dynamic_banded
 * function.
 *
 * This file requires some constant definition to compile. See naive_affine_full.c
 * and table.th for details of required definitions.
 */
#include <cstdio>
#include <cstdint>
#include <cstdlib>						/* for definition of the NULL */
#include <cctype>

using namespace std;

#define BIT_WIDTH 			8
#ifdef BW
#  define BAND_WIDTH		BW
#else
#  define BAND_WIDTH		32
#endif

#include "seqreader.h"
#include "ddiag.h"

#define ALG 				XSEA

/**
 * function declarations
 */
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

/**
 * @fn diag_affine_dynamic_banded
 *
 * @brief a affine-gap cost banded SIMD implementation.
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
diag_affine_dynamic_banded(
	struct sea_result *aln,
	struct sea_params param,
	void *mat)
{
	sea_int_t retval = SEA_ERROR;
	struct mpos o;

	/**
	 * fill in a matrix
	 */
	o = diag_affine_dynamic_banded_fill(
		aln, param,
		(char *)mat + 3 * BYTES_PER_LINE);		/** barriar */
	/**
	 * search maximum score position
	 */
	o = diag_affine_dynamic_banded_search(
		aln, param,
		(char *)mat + 3 * BYTES_PER_LINE,		/** barriar */
		o);

#ifdef DEBUG
	fprintf(stderr, "ldmat, %d, %d\n", aln->alen+BAND_WIDTH, aln->blen+BAND_WIDTH);
	for(int i = 0; i < aln->alen+1; i++) {
		for(int j = 0; j < aln->blen+1; j++) {
			fprintf(stderr, "%d, %d, %d, %d, %d, %d\n", i, j, COP(i, j), COQ(i, j), AADDR(COP(i, j), COQ(i, j)), ASCOREV((char *)mat + AADDR(COP(i, j), COQ(i, j))));
		}
	}
	fprintf(stderr, "\n");

	for(int i = 0; i < 3*e.p+4; i++) {
		for(int j = 0; j < BAND_WIDTH; j++) {
			fprintf(stderr, "%4d, ", ASCOREV((char *)mat + i*BYTES_PER_LINE + j*BYTES_PER_CELL));
		}
		fprintf(stderr, "\n");
	}
#endif

	/**
	 * if aln->aln is not NULL and score did not overflow, do traceback.
	 */
	if(o.m.p < 0) { return o.m.p; }
	if(aln->aln == NULL) { return SEA_SUCCESS; }
	retval = diag_affine_dynamic_banded_trace(
		aln, param,
		(char *)mat + 3 * BYTES_PER_LINE,		/** barriar */
		o);
	return(retval);
}

/**
 * @fn diag_affine_dynamic_banded_fill
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
diag_affine_dynamic_banded_fill(
	struct sea_result *aln,
	struct sea_params param,
	char *mat)
{
	sea_int_t const bw = BAND_WIDTH;
	sea_int_t i, j;
	sea_int_t m = param.m,
			  x = param.x,
			  gi = param.gi,
			  ge = param.ge,
			  xdrop = param.xdrop;
	sea_int_t apos = 0,				/** seq.a read position */
			  bpos = 0;				/** seq.b read position */
	sea_int_t alen = aln->alen,
			  blen = aln->blen;
	sea_int_t alim = alen + bw/2,
			  blim = blen + bw/2;
	sea_int_t dir = -1;				/** direction (DIR_H or DIR_V in diag_variant.h) */
	struct mpos o;
	
	/**
	 * declarations of SIMD registers.
	 */
	DECLARE_VEC_CELL_REG(mv);
	DECLARE_VEC_CELL_REG(xv);
	DECLARE_VEC_CELL_REG(giv);
	DECLARE_VEC_CELL_REG(gev);
	DECLARE_VEC_CHAR_REG(wq);
	DECLARE_VEC_CHAR_REG(wt);
	DECLARE_VEC_CELL_REG(v);
	DECLARE_VEC_CELL_REG(pv);
	DECLARE_VEC_CELL_REG(f);
	DECLARE_VEC_CELL_REG(e);
	DECLARE_VEC_CELL_REG(tmp1);
	DECLARE_VEC_CELL_REG(tmp2);
	DECLARE_VEC_CELL_REG(maxv);

	/**
	 * seqreaders
	 */
	DECLARE_SEQ(a);
	DECLARE_SEQ(b);

	CLEAR_SEQ(a, aln->a, aln->apos, aln->alen);
	CLEAR_SEQ(b, aln->b, aln->bpos, aln->blen);

	xdrop -= gi;					/** compensation */
	VEC_SET(mv, m);					/** (m, m, m, ..., m) */
	VEC_SET(xv, x);					/** (x, x, x, ..., x) */
	VEC_SET(giv, gi);				/** (gi, gi, gi, ..., gi) */
	VEC_SET(gev, ge);				/** (ge, ge, ge, ..., ge) */
	VEC_SET(maxv, CELL_MIN);		/** init max score vector with CELL_MIN */
	VEC_SET(pv, CELL_MIN);			/** phantom vector at p = -1 */

	/**
	 * initialize band at p = 0; the center cell is at (i, j) = (0, 0)
	 */
	VEC_INIT_PVN(v, m, gi, i); VEC_STORE(mat, v);
	VEC_SET(f, CELL_MIN); VEC_STORE(mat, f);
	VEC_SET(e, CELL_MIN); VEC_STORE(mat, e);
	VEC_MAX(maxv, maxv, v);

	/**
	 * pre feed of sequences
	 */
	VEC_CHAR_SETZERO(wq);			/** padding of seq.a: 0 */
	VEC_CHAR_SETONES(wt);			/** padding of seq.b: 0xff */
	for(apos = 0; apos < bw/2; apos++) {
		FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0,    wq);
	}
	for(bpos = 0; bpos < bw/2-1; bpos++) {
		FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt);
	}

	/**
	 * @macro UPDATE_FHALF
	 * @brief a set of the first half of the band advancing operations.
	 */
	#define UPDATE_FHALF() { \
		VEC_ADD(tmp1, v, giv); \
		VEC_ADD(tmp2, f, gev); \
		VEC_MAX(f, tmp1, tmp2); \
		VEC_ADD(tmp2, e, gev); \
		VEC_MAX(e, tmp1, tmp2); \
	}
	/**
	 * @macro UPDATE_LHALF
	 * @brief the last half of the operations.
	 */
	#define UPDATE_LHALF() { \
		VEC_COMPARE(tmp2, wq, wt); \
		VEC_SELECT(tmp2, xv, mv, tmp2); \
		VEC_ADD(tmp2, pv, tmp2); \
		VEC_MAX(tmp2, f, tmp2); \
		VEC_MAX(tmp2, e, tmp2); \
		VEC_ASSIGN(pv, v); \
		VEC_ASSIGN(v, tmp2); \
		VEC_MAX(maxv, maxv, v); \
	}

	i = 0; j = 0;							/** the center cell of the init vector */
	while(i < alim && j < blim) {
		UPDATE_FHALF();
		if(VEC_MSB(v) > VEC_LSB(v)) {
			if(dir == DIR_V) { VEC_SHIFT_R(pv); VEC_INSERT_MSB(pv, CELL_MIN); }
			VEC_SHIFT_R(e); VEC_INSERT_MSB(e, CELL_MIN);
			j++; dir = DIR_V;				/** go downward */
			FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt); bpos++;
		} else {
			if(dir == DIR_H) { VEC_SHIFT_L(pv); VEC_INSERT_LSB(pv, CELL_MIN); }
			VEC_SHIFT_L(f); VEC_INSERT_LSB(f, CELL_MIN);
			i++; dir = DIR_H;				/** go rightward */
			FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0, wq); apos++;
		}
		UPDATE_LHALF();
		VEC_STORE(mat, v); VEC_STORE(mat, f); VEC_STORE(mat, e);
		if(ALG == NW && COP(i, j) == COP(alen, blen)) { break; }
		if(ALG == XSEA && VEC_CENTER(v) + xdrop - VEC_CENTER(maxv) < 0) { break; }
	}
	#undef VEC_SW_MAX
	#undef VEC_SEA_MAX
	#undef UPDATE_FHALF
	#undef UPDATE_LHALF

	aln->len = COP(alim, blim);
	o.e.i = i; o.e.j = j;
	o.e.p = COP(i, j); o.e.q = COQ(i, j);		/** o.e.q == 0 */
	if(ALG != NW) {
		VEC_STORE(mat, maxv);					/** store max vector at the end of the memory */
		VEC_ASSIGN(tmp1, maxv);
		for(i = 1; i < bw; i++) {
			VEC_SHIFT_R(tmp1);
			VEC_MAX(maxv, tmp1, maxv);			/** extract maximum score in the maxv vector */ 
		}
		VEC_STORE(mat, maxv);					/** store max of the max vector */
	}
	return(o);
}

/**
 * @fn diag_affine_dynamic_banded_search
 *
 * @brief search a cell with maximal score.
 *
 * @return a struct mpos which contains maximal score position.
 */
struct mpos
diag_affine_dynamic_banded_search(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t const bw = BAND_WIDTH,
					bc = BYTES_PER_CELL,
					bl = BYTES_PER_LINE;

	if(ALG == NW) {
		o.m.i = aln->alen; o.m.j = aln->blen;
		o.m.p = o.e.p; o.m.q = COQ(o.m.i, o.m.j) - o.e.q;
		if(o.m.p != COP(aln->alen, aln->blen) || o.m.q > bw/2-1 || o.m.q < -bw/2) {
			o.m.p = SEA_ERROR_OUT_OF_BAND; return(o);
		}
	} else {	/** ALG == SEA or ALG == SW */
		sea_int_t i, j, p, q;
		char *lmat = mat + AADDR(o.e.p+1, -bw/2),
			 *tmat, *pl, *pu;
		sea_int_t max = ASCOREV(lmat + bl);

		o.m.i = o.m.j = o.m.p = o.m.q = 0;
		if(max == CELL_MAX) { o.m.p = SEA_ERROR_OVERFLOW; return(o); }
		if(max == 0) { return(o); } /** wind back to (0, 0) if no meaningful score was found */
		for(q = -bw/2; q < bw/2; q++, lmat += bc) {
			if(ASCOREV(lmat) == max) {
				i = o.e.i - q; j = o.e.j + q;
				pl = mat + AADDR(o.e.p-1, bw/2-1);
				pu = mat + AADDR(o.e.p-1, -bw/2);
				for(p = o.e.p, tmat = lmat-3*bl;
					p >= 0 && ASCOREV(tmat) != max;
					p--, tmat -= 3*bl) {
					if(ASCOREV(pl) > ASCOREV(pu)) { j--; } else { i--; }
					pl -= 3*bl; pu -= 3*bl;
				}
				if(p >= o.m.p) { o.m.i = i; o.m.j = j; o.m.p = p; o.m.q = q; }
			}
		}
		aln->alen = o.m.i;
		aln->blen = o.m.j;
	}
	aln->score = ASCOREV(mat + AADDR(o.m.p, o.m.q));
	return(o);
}

/**
 * @fn diag_affine_dynamic_banded_trace
 *
 * @brief traceback function.
 *
 * @param aln : a pointer to struct sea_result.
 * @param mat : a pointer to dynamic programming matrix.
 * @param matlen : the length of the p-coordinate of the matrix.
 * @param m : the start position of the trace.
 */
sea_int_t
diag_affine_dynamic_banded_trace(
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
	sea_int_t gi = param.gi;
	char *tmat = (char *)mat + AADDR(mp, mq),
		 *pu = (char *)mat + AADDR(mp-1, -bw/2),
		 *pl = (char *)mat + AADDR(mp-1, bw/2-1),
		 *p = (char *)aln->aln + aln->len - 1;
	sea_int_t score, dir;
	sea_int_t icnt = 0, ibases = 0, dcnt = 0, dbases = 0;

	#define DET_DIR(dir, pl, pu) { \
		(dir) = (ASCOREV(pl) > ASCOREV(pu)) ? DIR_V : DIR_H; \
		(pl) -= 3*bl; (pu) -= 3*bl; \
	}

	score = ASCOREV(tmat);
	dir = 0;
	while(mi > 0 || mj > 0) {
		DET_DIR(dir, pl, pu);
		if(score == ASCOREF(tmat)) {
			while(ASCOREF(tmat) != ASCOREV(tmat + DATOP(dir)) + gi) {
				tmat += DATOP(dir);
				mj--; ibases++;
				*p-- = 'I';
				DET_DIR(dir, pl, pu);
			}
			tmat += DATOP(dir);
			mj--; icnt++; ibases++;
			*p-- = 'I';
			score = ASCOREV(tmat);
		} else if(score == ASCOREE(tmat)) {
			while(ASCOREE(tmat) != ASCOREV(tmat + DALEFT(dir)) + gi) {
				tmat += DALEFT(dir);
				mi--; dbases++;
				*p-- = 'D';
				DET_DIR(dir, pl, pu);
			}
			tmat += DALEFT(dir);
			mi--; dcnt++; dbases++;
			*p-- = 'D';
			score = ASCOREV(tmat);
		} else {
			tmat += DALEFT(dir); DET_DIR(dir, pl, pu);
			tmat += DATOP(dir);
			mi--;
			mj--;
			*p-- = 'M';
			score = ASCOREV(tmat);
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
 * @fn diag_affine_dynamic_banded_matsize
 *
 * @brief returns the size of matrix for diag_affine_dynamic_banded, in bytes.
 *
 * @param[in] alen : the length of the input sequence a (in base pairs).
 * @param[in] blen : the length of the input sequence b (in base pairs).
 * @param[in] bandwidth : unused. give zero for practice.
 *
 * @return the size of a matrix in bytes.
 */
sea_int_t
diag_affine_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth)
{
	/** barriar (1) + matrix (alen+blen+bw) + max vectors (2) */
	return(3 * (1 +	alen + blen + bandwidth + 4) * BYTES_PER_LINE);
}

/**
 * @fn ddiag_affine
 * @brief wrapper to main.cc
 */
int ddiag_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt)
{

	struct sea_result *aln = (struct sea_result *)work;
	aln->score = 0;
	aln->ctx = NULL;

	void *base = (void *)(aln + 1);
	memset(base, 0, 32 * 6 * sizeof(int16_t));
	char *mat = (char *)base + 32 * 6 * sizeof(int16_t);

	aln->score = 0;
	aln->a = (void const *)a; aln->apos = 0; aln->alen = alen;
	aln->b = (void const *)b; aln->bpos = 0; aln->blen = blen;
	aln->len = 1024 * 1024 * 1024ULL;


	struct sea_params param;
	param.m = score_matrix[0];
	param.x = score_matrix[1];
	param.gi = gi;
	param.ge = ge;
	param.xdrop = xt;

	struct mpos o = diag_affine_dynamic_banded_fill(aln, param, mat);
	diag_affine_dynamic_banded_search(aln, param, mat, o);
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
	matlen = diag_affine_dynamic_banded_matsize(
		aln.alen, aln.blen, param.bandwidth);

	aln.aln = (void *)malloc(alnlen);
	mat = (void *)malloc(matlen);
	diag_affine_dynamic_banded(&aln, param, mat);

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
naive_affine_dynamic_banded_matsize(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);

extern sea_int_t
naive_affine_dynamic_banded(
	struct sea_result *aln,
	struct sea_params param,
	void *mat);

/**
 * @fn test_2_diag_affine_dynamic_banded
 *
 * @brief a unittest function of diag_affine_dynamic_banded.
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
test_2_diag_affine_dynamic_banded(
	void)
{
	sea_int_t m = 2,
			  x = -3,
			  gi = -4,
			  ge = -1;
	struct sea_context *ctx;

	ctx = sea_init_fp(
		SEA_BANDWIDTH_64,
		diag_affine_dynamic_banded,
		diag_affine_dynamic_banded_matsize,
		m, x, gi, ge,		/** the default blast scoring scheme */
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
		sea_assert_align(ctx, "A", 				"AA", 			gi+m,			"IM");
		sea_assert_align(ctx, "A", 				"AAA", 			gi+ge+m,		"IIM");
		sea_assert_align(ctx, "AAAA", 			"AA", 			gi+ge+m+m,		"DDMM");
		sea_assert_align(ctx, "TTTT", 			"TTTTTTTT", 	gi+3*ge+4*m,	"IIIIMMMM");
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
	sea_assert_align(ctx, "AAAAATTTT", 		"AAAAAGTTTT", 	9*m+gi,			"MMMMMIMMMM");
	sea_assert_align(ctx, "TTTTCCCCC", 		"TTTTACCCCC", 	9*m+gi,			"MMMMIMMMMM");
	sea_assert_align(ctx, "CCCGGGGGG", 		"CCCTGGGGGG", 	9*m+gi,			"MMMIMMMMMM");
	sea_assert_align(ctx, "GGGAATTT", 		"GGGCAAGTTT", 	8*m+2*gi,		"MMMIMMIMMM");

	/**
	 * when gaps with a base occurs on seq b (deletion).
	 */
	sea_assert_align(ctx, "AAAAAGTTTT", 	"AAAAATTTT", 	9*m+gi,			"MMMMMDMMMM");
	sea_assert_align(ctx, "TTTTACCCCC", 	"TTTTCCCCC", 	9*m+gi,			"MMMMDMMMMM");
	sea_assert_align(ctx, "CCCTGGGGGG", 	"CCCGGGGGG", 	9*m+gi,			"MMMDMMMMMM");
	sea_assert_align(ctx, "GGGCAAGTTT", 	"GGGAATTT", 	8*m+2*gi,		"MMMDMMDMMM");

	/**
	 * when a gap longer than two bases occurs on seq a.
	 */
	sea_assert_align(ctx, "AAAATTTT", 		"AAAAGGTTTT", 	8*m+gi+ge,		"MMMMIIMMMM");
	sea_assert_align(ctx, "GGGGCCCC", 		"GGGGTTTCCCC", 	8*m+gi+2*ge,	"MMMMIIIMMMM");
	sea_assert_align(ctx, "GGGGGCCCCC", 	"GGGGGTTTTCCCCC",10*m+gi+3*ge,	"MMMMMIIIIMMMMM");
	sea_assert_align(ctx, "TTTTAAGGGG", 	"TTTTCCAACCGGGG",10*m+2*gi+2*ge,"MMMMIIMMIIMMMM");

	/**
	 * when a gap longer than two bases occurs on seq b.
	 */
	sea_assert_align(ctx, "AAAAGGTTTT",	 	"AAAATTTT", 	8*m+gi+ge,		"MMMMDDMMMM");
	sea_assert_align(ctx, "GGGGTTTCCCC",	"GGGGCCCC", 	8*m+gi+2*ge,	"MMMMDDDMMMM");
	sea_assert_align(ctx, "GGGGGTTTTCCCCC",	"GGGGGCCCCC", 	10*m+gi+3*ge,	"MMMMMDDDDMMMMM");
	sea_assert_align(ctx, "TTTTCCAACCGGGG", "TTTTAAGGGG", 	10*m+2*gi+2*ge, "MMMMDDMMDDMMMM");

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
 * @fn test_8_cross_diag_affine_dynamic_banded
 *
 * @brief cross test between naive_affine_dynamic_banded and diag_affine_dynamic_banded
 */
#if HAVE_NAIVE_BANDED

void
test_8_cross_diag_affine_dynamic_banded(
	void)
{
	int i;
	int const cnt = 5;
	sea_int_t m = 2,
			  x = -3,
			  gi = -4,
			  ge = -1;
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
		naive_affine_dynamic_banded,
		naive_affine_dynamic_banded_matsize,
		m, x, gi, ge,
		10000);
	band = sea_init_fp(
		SEA_BANDWIDTH_64,
		diag_affine_dynamic_banded,
		diag_affine_dynamic_banded_matsize,
		m, x, gi, ge,
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
 * end of diag_affine_dynamic_banded.c
 */
