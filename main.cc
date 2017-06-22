
/**
 * @file main.cc
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/time.h>
#include "util.h"
#include "full.h"
#include "kvec.h"
#include "bench.h"
#include "ssw.h"

int blast_linear(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int blast_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);

int simdblast_linear(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int simdblast_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);

int ddiag_linear(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);
int ddiag_affine(
	void *work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen,
	int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt);
	// int8_t m, int8_t x, int8_t gi, int8_t ge, int16_t xt);



/* wrapper of Myers' wavefront algorithm */
extern "C" {
	#include "wave/align.h"

	int enlarge_vector(Work_Data *work, int newmax);
	int enlarge_points(Work_Data *work, int newmax);
	int forward_wave(Work_Data *work, Align_Spec *spec, Alignment *align, Path *bpath,
	                        int *mind, int maxd, int mida, int minp, int maxp);
	int reverse_wave(Work_Data *work, Align_Spec *spec, Alignment *align, Path *bpath,
	                        int mind, int maxd, int mida, int minp, int maxp);
}

struct wavefront_work_s {
	Work_Data *w;
	Align_Spec *s;
	Path apath, bpath;
};

struct wavefront_work_s *wavefront_init_work(void)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)malloc(sizeof(struct wavefront_work_s));

	Work_Data *w = New_Work_Data();
	enlarge_vector(w, 10000 * (6 * sizeof(int) + sizeof(uint64_t)));
	enlarge_points(w, 4 * 4096 / 100 * sizeof(uint16_t) + sizeof(Path));


	float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
	Align_Spec *s = New_Align_Spec(0.6, 100, freq);


	work->w = w;
	work->s = s;
	work->apath.trace = malloc(4096);
	work->bpath.trace = malloc(4096);

	return(work);
}


void wavefront_clean_work(struct wavefront_work_s *_work)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)_work;

	free(work->apath.trace);
	free(work->bpath.trace);

	Free_Align_Spec(work->s);
	Free_Work_Data(work->w);

	free(work);
	return;	
}

int 
wavefront(
	void *_work,
	char const *a,
	uint64_t alen,
	char const *b,
	uint64_t blen)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)_work;
	Work_Data *w = work->w;
	Align_Spec *s = work->s;

	Alignment aln;
	aln.path = &work->apath;
	aln.flags = 0;
	aln.aseq = (char *)a;
	aln.bseq = (char *)b;
	aln.alen = alen;
	aln.blen = blen;

	int low = 0;
	int high = 0;

	forward_wave(w, s, &aln, &work->bpath, &low, high, 0, -INT32_MAX, INT32_MAX);
	// fprintf(stderr, "(%u, %u), (%u, %u)\n", alen, blen, aln.path->aepos, aln.path->bepos);
	return(0);
}




/**
 * random sequence generator, modifier.
 * rseq generates random nucleotide sequence in ascii,
 * mseq takes ascii sequence, modifies the sequence in given rate.
 */
static char rbase(void)
{
	switch(rand() % 4) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'A';
	}
}

char *rseq(int len)
{
	int i;
	char *seq;
	seq = (char *)malloc(sizeof(char) * (len + 1));

	for(i = 0; i < len; i++) {
		seq[i] = rbase();
	}
	seq[len] = '\0';
	return(seq);
}

char *mseq(char const *seq, int x, int ins, int del)
{
	int i;
	int len = strlen(seq);
	char *mod, *ptr;
	mod = (char *)malloc(sizeof(char) * 2 * (len + 1));

	ptr = mod;
	for(i = 0; i < len; i++) {
		if(rand() % x == 0) { *ptr++ = rbase(); }
		else if(rand() % ins == 0) { *ptr++ = rbase(); i--; }
		else if(rand() % del == 0) { /* skip a base */ }
		else { *ptr++ = seq[i]; }
	}
	*ptr = '\0';
	return(mod);
}

int print_msg(int flag, char const *fmt, ...)
{
	int r = 0;
	if(flag == 0) {
		va_list l;
		va_start(l, fmt);
		r = vfprintf(stdout, fmt, l);
		va_end(l);
	}
	return(r);
}

int print_bench(int flag, char const *name, int64_t l, int64_t a, int64_t sl, int64_t sa)
{
	if(flag == 0) {
		return(printf("%s\t%ld\t%ld\t%ld\t%ld\n", name, l / 1000, a / 1000, sl, sa));
	} else if(flag == 1) {
		return(printf("%ld\t%ld\t", l / 1000, a / 1000));
	} else if(flag == 2) {
		return(printf("%ld\t", a / 1000));
	}
	return(0);
}

int main(int argc, char *argv[])
{
	int flag = 0;
	if(argc > 1 && strcmp(argv[1], "-s") == 0) {
		flag = 1;
		argc--; argv++;
	} else if(argc > 1 && strcmp(argv[1], "-a") == 0) {
		flag = 2;
		argc--; argv++;
	}

	uint64_t i;
	int const m = 1, x = -1, gi = -1, ge = -1;
	int const xt = 30;
	char *a, *b, *at, *bt;
	bench_t bl, ba, sl, sa, ddl, dda, wl, fa;
	volatile int64_t sbl = 0, sba = 0, ssl = 0, ssa = 0, sddl = 0, sdda = 0, sfa = 0;
	struct timeval tv;

	int8_t score_matrix[16] __attribute__(( aligned(16) ));
	build_score_matrix(score_matrix, m, x);

	gettimeofday(&tv, NULL);
	unsigned long s = (argc > 3) ? atoi(argv[3]) : tv.tv_usec;
	srand(s);
	print_msg(flag, "%lu\n", s);

	/* malloc work */
	void *work = aligned_malloc(1024 * 1024 * 1024, sizeof(__m128i));

	if(argc > 1 && strcmp(argv[1], "-") != 0) {
		uint64_t len = (argc > 1) ? atoi(argv[1]) : 1000;
		uint64_t cnt = (argc > 2) ? atoi(argv[2]) : 1000;
		a = rseq(len * 9 / 10);
		b = mseq(a, 10, 40, 40);
		at = rseq(len / 10);
		bt = rseq(len / 10);
		a = (char *)realloc(a, 2*len); strcat(a, at); free(at);
		b = (char *)realloc(b, 2*len); strcat(b, bt); free(bt);
		print_msg(flag, "%p, %p\n", a, b);

		// print_msg(flag, "%s\n%s\n", a, b);

		print_msg(flag, "len:\t%d\ncnt:\t%d\n", len, cnt);
		print_msg(flag, "m: %d\tx: %d\tgi: %d\tge: %d\n", m, x, gi, ge);
		print_msg(flag, "alg\tlinear\taffine\tsc(l)\tsc(a)\n");

		/* blast */
		bench_init(bl);
		bench_init(ba);
		bench_start(bl);
		for(i = 0; i < cnt; i++) {
			sbl += blast_linear(work, a, strlen(a), b, strlen(b), score_matrix, ge, xt);
		}
		bench_end(bl);
		bench_start(ba);
		for(i = 0; i < cnt; i++) {
			sba += blast_affine(work, a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
		}
		bench_end(ba);
		print_bench(flag, "blast", bench_get(bl), bench_get(ba), sbl, sba);

		/* simdblast */
		bench_init(sl);
		bench_init(sa);
		bench_start(sl);
		for(i = 0; i < cnt; i++) {
			ssl += simdblast_linear(work, a, strlen(a), b, strlen(b), score_matrix, ge, xt);
		}
		bench_end(sl);
		bench_start(sa);
		for(i = 0; i < cnt; i++) {
			ssa += simdblast_affine(work, a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
		}
		bench_end(sa);
		print_bench(flag, "simdblast", bench_get(sl), bench_get(sa), ssl, ssa);


		/* adaptive banded */
		bench_init(ddl);
		bench_init(dda);
		bench_start(ddl);
		for(i = 0; i < cnt; i++) {
			sddl += ddiag_linear(work, a, strlen(a), b, strlen(b), score_matrix, ge, xt);
		}
		bench_end(ddl);
		bench_start(dda);
		for(i = 0; i < cnt; i++) {
			sdda += ddiag_affine(work, a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
		}
		bench_end(dda);
		print_bench(flag, "aband", bench_get(ddl), bench_get(dda), sddl, sdda);

		/* wavefront */
		bench_init(wl);
		if(len >= 1000) {
			/* wavefront algorithm fails to align sequences shorter than 1000 bases */

			struct wavefront_work_s *wwork = wavefront_init_work();

			bench_start(wl);
			for(i = 0; i < cnt; i++) {
				wavefront(wwork, a, strlen(a), b, strlen(b));
			}
			bench_end(wl);
		}
		print_bench(flag, "wavefront", bench_get(wl), bench_get(wl), 0, 0);


		if(flag != 0) {
			printf("\n");
		}
		free(a);
		free(b);

	} else {
		int c;
		double frac = (argc > 2) ? atof(argv[2]) : 1.0;
		uint64_t const rl = 100;

		kvec_t(char) buf;
		kvec_t(char *) seq;
		kvec_t(uint64_t) len;
		kvec_t(uint32_t) lscore, ascore;

		kv_init(buf);
		kv_init(seq);
		kv_init(len);
		kv_init(lscore);
		kv_init(ascore);

		uint64_t base = 0;
		while((c = getchar()) != EOF) {
			if(c == '\n') {
				uint64_t l = (uint64_t)(frac * (kv_size(buf) - base));
				for(i = 0; i < rl; i++) {
					kv_push(buf, rbase());
				}
				kv_push(len, l + rl);
				kv_push(seq, (char *)base);
				buf.n = base + l + rl;
				kv_push(buf, '\0');
				base = kv_size(buf);
			} else {
				kv_push(buf, c);
			}
		}
		kv_push(len, kv_size(buf) - base);
		kv_push(seq, (char *)base);
		kv_push(buf, '\0');

		for(i = 0; i < kv_size(seq); i++) {
			kv_at(seq, i) += (ptrdiff_t)buf.a;
		}

		/* collect scores with full-sized dp */
		for(i = 0; i < kv_size(seq) / 2; i++) {
			sw_result_t l = sw_linear(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, ge);
			kv_push(lscore, l.score);
			free(l.path);

			sw_result_t a = sw_affine(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge);
			kv_push(ascore, a.score);
			free(a.path);
		}

		/* blast */
		bench_init(bl);
		bench_init(ba);
		bench_start(bl);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			uint32_t s = blast_linear(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, ge, xt);
			sbl += s > 0.8 * kv_at(lscore, i);
		}
		bench_end(bl);
		bench_start(ba);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			uint32_t s = blast_affine(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge, xt);
			sba += s > 0.8 * kv_at(ascore, i);
		}
		bench_end(ba);
		print_bench(flag, "blast", bench_get(bl), bench_get(ba), sbl, sba);

		/* simdblast */
		bench_init(sl);
		bench_init(sa);
		bench_start(sl);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			uint32_t s = simdblast_linear(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, ge, xt);
			ssl += s > 0.8 * kv_at(lscore, i);
		}
		bench_end(sl);
		bench_start(sa);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			uint32_t s = simdblast_affine(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge, xt);
			ssa += s > 0.8 * kv_at(ascore, i);
		}
		bench_end(sa);
		print_bench(flag, "simdblast", bench_get(sl), bench_get(sa), ssl, ssa);


		/* adaptive banded */
		bench_init(ddl);
		bench_init(dda);
		bench_start(ddl);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			uint32_t s = ddiag_linear(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, ge, xt);
			sddl += s > 0.8 * kv_at(lscore, i);
		}
		bench_end(ddl);
		bench_start(dda);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			uint32_t s = ddiag_affine(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge, xt);
			sdda += s > 0.8 * kv_at(ascore, i);
		}
		bench_end(dda);
		print_bench(flag, "aband", bench_get(ddl), bench_get(dda), sddl, sdda);

		/* wavefront */
		bench_init(wl);
		struct wavefront_work_s *wwork = wavefront_init_work();

		bench_start(wl);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			wavefront(wwork, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1));
		}
		bench_end(wl);
		print_bench(flag, "wavefront", bench_get(wl), bench_get(wl), 0, 0);

		/* SSW library */
		/* convert sequence to number string */
		for(i = 0; i < kv_size(buf); i++) {
			if(kv_at(buf, i) == '\0') { continue; }
			char c = kv_at(buf, i);
			kv_at(buf, i) = (((c>>2) ^ (c>>1)) & 0x03) + 1;
		}

		bench_init(fa);
		for(i = 0; i < kv_size(seq) / 2; i++) {
			s_profile *sp = ssw_init((int8_t *)kv_at(seq, i * 2), kv_at(len, i * 2), score_matrix, 4, 1);

			bench_start(fa);
			s_align *r = ssw_align(sp, (int8_t *)kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), -gi, -ge, 8, 0, 0, 30);
			bench_end(fa);

			sfa += r->score1;
			align_destroy(r);
			init_destroy(sp);
		}
		print_bench(flag, "farrar", bench_get(fa), bench_get(fa), sfa, 0);

		if(flag != 0) {
			printf("\n");
		}

		kv_destroy(buf);
		kv_destroy(seq);
		kv_destroy(len);
	}

	free(work);
	return 0;
}

/**
 * end of main.cc
 */
