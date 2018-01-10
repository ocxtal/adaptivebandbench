#define DEBUG

/**
 * @file main.cc
 */
#define __STDC_LIMIT_MACROS
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <sys/time.h>
#include "util.h"
#include "kvec.h"
#include "bench.h"
#include "parasail.h"
#include "ssw.h"
#include "full.h"

#define M 					( 1 )
#define X 					( 1 )
#define GI 					( 1 )
#define GE 					( 1 )
#define XDROP				( 70 )			// equal to the default of blastn (X = 100 (bit)) w/ (M, X, Gi, Ge) = (1, -1, 2, 1)

// #define OMIT_SCORE			1
#define PARASAIL_SCORE		1
#define RECALL_THRESH		0.99
// #define DEBUG_PATH
// #define DEBUG_BLAST			1

typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; void **a; } ptr_v;
typedef struct { size_t n, m; int32_t *a; } int32_v;

#define _base_signature			void *work, char const *a, uint64_t alen, char const *b, uint64_t blen, int8_t *score_matrix, int8_t gi, int8_t ge, int16_t xt, uint32_t bw
int scalar_affine(_base_signature);
int vertical_affine(_base_signature);
int diagonal_affine(_base_signature);
int striped_affine(_base_signature);
int blast_affine(_base_signature);
int simdblast_affine(_base_signature);
int adaptive_affine(_base_signature);

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

struct wavefront_work_s *wavefront_init_work(double id)
{
	struct wavefront_work_s *work = (struct wavefront_work_s *)malloc(sizeof(struct wavefront_work_s));

	Work_Data *w = New_Work_Data();
	enlarge_vector(w, 10000 * (6 * sizeof(int) + sizeof(uint64_t)));
	enlarge_points(w, 4 * 4096 / 100 * sizeof(uint16_t) + sizeof(Path));


	float freq[4] = { 0.25, 0.25, 0.25, 0.25 };
	Align_Spec *s = New_Align_Spec(id, 100, freq);


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
	return(aln.path->aepos + aln.path->bepos);
}

void print_alignment(char const *a, uint64_t alen, char const *b, uint64_t blen, char const *path)
{
	fprintf(stderr, "(%ld, %ld)\n", alen, blen);
	fprintf(stderr, "len(%lu)\n", strlen(path));

	for(uint64_t j = 0, k = 0; j < strlen(path); j++) {
		if(path[j] != 'I') {
			fprintf(stderr, "%c", a[k++]);
		} else {
			fprintf(stderr, "-");
		}
	}
	fprintf(stderr, "\n");
	for(uint64_t j = 0, k = 0; j < strlen(path); j++) {
		if(path[j] != 'D') {
			fprintf(stderr, "%c", b[k++]);
		} else {
			fprintf(stderr, "-");
		}
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "%s\n\n", path);
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

void revcomp(char *seq, int len)
{
	char map[32] = {
		0, 'T' & 0x1f, 0, 'G' & 0x1f, 0, 0, 0, 'C' & 0x1f, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 'A' & 0x1f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	};

	/* in-place revcomp */
	char *p = seq, *q = &seq[len - 1];
	while(p < q) {
		char t = *p;
		*p = map[*q & 0x1f] | (*q & ~0x1f);
		*q = map[t & 0x1f] | (t & ~0x1f);
		p++; q--;
	}

	if(p == q) { *p = map[*p]; }
	return;
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

int print_bench(int flag, char const *name, int64_t b, int64_t score)
{
	if(flag == 0) {
		return(printf("%s\t%ld\t%ld\n", name, b / 1000, score));
	} else if(flag == 1) {
		return(printf("%ld\n", b / 1000));
	} else if(flag == 2) {
		return(printf("%ld\t", b / 1000));
	}
	return(0);
}


struct params_s {
	int8_t score_matrix[16];
	int m, x, gi, ge, xt;
	uint32_t bw;
	uint64_t max_cnt, max_len;
	uint64_t flag, rdseed, pipe;
	char *list;

	uint8_v buf;
	ptr_v seq;
	uint64_v len;
	int32_v ascore;
	uint64_v apos;
	uint64_v bpos;

	void *work;
};

void init_args(struct params_s *p)
{
	build_score_matrix(p->score_matrix, M, -X);
	p->m = M; p->x = -X; p->gi = -GI; p->ge = -GE;
	p->xt = XDROP;
	p->bw = 32;
	p->max_len = 10000;
	p->max_cnt = 1000;
	p->flag = 0;
	p->rdseed = 0;
	p->pipe = 0;
	p->list = mm_strdup("scalar,vertical,diagonal,striped,adaptive,blast,simdblast");

	kv_init(p->buf);
	kv_init(p->seq);
	kv_init(p->len);
	kv_init(p->ascore);
	kv_init(p->apos);
	kv_init(p->bpos);

	/* malloc work */
	p->work = aligned_malloc(1024 * 1024 * 1024, sizeof(__m128i));

	struct timeval tv;
	gettimeofday(&tv, NULL);
	p->rdseed = tv.tv_usec;
	return;
}

void clean_args(struct params_s *p)
{
	free(p->list);
	free(p->buf.a);
	free(p->seq.a);
	free(p->len.a);
	free(p->ascore.a);
	free(p->apos.a);
	free(p->bpos.a);
	free(p->work);
	return;
}

int parse_args(struct params_s *p, int c, char *arg)
{
	switch(c) {
		case 'l': p->max_len = atoi(arg); break;
		case 'c': p->max_cnt = atoi(arg); break;
		case 's': p->flag = 1; break;
		case 'a': p->flag = 2; break;
		case 'n': free(p->list); p->list = mm_strdup(arg); break;
		case 'b': p->bw = atoi(arg); break;
		case 'x': p->xt = atoi(arg); break;
		case 'r': p->rdseed = atoi(arg); break;
		case 'i': p->pipe = 1; break;
	}
	return(0);
}

uint64_t read_seq(struct params_s *params)
{
	int c;
	uint64_t const rl = 100;

	uint64_t base = 0, i = 0;
	while((c = getchar()) != EOF) {
		if(c == '\n') {
			uint64_t l = MIN2(params->max_len, kv_size(params->buf) - base);
			if(kv_size(params->seq) & 0x02) {
				revcomp((char *)&kv_at(params->buf, base), kv_size(params->buf) - base);
			}

			for(i = 0; i < rl; i++) {
				kv_push(params->buf, rbase());
			}
			kv_push(params->len, l + rl);
			kv_push(params->seq, (void *)base);
			params->buf.n = base + l + rl;
			kv_push(params->buf, '\0');
			base = kv_size(params->buf);
		} else {
			kv_push(params->buf, c);
		}
		if(++i == 2 * params->max_cnt) { break; }
	}
	kv_push(params->len, kv_size(params->buf) - base);
	kv_push(params->seq, (char *)base);
	kv_push(params->buf, '\0');

	for(i = 0; i < kv_size(params->seq); i++) {
		kv_at(params->seq, i) = (uint8_t *)kv_at(params->seq, i) + (ptrdiff_t)params->buf.a;
	}
	return(i / 2);
}

uint64_t simulate_seq(struct params_s *params)
{
	for(uint64_t i = 0; i < params->max_cnt; i++) {
		char *a = rseq(params->max_len);
		char *b = mseq(a, 10, 40, 40);
		char *at = rseq(params->max_len / 10);
		char *bt = rseq(params->max_len / 10);

		a = (char *)realloc(a, 3 * params->max_len); strcat(a, at); free(at);
		b = (char *)realloc(b, 3 * params->max_len); strcat(b, bt); free(bt);

		kv_push(params->len, strlen(a));
		kv_push(params->len, strlen(b));

		kv_push(params->seq, (void *)kv_size(params->buf));
		kv_pushm(params->buf, a, strlen(a));
		kv_push(params->buf, '\0');

		kv_push(params->seq, (void *)kv_size(params->buf));
		kv_pushm(params->buf, b, strlen(b));
		kv_push(params->buf, '\0');

		free(a); free(b);
	}

	for(uint64_t i = 0; i < kv_size(params->seq); i++) {
		kv_at(params->seq, i) = (uint8_t *)kv_at(params->seq, i) + (ptrdiff_t)params->buf.a;
	}
	return(params->max_cnt);
}

void calc_score(struct params_s *params)
{
	kv_reserve(params->ascore, kv_size(params->seq) / 2);
	kv_reserve(params->apos, kv_size(params->seq) / 2);
	kv_reserve(params->bpos, kv_size(params->seq) / 2);
	#ifndef OMIT_SCORE
		parasail_matrix_t *_matrix = parasail_matrix_create("ACGT", params->m, params->x);
		#pragma omp parallel for
		for(uint64_t i = 0; i < kv_size(params->seq) / 2; i++) {
			parasail_result *r = parasail_sg_striped_sse41_128_16(
				(char const *)kv_at(params->seq, i * 2),     kv_at(params->len, i * 2),
				(char const *)kv_at(params->seq, i * 2 + 1), kv_at(params->len, i * 2 + 1),
				-params->gi - params->ge,
				-params->ge,
				_matrix
			);
			kv_at(params->ascore, i) = r->score;
			kv_at(params->apos, i) = r->end_query + 1;
			kv_at(params->bpos, i) = r->end_ref + 1;
			parasail_result_free(r);
		}
	#else
		for(i = 0; i < kv_size(params->seq) / 2; i++) {
			kv_at(params->ascore, i) = 0;
			kv_at(params->apos, i) = 0;
			kv_at(params->bpos, i) = 0;
		}
	#endif
	return;
}

struct mapping_s {
	char const *name;
	int (*fp)(_base_signature);
};
void bench_function(struct params_s *params, struct mapping_s *map, char const *name)
{
	uint32_t bw = params->bw, xt = params->xt;
	mm_split_foreach(name, ".", {
		switch(i) {
			case 1: bw = atoi(p); break;
			case 2: xt = atoi(p); break;
		}
	});

	int64_t score = 0;
	bench_t b;
	bench_init(b);
	for(uint64_t i = 0; i < kv_size(params->seq) / 2; i++) {
		bench_start(b);
		int32_t s = map->fp(params->work,
			(char const *)kv_at(params->seq, i * 2),     kv_at(params->len, i * 2),
			(char const *)kv_at(params->seq, i * 2 + 1), kv_at(params->len, i * 2 + 1),
			params->score_matrix,
			params->gi, params->ge,
			xt, bw
		);
		bench_end(b);
		score += s;

		debug("a(%s), b(%s)", kv_at(params->seq, i * 2), kv_at(params->seq, i * 2 + 1));

		sw_result_t a = sw_affine(
			(char const *)kv_at(params->seq, i * 2),     kv_at(params->len, i * 2),
			(char const *)kv_at(params->seq, i * 2 + 1), kv_at(params->len, i * 2 + 1),
			params->score_matrix, params->gi, params->ge
		);

		struct maxpos_s *mp = (struct maxpos_s *)params->work;
		fprintf(stderr, "i(%llu), score(%d, %d, %d), apos(%llu, %llu, %llu, %llu), bpos(%llu, %llu, %llu, %llu)\n",
			i, s, kv_at(params->ascore, i), a.score,
			mp->apos, kv_at(params->apos, i), a.apos, kv_at(params->len, i * 2),
			mp->bpos, kv_at(params->bpos, i), a.bpos, kv_at(params->len, i * 2 + 1));

		free(a.path);
	}
	bench_end(b);
	print_bench(params->flag, name, bench_get(b), score);
	return;
}

int main(int argc, char *argv[])
{
	/* name -> pointer mapping */
	#define fn(_name)	{ #_name, _name##_affine }
	struct mapping_s map[] = {
		/* static banded w/ standard matrix */
		fn(scalar), fn(vertical), fn(diagonal), fn(striped),
		/* non-standard banded */
		fn(blast), fn(simdblast), fn(adaptive)
	};
	#undef fn

	int i;
	struct params_s params __attribute__(( aligned(16) ));
	init_args(&params);
	while((i = getopt(argc, argv, "l:c:s:a:n:b:x:r:i:")) != -1) {
		if(parse_args(&params, i, optarg) != 0) { exit(1); }
	}

	srand(params.rdseed);
	print_msg(params.flag, "seed:%lu\tm: %d\tx: %d\tgi: %d\tge: %d\txdrop: %d\tbw: %d\tmax_len: %d\tmax_cnt: %d\n",
		params.rdseed,
		params.m, params.x, params.gi, params.ge,
		params.xt, params.bw,
		params.max_len, params.max_cnt
	);

	uint64_t cnt;
	if(argc > 1 && params.pipe != 0) {
		cnt = read_seq(&params);
	} else {
		cnt = simulate_seq(&params);
	}

	/* collect scores with full-sized dp */
	calc_score(&params);
	mm_split_foreach(params.list, ",", {
		for(uint64_t j = 0; j < sizeof(map) / sizeof(struct mapping_s); j++) {
			debug("%s, %s", p, map[j].name);
			if(strncmp(p, map[j].name, strlen(map[j].name)) == 0) {
				char name[l + 1];
				memcpy(name, p, l); name[l] = '\0';
				bench_function(&params, &map[j], name);
			}
		}
	});

	clean_args(&params);
	return(0);
}

#if 0
	if(argc > 1 && strcmp(argv[1], "-") != 0) {

		// print_msg(flag, "%s\n%s\n", a, b);

		print_msg(flag, "len:\t%d\ncnt:\t%d\n", len, cnt);
		print_msg(flag, "m: %d\tx: %d\tgi: %d\tge: %d\n", m, x, gi, ge);
		print_msg(flag, "alg\taffine\tsc(l)\tsc(a)\n");

		/* blast */
		bench_init(ba);
		bench_start(ba);
		for(i = 0; i < cnt; i++) {
			sba += blast_affine(work, a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
		}
		bench_end(ba);
		print_bench(flag, "blast", xt, bench_get(bl), bench_get(ba), sbl, sba);

		/* simdblast */
		bench_init(sa);
		bench_start(sa);
		for(i = 0; i < cnt; i++) {
			ssa += simdblast_affine(work, a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
		}
		bench_end(sa);
		print_bench(flag, "simdblast", xt, bench_get(sl), bench_get(sa), ssl, ssa);


		/* adaptive banded */
		bench_init(dda);
		bench_start(dda);
		for(i = 0; i < cnt; i++) {
			sdda += aband_affine_32(work, a, strlen(a), b, strlen(b), score_matrix, gi, ge, xt);
		}
		bench_end(dda);
		print_bench(flag, "aband", 32, bench_get(ddl), bench_get(dda), sddl, sdda);

		/* wavefront */
		bench_init(wl);
		if(len >= 1000) {
			/* wavefront algorithm fails to align sequences shorter than 1000 bases */

			struct wavefront_work_s *wwork = wavefront_init_work(0.6);

			bench_start(wl);
			for(i = 0; i < cnt; i++) {
				wavefront(wwork, a, strlen(a), b, strlen(b));
			}
			bench_end(wl);
		}
		print_bench(flag, "wavefront-0.6", 0, bench_get(wl), bench_get(wl), 0, 0);

		bench_init(wl);
		if(len >= 1000) {
			/* wavefront algorithm fails to align sequences shorter than 1000 bases */
			// default: an average correlation rate of -e (default 70%)
			struct wavefront_work_s *wwork = wavefront_init_work(0.7);

			bench_start(wl);
			for(i = 0; i < cnt; i++) {
				wavefront(wwork, a, strlen(a), b, strlen(b));
			}
			bench_end(wl);
		}
		print_bench(flag, "wavefront-0.7", 0, bench_get(wl), bench_get(wl), 0, 0);

		if(flag != 0) {
			printf("\n");
		}
		free(a);
		free(b);

	} else {


		/* blast */
		for(int xdrop = 40; xdrop <= 100; xdrop += 5) {
			if(flag > 0 && xdrop != 70) { continue; }

			bench_init(bl);
			bench_init(ba);
			bench_start(ba);
			sba = 0;
			for(i = 0; i < kv_size(seq) / 2; i++) {
				uint32_t s = blast_affine(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge, xdrop);
				sba += s >= RECALL_THRESH * kv_at(ascore, i);

				#if defined(DEBUG_PATH) && defined(DEBUG_BLAST)
				if(s < RECALL_THRESH * kv_at(ascore, i)) {
					sw_result_t a = sw_affine(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge);
					print_alignment(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), a.path);
					free(a.path);
				}
				#endif
			}
			bench_end(ba);
			print_bench(flag, "blast", xdrop, bench_get(bl), bench_get(ba), sbl, sba);
		}

		#ifndef DEBUG_PATH
			/* simdblast */
			for(int xdrop = 40; xdrop <= 100; xdrop += 5) {
				if(flag > 0 && xdrop != 70) { continue; }

				bench_init(sl);
				bench_init(sa);
				bench_start(sa);
				ssa = 0;
				for(i = 0; i < kv_size(seq) / 2; i++) {
					uint32_t s = simdblast_affine(work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge, xdrop);
					ssa += s >= RECALL_THRESH * kv_at(ascore, i);
				}
				bench_end(sa);
				print_bench(flag, "simdblast", xdrop, bench_get(sl), bench_get(sa), ssl, ssa);
			}
		#endif

		/* adaptive banded */
		int (*aband_affine_table[])(
			void *work,
			char const *a,
			uint64_t alen,
			char const *b,
			uint64_t blen,
			int8_t score_matrix[16], int8_t gi, int8_t ge, int16_t xt) = {

			aband_affine_32,
			aband_affine_40,
			aband_affine_48,
			aband_affine_56,
			aband_affine_64,
			aband_affine_72,
			aband_affine_80,
			aband_affine_88,
			aband_affine_96,
			aband_affine_104,
			aband_affine_112,
			aband_affine_120,
			aband_affine_128,
			aband_affine_136,
			aband_affine_144,
			aband_affine_152,
			aband_affine_160,
			aband_affine_168,
			aband_affine_176,
			aband_affine_184,
			aband_affine_192,
			aband_affine_200,
			aband_affine_208,
			aband_affine_216,
			aband_affine_224,
			aband_affine_232,
			aband_affine_240,
			aband_affine_248,
			aband_affine_256
		};

		for(int j = 0; j < sizeof(aband_affine_table) / sizeof(void *); j++) {
			if(flag > 0 && j != 8) { continue; }

			bench_init(ddl);
			bench_init(dda);
			bench_start(dda);
			sdda = 0;
			for(i = 0; i < kv_size(seq) / 2; i++) {
				uint32_t s = aband_affine_table[j](work, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge, xt);
				sdda += s >= RECALL_THRESH * kv_at(ascore, i);

				#if defined(DEBUG_PATH) && !defined(DEBUG_BLAST)
				if(s < RECALL_THRESH * kv_at(ascore, i)) {
					sw_result_t a = sw_affine(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), score_matrix, gi, ge);
					print_alignment(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), a.path);
					free(a.path);
				}
				#endif
			}
			bench_end(dda);
			print_bench(flag, "aband", j * 8 + 32, bench_get(ddl), bench_get(dda), sddl, sdda);
		}

		#ifndef DEBUG_PATH
			/* wavefront */
			bench_init(wl);
			struct wavefront_work_s *wwork = wavefront_init_work(0.7);
			if(max_len >= 1) {
				bench_start(wl);
				for(i = 0; i < kv_size(seq) / 2; i++) {
					uint32_t l = wavefront(wwork, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1));
					swl += l >= RECALL_THRESH * (kv_at(len, i * 2) + kv_at(len, i * 2 + 1) - 256);		// 128 bp margins for the both sequences
				}
				bench_end(wl);
			}
			wavefront_clean_work(wwork);
			print_bench(flag, "wavefront-0.7", 0, bench_get(bl), bench_get(wl), 0, swl);

			/* wavefront */
			bench_init(wl);
			wwork = wavefront_init_work(0.6);
			if(max_len >= 100) {
				bench_start(wl);
				for(i = 0; i < kv_size(seq) / 2; i++) {
					uint32_t l = wavefront(wwork, kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1));
					swl += l >= RECALL_THRESH * (kv_at(len, i * 2) + kv_at(len, i * 2 + 1) - 256);		// 128 bp margins for the both sequences
				}
				bench_end(wl);
			}
			wavefront_clean_work(wwork);
			print_bench(flag, "wavefront-0.6", 0, bench_get(bl), bench_get(wl), 0, swl);

			/* parasail */
			parasail_matrix_t *matrix = parasail_matrix_create("ACGT", M, -X);

			bench_init(pa);
			for(i = 0; i < kv_size(seq) / 2; i++) {
				bench_start(pa);
				parasail_result *r = parasail_sg_striped_sse41_128_16(kv_at(seq, i * 2), kv_at(len, i * 2), kv_at(seq, i * 2 + 1), kv_at(len, i * 2 + 1), -gi-ge, -ge, matrix);
				bench_end(pa);
				spa += r->score >= RECALL_THRESH * kv_at(ascore, i);
				parasail_result_free(r);
			}
			print_bench(flag, "parasail", 0, bench_get(bl), bench_get(pa), 0, spa);


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

				sfa += r->score1 >= RECALL_THRESH * kv_at(ascore, i);
				align_destroy(r);
				init_destroy(sp);
			}
			print_bench(flag, "farrar", 0, bench_get(bl), bench_get(fa), 0, sfa);
		#endif

		if(flag != 0) {
			printf("\n");
		}

		kv_destroy(buf);
		kv_destroy(seq);
		kv_destroy(len);
		kv_destroy(ascore);
	}
	return 0;
}
#endif
/**
 * end of main.cc
 */
