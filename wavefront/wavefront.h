

/* wrapper of Myers' wavefront algorithm */
#include "align.h"

int enlarge_vector(Work_Data *work, int newmax);
int enlarge_points(Work_Data *work, int newmax);
int forward_wave(Work_Data *work, Align_Spec *spec, Alignment *align, Path *bpath,
                        int *mind, int maxd, int mida, int minp, int maxp);
int reverse_wave(Work_Data *work, Align_Spec *spec, Alignment *align, Path *bpath,
                        int mind, int maxd, int mida, int minp, int maxp);

typedef struct {
	Work_Data *w;
	Align_Spec *s;
	Path apath, bpath;
} wavefront_work_t;

static inline
wavefront_work_t *wavefront_init_work(double id)
{
	wavefront_work_t *work = (wavefront_work_t *)malloc(sizeof(wavefront_work_t));

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

static inline
void wavefront_clean_work(wavefront_work_t *_work)
{
	wavefront_work_t *work = (wavefront_work_t *)_work;

	free(work->apath.trace);
	free(work->bpath.trace);

	Free_Align_Spec(work->s);
	Free_Work_Data(work->w);

	free(work);
	return;	
}

static inline
int wavefront(wavefront_work_t *work, char const *a, uint64_t alen, char const *b, uint64_t blen)
{
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
