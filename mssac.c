#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "libsais.h"
#include "libsais64.h"

#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define Malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define Calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define Realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define Grow(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = Realloc(type, (ptr), (__m)); \
		} \
	} while (0)

int ksa_sa(const unsigned char *T, int *SA, int n, int k);
int ksa_sa64(const unsigned char *T, int64_t *SA, int n, int k);

unsigned char seq_nt6_table[128];
void seq_char2nt6(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);
uint32_t SA_checksum(int64_t l, const int *s);
uint32_t SA_checksum64(int64_t l, const int64_t *s);
long peakrss(void);
double cputime(void);
double realtime(void);

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	kseq_t *seq;
	gzFile fp;
	int64_t l = 0, max = 0, n_sentinels = 0;
	int32_t c, algo, add_rev = 0, n_threads = 1;
	uint32_t checksum = 0;
	uint8_t *s = 0;
	double t_real, t_cpu;

	while ((c = ketopt(&o, argc, argv, 1, "a:rt:", 0)) >= 0) {
		if (c == 'r') add_rev = 1;
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'a') {
			if (strcmp(o.arg, "ksa64") == 0) algo = 1;
			else if (strcmp(o.arg, "ksa") == 0) algo = 2;
			else if (strcmp(o.arg, "sais64") == 0) algo = 3;
			else if (strcmp(o.arg, "sais") == 0) algo = 4;
			else {
				fprintf(stderr, "(EE) Unknown algorithm.\n");
				return 1;
			}
		}
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: mssa-bench [options] input.fasta\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -a STR    algorithm: ksa64, ksa, sais64 or sais [ksa64]\n");
#ifdef LIBSAIS_OPENMP
		fprintf(stderr, "  -t INT    number of threads for sais [%d]\n", n_threads);
#endif
		fprintf(stderr, "  -r        include reverse complement sequences\n");
		return 1;
	}

	// read FASTA/Q
	t_real = realtime();
	t_cpu = cputime();
	fp = gzopen(argv[o.ind], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		Grow(uint8_t, s, l + (seq->seq.l + 1), max);
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		memcpy(s + l, seq->seq.s, seq->seq.l + 1); // NB: we are copying 0
		l += seq->seq.l + 1;
		++n_sentinels;
		if (add_rev) {
			Grow(uint8_t, s, l + (seq->seq.l + 1), max);
			seq_revcomp6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
			++n_sentinels;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	printf("(MM) Read file in %.3f*%.3f sec (Peak RSS: %.3f MB)\n", realtime() - t_real, (cputime() - t_cpu) / (realtime() - t_real), peakrss() / 1024.0 / 1024.0);

	t_real = realtime();
	t_cpu = cputime();
	if (algo == 1) { // ksa64
		int64_t *SA = Malloc(int64_t, l);
		ksa_sa64(s, SA, l, 6);
		checksum = SA_checksum64(l, SA);
		free(SA); free(s);
	} else if (algo == 2) { // ksa
		int32_t *SA = Malloc(int32_t, l);
		ksa_sa(s, SA, l, 6);
		checksum = SA_checksum(l, SA);
		free(SA); free(s);
	} else if (algo == 3) { // libsais64
		int64_t i, k = 0, *tmp = Malloc(int64_t, l);
		for (i = 0; i < l; ++i)
			tmp[i] = s[i]? n_sentinels + s[i] : ++k;
		free(s);
		int64_t *SA = Malloc(int64_t, l + 10000);
#ifdef LIBSAIS_OPENMP
		if (n_threads > 1) {
			libsais64_long_omp(tmp, SA, l, n_sentinels + 6, 10000, n_threads);
		} else {
			libsais64_long(tmp, SA, l, n_sentinels + 6, 10000);
		}
#else
		libsais64_long(tmp, SA, l, n_sentinels + 6, 10000);
#endif
		checksum = SA_checksum64(l, SA);
		free(SA); free(tmp);
	} else if (algo == 4) { // libsais64
		int32_t i, k = 0, *tmp = Malloc(int32_t, l);
		for (i = 0; i < l; ++i)
			tmp[i] = s[i]? n_sentinels + s[i] : ++k;
		free(s);
		int32_t *SA = Malloc(int32_t, l + 10000);
#ifdef LIBSAIS_OPENMP
		if (n_threads > 1) {
			libsais_int_omp(tmp, SA, l, n_sentinels + 6, 10000, n_threads);
		} else {
			libsais_int(tmp, SA, l, n_sentinels + 6, 10000);
		}
#else
		libsais_int(tmp, SA, l, n_sentinels + 6, 10000);
#endif
		checksum = SA_checksum(l, SA);
		free(SA); free(tmp);
	} else {
		fprintf(stderr, "(EE) unknown algorithms\n");
		return 1;
	}
	printf("(MM) Generated SA in %.3f*%.3f sec (Peak RSS: %.3f MB; checksum: %x)\n", realtime() - t_real, (cputime() - t_cpu) / (realtime() - t_real), peakrss() / 1024.0 / 1024.0, checksum);
	return 0;
}

unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void seq_char2nt6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

uint32_t SA_checksum(int64_t len, const int32_t *s)
{
	uint32_t h = 2166136261U;
	int64_t i;
	for (i = 0; i < len; ++i)
		h ^= s[i], h *= 16777619;
	return h;
}

uint32_t SA_checksum64(int64_t len, const int64_t *s)
{
	uint32_t h = 2166136261U;
	int64_t i;
	for (i = 0; i < len; ++i)
		h ^= s[i], h *= 16777619;
	return h;
}

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#if defined(__linux__)
	return r.ru_maxrss * 1024;
#elif defined(__APPLE__)
	return r.ru_maxrss;
#endif
}

double realtime(void)
{
	static double realtime0 = -1.0;
	struct timeval tp;
	double t;
	gettimeofday(&tp, NULL);
	t = tp.tv_sec + tp.tv_usec * 1e-6;
	if (realtime0 < 0.0) realtime0 = t;
	return t - realtime0;
}
