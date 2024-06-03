#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "libsais64.h"

#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

int ksa_sa(const unsigned char *T, int *SA, int n, int k); // ksa, based on old sais
int ksa_sa64(const unsigned char *T, int64_t *SA, int n, int k); // ksa, based on old sais
int sais_int(const int *T, int *SA, int n, int k); // sais
int sais(const unsigned char *T, int *SA, int n); // sais
void suffixsort(int *x, int *p, int n, int k, int l); // qsufsort
int divsufsort(const unsigned char *T, int *SA, int n); // divsufsort
int ssort(int a[], int s[]); // ssort
void suffixArray(int* s, int* SA, int n, int K); // dc3
void SA_IS(unsigned char *s, int *SA, int n, int K, int cs); // the original SA-IS algorithm

unsigned char seq_nt6_table[128];
void seq_char2nt6(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);
uint32_t SA_checksum(int l, const int *s);
uint32_t SA_checksum64(int l, const int64_t *s);
long peakrss(void);
double cputime(void);
double realtime(void);

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	kseq_t *seq;
	gzFile fp;
	int *SA = 0, c, l = 0, max = 0, algo = 0, n_sentinels = 0, has_sentinel = 1, no_rev = 0;
	uint8_t *s = 0;
	double t_real, t_cpu;

	while ((c = ketopt(&o, argc, argv, 1, "a:xR", 0)) >= 0) {
		if (c == 'R') no_rev = 1;
		else if (c == 'x') has_sentinel = 0;
		else if (c == 'a') {
			if (strcmp(o.arg, "ksa") == 0) algo = 0;
			else if (strcmp(o.arg, "qsufsort") == 0) algo = 1;
			else if (strcmp(o.arg, "sais") == 0) algo = 2;
			else if (strcmp(o.arg, "divsufsort") == 0) algo = 3;
			else if (strcmp(o.arg, "ssort") == 0) algo = 4;
			else if (strcmp(o.arg, "dc3") == 0) algo = 5;
			else if (strcmp(o.arg, "is") == 0) algo = 6;
			else if (strcmp(o.arg, "ksa64") == 0) algo = 7;
			else if (strcmp(o.arg, "libsais64") == 0) algo = 8;
			else {
				fprintf(stderr, "(EE) Unknown algorithm.\n");
				return 1;
			}
		}
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: mssa-bench [options] input.fasta\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -a STR    algorithm: ksa, ksa64, sais, qsufsort, divsufsort, ssort, dc3, is [ksa]\n");
		fprintf(stderr, "  -x        do not regard a NULL as a sentinel\n");
		return 1;
	}

	t_real = realtime();
	t_cpu = cputime();
	fp = gzopen(argv[o.ind], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		if (l + (seq->seq.l + 1) * 2 + 1 >= max) {
			max = l + (seq->seq.l + 1) * 2 + 2;
			kroundup32(max);
			s = realloc(s, max);
		}
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		memcpy(s + l, seq->seq.s, seq->seq.l + 1);
		l += seq->seq.l + 1;
		++n_sentinels;
		if (!no_rev) {
			seq_revcomp6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
			++n_sentinels;
		}
	}
	s[l] = 0;
	kseq_destroy(seq);
	gzclose(fp);
	printf("(MM) Read file in %.3f*%.3f sec (Peak RSS: %.3f MB)\n", realtime() - t_real, (cputime() - t_cpu) / (realtime() - t_real), peakrss() / 1024.0 / 1024.0);

	t_real = realtime();
	t_cpu = cputime();
	if (has_sentinel) { // A NULL is regarded a sentinel
		if (algo == 0) { // ksa
			SA = (int*)malloc(sizeof(int) * l);
			ksa_sa(s, SA, l, 6);
			SA_checksum(l, SA);
			free(SA); free(s);
		} else if (algo == 7) { // ksa64
			int64_t *SA = (int64_t*)malloc(sizeof(int64_t) * l);
			ksa_sa64(s, SA, l, 6);
			SA_checksum64(l, SA); // TODO: not working for 64-bit arrays
			free(SA); free(s);
		} else if (algo == 8) { // libsais64
			int64_t i, *tmp, k = 0;
			tmp = (int64_t*)malloc(sizeof(int64_t) * (l + 3));
			for (i = 0; i < l; ++i)
				tmp[i] = s[i]? n_sentinels + s[i] : ++k;
			free(s);
			int64_t *SA = (int64_t*)malloc(sizeof(int64_t) * (l + 10000));
			libsais64_long(tmp, SA, l, n_sentinels + 6, 10000);
			SA_checksum64(l, SA);
			free(SA); free(tmp);
		} else if (algo == 1 || algo == 2 || algo == 4 || algo == 5) { // sais, qsufsort, ssort or dc3
			int i, *tmp, k = 0;
			tmp = (int*)malloc(sizeof(int) * (l + 3));
			for (i = 0; i < l; ++i)
				tmp[i] = s[i]? n_sentinels + s[i] : ++k;
			tmp[l] = tmp[l+1] = tmp[l+2] = 0; // required by dc3
			free(s);
			SA = (int*)malloc(sizeof(int) * (l + 1));
			if (algo == 1) { // qsufsort
				suffixsort(tmp, SA, l, n_sentinels + 6, 1);
				SA_checksum(l, SA + 1);
			} else if (algo == 2) { // sais
				sais_int(tmp, SA, l, n_sentinels + 6);
				SA_checksum(l, SA);
			} else if (algo == 4) { // ssort
				ssort(tmp, SA);
				SA_checksum(l, tmp + 1);
			} else if (algo == 5) { // dc3
				suffixArray(tmp, SA, l, n_sentinels + 5);
				SA_checksum(l, SA);
			}
			free(SA); free(tmp);
		} else {
			fprintf(stderr, "(EE) The selected algorithm cannot construct SA for strings with multiple sentinels.\n");
			return 1;
		}
	} else { // A NULL is regarded as an ordinary symbol
		if (algo == 0 || algo == 6) { // ksa or is
			int i;
			for (i = 0; i < l; ++i) ++s[i];
			SA = (int*)malloc(sizeof(int) * (l + 1));
			if (algo == 0) ksa_sa(s, SA, l + 1, 7);
			else SA_IS(s, SA, l + 1, 7, 1);
			SA_checksum(l, SA + 1);
			free(SA); free(s);
		} else if (algo == 1 || algo == 4 || algo == 5) { // qsufsort, ssort or dc3
			int i, *tmp;
			tmp = (int*)malloc(sizeof(int) * (l + 3));
			for (i = 0; i < l; ++i) tmp[i] = s[i] + 1;
			tmp[l] = tmp[l+1] = tmp[l+2] = 0; // required by dc3
			free(s);
			SA = (int*)malloc(sizeof(int) * (l + 1));
			if (algo == 1) {
				suffixsort(tmp, SA, l, 7, 1);
				SA_checksum(l, SA + 1);
			} else if (algo == 4) {
				ssort(tmp, SA);
				SA_checksum(l, tmp + 1);
			} else if (algo == 5) {
				suffixArray(tmp, SA, l, 6);
				SA_checksum(l, SA);
			}
			free(SA); free(tmp);
		} else if (algo == 2 || algo == 3) { // sais or divsufsort
			SA = (int*)malloc(sizeof(int) * (l + 1));
			if (algo == 2) sais(s, SA, l);
			else if (algo == 3) divsufsort(s, SA, l);
			SA_checksum(l, SA);
			free(SA); free(s);
		} else {
			fprintf(stderr, "(EE) not implemented yet\n");
			abort();
		}
	}
	printf("(MM) Generated SA in %.3f*%.3f sec (Peak RSS: %.3f MB)\n", realtime() - t_real, (cputime() - t_cpu) / (realtime() - t_real), peakrss() / 1024.0 / 1024.0);
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

uint32_t SA_checksum(int l, const int *s)
{
	uint32_t h = *s;
	const int *end = s + l;
	double t = cputime();
	for (; s < end; ++s) h = (h << 5) - h + *s;
	printf("(MM) Computed SA X31 checksum in %.3f seconds (checksum = %x)\n", cputime() - t, h);
	return h;
}

uint32_t SA_checksum64(int l, const int64_t *s)
{
	uint32_t h = *s;
	const int64_t *end = s + l;
	double t = cputime();
	for (; s < end; ++s) h = (h << 5) - h + *s;
	printf("(MM) Computed SA X31 checksum in %.3f seconds (checksum = %x)\n", cputime() - t, h);
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
