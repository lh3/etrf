#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8

#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key) \
	typedef struct { \
		rstype_t *b, *e; \
	} rsbucket_##name##_t; \
	void rs_insertsort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		rstype_t *i; \
		for (i = beg + 1; i < end; ++i) \
			if (rskey(*i) < rskey(*(i - 1))) { \
				rstype_t *j, tmp = *i; \
				for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j) \
					*j = *(j - 1); \
				*j = tmp; \
			} \
	} \
	void rs_sort_##name(rstype_t *beg, rstype_t *end, int n_bits, int s) \
	{ \
		rstype_t *i; \
		int size = 1<<n_bits, m = size - 1; \
		rsbucket_##name##_t *k, b[1<<RS_MAX_BITS], *be = b + size; \
		assert(n_bits <= RS_MAX_BITS); \
		for (k = b; k != be; ++k) k->b = k->e = beg; \
		for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; \
		for (k = b + 1; k != be; ++k) \
			k->e += (k-1)->e - beg, k->b = (k-1)->e; \
		for (k = b; k != be;) { \
			if (k->b != k->e) { \
				rsbucket_##name##_t *l; \
				if ((l = b + (rskey(*k->b)>>s&m)) != k) { \
					rstype_t tmp = *k->b, swap; \
					do { \
						swap = tmp; tmp = *l->b; *l->b++ = swap; \
						l = b + (rskey(tmp)>>s&m); \
					} while (l != k); \
					*k->b++ = tmp; \
				} else ++k->b; \
			} else ++k; \
		} \
		for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; \
		if (s) { \
			s = s > n_bits? s - n_bits : 0; \
			for (k = b; k != be; ++k) \
				if (k->e - k->b > RS_MIN_SIZE) rs_sort_##name(k->b, k->e, n_bits, s); \
				else if (k->e - k->b > 1) rs_insertsort_##name(k->b, k->e); \
		} \
	} \
	void radix_sort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		if (end - beg <= RS_MIN_SIZE) rs_insertsort_##name(beg, end); \
		else rs_sort_##name(beg, end, RS_MAX_BITS, (sizeof_key - 1) * RS_MAX_BITS); \
	}

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct {
	int st, en, k, keep;
} etrf_elem_t;

#define etrf_elem_key(x) ((x).st)
KRADIX_SORT_INIT(etrf, etrf_elem_t, etrf_elem_key, 4)

static etrf_elem_t *trf_k(int *n_, int *m_, etrf_elem_t *a, int min_reg_len, int k, int len, const uint8_t *str)
{
	int i, n = *n_, m = *m_, streak = 0;
	for (i = k; i <= len; ++i) {
		if (i < len && str[i] == str[i - k] && str[i] < 4) {
			++streak;
		} else {
			if (streak >= k && streak + k >= min_reg_len) {
				if (n == m) {
					m = m < 16? 16 : m + (m<<1);
					a = (etrf_elem_t*)realloc(a, sizeof(etrf_elem_t) * m);
				}
				a[n].st = i - streak - k, a[n].en = i, a[n].keep = 0, a[n++].k = k;
			}
			streak = 0;
		}
	}
	*n_ = n, *m_ = m;
	return a;
}

static void select_reg(int n, etrf_elem_t *a, int max)
{
	int i, i0, i00, j;
	for (i = 0; i < n; ++i)
		if (a[i].en <= max) break;
	if (i == n) return;
	for (i00 = i0 = i, i = i + 1; i <= n; ++i) {
		if (a[i].en > max) continue;
		if (i == n || a[i].st >= a[i0].en) {
			for (j = i00 + 1; j < i0; ++j)
				if (a[i00].en <= a[j].st) break;
			if (i0 - j > 0)
				select_reg(i0 - j, &a[j], a[i0].st);
			a[i0].keep = 1;
			i00 = i0, i0 = i;
		} else {
			if (a[i].en - a[i].st > a[i0].en - a[i0].st || (a[i].en - a[i].st == a[i0].en - a[i0].st && a[i].k < a[i0].k))
				i0 = i, a[i0].keep = 0;
			else a[i].keep = 0;
		}
	}
}

static void get_motif(const uint8_t *seq, const etrf_elem_t *a, char *motif)
{
	int i, min_i = 0;
	for (i = 1; i < a->k; ++i)
		if (memcmp(&seq[a->st + i], &seq[a->st + min_i], a->k) < 0)
			min_i = i;
	for (i = 0; i < a->k; ++i)
		motif[i] = "ACGTN"[seq[a->st + min_i + i]];
	motif[a->k] = 0;
}

static void process_seq(int max_motif_len, int min_len, const char *name, int len, uint8_t *seq)
{
	int i, k, n_a = 0, m_a = 0;
	char *motif;
	etrf_elem_t *a = 0;
	for (i = 0; i < len; ++i)
		seq[i] = seq_nt4_table[seq[i]];
	for (k = 1; k <= max_motif_len; ++k)
		a = trf_k(&n_a, &m_a, a, min_len, k, len, seq);
	radix_sort_etrf(a, a + n_a);
	select_reg(n_a, a, INT_MAX);
	for (i = k = 0; i < n_a; ++i)
		if (a[i].keep) a[k++] = a[i];
	n_a = k;
	motif = (char*)malloc(max_motif_len + 1);
	for (i = 0; i < n_a; ++i) {
		get_motif(seq, &a[i], motif);
		printf("%s\t%d\t%d\t%d\t%s\n", name, a[i].st, a[i].en, a[i].k, motif);
	}
	free(motif);
	free(a);
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	ketopt_t o = KETOPT_INIT;
	int c, ret, min_len = 13, max_motif_len = 100;

	while ((c = ketopt(&o, argc, argv, 1, "l:m:", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
		else if (c == 'm') max_motif_len = atoi(o.arg);
	}

	if (o.ind == argc) {
		fprintf(stderr, "Usage: etrf [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT    max motif length [%d]\n", max_motif_len);
		fprintf(stderr, "  -l INT    min region length [%d]\n", min_len);
		return 1;
	}
	fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
	if (fp == 0) {
		fprintf(stderr, "ERROR: failed to open file '%s'\n", argv[o.ind]);
		return 1;
	}
	ks = kseq_init(fp);
	while ((ret = kseq_read(ks)) >= 0)
		process_seq(max_motif_len, min_len, ks->name.s, ks->seq.l, (uint8_t*)ks->seq.s);
	if (ret != -1)
		fprintf(stderr, "WARNING: FASTA/Q parser returns %d\n", ret);
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
