#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <omp.h>

#define MAX_READ_LEN 50000
#define MAX_THREAD 100
#define MAX_EXTENSION_NAME 100
#define GRANULE 128
#define SMOOTHING_WINDOW 7
#define MEM_MARGIN 0.5
#define LINE_LENGTH 80 
#define MAX_LOAD_FACTOR 0.8
#define MAX_FILE_NUM 100
#define MAX_FILE_LEN 200
#define FASTQ_QUALITY_BASE 64
#define VARR_BUF_SIZE 100000ull
#define GIBIBYTE 1073741824ull

typedef enum {
	FASTA,
	FASTQ,
} seq_file_format_t;

typedef struct {
	FILE *fp;
	seq_file_format_t fmt;
} single_seq_file_t;

typedef struct {
    uint16_t len;
    uint8_t base[MAX_READ_LEN];
    uint16_t n_unknown;
    uint16_t unknown_pos[MAX_READ_LEN+1];
} seq_t;

typedef struct {
    uint64_t fwd;
    uint64_t rev;
} kmer_t;

typedef struct {
    int8_t *seq;
    int32_t len;
} lseq_t;

typedef struct {
    int32_t ofst;
    int32_t id;
} position_t;

typedef struct {
    FILE *seq_fp;
    FILE *cov_fp;
    lseq_t *seq;
    int64_t mem_usage;
    int64_t n_seq;
    int64_t seq_pool_size;
    int8_t *seq_pool;
    uint16_t *cov;
} contig_t;

typedef struct {
    uint64_t size;
    uint64_t capa;
    uint8_t *ele;
} varr8_t;

typedef struct {
    uint64_t size;
    uint64_t capa;
    uint64_t *ele;
} varr64_t;

extern uint64_t g_l_len;

void print_version(void);
void *my_malloc(size_t size);
void *my_calloc(size_t n, size_t size);
void *my_realloc(void *ptr, size_t size);
void my_free(void *ptr);
void all_free(void);
void my_exit(const uint16_t i);
void check_alloc(const void *p);
void check_file_open(const FILE *fp, const char *fname);
void check_tmpfile(const FILE *fp);
void check_mem_usage(const uint64_t usage, const uint64_t limit);
void put_command_log(int argc, char **argv, FILE *stream);
FILE *tmpfile_open(void);
void rev_comp(uint8_t *seq, const int64_t len);
int suffix_atoi(char *str);
uint64_t option_int(int argc, char **argv, uint64_t i, const int64_t min, const int64_t max, int64_t *val);
uint64_t option_multi_int(int argc, char **argv, uint64_t i, const int64_t min, const int64_t max, int64_t *val, int64_t *n_val);
uint64_t option_string(int argc, char **argv, uint64_t i, const int64_t max, char *val);
uint64_t option_float(int argc, char **argv, uint64_t i, const double min, const double max, double *val);
uint64_t option_multi_file(int argc, char **argv, uint64_t i, single_seq_file_t *val, int64_t *n_val);
uint64_t option_multi_file_name(int argc, char **argv, uint64_t i, FILE **fp_val, char **name_val, int64_t *n_val);
uint64_t option_multi_fasta_file(int argc, char **argv, uint64_t i, FILE **val, int64_t *n_val);
uint64_t get_kmer_cov_threshold(uint64_t *distr, uint64_t size);
seq_file_format_t check_seq_file_format(char *name);
void read_fasta(FILE *in, FILE *out);
void read_fastq(FILE *in, FILE *out);
double get_ave(const uint64_t *distr, const uint64_t min, const uint64_t max);
uint64_t get_left_minimal(const uint64_t *distr, const uint64_t max);
uint64_t get_left_minimal_smooth(const uint64_t *distr, const uint64_t max, const uint64_t wsize);
void read_fasta_mt(FILE *in, FILE **out, const uint64_t n_thread);
void read_fastq_mt(FILE *in, FILE **out, const uint64_t n_thread);
int cmp_int64(const void *a, const void *b);
int cmp_kmer_uniq(const void *a, const void *b);
int cmp_lmer_uniq(const void *a, const void *b);
int cmp_position(const void *a, const void *b);
void print_contig(FILE *contig_fp, char *out_name, const double occ_ratio);
void varr8_init(varr8_t *varr);
void varr8_destroy(varr8_t *varr);
void contig_init(contig_t *con);
void contig_destroy(contig_t *con);
void contig_read_fasta(contig_t *con, FILE *fp);
void contig_read_fasta_cov(contig_t *con, FILE *fp);
void contig_load_seq(contig_t *con);
void contig_load_cov(contig_t *con);

static inline int8_t ascii2num(const int8_t c)
{
    static const int8_t array[] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
    };
    return array[c];
}

static inline int8_t num2ascii(const int8_t c)
{
    static const int8_t array[] = {'A', 'C', 'G', 'T', 'N'};
    return array[c];
}

static inline int64_t max_64(const int64_t a, const int64_t b)
{
    return a > b ? a : b;
}

static inline int64_t min_64(const int64_t a, const int64_t b)
{
    return a < b ? a : b;
}

static inline uint64_t min_u64(const uint64_t a, const uint64_t b)
{
    return a < b ? a : b;
}

static inline int64_t abs64(const int64_t a)
{
    return a < 0 ? -a : a;
}

static inline void seq_write(seq_t *seq, FILE *seq_file)
{
    fwrite(&seq->n_unknown, sizeof(uint16_t), 1, seq_file); 
    fwrite(seq->unknown_pos, sizeof(uint16_t), seq->n_unknown, seq_file); 
    fwrite(&seq->len, sizeof(uint16_t), 1, seq_file); 
    fwrite(seq->base, sizeof(uint8_t), seq->len, seq_file); 
}

static inline int seq_read(seq_t *seq, FILE *seq_file)
{
    if (!fread(&seq->n_unknown, sizeof(uint16_t), 1, seq_file))
        return 0;
    fread(seq->unknown_pos, sizeof(uint16_t), seq->n_unknown, seq_file);
    fread(&seq->len, sizeof(uint16_t), 1, seq_file); 
    fread(seq->base, sizeof(uint8_t), seq->len, seq_file);
    return 1;
}

static inline uint64_t hash64(const uint64_t key, const uint64_t index_len)
{
    return (~key + (key >> index_len) + (~key >> 2*index_len));
}

static inline uint64_t rehash64(const uint64_t key, const uint64_t index_len)
{
    return ((~key + (key >> (index_len+1))) ^ (~key >> (2*index_len-1))) | 1ull;
}

static inline uint64_t flag_degree(const uint8_t flag)
{
    static const uint64_t array[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
    return array[flag];
}

static inline uint64_t flag_base(const uint8_t flag)
{
    static const uint64_t array[] = {4,0,1,4,2,4,4,4,3,4,4,4,4,4,4,4};
    return array[flag];
}

static inline void bseq_set(uint8_t *seq, const uint64_t pos, const uint8_t val)
{
    seq[pos >> 2] = (seq[pos >> 2] & ~(3 << ((pos&3) << 1))) | (val << ((pos&3) << 1));
}

static inline uint8_t bseq_get(uint8_t *seq, const uint64_t pos)
{
    return (seq[pos >> 2] >> ((pos&3) << 1)) & 3;
}

static inline void varr8_push(varr8_t *varr, const uint8_t val)
{
    if (varr->size == varr->capa) {
        varr->capa += VARR_BUF_SIZE;
        varr->ele = (uint8_t *)my_realloc(varr->ele, varr->capa * sizeof(uint8_t)); 
        check_alloc(varr->ele);
    }
    varr->ele[varr->size] = val;
    ++varr->size;
}

static inline uint8_t varr8_peek(const varr8_t *varr)
{
    return varr->ele[varr->size-1];
}

static inline void varr8_pop(varr8_t *varr)
{
    --varr->size;
}

static inline void varr8_clear(varr8_t *varr)
{
    varr->size = 0ull;
}

static inline void varr8_resize(varr8_t *varr, const uint64_t size)
{
    if (size > varr->capa) {
        varr->capa = size + VARR_BUF_SIZE;
        varr->ele = (uint8_t *)my_realloc(varr->ele, varr->capa * sizeof(uint8_t)); 
        check_alloc(varr->ele);
    }
    varr->size = size;
}

void varr64_init(varr64_t *varr);
void varr64_destroy(varr64_t *varr);

static inline void varr64_push(varr64_t *varr, const uint64_t val)
{
    if (varr->size == varr->capa) {
        varr->capa += VARR_BUF_SIZE;
        varr->ele = (uint64_t *)my_realloc(varr->ele, varr->capa * sizeof(uint64_t)); 
    }
    varr->ele[varr->size] = val;
    ++varr->size;
}

static inline uint64_t varr64_peek(const varr64_t *varr)
{
    return varr->ele[varr->size-1];
}

static inline void varr64_pop(varr64_t *varr)
{
    --varr->size;
}

static inline void varr64_clear(varr64_t *varr)
{
    varr->size = 0ull;
}

static inline void varr64_resize(varr64_t *varr, const uint64_t size)
{
    if (size > varr->capa) {
        varr->capa = size + VARR_BUF_SIZE;
        varr->ele = (uint64_t *)my_realloc(varr->ele, varr->capa * sizeof(uint64_t)); 
    }
    varr->size = size;
}

static inline int fskip(const int del, FILE *in)
{
	int c;
	while ((c = getc(in)) != del && c != EOF)
		;
	return c;
}

static inline int fastq_skip_qual(uint16_t len, FILE *in)
{
	while (getc(in) == '\n' || len-- > 0)
	   ;
	return getc(in);
}

static inline int fget_seq(seq_t *seq, const int del, FILE *in)
{
	int c;

	seq->len = seq->n_unknown = 0;
	while ((c = getc(in)) != del && c != EOF) {
		if (c == '\n')
			continue;
		if (seq->len >= MAX_READ_LEN) {
			fputs("error: too long sequence!\n", stderr);
			my_exit(1);
		}
		if (ascii2num(c) == 4) {
			seq->unknown_pos[seq->n_unknown] = seq->len;
			++seq->n_unknown;
		}
		else
			seq->base[seq->len] = ascii2num(c);
		++(seq->len);
	}
	return c;
}

static inline int fget_seq_simple(seq_t *seq, const int del, FILE *in)
{
	int c;

	seq->len = 0;
	while ((c = getc(in)) != del && c != EOF) {
		if (c == '\n')
			continue;
		seq->base[seq->len] = ascii2num(c);
		++(seq->len);
	}
	return c;
}

#endif
