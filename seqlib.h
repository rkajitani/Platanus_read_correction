#ifndef SEQLIB_H
#define SEQLIB_H

#include "common.h"
#include "seqlib.h"

#define INS_DISTR_TRUNC 0.01

typedef enum {
    PE1,
    PE2,
    MP1,
    MP2,
} pair_file_format_t;

typedef struct {
    FILE *fp[2];
    seq_file_format_t seq_fmt;
    pair_file_format_t pair_fmt;
    int64_t id;
} seqlib_input_t;

typedef struct {
    FILE *pair_fp;
    FILE *ins_len_fp;
    FILE *mapped_fp;
    int64_t n_pair;
    int64_t total_len;
    int64_t ave_len;
    int64_t ave_ins;
    int64_t sd_ins;
    double ave_cov;
} seqlib_t;

void seqlib_init(seqlib_t *lib);
void seqlib_destroy(seqlib_t *lib);
void seqlib_estimate_ins(seqlib_t *lib);
int cmp_seqlib(const void *a, const void *b);
int cmp_seqlib_ptr(const void *a, const void *b);
int cmp_seqlib_input(const void *a, const void *b);
void distr_truncate(int64_t *distr, const int64_t size, const double edge);
double distr_ave(const int64_t *distr, const int64_t size);
double distr_sd(const int64_t *distr, const int64_t size, const double ave);
void seqlib_read_fasta_pe1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread);
void seqlib_read_fasta_mp1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread);
void seqlib_read_fasta_pe2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread);
void seqlib_read_fasta_mp2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread);
void seqlib_read_fastq_pe1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread);
void seqlib_read_fastq_mp1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread);
void seqlib_read_fastq_pe2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread);
void seqlib_read_fastq_mp2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread);
int option_pair(int argc, char **argv, int i, const pair_file_format_t pair_fmt, seqlib_input_t *val, int64_t *n_val);
int option_indexed_int(int argc, char **argv, int i, int64_t *val);

#endif
