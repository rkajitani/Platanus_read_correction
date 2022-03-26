#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include "common.h"
#include "mcounter.h"
#include "mgraph.h"
#include "lcounter.h"
#include "lgraph.h"

typedef struct {
	char o[MAX_FILE_LEN + 1];
	int64_t t;
	int64_t m;
	int64_t k;
	int64_t s;
	int64_t n;
	int64_t c;
	double a;
	double u;
	double d;
	int64_t n_f;
	single_seq_file_t f[MAX_FILE_NUM];
} option_assemble_t;

void print_assemble_usage(void);
void option_assemble_init(option_assemble_t *opt);
void option_assemble_destroy(const option_assemble_t *opt);
uint8_t option_assemble_parse(option_assemble_t *opt, int argc, char **argv);
void show_distribution(const uint64_t *distr, const uint64_t min, const uint64_t max, const char *out_name);
void show_kmer_extension(const double min_log_p_join, const double ave_cov, const double ave_len, const uint64_t min_cov, uint64_t k, const uint64_t len_step, uint64_t cov_cut);
void assemble(const option_assemble_t *opt);
void assemble_mt(const option_assemble_t *opt);
double decrease_coverage(const double cov, const double ave_len, const uint64_t large_k, const uint64_t small_k);
double log_prob_join(const uint64_t cov_cut, const double ave_cov, const double ave_len, const uint64_t large_k, const uint64_t small_k);
uint64_t decrease_cov_cut(const uint64_t cov_cut, const double ave_cov, const double ave_len, const double min_log_p_join, const uint64_t large_k, const uint64_t small_k);
uint64_t max_kmer_len(const double min_log_p_join, const double ave_cov, const double ave_len, const uint64_t min_cov, uint64_t k_len, const uint64_t len_step);

#endif
