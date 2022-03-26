#ifndef GAP_CLOSE_H
#define GAP_CLOSE_H

#include "common.h"
#include "mcounter.h"
#include "mgraph.h"
#include "lcounter.h"
#include "lgraph.h"
#include "seqlib.h"
#include "mapper.h"

#define CLOSER_MIN_K 24
#define CLOSER_MAX_K 72
#define CLOSER_MIN_COV 3
#define QASSEMBLER_BRANCH_TH 0.5
#define QASSEMBLER_BUBBLE_TH 0.1

typedef struct {
	char o[MAX_FILE_LEN + 1];
	int64_t m;
	int64_t t;
	int64_t s;
	int64_t v;
	double e;
	int64_t n_pair;
	seqlib_input_t pair[MAX_FILE_NUM + 1];
	int64_t n_c;
	FILE *c[MAX_FILE_NUM];
} option_gap_close_t;

typedef struct {
	int32_t id;
	int32_t st;
	int32_t ed;
	int32_t len;
	int8_t *seq;
} gap_t;

typedef struct {
	uint64_t key;
	int64_t val;
} gap_rec_t;

typedef struct {
	int64_t mem_lim;
	int64_t mem_usage;
	int64_t min_ol;
	double max_mis_rate;
	int64_t n_sca;
	lseq_t *sca;
	int64_t n_gap;
	gap_t *gap;
	int64_t index_len;
	int64_t gap_table_size;
	gap_rec_t *gap_table;
} closer_t;

typedef struct {
	mcounter_t mc;
	mgraph_t mg;
	lcounter_t lc;
	lgraph_t lg;
	int64_t min_k;
	int64_t max_k;
	int64_t min_cov;
} qassembler_t;

void print_gap_close_usage(void);
void option_gap_close_init(option_gap_close_t *opt);
void option_gap_close_destroy(const option_gap_close_t *opt);
uint8_t option_gap_close_parse(option_gap_close_t *opt, int argc, char **argv);
void gap_close_mt(option_gap_close_t *opt);
void closer_init(closer_t *cl, const contig_t *con, const int64_t min_ol, const double max_mis_rate, const int64_t mem_lim);
void closer_destroy(closer_t *cl);
void closer_make_gap_table(closer_t *cl);
void closer_insert_gap(closer_t *cl, const int32_t id, const int32_t ofst, const int64_t val);
int64_t closer_get_gap_id(const closer_t *cl, const int32_t id, const int32_t ofst);
void closer_print(const closer_t *cl, FILE *close_seq_fp, const char *out_name);
void closer_print_seq(closer_t *cl, const char *out_name);
FILE *closer_save_gap_covering_reads_mt(const closer_t *cl, const seqlib_t *lib, const int64_t n_thread);
void closer_load_local_reads(closer_t *cl, FILE *gap_seq_fp);
void closer_local_assemble_mt(const closer_t *cl, const int64_t n_thread);
bool make_consensus_from_reads(const lseq_t *seq, const int64_t n_seq, const double threshold, lseq_t *cons);
void closer_save_unused_reads(closer_t *cl, FILE **unused_fpp);
void closer_load_unused_reads(closer_t *cl, FILE *unused_fp);
void qassembler_init(qassembler_t *qa, const int64_t min_k, const int64_t max_k, const int64_t min_cov);
void qassembler_destroy(qassembler_t *qa);
void qassembler_assemble(qassembler_t *qa, const int8_t *seq, const int64_t len);
void qassembler_init_lgraph(qassembler_t *qa);
int64_t qassembler_close_gap(const qassembler_t *qa, const lseq_t *sca, gap_t *gap, const int64_t min_ol, const double max_mis_rate);
void merge_mgraph_with_lgraph(const mgraph_t *mg, const lgraph_t *lg, lcounter_t *lc);
int64_t get_exact_overlap(const int8_t *lseq, const int64_t llen, const int8_t *rseq, const int64_t rlen, const int64_t min_ol);
int64_t get_similer_overlap(const int8_t *lseq, const int64_t llen, const int8_t *rseq, const int64_t rlen, const int64_t min_ol);
int cmp_gap_ptr(const void *a, const void *b);

static inline bool gap_rec_is_emp(gap_rec_t gap)
{
	return gap.key == 0ull;
}

#endif
