#ifndef CORRECT_H
#define CORRECT_H

#include "common.h"

#define MAX_NUM_K 100
#define MAX_CRR_K 32

typedef struct {
	char o[MAX_FILE_LEN + 1];
	int64_t t;
	int64_t m;
	int64_t n_fa;
	int64_t n_fq;
	char *fa_name[MAX_FILE_NUM];
	char *fq_name[MAX_FILE_NUM];
	FILE *fa_fp[MAX_FILE_NUM];
	FILE *fq_fp[MAX_FILE_NUM];
	int64_t n_k;
	int64_t k[MAX_NUM_K];
	int64_t q;
	double e;
} option_correct_t;

typedef struct {
	FILE *seq_file;
	uint64_t id;
	uint64_t mem_lim;
	uint64_t th;
	uint64_t n_corrected;
	uint64_t k_len;
	uint64_t k_mask;
	uint64_t index_len;
	uint64_t table_size;
	uint64_t *key_table;
	uint16_t *val_table;
	uint64_t n_file;
	uint64_t n_read[MAX_FILE_NUM * 2];
	double max_edit;
	char fname[MAX_FILE_NUM * 2][MAX_FILE_LEN * 2 + 1];
} corrector_t;

typedef struct {
	uint64_t n_thread;
	uint64_t k_len;
	uint64_t index_len;
	uint64_t table_size;
	uint64_t *key_table;
	uint16_t *val_table;
	uint64_t th;
	uint64_t n_corrected;
	uint64_t mem_lim;
	uint64_t fid;
	uint64_t n_file;
	uint64_t n_read[MAX_FILE_NUM * 2];
	char fname[MAX_FILE_NUM * 2][MAX_FILE_LEN + 5];
	double max_edit;
	corrector_t ct[MAX_THREAD];
	pthread_t thread[MAX_THREAD];
} corrector_threads_t;

void print_correct_usage(void);
void option_correct_init(option_correct_t *opt);
void option_correct_destroy(option_correct_t *opt);
int option_correct_parse(option_correct_t *opt, int argc, char **argv);
void correct(option_correct_t *opt);
void correct_mt(option_correct_t *opt);
void corrector_init(corrector_t *ct, uint64_t mem_lim);
void corrector_destroy_table(corrector_t *ct);
void corrector_destroy(corrector_t *ct);
void corrector_read_fasta(corrector_t *ct, FILE *fp);
void corrector_read_fastq(corrector_t *ct, FILE *fp, int8_t qual_cut);
void corrector_make_table(corrector_t *ct, uint64_t k_len);
void corrector_ungap_correct(corrector_t *ct);
void corrector_gapped_correct(corrector_t *ct);
void corrector_show_seq(corrector_t *ct);
void corrector_set_max_edit(corrector_t *ct, double max_edit);
void corrector_threads_init(corrector_threads_t *ctt, uint64_t n_thread, uint64_t mem_lim);
void corrector_threads_destroy_table(corrector_threads_t *ctt);
void corrector_threads_destroy(corrector_threads_t *ctt);
void corrector_threads_read_fasta(corrector_threads_t *ctt, FILE *fp);
void corrector_threads_read_fastq(corrector_threads_t *ctt, FILE *fp, int8_t qual_cut);
void corrector_threads_make_table(corrector_threads_t *ctt, uint64_t k_len);
void corrector_threads_count(corrector_t *ct);
void corrector_threads_recount(corrector_t *ct);
void corrector_threads_finish_table(corrector_t *ct);
void corrector_threads_show_seq(corrector_threads_t *ctt);
void corrector_threads_set_max_edit(corrector_threads_t *ctt, double max_edit);
void corrector_threads_ungap_correct(corrector_threads_t *ctt);
void corrector_threads_gapped_correct(corrector_threads_t *ctt);

static inline uint16_t corrector_occ(corrector_t *ct, uint64_t key) 
{
	uint64_t index;
	uint64_t step;

	index = hash64(key, ct->index_len) & (ct->table_size - 1);
	if (!(ct->val_table[index]))
		return 0;
	else if (key == ct->key_table[index])
		return ct->val_table[index];
	step = rehash64(key, ct->index_len);
	index = (index + step) & (ct->table_size-1);
	while (ct->val_table[index]) {
		if (key == ct->key_table[index])
			return ct->val_table[index];
		index = (index+step) & (ct->table_size-1);
	}
	return 0;
}

#endif
