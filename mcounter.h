#ifndef MCOUNTER_H
#define MCOUNTER_H

#include "common.h"

typedef struct {
    FILE *kmer_fp;
    FILE *table_fp;
    uint64_t mem_lim;
    uint64_t len_distr[MAX_READ_LEN + 1];
    uint64_t max_occ;
    uint64_t k_len;
    uint64_t k_mask;
    uint64_t index_len;
    uint64_t table_size;
    uint64_t *key_table;
    uint16_t *val_table;
    uint64_t *occ_distr;
} mcounter_t;

typedef struct {
    FILE *unstored_fp;
    FILE *in_kmer_fp;
    mcounter_t *mc;
    pthread_mutex_t *mutex;
} arg_mcounter_lock_and_count_t;

void mcounter_init(mcounter_t *mc, const uint64_t mem_lim);
void mcounter_destroy_table(mcounter_t *mc);
void mcounter_destroy(mcounter_t *mc);
void mcounter_save_kmer_read(mcounter_t *mc, const uint64_t k_len, FILE *read_fp);
void mcounter_save_additional_kmer_read(mcounter_t *mc, const uint64_t k_len, FILE *read_fp);
void mcounter_load_kmer(mcounter_t *mc, const uint64_t min_occ);
void mcounter_load_contig(mcounter_t *mc, const uint64_t k_len, const double occ_ratio, FILE *contig_fp);
uint64_t mcounter_get_med_len(mcounter_t *mc);
uint64_t mcounter_get_ave_len(mcounter_t *mc);
double mcounter_get_ave_cov(mcounter_t *mc);
void mcounter_save_kmer_read_mt(mcounter_t *mc, const uint64_t k_len, FILE **read_fp, const uint64_t n_thread);
void mcounter_save_additional_kmer_read_mt(mcounter_t *mc, const uint64_t k_len, FILE **read_fp, const uint64_t n_thread);
void mcounter_lock_and_count1(const arg_mcounter_lock_and_count_t *arg);
void mcounter_lock_and_count2(const arg_mcounter_lock_and_count_t *arg);
uint64_t mcounter_sorted_kmer_list(const mcounter_t *mc, uint64_t ** buf);
void mcounter_extract_read(const mcounter_t *mc, FILE **read_fpp);
double mcounter_estimate_error_rate(mcounter_t *mc, const uint64_t k_len, FILE *read_fp);
double mcounter_estimate_error_rate_mt(mcounter_t *mc, const uint64_t k_len, FILE **read_fp, const uint64_t n_thread);
double mcounter_calculate_error_rate(const mcounter_t *mc1, const mcounter_t *mc2, const double ave_len);
uint64_t mcounter_get_coverage_cutoff(const mcounter_t *mc, const double err_rate);
FILE *sorted_key_from_kmer_file(FILE *kmer_fp, const uint64_t min_occ);
FILE *sorted_key_from_contig_file(FILE *contig_fp, const uint64_t k_len);
void mcounter_count_quick(mcounter_t *mc, const uint8_t *seq, const uint64_t len, const uint64_t k_len);
void mcounter_save_table(mcounter_t *mc);

static inline uint16_t mcounter_occ(const mcounter_t *mc, const uint64_t key) 
{
    uint64_t index;
    uint64_t step;

    index = hash64(key, mc->index_len) & (mc->table_size - 1);
    if (!(mc->val_table[index]))
        return 0;
    else if (key == mc->key_table[index])
        return mc->val_table[index];
    step = rehash64(key, mc->index_len);
    index = (index + step) & (mc->table_size-1);
    while (mc->val_table[index]) {
        if (key == mc->key_table[index])
            return mc->val_table[index];
        index = (index+step) & (mc->table_size-1);
    }
    return 0;
}

static inline uint16_t *mcounter_ptr(const mcounter_t *mc, const uint64_t key) 
{
    uint64_t index;
    uint64_t step;

    index = hash64(key, mc->index_len) & (mc->table_size - 1);
    index = hash64(key, mc->index_len) & (mc->table_size - 1);
    if (!(mc->val_table[index]))
        return NULL;
    else if (key == mc->key_table[index])
        return &(mc->val_table[index]);
    step = rehash64(key, mc->index_len);
    index = (index + step) & (mc->table_size-1);
    while (mc->val_table[index]) {
        if (key == mc->key_table[index])
            return &(mc->val_table[index]);
        index = (index+step) & (mc->table_size-1);
    }
    return NULL;
}

#endif
