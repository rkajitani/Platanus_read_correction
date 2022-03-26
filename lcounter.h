#ifndef LCOUNTER_H
#define LCOUNTER_H

#include "common.h"
#include "lmer.h"

typedef struct {
    FILE *kmer_fp;
    uint64_t mem_lim;
    uint64_t len_distr[MAX_READ_LEN + 1];
    uint64_t l_len;
    uint64_t bl;
    uint64_t index_len;
    uint64_t table_size;
    uint64_t *key_table;
    uint16_t *val_table;
    uint64_t max_occ;
    uint64_t *occ_distr;
} lcounter_t;

typedef struct {
    FILE *unstored_fp;
    FILE *in_kmer_fp;
    lcounter_t *lc;
    pthread_mutex_t *mutex;
} arg_lcounter_lock_and_count_t;

void lcounter_init(lcounter_t *lc, const uint64_t mem_lim);
void lcounter_destroy_table(lcounter_t *lc);
void lcounter_destroy(lcounter_t *lc);
void lcounter_save_kmer_read(lcounter_t *lc, const uint64_t l_len, FILE *read_fp);
void lcounter_save_additional_kmer_read(lcounter_t *lc, uint64_t l_len, FILE *read_fp);
void lcounter_load_kmer(lcounter_t *lc, const uint64_t min_occ);
void lcounter_load_contig(lcounter_t *lc, const uint64_t l_len, const double occ_ratio, FILE *contig_fp);
void lcounter_save_kmer_read_mt(lcounter_t *lc, const uint64_t l_len, FILE **read_fp, const uint64_t n_thread);
void lcounter_save_additional_kmer_read_mt(lcounter_t *lc, const uint64_t l_len, FILE **read_fp, const uint64_t n_thread);
void lcounter_lock_and_count1(const arg_lcounter_lock_and_count_t *arg);
void lcounter_lock_and_count2(const arg_lcounter_lock_and_count_t *arg);
uint64_t lcounter_sorted_lmer_list(const lcounter_t *lc, uint64_t **buf);
void lcounter_extract_read(const lcounter_t *lc, FILE **read_fpp);
double lcounter_estimate_error_rate(lcounter_t *lc, const uint64_t l_len, FILE *read_fp);
double lcounter_estimate_error_rate_mt(lcounter_t *lc, const uint64_t l_len, FILE **read_fp, const uint64_t n_thread);
double lcounter_calculate_error_rate(const lcounter_t *lc1, const lcounter_t *lc2, const double ave_len);
uint64_t lcounter_get_coverage_cutoff(const lcounter_t *lc, const double err_rate);
FILE *sorted_lkey_from_lmer_file(FILE *lmer_fp, const uint64_t l_len, const uint64_t min_occ);
FILE *sorted_lkey_from_contig_file(FILE *contig_fp, const uint64_t l_len);
void lcounter_count_quick(lcounter_t *lc, const uint8_t *seq, const uint64_t len, const uint64_t l_len);

static inline uint16_t lcounter_occ(const lcounter_t *lc, uint64_t *key) 
{
    uint64_t index;
    uint64_t step;

    index = lmer_hash(key, lc->l_len, lc->index_len) & (lc->table_size-1);
    if (!lc->val_table[index])
        return 0;
    else if (lmer_cmp(key, &(lc->key_table[index * lc->bl]), lc->l_len) == 0)
        return lc->val_table[index];
    step = lmer_rehash(key, lc->l_len, lc->index_len);
    index = (index + step) & (lc->table_size-1);
    while (lc->val_table[index]) {
        if (lmer_cmp(key, &(lc->key_table[index * lc->bl]), lc->l_len) == 0)
            return lc->val_table[index];
        index = (index+step) & (lc->table_size-1);
    }
    return 0;
}

static inline uint16_t *lcounter_ptr(const lcounter_t *lc, uint64_t *key) 
{
    uint64_t index;
    uint64_t step;

    index = lmer_hash(key, lc->l_len, lc->index_len) & (lc->table_size-1);
    if (!lc->val_table[index])
        return NULL;
    else if (lmer_cmp(key, &(lc->key_table[index * lc->bl]), lc->l_len) == 0)
        return &(lc->val_table[index]);
    step = lmer_rehash(key, lc->l_len, lc->index_len);
    index = (index + step) & (lc->table_size-1);
    while (lc->val_table[index]) {
        if (lmer_cmp(key, &(lc->key_table[index * lc->bl]), lc->l_len) == 0)
            return &(lc->val_table[index]);
        index = (index+step) & (lc->table_size-1);
    }
    return NULL;
}

#endif
