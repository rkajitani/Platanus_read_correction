#ifndef LGRAPH_H
#define LGRAPH_H

#include "common.h"
#include "lmer.h"
#include "lcounter.h"

typedef struct {
    uint16_t cov;
    uint8_t out;
} ljunc_t;

typedef struct {
    uint64_t len;
    uint16_t cov;
    uint8_t out;
    uint8_t seq[0];
} lstra_t;

typedef struct {
    uint64_t mem_lim;
    uint64_t l_len;
    uint64_t bl;
    uint64_t n_junc;
    uint64_t n_stra;
    uint64_t jindex_len;
    uint64_t sindex_len;
    uint64_t jtable_size;
    uint64_t stable_size;
    uint64_t *jkey_table;
    ljunc_t *jval_table;
    uint64_t *lskey_table;
    uint64_t *rskey_table;
    lstra_t **lsval_table;
    lstra_t **rsval_table;
    double bubble_th;
    double branch_th;
    FILE *junc_file;
    FILE *stra_file;
} lgraph_t;

extern lstra_t g_lstra_deleted;

void lgraph_init(lgraph_t *lg, const double branch_th, const double bubble_th);
void lgraph_destroy_table(lgraph_t *lg);
void lgraph_destroy(lgraph_t *lg);
void lgraph_save(lgraph_t *lg, const lcounter_t *lc, FILE *sorted_key_fp);
void lgraph_load(lgraph_t *lg);
void lgraph_extract_read(const lgraph_t *lg, FILE **read_fpp);
void lgraph_show_contig(const lgraph_t *lg, const double cov_ratio, const char *out_name);
void lgraph_delete_clear_junc(const lgraph_t *lg, const uint64_t len);
void lgraph_lstra_insert(lgraph_t *lg, lstra_t *sp);
void lgraph_lstra_delete(lgraph_t *lg, lstra_t *sp);
void lgraph_join_nodes(lgraph_t *lg);
void lgraph_simplify(lgraph_t *lg);
uint64_t lgraph_crush_bubble(lgraph_t *lg, const double ave_cov, FILE **bubble_fpp);
uint64_t lgraph_cut_branch(lgraph_t *lg);
FILE *lgraph_save_contig_simple(const lgraph_t *lg);
FILE *lgraph_save_contig(const lgraph_t *lg, const uint64_t next_l);
FILE *lgraph_merge(lgraph_t *lg, lgraph_t *ext_lg);
lstra_t *lgraph_join_ss(lgraph_t *lg, lstra_t *lsp, lstra_t *rsp);
lstra_t *lgraph_join_sj(lgraph_t *lg, lstra_t *sp, ljunc_t *jp);
lstra_t *lgraph_join_js(lgraph_t *lg, ljunc_t *jp, lstra_t *sp);
lstra_t *lgraph_join_jj(lgraph_t *lg, ljunc_t *ljp, ljunc_t *rjp);
lstra_t *lgraph_change_j2s(lgraph_t *lg, ljunc_t *jp);
void lgraph_cut_branch_iterative(lgraph_t *lg);
FILE *lgraph_crush_bubble_iterative(lgraph_t *lg, const double ave_cov);
uint64_t lgraph_sorted_prefix_list(const lgraph_t *lg, uint64_t **buf);
uint64_t lgraph_sorted_junction_list(const lgraph_t *lg, uint64_t **buf);
void lgraph_save_edge_kmer(const lgraph_t *lg, lcounter_t *lc, const uint64_t next_l);
void lgraph_print_dot(const lgraph_t *lg, const double cov_ratio, const char *out_name);
void lgraph_make_graph_quick(lgraph_t *lg, const lcounter_t *lc, const uint64_t min_cov);
void lgraph_straight_stats(const lgraph_t *lg, uint64_t *len_cutoff, uint64_t *cov_cutoff, double *ave_cov);
void lgraph_delete_erroneous_straight_node(lgraph_t *lg, const uint64_t len_cutoff, const uint64_t cov_cutoff);

static inline bool ljunc_is_del(const ljunc_t lj) 
{
    return lj.cov == UINT16_MAX;
}

static inline bool ljunc_is_emp(const ljunc_t lj) 
{
    return lj.cov == 0;
}

static inline ljunc_t *lgraph_ljunc_find(const lgraph_t *lg, uint64_t *key)
{
    uint64_t index;
    uint64_t step;

    index = lmer_hash(key, lg->l_len, lg->jindex_len) & (lg->jtable_size-1);
    if (ljunc_is_emp(lg->jval_table[index]))
        return NULL;
    else if (!ljunc_is_del(lg->jval_table[index]) && lmer_cmp(key, &(lg->jkey_table[index * lg->bl]), lg->l_len) == 0) {
        return &(lg->jval_table[index]);
    }
    step = lmer_rehash(key, lg->l_len, lg->jindex_len);
    index = (index + step) & (lg->jtable_size-1);
    while (!ljunc_is_emp(lg->jval_table[index])) {
        if (!ljunc_is_del(lg->jval_table[index]) && lmer_cmp(key, &(lg->jkey_table[index * lg->bl]), lg->l_len) == 0) {
            return &(lg->jval_table[index]);
        }
        index = (index+step) & (lg->jtable_size-1);
    }
    return NULL;
}

static inline void ljunc_delete(ljunc_t *lj_p)
{
    lj_p->cov = UINT16_MAX;
}

static inline void lgraph_ljunc_insert(lgraph_t *lg, const ljunc_t *jp, uint64_t *key)
{
    uint64_t index;
    uint64_t step;

    index = lmer_hash(key, lg->l_len, lg->jindex_len) & (lg->jtable_size-1);
    if (ljunc_is_emp(lg->jval_table[index]) || ljunc_is_del(lg->jval_table[index])) {
        lg->jval_table[index] = *jp;
        lmer_cpy(&(lg->jkey_table[index * lg->bl]), key, lg->l_len);
    }
    else {
        step = lmer_rehash(key, lg->l_len, lg->jindex_len);
        index = (index+step) & (lg->jtable_size-1);
        while (!(ljunc_is_emp(lg->jval_table[index]) || ljunc_is_del(lg->jval_table[index])))
            index = (index+step) & (lg->jtable_size-1);
        lg->jval_table[index] = *jp;
        lmer_cpy(&(lg->jkey_table[index * lg->bl]), key, lg->l_len);
    }
}

static inline bool lstra_is_del(const lstra_t *ls_p)
{
    return ls_p == &g_lstra_deleted;
}

static inline bool lstra_is_emp(const lstra_t *ls_p)
{
    return ls_p == NULL;
}

static inline lstra_t *lgraph_lstra_find_pre(const lgraph_t *lg, uint64_t *key)
{
    uint64_t index;
    uint64_t step;

    index = lmer_hash(key, lg->l_len, lg->sindex_len) & (lg->stable_size-1);
    if (lstra_is_emp(lg->lsval_table[index]))
        return NULL;
    else if (!lstra_is_del(lg->lsval_table[index]) && lmer_cmp(key, &(lg->lskey_table[index * lg->bl]), lg->l_len) == 0)
        return lg->lsval_table[index];
    step = lmer_rehash(key, lg->l_len, lg->sindex_len);
    index = (index + step) & (lg->stable_size-1);
    while (!lstra_is_emp(lg->lsval_table[index])) {
        if (!lstra_is_del(lg->lsval_table[index]) && lmer_cmp(key, &(lg->lskey_table[index * lg->bl]), lg->l_len) == 0)
            return lg->lsval_table[index];
        index = (index+step) & (lg->stable_size-1);
    }
    return NULL;
}

static inline lstra_t *lgraph_lstra_find_suf(const lgraph_t *lg, uint64_t *key)
{
    uint64_t index;
    uint64_t step;

    index = lmer_hash(key, lg->l_len, lg->sindex_len) & (lg->stable_size-1);
    if (lstra_is_emp(lg->rsval_table[index]))
        return NULL;
    else if (!lstra_is_del(lg->rsval_table[index]) && lmer_cmp(key, &(lg->rskey_table[index * lg->bl]), lg->l_len) == 0)
        return lg->rsval_table[index];
    step = lmer_rehash(key, lg->l_len, lg->sindex_len);
    index = (index + step) & (lg->stable_size-1);
    while (!lstra_is_emp(lg->rsval_table[index])) {
        if (!lstra_is_del(lg->rsval_table[index]) && lmer_cmp(key, &(lg->rskey_table[index * lg->bl]), lg->l_len) == 0)
            return lg->rsval_table[index];
        index = (index+step) & (lg->stable_size-1);
    }
    return NULL;
}

#endif
