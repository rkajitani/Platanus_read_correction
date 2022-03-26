#ifndef MGRAPH_H
#define MGRAPH_H

#include "common.h"
#include "mcounter.h"
#include "lcounter.h"

typedef struct {
    uint64_t key;
    uint16_t cov;
    uint8_t out;
} mjunc_t;

typedef struct {
    uint64_t lkey;
    uint64_t rkey;
    uint64_t len;
    uint16_t cov;
    uint8_t out;
    uint8_t seq[0];
} mstra_t;

typedef struct {
    uint64_t mem_lim;
    uint64_t k_len;
    uint64_t n_junc;
    uint64_t n_stra;
    uint64_t jindex_len;
    uint64_t sindex_len;
    uint64_t jtable_size;
    uint64_t stable_size;
    double bubble_th;
    double branch_th;
    mjunc_t *jtable;
    mstra_t **lstable;
    mstra_t **rstable;
    FILE *junc_file;
    FILE *stra_file;
} mgraph_t;

extern mstra_t g_mstra_deleted;

void mgraph_init(mgraph_t *mg, const double branch_th, const double bubble_th);
void mgraph_destroy_table(mgraph_t *mg);
void mgraph_destroy(mgraph_t *mg);
void mgraph_save(mgraph_t *mg, const mcounter_t *mc, FILE *sorted_key_fp);
void mgraph_load(mgraph_t *mg);
void mgraph_extract_read(const mgraph_t *mg, FILE **read_fpp);
void mgraph_show_contig(const mgraph_t *mg, const double cov_ratio, const char *out_name);
void mgraph_join_nodes(mgraph_t *mg);
void mgraph_simplify(mgraph_t *mg);
uint64_t mgraph_crush_bubble(mgraph_t *mg, const double ave_cov, FILE **bubble_fpp);
FILE *mgraph_save_contig(const mgraph_t *mg, const uint64_t next_k);
FILE *mgraph_merge(mgraph_t *mg, mgraph_t *ext_mg);
mstra_t *mgraph_join_ss(mgraph_t *mg, mstra_t *lsp, mstra_t *rsp);
mstra_t *mgraph_join_sj(mgraph_t *mg, mstra_t *sp, mjunc_t *jp);
mstra_t *mgraph_join_js(mgraph_t *mg, mjunc_t *jp, mstra_t *sp);
mstra_t *mgraph_join_jj(mgraph_t *mg, mjunc_t *ljp, mjunc_t *rjp);
mstra_t *mgraph_change_j2s(mgraph_t *mg, mjunc_t *jp);
uint64_t mgraph_cut_branch(mgraph_t *mg);
void mgraph_cut_branch_iterative(mgraph_t *mg);
FILE *mgraph_crush_bubble_iterative(mgraph_t *mg, const double ave_cov);
int cmp_kmer_uniq(const void *a, const void *b);
uint64_t mgraph_sorted_prefix_list(const mgraph_t *mg, uint64_t **buf);
uint64_t mgraph_sorted_junction_list(const mgraph_t *mg, uint64_t **buf);
void mgraph_save_edge_kmer(const mgraph_t *mg, mcounter_t *mc, const uint64_t next_k);
void mgraph_make_graph_quick(mgraph_t *mg, const mcounter_t *mc, const uint64_t min_cov);
void mgraph_straight_stats(const mgraph_t *mg, uint64_t *len_cutoff, uint64_t *cov_cutoff, double *ave_cov);
void mgraph_delete_erroneous_straight_node(mgraph_t *mg, const uint64_t len_cutoff, const uint64_t cov_cutoff);

static inline bool mjunc_is_del(const mjunc_t mj) 
{
    return mj.cov == UINT16_MAX;
}

static inline bool mjunc_is_emp(const mjunc_t mj) 
{
    return mj.cov == 0;
}

static inline mjunc_t *mgraph_mjunc_find(const mgraph_t *mg, const uint64_t key)
{
    uint64_t index;
    uint64_t step;

    if (mjunc_is_emp(mg->jtable[(index = hash64(key, mg->jindex_len)&(mg->jtable_size-1))]))
        return NULL;
    else if (!mjunc_is_del(mg->jtable[index]) && key == mg->jtable[index].key)
        return &(mg->jtable[index]);
    index = (index + (step = rehash64(key, mg->jindex_len))) & (mg->jtable_size-1);
    while (!mjunc_is_emp(mg->jtable[index])) {
        if (!mjunc_is_del(mg->jtable[index]) && key == mg->jtable[index].key)
            return &(mg->jtable[index]);
        index = (index+step) & (mg->jtable_size-1);
    }
    return NULL;
}

static inline void mjunc_delete(mjunc_t *mj_p)
{
    mj_p->cov = UINT16_MAX;
}

static inline void mgraph_mjunc_insert(mgraph_t *mg, const mjunc_t *jp)
{
    uint64_t index;
    uint64_t step;

    index = hash64(jp->key, mg->jindex_len) & (mg->jtable_size-1);
    if (mjunc_is_emp(mg->jtable[index]) || mjunc_is_del(mg->jtable[index]))
        mg->jtable[index] = *jp;
    else {
        step = rehash64(jp->key, mg->jindex_len);
        index = (index+step) & (mg->jtable_size-1);
        while (!(mjunc_is_emp(mg->jtable[index]) || mjunc_is_del(mg->jtable[index])))
            index = (index+step) & (mg->jtable_size-1);
        mg->jtable[index] = *jp;
    }
}

static inline bool mstra_is_del(const mstra_t *ms_p)
{
    return ms_p == &g_mstra_deleted;
}

static inline bool mstra_is_emp(const mstra_t *ms_p)
{
    return ms_p == NULL;
}

static inline mstra_t *mgraph_mstra_find_pre(const mgraph_t *mg, const uint64_t key)
{
    uint64_t index;
    uint64_t step;

    if (mstra_is_emp(mg->lstable[(index = hash64(key, mg->sindex_len)&(mg->stable_size-1))]))
        return NULL;
    else if (!mstra_is_del(mg->lstable[index]) && key == mg->lstable[index]->lkey)
        return mg->lstable[index];
    index = (index + (step = rehash64(key, mg->sindex_len))) & (mg->stable_size-1);
    while (!mstra_is_emp(mg->lstable[index])) {
        if (!mstra_is_del(mg->lstable[index]) && key == mg->lstable[index]->lkey)
            return mg->lstable[index];
        index = (index+step) & (mg->stable_size-1);
    }
    return NULL;
}

static inline mstra_t *mgraph_mstra_find_suf(const mgraph_t *mg, const uint64_t key)
{
    uint64_t index;
    uint64_t step;

    if (mstra_is_emp(mg->rstable[(index = hash64(key, mg->sindex_len)&(mg->stable_size-1))]))
        return NULL;
    else if (!mstra_is_del(mg->rstable[index]) && key == mg->rstable[index]->rkey)
        return mg->rstable[index];
    index = (index + (step = rehash64(key, mg->sindex_len))) & (mg->stable_size-1);
    while (!mstra_is_emp(mg->rstable[index])) {
        if (!mstra_is_del(mg->rstable[index]) && key == mg->rstable[index]->rkey)
            return mg->rstable[index];
        index = (index+step) & (mg->stable_size-1);
    }
    return NULL;
}

static inline void mgraph_mstra_delete(mgraph_t *mg, mstra_t *sp)
{
    uint64_t index;
    uint64_t step;

    index = hash64(sp->rkey, mg->sindex_len) & (mg->stable_size-1);
    if (!mstra_is_del(mg->rstable[index]) && sp->rkey == mg->rstable[index]->rkey)
        mg->rstable[index] = &g_mstra_deleted;
    else {
        index = (index + (step = rehash64(sp->rkey, mg->sindex_len))) & (mg->stable_size-1);
        while (mstra_is_del(mg->rstable[index]) || sp->rkey != mg->rstable[index]->rkey)
            index = (index+step) & (mg->stable_size-1);
        mg->rstable[index] = &g_mstra_deleted;
    }

    index = hash64(sp->lkey, mg->sindex_len) & (mg->stable_size-1);
    if (!mstra_is_del(mg->lstable[index]) && sp->lkey == mg->lstable[index]->lkey)
        mg->lstable[index] = &g_mstra_deleted;
    else {
        index = (index + (step = rehash64(sp->lkey, mg->sindex_len))) & (mg->stable_size-1);
        while (mstra_is_del(mg->lstable[index]) || sp->lkey != mg->lstable[index]->lkey)
            index = (index+step) & (mg->stable_size-1);
        mg->lstable[index] = &g_mstra_deleted;
    }

    my_free(sp);
}

static inline void mgraph_mstra_insert(mgraph_t *mg, mstra_t *sp)
{
    uint64_t index;
    uint64_t step;

    index = hash64(sp->lkey, mg->sindex_len) & (mg->stable_size-1);
    if (mstra_is_emp(mg->lstable[index]) || mstra_is_del(mg->lstable[index]))
        mg->lstable[index] = sp;
    else {
        step = rehash64(sp->lkey, mg->sindex_len);
        index = (index+step) & (mg->stable_size-1);
        while (!(mstra_is_emp(mg->lstable[index]) || mstra_is_del(mg->lstable[index])))
            index = (index+step) & (mg->stable_size-1);
        mg->lstable[index] = sp;
    }

    index = hash64(sp->rkey, mg->sindex_len) & (mg->stable_size-1);
    if (mstra_is_emp(mg->rstable[index]) || mstra_is_del(mg->rstable[index]))
        mg->rstable[index] = sp;
    else {
        step = rehash64(sp->rkey, mg->sindex_len);
        index = (index+step) & (mg->stable_size-1);
        while (!(mstra_is_emp(mg->rstable[index]) || mstra_is_del(mg->rstable[index])))
            index = (index+step) & (mg->stable_size-1);
        mg->rstable[index] = sp;
    }
}

#endif
