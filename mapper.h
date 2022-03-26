#ifndef MAPPER_H
#define MAPPER_H

#include "common.h"
#include "seqlib.h"

typedef struct {
    uint64_t key;
    union {
        position_t pos;
        int64_t val;
    };
} map_rec_t;

typedef struct {
    int64_t mem_lim;
    int64_t mem_usage;

    int64_t key_len;
    uint64_t mask;
    int64_t seed_len;
    int64_t index_len;

    map_rec_t *table;
    uint64_t table_size;
    int64_t pos_pool_size;
    position_t *pos_pool;
    int64_t max_occ;

    FILE *seq_fp;
    lseq_t *seq;
    int64_t n_seq;
    int64_t seq_pool_size;
    int8_t *seq_pool;
} mapper_t;

typedef struct {
    int64_t key_len;
    int64_t seed_len;
    position_t *bubble_map;
    mapper_t cmp;
    mapper_t bmp;
} hetero_mapper_t;


static inline bool map_rec_is_emp(const map_rec_t mr)
{
    return mr.val == 0;
}

void mapper_init(mapper_t *mp, const int64_t seed_len, const int64_t key_len, const int64_t mem_lim);
void mapper_destroy_table(mapper_t *mp);
void mapper_destroy(mapper_t *sc);
void mapper_set_contig(mapper_t *mp, const contig_t *con);
void mapper_make_table(mapper_t *sc);
int64_t mapper_map_seed(const mapper_t *sc, const kmer_t *kmer, position_t *buf);
position_t mapper_map_read(const mapper_t *sc, const seq_t *read, position_t *buf);
void mapper_map_pair_mt(const mapper_t *mp, seqlib_t *lib, const int64_t min_ins, const int64_t n_thread);
void mapper_map_small_gap_mt(const mapper_t *mp, seqlib_t *lib, FILE **gap_seq_fp, const int64_t n_thread);
void mapper_map_pair_nolink_mt(const mapper_t *mp, seqlib_t *lib, const int64_t n_thread);
void hetero_mapper_init(hetero_mapper_t *mp, const int64_t seed_len, const int64_t key_len, const int64_t mem_lim);
void hetero_mapper_destroy(hetero_mapper_t *mp);
void hetero_mapper_merge_mt(hetero_mapper_t *mp, const int64_t n_thread);
position_t hetero_mapper_map_read(const hetero_mapper_t *mp, const seq_t *read, position_t *buf);
void hetero_mapper_map_pair_mt(const hetero_mapper_t *mp, seqlib_t *lib, const int64_t min_ins, const int64_t n_thread);
void mapper_map_bubble(const mapper_t *mp, const contig_t *bubble, const int64_t min_match_len, FILE **pos_fp, const int64_t n_thread);

#endif
