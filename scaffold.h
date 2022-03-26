#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "common.h"
#include "mapper.h"

#define SCA_HASH_OVERLAP 32
#define MIN_SCA_LEN 100
#define SC_EDGE_NUM_RATE_TH 0.5
#define SC_EDGE_EXPECTED_RATE_TH 0.5
#define MIN_TOL_FACTOR 2
#define MAX_TOL_FACTOR 3
#define SC_REP 0x1 
#define SC_INC 0x2
#define SC_DEL 0x4 

typedef struct {
    char o[MAX_FILE_LEN + 1];
    int64_t m;
    int64_t t;
    int64_t s;
    int64_t v;
    int64_t l;
    double u;
    int64_t n_c;
    FILE *c[MAX_FILE_NUM];
    int64_t n_b;
    FILE *b[MAX_FILE_NUM];
    int64_t n_pair;
    seqlib_input_t pair[MAX_FILE_NUM + 1];
    int64_t n[MAX_FILE_NUM];
    int64_t a[MAX_FILE_NUM];
    int64_t d[MAX_FILE_NUM];
} option_scaffold_t;

typedef struct {
    int64_t id1;
    int64_t id2;
    int64_t ofst1;
    int64_t ofst2;
    int64_t gap;
} link_t;

typedef struct {
    int32_t id1;
    int32_t id2;
    int32_t len;
} overlap_rec_t;

typedef struct {
    int8_t dir;
    int64_t ed;
    int64_t len;
    int64_t n_link;
    int64_t *bkdn;
} sc_edge_t;

typedef struct {
    int64_t id;
    int64_t st;
    int64_t ed;
} scaffold_part_t;

typedef struct {
    int8_t state;
    int64_t len;
    int64_t n_edge;
    sc_edge_t *edge;
    int64_t n_con;
    scaffold_part_t *con;
} sc_node_t;

typedef struct {
    int64_t id;
    int64_t st;
    int64_t ed;
    int64_t dist;
    int64_t n_link;
} sc_layout_part_t;

typedef struct {
    int64_t size;
    int64_t capa;
    sc_layout_part_t *part;
} sc_layout_t;

typedef struct {
    int64_t mem_lim;
    int64_t mem_usage;

    int64_t seed_len;
    int64_t min_overlap;
    int64_t hash_overlap;
    int64_t index_len;

    FILE *con_fp;
    lseq_t *con;
    int64_t n_con;
    position_t *con_inc;
    int64_t con_pool_size;
    scaffold_part_t *con_pool;
    int64_t seq_pool_size;
    int8_t *seq_pool;
    uint16_t *cov;
	double ave_cov;

    seqlib_t *lib;
    int64_t min_link;
    FILE *link_fp;
    int64_t tol;
    sc_node_t *node;
    int64_t n_node;
    int64_t edge_pool_size;
    sc_edge_t *edge_pool;
    int64_t bkdn_pool_size;
    int64_t *bkdn_pool;

    FILE *ol_fp;
    int64_t ol_table_size;
    int64_t ol_index_len;
    overlap_rec_t *ol_table;
} scaffolder_t;

void print_scaffold_usage();
void scaffold_mt(option_scaffold_t *opt);
void option_scaffold_init(option_scaffold_t *opt);
void option_scaffold_destroy(const option_scaffold_t *opt);
uint8_t option_scaffold_parse(option_scaffold_t *opt, int argc, char **argv);
int64_t contig_median_len(const contig_t *con);
int64_t contig_n50_len(const contig_t *con);
double contig_ave_cov(const contig_t *con, int64_t min_len);
void scaffolder_init(scaffolder_t *sc, const int64_t s, const int64_t o, const int64_t mem_lim);
void scaffolder_destroy(scaffolder_t *sc);
void scaffolder_destroy_graph(scaffolder_t *sc);
void scaffolder_save_overlap(scaffolder_t *sc, const mapper_t *mp, const int64_t hash_overlap, const int64_t len_cutoff);
void scaffolder_load_overlap(scaffolder_t *sc);
void scaffolder_insert_overlap(scaffolder_t *sc, overlap_rec_t ol);
int64_t scaffolder_get_overlap(scaffolder_t *sc, int64_t id1, int64_t id2);
int64_t scaffolder_get_short_overlap(scaffolder_t *sc, int64_t id1, int64_t id2);
int64_t scaffolder_get_sca_overlap(scaffolder_t *sc, int64_t id1, int64_t id2);
int64_t scaffolder_estimate_link(scaffolder_t *sc);
void scaffolder_init_sca(scaffolder_t *sc, contig_t *con, const int64_t seed_len, const double ave_cov);
void scaffolder_set_seqlib(scaffolder_t *sc, seqlib_t *lib);
void scaffolder_set_tol(scaffolder_t *sc, int64_t tol);
void scaffolder_link(scaffolder_t *sc, int64_t min_link);
void scaffolder_make_graph(scaffolder_t *sc);
int64_t scaffolder_delete_erroneous_edge(scaffolder_t *sc);
int64_t scaffolder_crush_bubble(scaffolder_t *sc, double bubble_th, FILE **bubble_fpp);
void scaffolder_layout_nodes(scaffolder_t *sc, sc_node_t *new_node, sc_layout_t *ret, sc_layout_t *work);
void scaffolder_crush_bubble_iterative(scaffolder_t *sc, double bubble_th, FILE **bubble_fpp);
void scaffolder_delete_erroneous_edge_iterative(scaffolder_t *sc);
void scaffolder_delete_repeat_edge(scaffolder_t *sc);
void scaffolder_detect_repeat(scaffolder_t *sc);
void scaffolder_delete_edges(scaffolder_t *sc, varr64_t *ids);
void scaffolder_split(scaffolder_t *sc);
void scaffolder_scaffold(scaffolder_t *sc);
void scaffolder_remake(scaffolder_t *sc, int64_t new_n_node, int64_t new_con_pool_size, FILE *sca_fp);
void scaffolder_print_layout(scaffolder_t *sc, char *out_name);
void scaffolder_print_seq(scaffolder_t *sc, char *out_name);
void scaffolder_cut_and_print_seq(scaffolder_t *sc, int64_t min_len, char *out_name);
double scaffolder_node_cov(scaffolder_t *sc, sc_node_t *node);
double expected_link(seqlib_t *lib, double l1, double l2, double g);
/* int cmp_position(const void *a, const void *b); */
int cmp_link(const void *a, const void *b);
int cmp_sc_layout_part(const void *a, const void *b);
void sc_layout_init(sc_layout_t *lo);
void sc_layout_destroy(sc_layout_t *lo);
void sc_layout_resize(sc_layout_t *lo, int64_t size);
void sc_layout_delete(sc_layout_t *lo, int64_t i);
void scaffolder_layout2seq(scaffolder_t *sc, sc_layout_part_t *lo, int64_t lo_size, varr8_t *ret);
int64_t align_scaffold(varr8_t *sca1, varr8_t *sca2, varr64_t *work);
double scaffolder_layout_ave_cov(scaffolder_t *sc, sc_layout_part_t *lo, int64_t lo_size);
void print_scaffold_bubble(FILE *bubble_fp, char *out_name);
int64_t scaffolder_delete_hetero_edge(scaffolder_t *sc);
int64_t scaffolder_crush_hetero_bubble(scaffolder_t *sc, FILE **bubble_fpp);
double scaffolder_ave_len(scaffolder_t *sc);
int64_t scaffolder_get_similer_sca_overlap(scaffolder_t *sc, int64_t id1, int64_t id2);
int64_t scaffolder_get_similer_overlap(const scaffolder_t *sc, int64_t id1, int64_t id2);

static inline bool overlap_rec_is_emp(const overlap_rec_t ol)
{
    return ol.id1 == 0ull;
}

#endif
