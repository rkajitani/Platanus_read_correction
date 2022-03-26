#ifndef BUBBLE_MAP_H
#define BUBBLE_MAP_H

#include "common.h"
#include "mapper.h"

enum bubble_mutation_type {
	BUB_MIS,
	BUB_DEL,
	BUB_INS,
};

typedef struct {
	char o[MAX_FILE_LEN + 1];
	int64_t m;
	int64_t t;
	int64_t k;
	int64_t n_c;
	FILE *c[MAX_FILE_NUM];
	int64_t n_b;
	FILE *b[MAX_FILE_NUM];
	double u;
} option_bubble_map_t;

typedef struct {
	enum bubble_mutation_type type;
	int32_t bub_id;
	int32_t con_id;
	int32_t ofst;
	int32_t len;
	uint8_t seq[0];
} bubble_mutation_t;

typedef struct {
	int8_t *pool;
	int64_t pool_size;
	int8_t **name;
	int64_t n_name;
	FILE *name_fp;
} name_list_t;

void print_bubble_map_usage(void);
void option_bubble_map_init(option_bubble_map_t *opt);
void option_bubble_map_destroy(const option_bubble_map_t *opt); 
uint8_t option_bubble_map_parse(option_bubble_map_t *opt, int argc, char **argv);
void bubble_map(const option_bubble_map_t *opt);
void print_bubble(FILE **pos_fp, const name_list_t *cname, const name_list_t *bname, const int64_t n_thread, char *out_name);
void detect_mutation(FILE **pos_fp, const name_list_t *cname, const name_list_t *bname, const int64_t min_match_len, const double bubble_th, const int64_t n_thread, char *out_name);
int64_t align_bubble(const varr8_t *bub, const varr8_t *con, varr64_t *score, varr8_t *trace);
int64_t trace_bubble_alignment(const varr8_t *bub, const varr8_t *con, const varr8_t *trace, const int32_t bub_id, const position_t *pos, FILE **mut_fpp);
int cmp_bubble_mutation_p(const void *a, const void *b);
void name_list_init(name_list_t *list);
void name_list_destroy(name_list_t *list);
void name_list_read_fasta(name_list_t *list, FILE *fp);
void name_list_load(name_list_t *list);

#endif
