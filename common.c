#include "common.h"
#include <omp.h>

uint64_t g_l_len;

void *my_malloc(size_t size)
{
    void *ptr = (void *)malloc(size);
    check_alloc(ptr);
    return ptr;
}

void *my_calloc(size_t n, size_t size)
{
    void *ptr = (void *)calloc(n, size);
    check_alloc(ptr);
    return ptr;
}

void *my_realloc(void *ptr, size_t size)
{
    void *new_ptr = (void *)realloc(ptr, size);
    check_alloc(new_ptr);
    return new_ptr;
}

void my_free(void *ptr)
{
    free(ptr);
}

void my_exit(const uint16_t i)
{
    exit(i);
}

void check_alloc(const void *p)
{
    if (p == NULL) {
        fputs("error: cannot allocate memory!\n", stderr);
        my_exit(1);
    }
}

void check_file_open(const FILE *fp, const char *fname)
{
    if (fp == NULL) {
        fprintf(stderr, "error: cannot open file %s\n",  fname);
        my_exit(1);
    }
}

void check_tmpfile(const FILE *fp)
{
    if (fp == NULL) {
        fputs("error: cannot create tmpfile_open!\n", stderr);
        my_exit(1);
    }
}

void check_mem_usage(const uint64_t usage, const uint64_t limit)
{
    if (usage > limit)
        fprintf(stderr, "WARNING: memory usage(%luB) exceeds the limit(%luB), sorry!\n", usage, limit);
}

void put_command_log(int argc, char **argv, FILE *stream)
{
    int i;

    fprintf(stream, "%s", argv[0]);
    for (i = 1; i < argc; ++i)
        fprintf(stream, " %s", argv[i]);
    fputs("\n\n", stream);
}

FILE *tmpfile_open(void)
{
    int fd;
    char template[] = TMPFILE_DIR "/XXXXXX";
    FILE *fp;

    fd = mkstemp(template);
    if (fd == -1) {
        fputs("error: cannot create tmpfile_open!\n", stderr);
        my_exit(1);
    }

    fp = fdopen(fd, "wb+");
    unlink(template);
    return fp;
}

void rev_comp(uint8_t *seq, const int64_t len)
{
    int64_t i;
    int8_t tmp;

    for (i = 0; i < len / 2; ++i) {
        tmp = seq[i];
        seq[i] = seq[len - i - 1] != 4 ? 3^seq[len - i - 1] : 4;
        seq[len - i - 1] = tmp != 4 ? 3^tmp : 4;
    }
    if (len % 2 == 1)
        seq[i] = seq[i] != 4 ? 3^seq[i] : 4;
}

int suffix_atoi(char *str)
{
    for (; *str != '\0' && !isdigit(*str); ++str)
        ;
    return atoi(str);
}

uint64_t option_int(int argc, char **argv, uint64_t i, const int64_t min, const int64_t max, int64_t *val)
{
    char *name;

    name = argv[i];

    ++i;
    if (i < argc)
        *val = atoi(argv[i]);
    if (i >= argc || *val < min || *val > max) {
        fprintf(stderr, "error: %s must specify integer[", name);
        if (min != INT64_MIN)
            fprintf(stderr, "%ld", min);
        fputs(", ", stderr);
        if (min != INT64_MAX)
            fprintf(stderr, "%ld", max);
        fputs("]\n\n", stderr);
        my_exit(1);
    }
    ++i;
    return i;
}

uint64_t option_multi_int(int argc, char **argv, uint64_t i, const int64_t min, const int64_t max, int64_t *val, int64_t *n_val)
{
    char *name;

    name = argv[i];

    for (++i; i < argc && argv[i][0] != '-'; ++i) {
        if (i < argc) {
            val[*n_val] = atoi(argv[i]);
            ++(*n_val);
        }
        if (atoi(argv[i]) > max || atoi(argv[i]) < min) {
            fprintf(stderr, "error: %s must specify integer[", name);
            if (min != INT64_MIN)
                fprintf(stderr, "%ld", min);
            fputs(", ", stderr);
            if (min != INT64_MAX)
                fprintf(stderr, "%ld", max);
            fputs("]\n\n", stderr);
            my_exit(1);
        }
    }
    return i;
}

uint64_t option_string(int argc, char **argv, uint64_t i, const int64_t max, char *val)
{
    char *name;

    name = argv[i];

    ++i;
    if (i < argc)
        strcpy(val, argv[i]);
    if (i >= argc || strlen(argv[i]) > max) {
        fprintf(stderr, "error: %s must specify string(length <= %ld)\n", name, max);    
        my_exit(1);
    }
    ++i;
    return i;
}

uint64_t option_float(int argc, char **argv, uint64_t i, const double min, const double max, double *val)
{
    char *name;

    name = argv[i];

    ++i;
    if (i < argc)
        *val = atof(argv[i]);
    if (i >= argc || *val < min || *val > max) {
        fprintf(stderr, "error: %s must specify integer[", name);
        if (min != -DBL_MAX)
            fprintf(stderr, "%f", min);
        fputs(", ", stderr);
        if (min != DBL_MAX)
            fprintf(stderr, "%f", max);
        fputs("]\n\n", stderr);
        my_exit(1);
    }
    ++i;
    return i;
}

uint64_t option_multi_file(int argc, char **argv, uint64_t i, single_seq_file_t *val, int64_t *n_val)
{
    ++i;
    if (i >= argc) {
        fprintf(stderr, "error: %s must specify file(number <= %u, name_length <= %u)\n\n", argv[i], MAX_FILE_NUM, MAX_FILE_LEN);
        my_exit(1);
    }
    for (; i < argc && argv[i][0] != '-'; ++i) {
        val[*n_val].fmt = check_seq_file_format(argv[i]);
        val[*n_val].fp = fopen(argv[i], "r");
        ++(*n_val);
    }
    return i;
}

uint64_t option_multi_file_name(int argc, char **argv, uint64_t i, FILE **fp_val, char **name_val, int64_t *n_val)
{
    char *name;

    name = argv[i];

    ++i;
    if (i >= argc)
        fprintf(stderr, "error: %s must specify file(number <= %u, name_length <= %u)\n\n", name, MAX_FILE_NUM, MAX_FILE_LEN);
    for (; i < argc && argv[i][0] != '-'; ++i) {
        if((fp_val[*n_val] = fopen(argv[i], "r")) == NULL) {
            fprintf(stderr, "error: cannot open \"%s\"\n", argv[i]);
            my_exit(1);
        }
        name_val[*n_val] = (char *)my_malloc((strlen(argv[i]) + 1) * sizeof(char));
        check_alloc(name_val[*n_val]);
        strcpy(name_val[*n_val], argv[i]);
        ++(*n_val);
    }
    return i;
}


uint64_t option_multi_fasta_file(int argc, char **argv, uint64_t i, FILE **val, int64_t *n_val)
{
    ++i;
    if (i >= argc) {
        fprintf(stderr, "error: %s must specify file(number <= %u, name_length <= %u)\n\n", argv[i], MAX_FILE_NUM, MAX_FILE_LEN);
        my_exit(1);
    }
    for (; i < argc && argv[i][0] != '-'; ++i) {
        if (check_seq_file_format(argv[i]) != FASTA) {
            fprintf(stderr, "error: invalid file format %s\nfasta format expected\n", argv[i]);
            my_exit(1);
        }
        val[*n_val] = fopen(argv[i], "r");
        ++(*n_val);
    }
    return i;
}

seq_file_format_t check_seq_file_format(char *name)
{
    char c;
    FILE *fp;

    fp = fopen(name, "r");
    if (fp == NULL) {
        fprintf(stderr, "error: cannot open file %s\n", name);
        my_exit(1);
    }

    c = getc(fp);
    while (c != EOF) {
        if (c == '>')
            return FASTA;
        else if (c == '@')
            return FASTQ;
        c = fskip('\n', fp);
        c = getc(fp);
    }
    fclose(fp);

    fprintf(stderr, "error: invalid file format %s\nfasta or fastq format expected\n", name);
    my_exit(1);
    return FASTA;
}

void read_fasta(FILE *in, FILE *out)
{
    int c;
    seq_t seq = {0};

    fseek(out, 0, SEEK_END);
    rewind(in);
    c = fskip('>', in);
    while (c != EOF) {
        fskip('\n', in);
        c = fget_seq(&seq, '>', in);

        seq_write(&seq, out);
    }
}

void read_fastq(FILE *in, FILE *out)
{
    int c;
    seq_t seq = {0};

    fseek(out, 0, SEEK_END);
    rewind(in);
    c = fskip('@', in);
    while (c != EOF) {
        fskip('\n', in);
        fget_seq(&seq, '+', in);
        fskip('\n', in);
        c = fastq_skip_qual(seq.len, in);

        seq_write(&seq, out);
    }
}

double get_ave(const uint64_t *distr, const uint64_t min, const uint64_t max)
{
    uint64_t i;
    uint64_t sum;
    uint64_t num;

    sum = num = 0;
    for (i = min; i <= max; ++i) {
        sum += i * distr[i];
        num += distr[i];
    }
    if (num > 0)
        return (double)sum / num;
    else
        return 0.0;
}

uint64_t get_left_minimal(const uint64_t *distr, const uint64_t max)
{
    uint64_t i;
    uint64_t pre;

    if (max <= 0)
        return 0;

    for (i = 0; i <= max; ++i) {
        if (distr[i] > 0)
            break;
    }

    pre = distr[i];

    for (++i; i <= max; ++i) {
        if (distr[i] >= pre)
            break;
        pre = distr[i];
    }
    if (i <= max)
        return i - 1;
    else
        return 1;
}

uint64_t get_left_minimal_smooth(const uint64_t *distr, const uint64_t max, const uint64_t wsize)
{
    uint64_t i;
    uint64_t pre;
    uint64_t *window;

    if (max <= wsize)
        return 0;
    
    window = (uint64_t *)my_calloc(max - wsize + 2, sizeof(uint64_t));
    for (i = 0; i < wsize; ++i)
        window[1] += distr[1 + i];

    pre = window[1];
    for (i = 2; i <= (max - wsize + 1); ++i) {
        window[i] = window[i - 1] - distr[i - 1]  + distr[i + wsize - 1];
        if (window[i] >= pre)
            break;
        pre = window[i];
    }

    my_free(window);

    if (i <= max)
        return (i - 1 + wsize / 2);
    else
        return (1 + wsize / 2);
}

void read_fasta_mt(FILE *in, FILE **out, const uint64_t n_thread)
{
    uint64_t i;
    int c;
    seq_t seq = {0};

    for (i = 0; i < n_thread; ++i) 
        fseek(out[i], 0, SEEK_END);

    rewind(in);
    c = fskip('>', in);
    i = 0;
    while (c != EOF) {
        fskip('\n', in);
        c = fget_seq(&seq, '>', in);

        seq_write(&seq, out[i]);
        i = (i + 1) % n_thread;
    }
}

void read_fastq_mt(FILE *in, FILE **out, uint64_t const n_thread)
{
    uint64_t i;
    int8_t c;
    seq_t seq = {0};

    for (i = 0; i < n_thread; ++i) 
        fseek(out[i], 0, SEEK_END);

    rewind(in);
    c = fskip('@', in);
    i = 0;
    while (c != EOF) {
        fskip('\n', in);
        fget_seq(&seq, '+', in);
        fskip('\n', in);
        c = fastq_skip_qual(seq.len, in);

        seq_write(&seq, out[i]);
        i = (i + 1) % n_thread;
    }
}

int cmp_int64(const void *a, const void *b)
{
    return *(int64_t *)a - *(int64_t *)b;
}

int cmp_kmer_uniq(const void *a, const void *b)
{
/*
    if (*(uint64_t *)a < *(uint64_t *)b)
        return -1;
    else
        return 1;
*/
    if (*(uint64_t *)a < *(uint64_t *)b)
        return -1;
    else if (*(uint64_t *)a > *(uint64_t *)b)
        return 1;
    else
        return 0;
}

int cmp_lmer_uniq(const void *a, const void *b)
{
    uint64_t len;

    len = g_l_len;

    if (len % 32) {
        ((uint64_t *)a)[(len-1)/32] &= ~(~0ull << (len%32*2));
        ((uint64_t *)b)[(len-1)/32] &= ~(~0ull << (len%32*2));
    }
    for (len = (len-1)/32; len > 0; --len)
        if (((uint64_t *)a)[len] != ((uint64_t *)b)[len])
            return ((int64_t)((uint64_t *)a)[len] < (int64_t)((uint64_t *)b)[len] ? -1 : 1);
    return ((int64_t)((uint64_t *)a)[len] < (int64_t)((uint64_t *)b)[len] ? -1 : 1);
}

int cmp_position(const void *a, const void *b)
{
    return ((position_t *)a)->id - ((position_t *)b)->id;
}


void print_contig(FILE *contig_fp, char *out_name, const double occ_ratio)
{
    uint8_t *seq;
    uint16_t occ;
    uint64_t i;
    uint64_t len;
    uint64_t max;
    uint64_t n_seq;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    max = 0;
    rewind(contig_fp);
    while (fread(&len, sizeof(uint64_t), 1, contig_fp)) {
        fread(&occ, sizeof(uint16_t), 1, contig_fp);
        fseek(contig_fp, (len+3)/4 * sizeof(uint8_t), SEEK_CUR);
        if ((len+3)/4 > max)
            max = len;
    }

    seq = (uint8_t *)my_malloc(max * sizeof(uint8_t));

    n_seq = 0;
    rewind(contig_fp);
    while (fread(&len, sizeof(uint64_t), 1, contig_fp)) {
        fread(&occ, sizeof(uint16_t), 1, contig_fp);
        occ = occ * occ_ratio + 0.5;
        fread(seq, sizeof(uint8_t), (len+3)/4, contig_fp);

        ++n_seq;
        fprintf(out, ">seq%lu_len%lu_cov%u\n", n_seq, len, occ);
        for (i = 0; i < len; ++i) {
            putc(num2ascii(bseq_get(seq, i)), out);
            if ((i+1)%LINE_LENGTH == 0)
                putc('\n', out);
        }
        if (i % LINE_LENGTH != 0)
            putc('\n', out);
    }

    my_free(seq);
    fclose(out);
}

void contig_init(contig_t *con)
{
    memset(con, 0, sizeof(contig_t));
}

void contig_destroy(contig_t *con)
{
    if (con->seq != NULL) {
        my_free(con->seq);
        my_free(con->seq_pool);
    }
    if (con->cov != NULL)
        my_free(con->cov);
    if (con->seq_fp != NULL)
        fclose(con->seq_fp);
    if (con->cov_fp != NULL)
        fclose(con->cov_fp);
    memset(con, 0, sizeof(contig_t));
}

void contig_read_fasta(contig_t *con, FILE *fp)
{
    int8_t c;
    varr8_t seq;

    varr8_init(&seq);

    if (con->seq_fp == NULL && (con->seq_fp = tmpfile_open()) == NULL) {
        fputs("error: cannot create tmpfile_open!\n", stderr);
        my_exit(1);
    }
    fseek(con->seq_fp, 0, SEEK_END);

    rewind(fp);
    c = fskip('>', fp);
    while (c != EOF) {
        fskip('\n', fp);
        varr8_clear(&seq);
        while ((c = getc(fp)) != '>' && c != EOF) {
            if (c == '\n')
                continue;
            varr8_push(&seq, ascii2num(c));
        }

        fwrite(&(seq.size), sizeof(int64_t), 1, con->seq_fp);
        fwrite(seq.ele, sizeof(int8_t), seq.size, con->seq_fp);
        ++con->n_seq;
        con->seq_pool_size += seq.size;
    }

    varr8_destroy(&seq);
}

void contig_read_fasta_cov(contig_t *con, FILE *fp)
{
    const char cov_symbol[] = "cov";
    uint64_t i;
    int8_t c;
    uint64_t cov;
    uint16_t tmp_u16;
    varr8_t seq;

    varr8_init(&seq);

    if (con->seq_fp == NULL) {
        con->seq_fp = tmpfile_open();
        check_tmpfile(con->seq_fp);
    }

    fseek(con->seq_fp, 0, SEEK_END);
    if (con->cov_fp == NULL) {
        con->cov_fp = tmpfile_open();
        check_tmpfile(con->cov_fp);
    }
    fseek(con->cov_fp, 0, SEEK_END);

    rewind(fp);
    c = fskip('>', fp);
    while (c != EOF) {
        cov = i = 0;
        while ((c = getc(fp)) != '\n' && c != EOF) {
            if (cov_symbol[i] == '\0') {
                if (isdigit(c))
                    cov = cov * 10 + c - 48;
                else {
                    c = fskip('\n', fp);
                    break;
                }
            }
            else if (c == cov_symbol[i])
                ++i;
            else
                i = 0;
        }
        tmp_u16 = (uint16_t)min_u64(cov, UINT16_MAX);
        fwrite(&tmp_u16, sizeof(uint16_t), 1, con->cov_fp);

        varr8_clear(&seq);
        while ((c = getc(fp)) != '>' && c != EOF) {
            if (c == '\n')
                continue;
            varr8_push(&seq, ascii2num(c));
        }

        fwrite(&(seq.size), sizeof(int64_t), 1, con->seq_fp);
        fwrite(seq.ele, sizeof(int8_t), seq.size, con->seq_fp);
        ++con->n_seq;
        con->seq_pool_size += seq.size;
    }

    varr8_destroy(&seq);
}

void contig_load_seq(contig_t *con)
{
    int64_t i;
    int64_t j;

    con->seq = (lseq_t *)my_calloc(con->n_seq, sizeof(lseq_t));
    con->mem_usage += con->n_seq * sizeof(lseq_t);

    con->seq_pool = (int8_t *)my_malloc(con->seq_pool_size * sizeof(int8_t));
    con->mem_usage += con->seq_pool_size * sizeof(int8_t);

    rewind(con->seq_fp);
    j = 0;
    for (i = 0; i < con->n_seq; ++i) {
        fread(&(con->seq[i].len), sizeof(int64_t), 1, con->seq_fp);
        fread(&(con->seq_pool[j]), sizeof(int8_t), con->seq[i].len, con->seq_fp);
        con->seq[i].seq = &(con->seq_pool[j]);
        j += con->seq[i].len;
    }
}

void contig_load_cov(contig_t *con)
{
    int64_t i;

    con->cov = (uint16_t *)my_malloc(con->n_seq * sizeof(uint16_t));
    con->mem_usage += con->n_seq * sizeof(uint16_t);

    rewind(con->cov_fp);
    for (i = 0; i < con->n_seq; ++i)
        fread(&(con->cov[i]), sizeof(int16_t), 1, con->cov_fp);
}

void varr8_init(varr8_t *varr)
{
    varr->ele = (uint8_t *)my_malloc(VARR_BUF_SIZE*sizeof(uint8_t));
    varr->size = 0ull;
    varr->capa = VARR_BUF_SIZE;
}

void varr8_destroy(varr8_t *varr)
{
    if (varr->ele != NULL)
        my_free(varr->ele);
    varr->size = 0ull;
    varr->capa = 0ull;
}

void varr64_init(varr64_t *varr)
{
    varr->ele = (uint64_t *)my_malloc(VARR_BUF_SIZE*sizeof(uint64_t));
    varr->size = 0ull;
    varr->capa = VARR_BUF_SIZE;
}

void varr64_destroy(varr64_t *varr)
{
    if (varr->ele != NULL)
        my_free(varr->ele);
    varr->size = 0ull;
    varr->capa = 0ull;
}
