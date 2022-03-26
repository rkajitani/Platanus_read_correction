#include "seqlib.h"

void seqlib_init(seqlib_t *lib)
{
    memset(lib, 0, sizeof(seqlib_t));
}

void seqlib_destroy(seqlib_t *lib)
{
    if (lib->pair_fp != NULL)
        fclose(lib->pair_fp);
    if (lib->ins_len_fp != NULL)
        fclose(lib->ins_len_fp);
    if (lib->mapped_fp != NULL)
        fclose(lib->mapped_fp);
    memset(lib, 0, sizeof(seqlib_t));
}

void seqlib_estimate_ins(seqlib_t *lib)
{
    int64_t ins_len;
    int64_t tmp;
    varr64_t distr;

    varr64_init(&distr);

    varr64_clear(&distr);
    rewind(lib->ins_len_fp);
    while (fread(&ins_len, sizeof(int64_t), 1, lib->ins_len_fp)) {
        if (ins_len >= distr.size) {
            tmp = distr.size;
            varr64_resize(&distr, ins_len + 1);
            memset(&(distr.ele[tmp]), 0, (ins_len - tmp + 1) * sizeof(int64_t));
        }
        ++distr.ele[ins_len];
    }

/*
	int64_t i;
	for (i = 0; i < distr.size; ++i)
		printf("%ld\t%ld\n", i, distr.ele[i]);
	exit(0);
*/

    if (distr.size == 0) {
		lib->ave_ins = 0;
		lib->sd_ins = 0;
		return;
    }

    distr_truncate((int64_t *)distr.ele, distr.size, INS_DISTR_TRUNC);

    lib->ave_ins = distr_ave((int64_t *)distr.ele, distr.size) + 0.5;
    lib->sd_ins = distr_sd((int64_t *)distr.ele, distr.size, lib->ave_ins) + 0.5;

    varr64_destroy(&distr);
}

int cmp_seqlib(const void *a, const void *b)
{
    if (((seqlib_t *)a)->ave_ins < ((seqlib_t *)b)->ave_ins)
        return -1;
    if (((seqlib_t *)a)->ave_ins > ((seqlib_t *)b)->ave_ins)
        return 1;
    else 
        return 0;
}

int cmp_seqlib_ptr(const void *a, const void *b)
{
    if ((*(seqlib_t **)a)->ave_ins < (*(seqlib_t **)b)->ave_ins)
        return -1;
    if ((*(seqlib_t **)a)->ave_ins > (*(seqlib_t **)b)->ave_ins)
        return 1;
    else 
        return 0;
}

int cmp_seqlib_input(const void *a, const void *b)
{
    if (((seqlib_input_t *)a)->id < ((seqlib_input_t *)b)->id)
        return -1;
    if (((seqlib_input_t *)a)->id > ((seqlib_input_t *)b)->id)
        return 1;
    else 
        return 0;
}

void distr_truncate(int64_t *distr, const int64_t size, const double edge)
{
    int64_t i;
    int64_t sum;
    int64_t cum;

    if (size <= 1)
        return;

    sum = 0;
    for (i = 1; i < size; ++i)
        sum += i * distr[i];
    cum = 0;
    for (i = 1; (double)cum < sum * edge; ++i) {
        cum += i * distr[i];
        distr[i - 1] = 0;
    }
    cum = (size - 1) * distr[size - 1];
    for (i = size - 2; (double)cum < sum * edge; --i) {
        cum += i * distr[i];
        distr[i + 1] = 0;
    }
}

double distr_ave(const int64_t *distr, const int64_t size)
{
    int64_t i;
    int64_t sum;
    int64_t num;

    sum = num = 0;
    for (i = 0; i < size; ++i) {
        sum += i * distr[i];
        num += distr[i];
    }
    return (double)sum / (double)num;
}

double distr_sd(const int64_t *distr, const int64_t size, const double ave)
{
    int64_t i;
    int64_t num;
    double sum;

    sum = 0.0;
    num = 0;
    for (i = 0; i < size; ++i) {
        sum += ((double)i - ave) * ((double)i - ave) * (double)distr[i];
        num += distr[i];
    }
    if (num > 1)
        return sqrt(sum / (double)(num - 1));
    else
        return 0.0;
}

void seqlib_read_fasta_pe1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '>';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(fp);
	fc = fskip('>', fp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', fp);
		fc = fget_seq_simple(&fseq, '>', fp);

		fskip('\n', fp);
		rc = fget_seq_simple(&rseq, '>', fp);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++(lib->n_pair);
    }

    if (rc == '>') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

void seqlib_read_fasta_pe2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '>';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(ffp);
    rewind(rfp);
	fc = fskip('>', ffp);
	rc = fskip('>', rfp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', ffp);
		fc = fget_seq_simple(&fseq, '>', ffp);

		fskip('\n', rfp);
		rc = fget_seq_simple(&rseq, '>', rfp);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++lib->n_pair;
    }

    if (fc == '>' || rc == '>') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }

}

void seqlib_read_fasta_mp1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '>';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(fp);
	fc = fskip('>', fp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', fp);
		fc = fget_seq_simple(&fseq, '>', fp);
		rev_comp(fseq.base, fseq.len);

		fskip('\n', fp);
		rc = fget_seq_simple(&rseq, '>', fp);
		rev_comp(rseq.base, rseq.len);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++lib->n_pair;
    }

    if (rc == '>') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

void seqlib_read_fasta_mp2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '>';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(ffp);
    rewind(rfp);
	fc = fskip('>', ffp);
	rc = fskip('>', rfp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', ffp);
		fc = fget_seq_simple(&fseq, '>', ffp);
		rev_comp(fseq.base, fseq.len);

		fskip('\n', rfp);
		rc = fget_seq_simple(&rseq, '>', rfp);
		rev_comp(rseq.base, rseq.len);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++lib->n_pair;
    }

    if (fc == '>' || rc == '>') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

void seqlib_read_fastq_pe1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '@';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(fp);
	fc = fskip('@', fp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', fp);
		fget_seq_simple(&fseq, '+', fp);
		fskip('\n', fp);
		fc = fastq_skip_qual(fseq.len, fp);

		fskip('\n', fp);
		fget_seq_simple(&rseq, '+', fp);
		fskip('\n', fp);
		rc = fastq_skip_qual(rseq.len, fp);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++(lib->n_pair);
    }

    if (rc == '@') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

void seqlib_read_fastq_pe2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '@';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(ffp);
    rewind(rfp);
	fc = fskip('@', ffp);
	rc = fskip('@', rfp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', ffp);
		fget_seq_simple(&fseq, '+', ffp);
		fskip('\n', ffp);
		fc = fastq_skip_qual(fseq.len, ffp);

		fskip('\n', rfp);
		fget_seq_simple(&rseq, '+', rfp);
		fskip('\n', rfp);
		rc = fastq_skip_qual(rseq.len, rfp);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++lib->n_pair;
    }

    if (fc == '@' || rc == '@') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

void seqlib_read_fastq_mp1_mt(seqlib_t *lib, FILE *fp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '@';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(fp);
	fc = fskip('@', fp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', fp);
		fget_seq_simple(&fseq, '+', fp);
		fskip('\n', fp);
		fc = fastq_skip_qual(fseq.len, fp);
		rev_comp(fseq.base, fseq.len);

		fskip('\n', fp);
		fget_seq_simple(&rseq, '+', fp);
		fskip('\n', fp);
		rc = fastq_skip_qual(rseq.len, fp);
		rev_comp(rseq.base, rseq.len);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++lib->n_pair;
    }

    if (rc == '@') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

void seqlib_read_fastq_mp2_mt(seqlib_t *lib, FILE *ffp, FILE *rfp, const int64_t n_thread)
{
    int64_t i;
    int fc;
    int rc = '@';
    seq_t fseq = {0};
    seq_t rseq = {0};

    for (i = 0; i < n_thread; ++i)
        fseek(lib[i].pair_fp, 0, SEEK_END);

    rewind(ffp);
    rewind(rfp);
	fc = fskip('@', ffp);
	rc = fskip('@', rfp);
    i = 0;
    while (fc != EOF && rc != EOF) {
		fskip('\n', ffp);
		fget_seq_simple(&fseq, '+', ffp);
		fskip('\n', ffp);
		fc = fastq_skip_qual(fseq.len, ffp);
		rev_comp(fseq.base, fseq.len);

		fskip('\n', rfp);
		fget_seq_simple(&rseq, '+', rfp);
		fskip('\n', rfp);
		rc = fastq_skip_qual(rseq.len, rfp);
		rev_comp(rseq.base, rseq.len);

        seq_write(&fseq, lib[i].pair_fp);
        seq_write(&rseq, lib[i].pair_fp);
        i = (i + 1) % n_thread;

        lib->total_len += fseq.len + rseq.len;
        ++lib->n_pair;
    }

    if (fc == '@' || rc == '@') {
        fputs("error: the number of reads differ from the other!\n", stderr);
        my_exit(1);
    }
}

int option_pair(int argc, char **argv, int i, const pair_file_format_t pair_fmt, seqlib_input_t *val, int64_t *n_val)
{
    char *name;
    int64_t lib_id;

    name = argv[i];
    lib_id = suffix_atoi(argv[i]);

    if (pair_fmt == PE1 || pair_fmt == MP1) {
        for (++i; i < argc && argv[i][0] != '-'; ++i) {
            val[*n_val].id = lib_id;
            val[*n_val].pair_fmt = pair_fmt;
			val[*n_val].seq_fmt = check_seq_file_format(argv[i]);
            val[*n_val].fp[0] = fopen(argv[i], "r");
            ++(*n_val);
        }
    }
    else {
        for (++i; i < argc && argv[i][0] != '-'; i += 2) {
            val[*n_val].id = lib_id;
            val[*n_val].pair_fmt = pair_fmt;
			val[*n_val].seq_fmt = check_seq_file_format(argv[i]);
            val[*n_val].fp[0] = fopen(argv[i], "r");
            if (i+1 == argc || argv[i+1][0] == '-') {
                fprintf(stderr, "error: %s must specify lib_id fwd1 rev1 [fwd2 rev2 ...]\n\n", name);    
                my_exit(1);
            }
			if (val[*n_val].seq_fmt != check_seq_file_format(argv[i+1])) {
                fprintf(stderr, "error: pair-file %s and %s must have common file format (fasta or fastq)\n\n", argv[i], argv[i+1]);
				my_exit(1);
			}
            val[*n_val].fp[1] = fopen(argv[i+1], "r");
            ++(*n_val);
        }
    }

    return i;
}

int option_indexed_int(int argc, char **argv, int i, int64_t *val)
{
    char *name;
	int64_t lib_id;

    name = argv[i];
	lib_id = suffix_atoi(argv[i]);

    ++i;
    if (i >= argc || lib_id < 0 || lib_id > MAX_FILE_NUM || atoi(argv[i]) < 0) {
        fprintf(stderr, "error: %s must specify INT (min_ins[0, ])\n\n", name);
        my_exit(1);
    }
    val[lib_id - 1] = atoi(argv[i]);
    ++i;

    return i;
}
