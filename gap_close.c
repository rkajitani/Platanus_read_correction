#include "gap_close.h"

void print_gap_close_usage(void)
{
    option_gap_close_t def = {};

    option_gap_close_init(&def);

    fprintf(stderr, "\nUsage: platanus gap_close [Options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR                             : prefix of output file (default %s, length <= %u)\n", def.o, MAX_FILE_LEN);
    fprintf(stderr, "    -c FILE1 [FILE2 ...]               : scaffold_file (fasta format)\n");
    fprintf(stderr, "    -ip{INT} PAIR1 [PAIR2 ...]         : INT = lib_id, PAIR1 = inward_pair_file (reads in 1 file, fasta or fastq)\n");
    fprintf(stderr, "    -op{INT} PAIR1 [PAIR2 ...]         : INT = lib_id, PAIR1 = outward_pair_file (reads in 1 file, fasta or fastq)\n");
    fprintf(stderr, "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : INT = lib_id, (FWD1, REV1) = inward_pair_files (reads in 2 files, fasta or fastq)\n");
    fprintf(stderr, "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : INT = lib_id, (FWD1, REV1) = outward_pair_files (reads in 2 files, fasta or fastq)\n");
    fprintf(stderr, "    -s INT                             : mapping seed length (default %ld)\n", def.s);
    fprintf(stderr, "    -v INT                             : minimum overlap length (default %ld)\n", def.v);
    fprintf(stderr, "    -e FLOAT                           : maximum error rate of overlap (identity, default %.2f)\n", def.e);
    fprintf(stderr, "    -m INT                             : memory limit (GB, >= 1, default %lu)\n", def.m / 1073741824ul);
    fprintf(stderr, "    -t INT                             : number of threads (<= %ld, default 1)\n", def.t);

    fprintf(stderr, "Output:\n");
    fprintf(stderr, "    PREFIX_gapClosed.fa\n");

    option_gap_close_destroy(&def);
}

void option_gap_close_init(option_gap_close_t *opt)
{
    strcpy(opt->o, "out");
    opt->s = 32;
    opt->v = 32;
    opt->e = 0.05;
    opt->m = 4 * GIBIBYTE;
    opt->t = 1;
    opt->n_c = 0;
    opt->n_pair = 0;
}

void option_gap_close_destroy(const option_gap_close_t *opt)
{
    uint64_t i;

    for (i = 0; i < opt->n_c; ++i)
        fclose(opt->c[i]);
    for (i = 0; i < opt->n_pair; ++i) {
        switch (opt->pair[i].pair_fmt) {
            case PE1:
            case MP1:
                fclose(opt->pair[i].fp[0]);
                break;
            case PE2:
            case MP2:
                fclose(opt->pair[i].fp[0]);
                fclose(opt->pair[i].fp[1]);
        }
    }
}

uint8_t option_gap_close_parse(option_gap_close_t *opt, int argc, char **argv)
{
    uint64_t i;

    for (i = 2; i < argc;) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help") || !strcmp(argv[i], "--help"))
            return -1;
        else if (!strcmp(argv[i], "-o"))
            i = option_string(argc, argv, i, MAX_FILE_LEN, opt->o);
        else if (!strcmp(argv[i], "-c"))
            i = option_multi_fasta_file(argc, argv, i, opt->c, &(opt->n_c));
        else if (!strcmp(argv[i], "-s"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->s));
        else if (!strcmp(argv[i], "-v"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->v));
        else if (!strcmp(argv[i], "-e"))
            i = option_float(argc, argv, i, 0.0, 1.0, &(opt->e));
        else if (!strcmp(argv[i], "-t"))
            i = option_int(argc, argv, i, 1, MAX_THREAD, &(opt->t));
        else if (!strcmp(argv[i], "-m")) {
            i = option_int(argc, argv, i, 1, INT64_MAX, &(opt->m));
            opt->m *= GIBIBYTE;
        }
        else if (strstr(argv[i], "-ip") == argv[i])
            i = option_pair(argc, argv, i, PE1, opt->pair, &(opt->n_pair));
        else if (strstr(argv[i], "-IP") == argv[i])
            i = option_pair(argc, argv, i, PE2, opt->pair, &(opt->n_pair));
        else if (strstr(argv[i], "-op") == argv[i])
            i = option_pair(argc, argv, i, MP1, opt->pair, &(opt->n_pair));
        else if (strstr(argv[i], "-OP") == argv[i])
            i = option_pair(argc, argv, i, MP2, opt->pair, &(opt->n_pair));
        else {
            fprintf(stderr, "error: wrong option \"%s\"\n\n", argv[i]);
            return -1;
        }
    }

    if (opt->n_pair == 0) {
        fputs("error: no pair_read file\n", stderr);
        return -1;
    }

    if (opt->n_c == 0) {
        fputs("error: no contig file\n", stderr);
        return -1;
    }
    return 0;
}

void gap_close_mt(option_gap_close_t *opt)
{
    char out_name[MAX_FILE_LEN + MAX_EXTENSION_NAME + 1];
    const int64_t n_thread = opt->t;
    int64_t i;
    int64_t j;
    const int64_t seed_len = opt->s;
    const int64_t min_overlap = opt->v;
    int64_t key_len = seed_len;
    int64_t n_lib;
    int64_t n_pair[MAX_FILE_NUM] = {0};
    contig_t con;
    mapper_t mp;
    closer_t cl;
    seqlib_t *lib[MAX_FILE_NUM];
    seqlib_input_t *pair[MAX_FILE_NUM];
    double max_mis_rate = opt->e;
    FILE *gap_seq_fp;
    FILE *unused_fp;

    if (key_len > 32)
        key_len = 32;
    contig_init(&con);
    mapper_init(&mp, seed_len, key_len, opt->m);

    qsort(opt->pair, opt->n_pair, sizeof(seqlib_input_t), cmp_seqlib_input);

    opt->pair[opt->n_pair].id = 0;
    n_lib = 0;
    for (i = 0; i < opt->n_pair; ++i) {
        ++(n_pair[n_lib]);
        if (opt->pair[i].id != opt->pair[i + 1].id) {
            pair[n_lib] = &(opt->pair[i - n_pair[n_lib] + 1]);
            ++n_lib;
        }
    }

    omp_set_num_threads(n_thread);

#    pragma omp parallel for private(j) schedule(static, 1)
    for (i = -1; i < n_lib; ++i) {
        if (i == -1) {
            for (j = 0; j < opt->n_c; ++j)
                contig_read_fasta(&con, opt->c[j]);
            contig_load_seq(&con);
            mapper_set_contig(&mp, &con);
            mapper_make_table(&mp);
        }
        else {
            lib[i] = (seqlib_t *)my_malloc(n_thread * sizeof(seqlib_t));
            for (j = 0; j < n_thread; ++j) {
                seqlib_init(&(lib[i][j]));    
                lib[i][j].pair_fp = tmpfile_open();
                check_tmpfile(lib[i][j].pair_fp);
            }
            for (j = 0; j < n_pair[i]; ++j) {
                switch (pair[i][j].pair_fmt) {
                    case PE1:
						if (pair[i][j].seq_fmt == FASTA)
							seqlib_read_fasta_pe1_mt(lib[i], pair[i][j].fp[0], n_thread);
						else
							seqlib_read_fastq_pe1_mt(lib[i], pair[i][j].fp[0], n_thread);
                        break;
                    case PE2:
						if (pair[i][j].seq_fmt == FASTA)
							seqlib_read_fasta_pe2_mt(lib[i], pair[i][j].fp[0], pair[i][j].fp[1], n_thread);
						else
							seqlib_read_fastq_pe2_mt(lib[i], pair[i][j].fp[0], pair[i][j].fp[1], n_thread);
                        break;
                    case MP1:
						if (pair[i][j].seq_fmt == FASTA)
							seqlib_read_fasta_mp1_mt(lib[i], pair[i][j].fp[0], n_thread);
						else
							seqlib_read_fastq_mp1_mt(lib[i], pair[i][j].fp[0], n_thread);
                        break;
                    case MP2:
						if (pair[i][j].seq_fmt == FASTA)
							seqlib_read_fasta_mp2_mt(lib[i], pair[i][j].fp[0], pair[i][j].fp[1], n_thread);
						else
							seqlib_read_fastq_mp2_mt(lib[i], pair[i][j].fp[0], pair[i][j].fp[1], n_thread);
                        break;
                }
            }
        }
    }

    for (i = 0; i < n_lib; ++i) {
        fprintf(stderr, "[LIBRARY %ld]\n", i + 1);
        lib[i]->ave_len = (double)(lib[i]->total_len) / (2 * lib[i]->n_pair) + 0.5;
        mapper_map_pair_nolink_mt(&mp, lib[i], n_thread);
        seqlib_estimate_ins(lib[i]);
        fprintf(stderr, "AVE_INS = %ld\nSD_INS = %ld\n", lib[i]->ave_ins, lib[i]->sd_ins);
    }

    qsort(lib, n_lib, sizeof(seqlib_t *), cmp_seqlib_ptr);

    closer_init(&cl, &con, min_overlap, max_mis_rate, opt->m);
    mapper_destroy(&mp);

    closer_make_gap_table(&cl);

    unused_fp = tmpfile_open();
    check_tmpfile(unused_fp);

    for (i = 0; i < n_lib; ++i) {
        fprintf(stderr, "[LIBRARY %ld]\n", i + 1);
        gap_seq_fp = closer_save_gap_covering_reads_mt(&cl, lib[i], n_thread);

        for (j = 0; j < n_thread; ++j) {
            fclose(lib[i][j].mapped_fp);
            lib[i][j].mapped_fp = NULL;
        }
        seqlib_destroy(lib[i]);
        my_free(lib[i]);

        closer_load_local_reads(&cl, gap_seq_fp);
        fclose(gap_seq_fp);
        closer_local_assemble_mt(&cl, n_thread);
        if (n_lib > 1)
            closer_save_unused_reads(&cl, &unused_fp);
    }

    if (n_lib > 1) {
        fprintf(stderr, "[ALL LIBRARY]\n");
        closer_load_unused_reads(&cl, unused_fp);
        closer_local_assemble_mt(&cl, n_thread);
    }

    fclose(unused_fp);

    strcpy(out_name, opt->o);
    strcat(out_name, "_gapClosed.fa");
    closer_print_seq(&cl, out_name);

    closer_destroy(&cl);
    contig_destroy(&con);

    fputs("gap_close completed!\n", stderr);
}

void closer_init(closer_t *cl, const contig_t *con, const int64_t min_ol, const double max_mis_rate, const int64_t mem_lim)
{
    memset(cl, 0, sizeof(closer_t));
    cl->sca = con->seq;
    cl->n_sca = con->n_seq;
    cl->min_ol = min_ol;
    cl->max_mis_rate = max_mis_rate;
    cl->mem_lim = mem_lim;
    cl->mem_usage += con->mem_usage;
}

void closer_destroy(closer_t *cl)
{
    int64_t i;

    if (cl->gap_table != NULL)
        my_free(cl->gap_table);
    if (cl->gap != NULL) {
        for (i = 0; i < cl->n_gap; ++i) {
            if (cl->gap[i].len > 0)
                my_free(cl->gap[i].seq);
        }
        my_free(cl->gap);
    }
        
    memset(cl, 0, sizeof(closer_t));
}

void closer_make_gap_table(closer_t *cl)
{
    int64_t total_gap_len;
    int32_t i;
    int32_t j;

    fputs("making hash table of gaps...\n", stderr);

    cl->n_gap = total_gap_len = 0;
    for (i = 0; i < cl->n_sca; ++i) {
        j = 0;
        while (j < cl->sca[i].len) {
            if (cl->sca[i].seq[j] == 4) {
                while (cl->sca[i].seq[j] == 4 && j < cl->sca[i].len) {
                    ++total_gap_len;
                    ++j;
                }
                ++cl->n_gap;
            }
            ++j;
        }
    }

    cl->gap = (gap_t *)my_calloc(cl->n_gap, sizeof(gap_t));
    cl->mem_usage += cl->n_gap * sizeof(gap_t);

    cl->gap_table_size = 1;
    cl->index_len = 0;
    for (i = 0; i < 64; ++i) {
        cl->gap_table_size *= 2;
        ++cl->index_len;
        if (cl->gap_table_size * MAX_LOAD_FACTOR > total_gap_len)
            break;
    }
    cl->mem_usage += cl->gap_table_size * sizeof(gap_rec_t);
    check_mem_usage(cl->mem_usage, cl->mem_lim);
    fprintf(stderr, "MEM_USAGE=%ld\n", cl->mem_usage);

    cl->gap_table = my_calloc(cl->gap_table_size, sizeof(gap_rec_t));

    cl->n_gap = 0;
    for (i = 0; i < cl->n_sca; ++i) {
        j = 0;
        while (j < cl->sca[i].len) {
            if (cl->sca[i].seq[j] == 4) {
                cl->gap[cl->n_gap].id = i + 1;
                cl->gap[cl->n_gap].st = j;
                while (cl->sca[i].seq[j] == 4 && j < cl->sca[i].len) {
                    closer_insert_gap(cl, i + 1, j, cl->n_gap);
                    ++j;
                }
                cl->gap[cl->n_gap].ed = j;
                ++cl->n_gap;
            }
            ++j;
        }
    }
}

void closer_insert_gap(closer_t *cl, const int32_t id, const int32_t ofst, const int64_t val)
{
    uint64_t key;
    int64_t index;
    int64_t step;

    key = ((uint64_t)id << 32) | (uint64_t)ofst;
    index = hash64(key, cl->index_len) & (cl->gap_table_size-1);
    if (gap_rec_is_emp(cl->gap_table[index])) {
        cl->gap_table[index].key = key;
        cl->gap_table[index].val = val;
    }
    else {
        step = rehash64(key, cl->index_len);
        index = (index+step) & (cl->gap_table_size-1);
        while (!(gap_rec_is_emp(cl->gap_table[index])))
            index = (index+step) & (cl->gap_table_size-1);
        cl->gap_table[index].key = key;
        cl->gap_table[index].val = val;
    }
}

int64_t closer_get_gap_id(const closer_t *cl, const int32_t id, const int32_t ofst)
{
    uint64_t key;
    int64_t index;
    int64_t step;

    key = ((uint64_t)id << 32) | (uint64_t)ofst;
    index = hash64(key, cl->index_len) & (cl->gap_table_size-1);
    if (!(gap_rec_is_emp(cl->gap_table[index])) && cl->gap_table[index].key == key)
        return cl->gap_table[index].val;
    else {
        step = rehash64(key, cl->index_len);
        index = (index+step) & (cl->gap_table_size-1);
        while (!(gap_rec_is_emp(cl->gap_table[index]))) {
            if (cl->gap_table[index].key == key)
                return cl->gap_table[index].val;
            index = (index+step) & (cl->gap_table_size-1);
        }
        return -1;
    }
}

bool make_consensus_from_reads(const lseq_t *seq, const int64_t n_seq, const double threshold, lseq_t *cons)
{
    int32_t i;
    int32_t j;
    int32_t tmp;
    int32_t most;
    int32_t most_len = 0;
    int32_t len_count[MAX_READ_LEN + 1] = {0};
    int32_t base_count[5];

    for (i = 0; i < n_seq; ++i)
        ++len_count[seq[i].len];

    most = 0;
    for (i = 0; i <= MAX_READ_LEN; ++i) {
        if (len_count[i] > most) {
            most = len_count[i];
            most_len = i;
        }
    }
    if ((double)most / n_seq < threshold)
        return false;

    if (most_len ==  0) {
        cons->len = 0;
        return true;
    }
    
    most = 0;
    for (i = 0; i < most_len; ++i) {
        for (j = 0; j < 4; ++j)
            base_count[j] = 0;
        for (j = 0; j < n_seq; ++j) {
            if (seq[j].len == most_len)
                ++base_count[seq[j].seq[i]];
        }
        tmp = 0;
        for (j = 0; j < 4; ++j) {
            if (base_count[j] > tmp) {
                tmp = base_count[j];
                (cons->seq)[i] = j;
            }
        }
        most += tmp;
    }

    if ((double)most / (most_len * len_count[most_len]) < threshold)
        return false;
    
    cons->len = most_len;
    return true;
}

void closer_print(const closer_t *cl, FILE *close_seq_fp, const char *out_name)
{
    int8_t *base_pool;
    int32_t i;
    int32_t j;
    int32_t k;
    int32_t l;
    int64_t base_pool_size;
    int64_t gap_id;
    FILE *out;
    lseq_t *cons;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    cons = (lseq_t *)my_malloc(cl->n_gap * sizeof(lseq_t));
    for (i = 0; i < cl->n_gap; ++i)
        cons[i].len = -1;

    base_pool_size = 0;
    rewind(close_seq_fp);
    while (fread(&gap_id, sizeof(int64_t), 1, close_seq_fp)) {
        fread(&(cons[gap_id].len), sizeof(int32_t), 1, close_seq_fp);
        fseek(close_seq_fp, cons[gap_id].len * sizeof(int8_t), SEEK_CUR);
        base_pool_size += cons[gap_id].len;
    }

    base_pool = (int8_t *)my_malloc(base_pool_size * sizeof(int8_t));

    base_pool_size = 0;
    rewind(close_seq_fp);
    while (fread(&gap_id, sizeof(int64_t), 1, close_seq_fp)) {
        fread(&(cons[gap_id].len), sizeof(int32_t), 1, close_seq_fp);
        fread(&(base_pool[base_pool_size]), sizeof(int8_t), cons[gap_id].len, close_seq_fp);
        cons[gap_id].seq = &(base_pool[base_pool_size]);
        base_pool_size += cons[gap_id].len;
    }

    for (i = 0; i < cl->n_sca; ++i) {
        fprintf(out, ">scaffold%d\n", i + 1);
        l = 0;
        j = 0;
        while (j < cl->sca[i].len) {
            if (cl->sca[i].seq[j] == 4) {
                gap_id = closer_get_gap_id(cl, (int32_t)(i + 1), (int32_t)j);
                if (cons[gap_id].len >= 0) {
                    for (k = 0; k < cons[gap_id].len; ++k) {
                        putc(num2ascii(cons[gap_id].seq[k]), out);
                        ++l;
                        if (l % LINE_LENGTH == 0)
                            putc('\n', out);
                    }
                    while (cl->sca[i].seq[j] == 4 && j < cl->sca[i].len)
                        ++j;
                }
                else {
                    while (cl->sca[i].seq[j] == 4 && j < cl->sca[i].len) {
                        ++j;
                        putc('N', out);
                        ++l;
                        if (l % LINE_LENGTH == 0)
                            putc('\n', out);
                    }
                }
            }
            else {
                putc(num2ascii(cl->sca[i].seq[j]), out);
                ++l;
                if (l % LINE_LENGTH == 0)
                    putc('\n', out);
                ++j;
            }
        }
        if (l % LINE_LENGTH != 0)
            putc('\n', out);
    }

    fclose(out);
    my_free(base_pool);
}

FILE *closer_save_gap_covering_reads_mt(const closer_t *cl, const seqlib_t *lib, const int64_t n_thread)
{
    int8_t c;
    int64_t i;
    int64_t j;
    int64_t ed;
    int64_t tol;
    int64_t gap_id;
    position_t fpos;
    position_t rpos;
    seq_t fwd;
    seq_t rev;
    FILE *tmp_fp[MAX_THREAD];

    fputs("saving reads covering gaps...\n", stderr);

    tol = 3 * lib->sd_ins;

    omp_set_num_threads(n_thread);

#    pragma omp parallel for  schedule(static, 1) private(j, ed, gap_id, fwd, rev, fpos, rpos)
    for (i = 0; i < n_thread; ++i) {
        tmp_fp[i] = tmpfile_open();
        check_tmpfile(tmp_fp[i]);

        rewind (lib[i].mapped_fp);
        while (fread(&fpos, sizeof(position_t), 1, lib[i].mapped_fp)) {
            fread(&(fwd.len), sizeof(uint16_t), 1, lib[i].mapped_fp);
            fread(fwd.base, sizeof(uint8_t), fwd.len, lib[i].mapped_fp);

            fread(&rpos, sizeof(position_t), 1, lib[i].mapped_fp);
            fread(&(rev.len), sizeof(uint16_t), 1, lib[i].mapped_fp);
            fread(rev.base, sizeof(uint8_t), rev.len, lib[i].mapped_fp);
            
            if (fpos.id != 0) {
                if (fpos.id > 0) {
                    j = max_64(fpos.ofst, fpos.ofst + lib->ave_ins - tol - rev.len);
                    j = min_64(j, cl->sca[fpos.id - 1].len - 1);
                    j = max_64(j, 0);
                    ed = min_64(fpos.ofst + lib->ave_ins + tol, cl->sca[fpos.id - 1].len);
                    rev_comp(rev.base, rev.len);
                }
                else {
                    fpos.id *= -1;
                    j = min_64(fpos.ofst, fpos.ofst - lib->ave_ins - tol);
                    j = max_64(j, 0);
                    ed = min_64(fpos.ofst - lib->ave_ins + tol + rev.len, cl->sca[fpos.id - 1].len);
                }
                while (j < ed) {
                    if (cl->sca[fpos.id - 1].seq[j] == 4) {
                        gap_id = closer_get_gap_id(cl, fpos.id, j);
                        if (cl->gap[gap_id].len == 0) {
                            fwrite(&gap_id, sizeof(int64_t), 1, tmp_fp[i]);
                            fwrite(&(rev.len), sizeof(uint16_t), 1, tmp_fp[i]);
                            fwrite(rev.base, sizeof(uint8_t), rev.len, tmp_fp[i]);
                        }
                        while (cl->sca[fpos.id - 1].seq[j] == 4 && j < ed)
                            ++j;
                    }
                    ++j;
                }
            }

            if (rpos.id != 0) {
                if (rpos.id > 0) {
                    j = max_64(rpos.ofst, rpos.ofst + lib->ave_ins - tol - fwd.len);
                    j = min_64(j, cl->sca[rpos.id - 1].len - 1);
                    j = max_64(j, 0);
                    ed = min_64(rpos.ofst + lib->ave_ins + tol, cl->sca[rpos.id - 1].len);
                    rev_comp(fwd.base, fwd.len);
                }
                else {
                    rpos.id *= -1;
                    j = min_64(rpos.ofst, rpos.ofst - lib->ave_ins - tol);
                    j = max_64(j, 0);
                    ed = min_64(rpos.ofst - lib->ave_ins + tol + fwd.len, cl->sca[rpos.id - 1].len);
                }
                while (j < ed) {
                    if (cl->sca[rpos.id - 1].seq[j] == 4) {
                        gap_id = closer_get_gap_id(cl, rpos.id, j);
                        if (cl->gap[gap_id].len == 0) {
                            fwrite(&gap_id, sizeof(int64_t), 1, tmp_fp[i]);
                            fwrite(&(fwd.len), sizeof(uint16_t), 1, tmp_fp[i]);
                            fwrite(fwd.base, sizeof(uint8_t), fwd.len, tmp_fp[i]);
                        }
                        while (cl->sca[rpos.id - 1].seq[j] == 4 && j < ed)
                            ++j;
                    }
                    ++j;
                }
            }
        }
    }

    for (i = 1; i < n_thread; ++i) {
        rewind(tmp_fp[i]);
        while (fread(&c, sizeof(int8_t), 1, tmp_fp[i]))
            putc(c, tmp_fp[0]);
        fclose(tmp_fp[i]);
    }

    return tmp_fp[0];
}
            
void closer_load_local_reads(closer_t *cl, FILE *gap_seq_fp)
{
    int64_t i;
    int64_t gap_id;
    uint16_t len;

    fputs("loading reads covering gaps...\n", stderr);

    rewind(gap_seq_fp);
    while (fread(&gap_id, sizeof(int64_t), 1, gap_seq_fp)) {
        fread(&len, sizeof(uint16_t), 1, gap_seq_fp);
        cl->gap[gap_id].len += len + 1;
        fseek(gap_seq_fp, len * sizeof(int8_t), SEEK_CUR);
    }

    for (i = 0; i < cl->n_gap; ++i) {
        if (cl->gap[i].len <= 0 || cl->gap[i].seq != NULL)
            continue;
        cl->gap[i].seq = (int8_t *)my_malloc(cl->gap[i].len * sizeof(int8_t));
        cl->gap[i].len = 0;
    }

    rewind(gap_seq_fp);
    while (fread(&gap_id, sizeof(int64_t), 1, gap_seq_fp)) {
        fread(&len, sizeof(uint16_t), 1, gap_seq_fp);
        cl->gap[gap_id].seq[cl->gap[gap_id].len] = 4;
        fread(&(cl->gap[gap_id].seq[cl->gap[gap_id].len + 1]), sizeof(int8_t), len, gap_seq_fp);
        cl->gap[gap_id].len += len + 1;
    }
}

void closer_local_assemble_mt(const closer_t *cl, const int64_t n_thread)
{
    int64_t i;
    int64_t j;
    int64_t n_closed;
    qassembler_t qa[MAX_THREAD];
    gap_t *gap_p;
    gap_t **buf;
    int64_t buf_size;

    fputs("assembling localized reads...\n", stderr);

    omp_set_num_threads(n_thread);

    n_closed = 0;
#    pragma omp parallel for  schedule(static, 1) private(j, buf, buf_size, gap_p) reduction(+: n_closed)
    for (i = 0; i < n_thread; ++i) {
        buf = (gap_t **)my_malloc((cl->n_gap / n_thread + 1) * sizeof(gap_t *));
        buf_size = 0;
        for (j = i; j < cl->n_gap; j += n_thread) {
            buf[buf_size] = (&cl->gap[j]);
            ++buf_size;
        }
        qsort(buf, buf_size, sizeof(gap_t *), cmp_gap_ptr);
        qassembler_init(&(qa[i]), CLOSER_MIN_K, CLOSER_MAX_K, CLOSER_MIN_COV);
        for (j = 0; j < buf_size; ++j) {
            gap_p = buf[j];
            if (gap_p->len <= 0 || gap_p->seq[0] != 4)
                continue;
            qassembler_assemble(&(qa[i]), gap_p->seq, gap_p->len);
            n_closed += qassembler_close_gap(&(qa[i]), &(cl->sca[gap_p->id - 1]), gap_p, cl->min_ol, cl->max_mis_rate);
            qassembler_init_lgraph(&(qa[i]));
        }
        my_free(buf);
        qassembler_destroy(&(qa[i]));
    }

    fprintf(stderr, "NUM_CLOSED_GAPS = %ld\n", n_closed);
}

void closer_save_unused_reads(closer_t *cl, FILE **unused_fpp)
{
    int64_t i;

    fseek(*unused_fpp, 0, SEEK_END);
    for (i = 0; i < cl->n_gap; ++i) {
        if (cl->gap[i].len <= 0 || cl->gap[i].seq[0] != 4)
            continue;
        fwrite(&i, sizeof(int64_t), 1, *unused_fpp);
        fwrite(&(cl->gap[i].len), sizeof(int32_t), 1, *unused_fpp);
        fwrite(cl->gap[i].seq, sizeof(int8_t), cl->gap[i].len, *unused_fpp);
        my_free(cl->gap[i].seq);
        cl->gap[i].seq = NULL;
        cl->gap[i].len = 0;
    }
}

void closer_load_unused_reads(closer_t *cl, FILE *unused_fp)
{
    int64_t i;
    int64_t gap_id;
    int32_t len;
    int8_t *is_closed;

    is_closed = (int8_t *)my_calloc(cl->n_gap, sizeof(int8_t));

    for (i = 0; i < cl->n_gap; ++i) {
        if (cl->gap[i].len != 0)
            is_closed[i] = 1;
    }

    rewind(unused_fp);
    while (fread(&gap_id, sizeof(int64_t), 1, unused_fp)) {
        fread(&len, sizeof(int32_t), 1, unused_fp);
        if (!(is_closed[gap_id]))
            cl->gap[gap_id].len += len;
        fseek(unused_fp, len * sizeof(int8_t), SEEK_CUR);
    }

    for (i = 0; i < cl->n_gap; ++i) {
        if (!(is_closed[i]) && cl->gap[i].len > 0) {
            cl->gap[i].seq = (int8_t *)my_malloc(cl->gap[i].len * sizeof(int8_t));;
            cl->gap[i].len = 0;
        }
    }

    rewind(unused_fp);
    while (fread(&gap_id, sizeof(int64_t), 1, unused_fp)) {
        fread(&len, sizeof(int32_t), 1, unused_fp);
        if (is_closed[gap_id])
            fseek(unused_fp, len * sizeof(int8_t), SEEK_CUR);
        else {
            fread(&(cl->gap[gap_id].seq[cl->gap[gap_id].len]), sizeof(int8_t), len, unused_fp);
            cl->gap[gap_id].len += len;
        }
    }

    my_free(is_closed);
}

void closer_print_seq(closer_t *cl, const char *out_name)
{
    int32_t i;
    int32_t j;
    int32_t k;
    int32_t l;
    int64_t n_closed;
    int64_t gap_id;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    gap_id = n_closed = 0;
    for (i = 0; i < cl->n_sca; ++i) {
        fprintf(out, ">scaffold%d\n", i + 1);
        j = l = 0;
        while (j < cl->sca[i].len) {
            if (cl->sca[i].seq[j] != 4) {
                putc(num2ascii(cl->sca[i].seq[j]), out);
                ++l;
                if (l % LINE_LENGTH == 0)
                    putc('\n', out);
                ++j;
                continue;
            }

            if (cl->gap[gap_id].len > 0) {
                if (cl->gap[gap_id].seq[0] != 4) {
                    for (k = 0; k < cl->gap[gap_id].len; ++k) {
                        putc(num2ascii(cl->gap[gap_id].seq[k]), out);
                        ++l;
                        if (l % LINE_LENGTH == 0)
                            putc('\n', out);
                    }
                    j += cl->gap[gap_id].ed - cl->gap[gap_id].st;
                    ++n_closed;
                }
                else {
                    for (; j < cl->gap[gap_id].ed; ++j) {
                        putc('N', out);
                        ++l;
                        if (l % LINE_LENGTH == 0)
                            putc('\n', out);
                    }
                }
            }
            else if (cl->gap[gap_id].len < 0) {
                j += cl->gap[gap_id].ed - cl->gap[gap_id].st - (cl->gap[gap_id].len + 1);
                ++n_closed;
            }
            else {
                for (; j < cl->gap[gap_id].ed; ++j) {
                    putc('N', out);
                    ++l;
                    if (l % LINE_LENGTH == 0)
                        putc('\n', out);
                }
            }

            ++gap_id;
        }
        if (l % LINE_LENGTH != 0)
            putc('\n', out);
    }

    fclose(out);

    fprintf(stderr, "TOTAL_NUM_CLOSED_GAPS = %ld\n", n_closed);
}

int64_t get_exact_overlap(const int8_t *lseq, const int64_t llen, const int8_t *rseq, const int64_t rlen, const int64_t min_ol)
{
    int64_t i;
    int64_t j;
    int64_t max_ol;

    max_ol = min_64(llen, rlen);

    for (i = max_ol; i >= min_ol; --i) {
        for (j = 0; j < i; ++j) {
            if (lseq[llen - j - 1] != rseq[i - j - 1])
                break;
        }
        if (j == i)
            return i;
    }
    return 0;
}

int64_t get_similer_overlap(const int8_t *lseq, const int64_t llen, const int8_t *rseq, const int64_t rlen, const int64_t min_ol)
{
    int64_t i;
    int64_t j;
    int64_t max_ol;
    int64_t tol_mis;
    int64_t n_mis;

    max_ol = min_64(llen, rlen);

    for (i = max_ol; i >= min_ol; --i) {
        tol_mis = i * 0.1 + 0.5;
        n_mis = 0;
        for (j = 0; j < i; ++j) {
            if (lseq[llen - j - 1] != rseq[i - j - 1]) {
                ++n_mis;
                if (n_mis > tol_mis)
                    break;
            }
        }
        if (j == i)
            return i;
    }
    return 0;
}

int64_t qassembler_close_gap(const qassembler_t *qa, const lseq_t *sca, gap_t *gap, const int64_t min_ol, const double max_mis_rate)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t clen;
    int64_t l_ol;
    int64_t r_ol;
    int64_t max_ol;
    int64_t max_l_ol;
    int64_t max_r_ol;
    lstra_t *sp;
    lstra_t *max_sp = NULL;
    int64_t tol_mis;
    int64_t n_mis;
    double min_rate;
    double l_rate;
    double r_rate;

    max_l_ol = max_r_ol = 0;
    min_rate = 1.0;
    for (i = 0; i < qa->lg.stable_size; ++i) {
        if (lstra_is_emp(qa->lg.lsval_table[i]) || lstra_is_del(qa->lg.lsval_table[i]))
            continue;
        sp = qa->lg.lsval_table[i];
        clen = sp->len + qa->lg.l_len - 1;

        max_ol = min_64(gap->st, clen);
        l_rate = 1.0;
        l_ol = 0;
        for (j = max_ol; j >= min_ol; --j) {
            tol_mis = j * max_mis_rate + 0.5;
            n_mis = 0;
            for (k = 0; k < min_ol; ++k) {
                if (sca->seq[gap->st - k - 1] == 4)
                    break;
                else if (sca->seq[gap->st - k - 1] != bseq_get(sp->seq, j - k - 1)) {
                    ++n_mis;
                    if (n_mis > tol_mis)
                        break;
                }
            }
            if (k < min_ol)
                continue;

            for (; k < j; ++k) {
                if (sca->seq[gap->st - k - 1] == 4) {
                    --k;
                    tol_mis = k * max_mis_rate + 0.5;
                    break;
                }
                else if (sca->seq[gap->st - k - 1] != bseq_get(sp->seq, j - k - 1)) {
                    ++n_mis;
                    if (n_mis > tol_mis)
                        break;
                }
            }
            if (n_mis <= tol_mis && (double)n_mis / k < l_rate) {
                l_rate = (double)n_mis / k;
                l_ol = j;
            }
        }
        if (l_ol == 0)
            continue;

        max_ol = min_64(clen, sca->len - gap->ed);
        r_rate = 1.0;
        r_ol = 0;
        for (j = max_ol; j >= min_ol; --j) {
            tol_mis = j * max_mis_rate + 0.5;
            n_mis = 0;
            for (k = 0; k < min_ol; ++k) {
                if (sca->seq[gap->ed + j - k - 1] == 4)
                    break;
                else if (bseq_get(sp->seq, clen - k - 1) != sca->seq[gap->ed + j - k - 1]) {
                    ++n_mis;
                    if (n_mis > tol_mis)
                        break;
                }
            }
            if (k < min_ol)
                continue;

            for (; k < j; ++k) {
                if (sca->seq[gap->ed + j - k - 1] == 4) {
                    --k;
                    tol_mis = k * max_mis_rate + 0.5;
                    break;
                }
                else if (bseq_get(sp->seq, clen - k - 1) != sca->seq[gap->ed + j - k - 1]) {
                    ++n_mis;
                    if (n_mis > tol_mis)
                        break;
                }
            }
            if (n_mis <= tol_mis && (double)n_mis / j < r_rate) {
                r_rate = (double)n_mis / j;
                r_ol = j;
            }
        }
        if (r_ol == 0)
            continue;

        if ((l_rate * l_ol + r_rate * r_ol) / (l_ol + r_ol) < min_rate) {
            max_l_ol = l_ol;
            max_r_ol = r_ol;
            min_rate = (l_rate * l_ol + r_rate * r_ol) / (l_ol + r_ol);
            max_sp = sp;
        }
    }

    if (max_l_ol == 0)
        return 0;

    clen = max_sp->len + qa->lg.l_len - 1;
    if (clen - max_l_ol - max_r_ol > 0) {
        gap->len = clen - max_l_ol - max_r_ol;
        gap->seq = (int8_t *)my_realloc(gap->seq, gap->len * sizeof(int8_t));
        for (i = 0; i < gap->len; ++i)
            gap->seq[i] = bseq_get(max_sp->seq, max_l_ol + i);
    }
    else {
        i = max_l_ol + max_r_ol - clen;
        for (j = 0; j < i; ++j) {
            if (sca->seq[gap->st - j - 1] != sca->seq[gap->ed + i - j - 1] ||
                sca->seq[gap->st - j - 1] == 4 ||
                sca->seq[gap->ed + i - j - 1] == 4)
                break;
        }
        if (j != i)
            return 0;
        my_free(gap->seq);
        gap->seq = NULL;
        gap->len = clen - max_l_ol - max_r_ol - 1;
    }

    return 1;
}

void qassembler_init(qassembler_t *qa, const int64_t min_k, const int64_t max_k, const int64_t min_cov)
{
    memset(qa, 0, sizeof(qassembler_t));

    mcounter_init(&(qa->mc), 0);
    mgraph_init(&(qa->mg), QASSEMBLER_BRANCH_TH, QASSEMBLER_BUBBLE_TH);
    lcounter_init(&(qa->lc), 0);
    lgraph_init(&(qa->lg), QASSEMBLER_BRANCH_TH, QASSEMBLER_BUBBLE_TH);

    qa->min_k = min_k;
    qa->max_k = max_k;
    qa->min_cov = min_cov;
}

void qassembler_destroy(qassembler_t *qa)
{
    mcounter_destroy(&(qa->mc));
    mgraph_destroy(&(qa->mg));
    lcounter_destroy(&(qa->lc));
    lgraph_destroy(&(qa->lg));
}

void qassembler_assemble(qassembler_t *qa, const int8_t *seq, const int64_t len)
{
    int64_t i;
    uint64_t n_delete;

    mcounter_count_quick(&(qa->mc), (uint8_t *)seq, len, qa->min_k);
    mgraph_make_graph_quick(&(qa->mg), &(qa->mc), qa->min_cov);
    n_delete = 0;
    do {
        n_delete = mgraph_cut_branch(&(qa->mg));
        mgraph_join_nodes(&(qa->mg));
    } while (n_delete > 0);

    lcounter_count_quick(&(qa->lc), (uint8_t *)seq, len, qa->max_k);
    lgraph_make_graph_quick(&(qa->lg), &(qa->lc), qa->min_cov);
    n_delete = 0;
    do {
        n_delete = lgraph_cut_branch(&(qa->lg));
        lgraph_join_nodes(&(qa->lg));
    } while (n_delete > 0);

    merge_mgraph_with_lgraph(&(qa->mg), &(qa->lg), &(qa->lc));

    for (i = 0; i < qa->mg.stable_size; ++i) {
        if (mstra_is_emp(qa->mg.lstable[i]) || mstra_is_del(qa->mg.lstable[i]))
            continue;
        my_free(qa->mg.lstable[i]);
        qa->mg.lstable[i] = NULL;
    }

    for (i = 0; i < qa->lg.stable_size; ++i) {
        if (lstra_is_emp(qa->lg.lsval_table[i]) || lstra_is_del(qa->lg.lsval_table[i]))
            continue;
        my_free(qa->lg.lsval_table[i]);
        qa->lg.lsval_table[i] = NULL;
    }

    lgraph_make_graph_quick(&(qa->lg), &(qa->lc), qa->min_cov);
    n_delete = 0;
    do {
        n_delete = lgraph_cut_branch(&(qa->lg));
        lgraph_join_nodes(&(qa->lg));
    } while (n_delete > 0);
    n_delete = 0;
    do {
        n_delete = lgraph_crush_bubble(&(qa->lg), (double)UINT16_MAX, NULL);
        lgraph_join_nodes(&(qa->lg));
    } while (n_delete > 0);
}

void qassembler_init_lgraph(qassembler_t *qa)
{
    int64_t i;

    for (i = 0; i < qa->lg.stable_size; ++i) {
        if (lstra_is_emp(qa->lg.lsval_table[i]) || lstra_is_del(qa->lg.lsval_table[i]))
            continue;
        my_free(qa->lg.lsval_table[i]);
        qa->lg.lsval_table[i] = NULL;
    }
}

void merge_mgraph_with_lgraph(const mgraph_t *mg, const lgraph_t *lg, lcounter_t *lc)
{
    uint64_t i;
    uint64_t j;
    uint64_t n_key;
    uint64_t l_len;
    uint64_t bl;
    uint64_t index;
    uint64_t step;
    lmer_t key;
    mstra_t *msp;
    lstra_t *lsp;

    l_len = lc->l_len;
    bl = lc->bl;

    n_key = 0;
    for (i = 0; i < mg->stable_size; ++i) {
        msp = mg->lstable[i];
        if (mstra_is_emp(msp) || mstra_is_del(msp) || msp->len + mg->k_len <= lg->l_len)
            continue;
        n_key += msp->len + mg->k_len - lg->l_len;
    }
    for (i = 0; i < lg->stable_size; ++i) {
        lsp = lg->lsval_table[i];
        if (lstra_is_emp(lsp) || lstra_is_del(lsp))
            continue;
        n_key += lsp->len;
    }

    if ((double)lc->table_size * MAX_LOAD_FACTOR <= n_key) {
        lc->table_size = 2;
        lc->index_len = 1;
        while ((double)lc->table_size * MAX_LOAD_FACTOR <= n_key) {
            lc->table_size *= 2;
            ++lc->index_len;
        }
        lc->key_table = my_realloc(lc->key_table, lc->table_size * bl * sizeof(uint64_t));
        lc->val_table = my_realloc(lc->val_table, lc->table_size*sizeof(uint16_t));
    }

    memset(lc->val_table, 0, lc->table_size*sizeof(uint16_t));

    for (i = 0; i < mg->stable_size; ++i) {
        msp = mg->lstable[i];
        if (mstra_is_emp(msp) || mstra_is_del(msp))
            continue;

        if (msp->len + mg->k_len - 1 < l_len)
            continue;
        for (j = 0; j < l_len - 1; ++j)
            lmer_set(key, j+1, bseq_get(msp->seq, j));
        for (j = 0; j < msp->len + mg->k_len - l_len; ++j) {
            lmer_lshift(key, l_len);
            lmer_set(key, l_len - 1, bseq_get(msp->seq, j + l_len - 1));
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
            if (!lc->val_table[index]) {
                lc->val_table[index] = msp->cov;
                lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
            }
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                lc->val_table[index] = msp->cov;
            }
            else {
                step = lmer_rehash(key, l_len, lc->index_len);
                index = (index + step) & (lc->table_size-1);
                while (lc->val_table[index])
                    index = (index+step) & (lc->table_size-1);
                lc->val_table[index] = msp->cov;
                lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
            }
        }
    }

    for (i = 0; i < lg->stable_size; ++i) {
        lsp = lg->lsval_table[i];
        if (lstra_is_emp(lsp) || lstra_is_del(lsp))
            continue;

        for (j = 0; j < l_len - 1; ++j)
            lmer_set(key, j+1, bseq_get(lsp->seq, j));
        for (j = 0; j < lsp->len; ++j) {
            lmer_lshift(key, l_len);
            lmer_set(key, l_len - 1, bseq_get(lsp->seq, j + l_len - 1));
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
            if (!lc->val_table[index]) {
                lc->val_table[index] = lsp->cov;
                lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
            }
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) != 0) {
                step = lmer_rehash(key, l_len, lc->index_len);
                index = (index + step) & (lc->table_size-1);

                while (lc->val_table[index] && lmer_cmp(&(lc->key_table[index * bl]), key, l_len) != 0)
                    index = (index+step) & (lc->table_size-1);
                if (!(lc->val_table[index])) {
                    lc->val_table[index] = lsp->cov;
                    lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
                }
            }
        }
    }
}

int cmp_gap_ptr(const void *a, const void *b)
{
    if ((*(gap_t **)a)->len < (*(gap_t **)b)->len)
        return -1;
    else if ((*(gap_t **)a)->len > (*(gap_t **)b)->len)
        return 1;
    else
        return 0;
}
