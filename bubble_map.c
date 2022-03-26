#include "bubble_map.h"

void print_bubble_map_usage(void)
{
    option_bubble_map_t def = {};

    option_bubble_map_init(&def);

    fprintf(stderr, "\nUsage: platanus bubble_map [Options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR              : prefix of output file (default out, length <= %u)\n", MAX_FILE_LEN);
    fprintf(stderr, "    -b FILE1 [FILE2 ...]: bubble_file (fasta format)\n");
    fprintf(stderr, "    -c FILE1 [FILE2 ...]: contig(scaffold)_file (fasta format)\n");
    fprintf(stderr, "    -k INT              : maximum k of contig assembly (default %ld)\n", def.k);
    fprintf(stderr, "    -u FLOAT            : maximum difference for bubble crush (identity, default %.1f)\n", def.u);
    fprintf(stderr, "    -m INT              : memory limit (GB, >= 1, default %llu)\n", def.m / GIBIBYTE);
    fprintf(stderr, "    -t INT              : number of threads (<= %u, default 1)\n", MAX_THREAD);

    fprintf(stderr, "Outputs:\n");
    fprintf(stderr, "    PREFIX_mutation.tsv\n");
    fprintf(stderr, "    PREFIX_bubble.tsv\n");

    option_bubble_map_destroy(&def);
}

void option_bubble_map_init(option_bubble_map_t *opt)
{
    strcpy(opt->o, "out");
    opt->k = 50;
    opt->u = 0.1;
    opt->m = 4 * GIBIBYTE;
    opt->t = 1;
    opt->n_c = 0;
    opt->n_b = 0;
}

void option_bubble_map_destroy(const option_bubble_map_t *opt)
{
    uint64_t i;

    for (i = 0; i < opt->n_c; ++i)
        fclose(opt->c[i]);
    for (i = 0; i < opt->n_b; ++i)
        fclose(opt->b[i]);
}

uint8_t option_bubble_map_parse(option_bubble_map_t *opt, int argc, char **argv)
{
    uint64_t i;

    for (i = 2; i < argc;) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help") || !strcmp(argv[i], "--help"))
            return -1;
        else if (!strcmp(argv[i], "-o"))
            i = option_string(argc, argv, i, MAX_FILE_LEN, opt->o);
        else if (!strcmp(argv[i], "-c"))
            i = option_multi_fasta_file(argc, argv, i, opt->c, &(opt->n_c));
        else if (!strcmp(argv[i], "-b"))
            i = option_multi_fasta_file(argc, argv, i, opt->b, &(opt->n_b));
        else if (!strcmp(argv[i], "-k"))
            i = option_int(argc, argv, i, 1, INT64_MAX, &(opt->k));
        else if (!strcmp(argv[i], "-u"))
            i = option_float(argc, argv, i, 0.0, 1.0, &(opt->u));
        else if (!strcmp(argv[i], "-t"))
            i = option_int(argc, argv, i, 1, MAX_THREAD, &(opt->t));
        else if (!strcmp(argv[i], "-m")) {
            i = option_int(argc, argv, i, 1, INT64_MAX, &(opt->m));
            opt->m *= GIBIBYTE;
        }
        else {
            fprintf(stderr, "error: wrong option \"%s\"\n\n", argv[i]);
            return -1;
        }
    }

    if (opt->n_c == 0) {
        fputs("error: no contig file\n", stderr);
        return -1;
    }

    if (opt->n_b == 0) {
        fputs("error: no bubble file\n", stderr);
        return -1;
    }
    return 0;
}

void bubble_map(const option_bubble_map_t *opt)
{
    int64_t i;
    char out_name[MAX_FILE_LEN + MAX_EXTENSION_NAME + 1];
    const int64_t n_thread = opt->t;
    const int64_t min_match_len = opt->k - 1;
    const double bubble_th = opt->u;
    int64_t key_len = min_64(min_match_len, 32);
    contig_t con;
    contig_t bub;
    mapper_t mp;
    FILE *pos_fp[MAX_THREAD];
    name_list_t cname;
    name_list_t bname;

    if (key_len > 32)
        key_len = 32;

    contig_init(&con);
    contig_init(&bub);

    mapper_init(&mp, min_match_len, key_len, opt->m);

    omp_set_num_threads(n_thread);

#    pragma omp parallel sections private(i)
    {
#        pragma omp section
        {
            for (i = 0; i < opt->n_c; ++i)
                contig_read_fasta(&con, opt->c[i]);
            contig_load_seq(&con);
            mapper_set_contig(&mp, &con);
            mapper_make_table(&mp);
        }

#        pragma omp section
        {
            for (i = 0; i < opt->n_b; ++i)
                contig_read_fasta(&bub, opt->b[i]);
            contig_load_seq(&bub);
        }
    }
    
    mapper_map_bubble(&mp, &bub, min_match_len, pos_fp, n_thread);

    mapper_destroy(&mp);
    contig_destroy(&con);
    contig_destroy(&bub);

#    pragma omp parallel sections private(i)
    {
#        pragma omp section
        {
            name_list_init(&cname);
            for (i = 0; i < opt->n_c; ++i)
                name_list_read_fasta(&cname, opt->c[i]);
            name_list_load(&cname);
        }

#        pragma omp section
        {
            name_list_init(&bname);
            for (i = 0; i < opt->n_b; ++i)
                name_list_read_fasta(&bname, opt->b[i]);
            name_list_load(&bname);
        }
    }

    strcpy(out_name, opt->o);
    strcat(out_name, "_mutation.tsv");
    detect_mutation(pos_fp, &cname, &bname, min_match_len, bubble_th, n_thread, out_name);

    strcpy(out_name, opt->o);
    strcat(out_name, "_bubble.tsv");
    print_bubble(pos_fp, &cname, &bname, n_thread, out_name);

    name_list_destroy(&cname);
    name_list_destroy(&bname);
    for (i = 0; i < n_thread; ++i)
        fclose(pos_fp[i]);
    
    fputs("bubble_map completed!\n", stderr);
}

void print_bubble(FILE **pos_fp, const name_list_t *cname, const name_list_t *bname, const int64_t n_thread, char *out_name)
{
    int64_t i;
    int64_t t;
    uint64_t tmp;
    int32_t bub_id;
    bool is_mapped;
    position_t pos;
    varr8_t bub;
    varr8_t con;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    fputs("#contig_name\tcontig_pos\tbubble_name\tcontig_seq\tbubble_seq\n", out);

    varr8_init(&bub);
    varr8_init(&con);

    for (t = 0; t < n_thread; ++t)
        rewind(pos_fp[t]);

    t = bub_id = 0;
    while (fread(&is_mapped, sizeof(bool), 1, pos_fp[t])) {
        ++bub_id;
        if (is_mapped) {
            fread(&tmp, sizeof(uint64_t), 1, pos_fp[t]);
            varr8_resize(&bub, tmp);
            fread(bub.ele, sizeof(uint8_t), bub.size, pos_fp[t]);
            fread(&tmp, sizeof(uint64_t), 1, pos_fp[t]);
            varr8_resize(&con, tmp);
            fread(con.ele, sizeof(uint8_t), con.size, pos_fp[t]);
            fread(&pos, sizeof(position_t), 1, pos_fp[t]);

            fprintf(out, "%s\t%d\t%s\t", cname->name[pos.id - 1], pos.ofst, bname->name[bub_id - 1]);
            if (con.size + bub.size > 0) {
                if (con.size > 0) {
                    for (i = 0; i < con.size; ++i)
                        putc(num2ascii(con.ele[i]), out);
                }
                else
                    putc('-', out);
                putc('\t', out);

                if (bub.size > 0) {
                    for (i = 0; i < bub.size; ++i)
                        putc(num2ascii(bub.ele[i]), out);
                }
                else
                    putc('-', out);
                putc('\n', out);
            }
            else
                fputs("perfect\tperfect\n", out);
        }
        else
            fprintf(out, "none\tnone\t%s\tnone\tnone\n", bname->name[bub_id - 1]);

        t = (t + 1) % n_thread;
    }
}

void detect_mutation(FILE **pos_fp, const name_list_t *cname, const name_list_t *bname, const int64_t min_match_len, const double bubble_th, const int64_t n_thread, char *out_name)
{
    int64_t i;
    int64_t j;
    int64_t t;
    int64_t tmp;
    int64_t n_mut = 0;
    int32_t bub_id;
    bool is_mapped;
    position_t pos;
    varr8_t bub;
    varr8_t con;
    varr8_t trace;
    varr64_t score;
    bubble_mutation_t mut;
    bubble_mutation_t **buf;
    FILE *mut_fp[MAX_THREAD];
    FILE *tmp_fp[MAX_THREAD];
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    omp_set_num_threads(n_thread);

    fputs("#contig_name\tcontig_pos\tbubble_name\tcontig_seq\tbubble_seq\n", out);

#    pragma omp parallel for schedule(static, 1) private(tmp, bub_id, is_mapped, pos, bub, con, trace, score) reduction(+: n_mut)
    for (t = 0; t < n_thread; ++t) {
        varr8_init(&bub);
        varr8_init(&con);
        varr8_init(&trace);
        varr64_init(&score);

        bub_id = t + 1;
        mut_fp[t] = tmpfile_open();
        tmp_fp[t] = tmpfile_open();
        rewind(pos_fp[t]);
        while (fread(&is_mapped, sizeof(bool), 1, pos_fp[t])) {
            bub_id += n_thread;
            if (!is_mapped) {
                fwrite(&is_mapped, sizeof(bool), 1, tmp_fp[t]);
                continue;
            }
            fread(&tmp, sizeof(int64_t), 1, pos_fp[t]);
            varr8_resize(&bub, tmp);
            fread(bub.ele, sizeof(uint8_t), bub.size, pos_fp[t]);
            fread(&tmp, sizeof(int64_t), 1, pos_fp[t]);
            varr8_resize(&con, tmp);
            fread(con.ele, sizeof(uint8_t), con.size, pos_fp[t]);
            fread(&pos, sizeof(position_t), 1, pos_fp[t]);

            if (abs64((int64_t)bub.size - (int64_t)con.size) > (max_64(bub.size, con.size) + 2 *  min_match_len) * bubble_th
                || align_bubble(&bub, &con, &score, &trace) > (max_64(bub.size, con.size) + 2 *  min_match_len) * bubble_th) {
                is_mapped = false;
                fwrite(&is_mapped, sizeof(bool), 1, tmp_fp[t]);
                continue;
            }

            n_mut += trace_bubble_alignment(&bub, &con, &trace, bub_id - n_thread, &pos, &(mut_fp[t]));

            fwrite(&is_mapped, sizeof(bool), 1, tmp_fp[t]);
            fwrite(&(bub.size), sizeof(int64_t), 1, tmp_fp[t]);
            fwrite(bub.ele, sizeof(uint8_t), bub.size, tmp_fp[t]);
            fwrite(&(con.size), sizeof(int64_t), 1, tmp_fp[t]);
            fwrite(con.ele, sizeof(uint8_t), con.size, tmp_fp[t]);
            fwrite(&pos, sizeof(position_t), 1, tmp_fp[t]);
        }

        varr8_destroy(&bub);
        varr8_destroy(&con);
        varr8_destroy(&trace);
        varr64_destroy(&score);

        fclose(pos_fp[t]);
        pos_fp[t] = tmp_fp[t];
    }

    buf = (bubble_mutation_t **)my_malloc(n_mut * sizeof(bubble_mutation_t *));

    n_mut = 0;
    for (t = 0; t < n_thread; ++t) {
        rewind(mut_fp[t]);
        while (fread(&mut, sizeof(bubble_mutation_t), 1, mut_fp[t])) {
            if (mut.type == BUB_MIS) {
                buf[n_mut] = (bubble_mutation_t *)my_malloc(sizeof(bubble_mutation_t) + sizeof(uint8_t) * mut.len * 2);
                *(buf[n_mut]) = mut;
                fread(buf[n_mut]->seq, sizeof(uint8_t), mut.len, mut_fp[t]);
                fread(&(buf[n_mut]->seq[mut.len]), sizeof(uint8_t), mut.len, mut_fp[t]);
            }
            else if (mut.type == BUB_INS) {
                buf[n_mut] = (bubble_mutation_t *)my_malloc(sizeof(bubble_mutation_t) + sizeof(uint8_t) * mut.len);
                *(buf[n_mut]) = mut;
                fread(buf[n_mut]->seq, sizeof(uint8_t), mut.len, mut_fp[t]);
            }
            else {
                buf[n_mut] = (bubble_mutation_t *)my_malloc(sizeof(bubble_mutation_t) + sizeof(uint8_t) * mut.len);
                *(buf[n_mut]) = mut;
                fread(buf[n_mut]->seq, sizeof(uint8_t), mut.len, mut_fp[t]);
            }
            ++n_mut;
        }
        fclose(mut_fp[t]);
    }

    qsort(buf, n_mut, sizeof(bubble_mutation_t *), cmp_bubble_mutation_p);

    for (i = 0; i < n_mut; ++i) {
        fprintf(out, "%s\t%d\t%s\t", cname->name[buf[i]->con_id - 1], buf[i]->ofst, bname->name[buf[i]->bub_id - 1]);
        if (buf[i]->type == BUB_MIS) {
            for (j = 0; j < buf[i]->len; ++j)
                putc(num2ascii(buf[i]->seq[j]), out);
            putc('\t', out);
            for (; j < 2 * buf[i]->len; ++j)
                putc(num2ascii(buf[i]->seq[j]), out);
        }
        else if (buf[i]->type == BUB_INS) {
            fprintf(out, "-\t");
            for (j = 0; j < buf[i]->len; ++j)
                putc(num2ascii(buf[i]->seq[j]), out);
        }
        else {
            for (j = 0; j < buf[i]->len; ++j)
                putc(num2ascii(buf[i]->seq[j]), out);
            fprintf(out, "\t-");
        }
        putc('\n', out);

        my_free(buf[i]);
    }

    my_free(buf);
    fclose(out);
}

int64_t align_bubble(const varr8_t *bub, const varr8_t *con, varr64_t *score, varr8_t *trace)
{
    int64_t m;
    int64_t n;
    int64_t n_col;
    
    n_col = con->size + 1;
    varr64_resize(score, n_col * 2);
    varr8_resize(trace, ((bub->size + 1) * (con->size + 1) + 3) / 4);

    score->ele[0] = 0;
    bseq_set(trace->ele, 0, 0);
    for (n = 0; n < con->size; ++n) {
        score->ele[n + 1] = n + 1;
        bseq_set(trace->ele, n + 1, 3);
    }
    for (m = 0; m < bub->size; ++m) {
        score->ele[(m&1)*n_col] = m;
        score->ele[(~m&1)*n_col] = m+1;
        bseq_set(trace->ele, (m + 1) * n_col, 2);
        for (n = 0; n < con->size; ++n) {
            if (bub->ele[m] == con->ele[n]) {
                score->ele[(~m&1)*n_col + n+1] = score->ele[(m&1)*n_col + n];
                bseq_set(trace->ele, (m + 1) * n_col + (n + 1), 0);
                continue;
            }
            score->ele[(~m&1)*n_col + n+1] = score->ele[(m&1)*n_col + n] + 1;
            bseq_set(trace->ele, (m + 1) * n_col + (n + 1), 1);
            if (score->ele[(~m&1)*n_col + n+1] >= score->ele[(m&1)*n_col + n+1] + 1) {
                if (score->ele[(~m&1)*n_col + n+1] > score->ele[(m&1)*n_col + n+1] + 1
                    || bseq_get(trace->ele, m * n_col + (n + 1)) == 2) {

                    score->ele[(~m&1)*n_col + n+1] = score->ele[(m&1)*n_col + n+1] + 1;
                    bseq_set(trace->ele, (m + 1) * n_col + (n + 1), 2);
                }
            }
            if (score->ele[(~m&1)*n_col + n+1] >= score->ele[(~m&1)*n_col + n] + 1) {
                if (score->ele[(~m&1)*n_col + n+1] > score->ele[(~m&1)*n_col + n] + 1
                    || bseq_get(trace->ele, (m + 1) * n_col + n) == 3) {

                    score->ele[(~m&1)*n_col + n+1] = score->ele[(~m&1)*n_col + n] + 1;
                    bseq_set(trace->ele, (m + 1) * n_col + (n + 1), 3);
                }
            }
        }
    }

    return score->ele[(m&1)*n_col + n];
}

int64_t trace_bubble_alignment(const varr8_t *bub, const varr8_t *con, const varr8_t *trace, const int32_t bub_id, const position_t *pos, FILE **mut_fpp)
{
    int64_t m;
    int64_t n;
    int64_t n_col;
    int64_t tmp;
    int64_t n_mut;
    bubble_mutation_t mut;

    n_mut = 0;
    mut.bub_id = bub_id;
    mut.con_id = pos->id;
    mut.ofst = pos->ofst;
    n_col = con->size + 1;

    m = bub->size;
    n = con->size;
    while (m > 0 || n > 0) {
        switch (bseq_get(trace->ele, m * n_col + n)) {
            case 0:
                --m;
                --n;
                break;
            case 1:
                tmp = m;
                --m;
                --n;
                while (bseq_get(trace->ele, m * n_col + n) == 1) {
                    --m;
                    --n;
                }
                mut.len = tmp - m;
                mut.type = BUB_MIS;
                mut.ofst = pos->ofst + n;
                fwrite(&mut, sizeof(bubble_mutation_t), 1, *mut_fpp);
                fwrite(&(con->ele[n]), sizeof(uint8_t), mut.len, *mut_fpp);
                fwrite(&(bub->ele[m]), sizeof(uint8_t), mut.len, *mut_fpp);
                ++n_mut;
                break;
            case 2:
                tmp = m;
                --m;
                while (bseq_get(trace->ele, m * n_col + n) == 2)
                    --m;
                mut.len = tmp - m;
                mut.type = BUB_INS;
                mut.ofst = pos->ofst + n;
                fwrite(&mut, sizeof(bubble_mutation_t), 1, *mut_fpp);
                fwrite(&(bub->ele[m]), sizeof(uint8_t), mut.len, *mut_fpp);
                ++n_mut;
                break;
            default:
                tmp = n;
                --n;
                while (bseq_get(trace->ele, m * n_col + n) == 3)
                    --n;
                mut.len = tmp - n;
                mut.type = BUB_DEL;
                mut.ofst = pos->ofst + n;
                fwrite(&mut, sizeof(bubble_mutation_t), 1, *mut_fpp);
                fwrite(&(con->ele[n]), sizeof(uint8_t), mut.len, *mut_fpp);
                ++n_mut;
        }
    }

    return n_mut;
}

int cmp_bubble_mutation_p(const void *a, const void *b)
{
    if ((*((bubble_mutation_t **)a))->con_id != (*((bubble_mutation_t **)b))->con_id)
        return (*((bubble_mutation_t **)a))->con_id - (*((bubble_mutation_t **)b))->con_id;
    else if ((*((bubble_mutation_t **)a))->ofst != (*((bubble_mutation_t **)b))->ofst)
        return (*((bubble_mutation_t **)a))->ofst - (*((bubble_mutation_t **)b))->ofst;
    else if ((*((bubble_mutation_t **)a))->bub_id != (*((bubble_mutation_t **)b))->bub_id)
        return (*((bubble_mutation_t **)a))->bub_id - (*((bubble_mutation_t **)b))->bub_id;
    else if ((*((bubble_mutation_t **)a))->type == BUB_INS)
        return -1;
    else if ((*((bubble_mutation_t **)b))->type == BUB_INS)
        return 1;
    else
        return 0;
}

void name_list_init(name_list_t *list)
{
    memset(list, 0, sizeof(name_list_t));
}

void name_list_destroy(name_list_t *list)
{
    if (list->name_fp != NULL)
        fclose(list->name_fp);
    if (list->pool != NULL) {
        my_free(list->pool);
        my_free(list->name);
    }
    memset(list, 0, sizeof(name_list_t));
}

void name_list_read_fasta(name_list_t *list, FILE *fp)
{
    int8_t c;

    if (list->name_fp == NULL) {
        list->name_fp = tmpfile_open();
        check_tmpfile(list->name_fp);
    }
    rewind(fp);
    fseek(list->name_fp, 0, SEEK_END);
    while ((c = getc(fp)) != '>') {
        if (c == EOF) {
            fputs("error: invalid file format!\n", stderr);
            my_exit(1);
        }
    }
    while (c != EOF) {
        while ((c = getc(fp)) != '\n') {
            putc(c, list->name_fp);
            ++(list->pool_size);
        }
        while ((c = getc(fp)) != '>' && c != EOF)
            ;
        putc('\0', list->name_fp);
        ++(list->pool_size);;
    }
}

void name_list_load(name_list_t *list)
{
    int64_t i;
    int64_t j;

    list->pool = (int8_t *)my_malloc(list->pool_size * sizeof(int8_t));
    
    rewind(list->name_fp);
    fread(list->pool, sizeof(int8_t), list->pool_size, list->name_fp);

    list->n_name = 0;
    for (i = 0; i < list->pool_size; ++i) {
        if (list->pool[i] == '\0')
            ++(list->n_name);
    }

    list->name = (int8_t **)my_malloc(list->n_name * sizeof(int8_t *));

    i = j = 0;
    while (i < list->pool_size) {
        list->name[j] = &(list->pool[i]);
        ++j;
        while (list->pool[i] != '\0')
            ++i;
        ++i;
    }

    fclose(list->name_fp);
    list->name_fp = NULL;
}
