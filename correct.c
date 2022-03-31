#include "correct.h"

static pthread_mutex_t *g_ctt_table_mutex;
static FILE *g_ctt_unstored_fp[MAX_THREAD];
static FILE *g_ctt_stored_fp[MAX_THREAD];

void print_correct_usage(void)
{
    option_correct_t def;

    option_correct_init(&def);

    fprintf(stderr, "Usage: platanus correct [Options]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o STR                       : suffix of output file (def %s, length <= %u)\n", ".corrected", MAX_FILE_LEN);
    fprintf(stderr, "  -count_fa FILE1 [FILE2 ...]  : reads file for k-mer count (fasta, number <= %u)\n", MAX_FILE_NUM);
    fprintf(stderr, "  -count_fq FILE1 [FILE2 ...]  : reads file for k-mer count (fastq, number <= %u)\n", MAX_FILE_NUM);
    fprintf(stderr, "  -target_fa FILE1 [FILE2 ...] : target reads file (fasta, number <= %u)\n", MAX_FILE_NUM);
    fprintf(stderr, "  -target_fq FILE1 [FILE2 ...] : target reads file (fastq, number <= %u)\n", MAX_FILE_NUM);
    fprintf(stderr, "  -k INT [INT ...]             : k values (<= %u, def %lu)\n", MAX_CRR_K, def.k[0]);
    fprintf(stderr, "  -c INT                       : lower occurrence number threshold (-1 means auto, def %ld)\n", def.c);
    fprintf(stderr, "  -d INT                       : upper occurrence number threshold (-1 means auto, def %ld)\n", def.d);
    fprintf(stderr, "  -e FLOAT                     : max_edit_distance / read_length (<= 1, def %.2f)\n", def.e);
    fprintf(stderr, "  -t INT                       : number of threads (<= %u, def %lu)\n", MAX_THREAD, def.t);
    fprintf(stderr, "  -m INT                       : memory limit (GB, >= 1, def %llu)\n", def.m / GIBIBYTE);
    fprintf(stderr, "  -no_ungap                    : not perform ungap correction (def false)\n");
    fprintf(stderr, "  -no_gapped                   : not perform gapped correction (def false)\n");
    fprintf(stderr, "  -h, -help, --help            : print usage\n");

    fprintf(stderr, "\nOutputs:\n");
    fprintf(stderr, "    INPUT_FILE1.corrected [INPUT_FILE2.corrected ...] \n");

    option_correct_destroy(&def);
}

void option_correct_init(option_correct_t *opt)
{
    memset(opt, 0, sizeof(option_correct_t));

    strcpy(opt->o, ".cor");

	opt->n_count_fa = 0;
	opt->n_count_fq = 0;
	opt->n_target_fa = 0;
	opt->n_target_fq = 0;

    opt->n_k = 1;
    opt->k[0] = 21;

    opt->c = -1;
    opt->d = INT64_MAX;
    opt->e = 0.03;
    opt->t = 1;
    opt->m = 4 * GIBIBYTE;

    opt->no_ungap = false;
    opt->no_gapped = false;
}

void option_correct_destroy(option_correct_t *opt)
{
    uint64_t i;

    for (i = 0; i < opt->n_count_fa; ++i) {
        fclose(opt->count_fa_fp[i]);
        my_free(opt->count_fa_name[i]);
    }
    for (i = 0; i < opt->n_count_fq; ++i) {
        fclose(opt->count_fq_fp[i]);
        my_free(opt->count_fq_name[i]);
    }
    for (i = 0; i < opt->n_target_fa; ++i) {
        fclose(opt->target_fa_fp[i]);
        my_free(opt->target_fa_name[i]);
    }
    for (i = 0; i < opt->n_target_fq; ++i) {
        fclose(opt->target_fq_fp[i]);
        my_free(opt->target_fq_name[i]);
    }
}

int option_correct_parse(option_correct_t *opt, int argc, char **argv)
{
    uint64_t i;
	bool is_k_set = false;

    for (i = 2; i < argc;) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help") || !strcmp(argv[i], "--help"))
            return -1;
        else if (!strcmp(argv[i], "-k")) {
			if (!is_k_set) {
				opt->n_k = 0;
				is_k_set = true;
			}
			i = option_multi_int(argc, argv, i, 1, MAX_CRR_K, opt->k, &(opt->n_k));
		}
        else if (!strcmp(argv[i], "-o"))
            i = option_string(argc, argv, i, MAX_FILE_LEN, opt->o);
        else if (!strcmp(argv[i], "-count_fa"))
            i = option_multi_file_name(argc, argv, i, opt->count_fa_fp, opt->count_fa_name, &(opt->n_count_fa));
        else if (!strcmp(argv[i], "-count_fq"))
            i = option_multi_file_name(argc, argv, i, opt->count_fq_fp, opt->count_fq_name, &(opt->n_count_fq));
        else if (!strcmp(argv[i], "-target_fa"))
            i = option_multi_file_name(argc, argv, i, opt->target_fa_fp, opt->target_fa_name, &(opt->n_target_fa));
        else if (!strcmp(argv[i], "-target_fq"))
            i = option_multi_file_name(argc, argv, i, opt->target_fq_fp, opt->target_fq_name, &(opt->n_target_fq));
        else if (!strcmp(argv[i], "-c"))
            i = option_int(argc, argv, i, INT64_MIN, INT64_MAX, &(opt->c));
        else if (!strcmp(argv[i], "-d"))
            i = option_int(argc, argv, i, INT64_MIN, INT64_MAX, &(opt->d));
        else if (!strcmp(argv[i], "-e"))
            i = option_float(argc, argv, i, 0.0, 1.0, &(opt->e));
        else if (!strcmp(argv[i], "-t"))
            i = option_int(argc, argv, i, 1, MAX_THREAD, &(opt->t));
        else if (!strcmp(argv[i], "-m")) {
            i = option_int(argc, argv, i, 1, INT64_MAX, &(opt->m));
            opt->m *= GIBIBYTE;
        }
        else if (!strcmp(argv[i], "-no_ungap")) {
            opt->no_ungap = true;
			++i;
		}
        else if (!strcmp(argv[i], "-no_gapped")) {
            opt->no_gapped = true;
			++i;
		}
        else {
            fprintf(stderr, "error: wrong option \"%s\"\n\n", argv[i]);
            return -1;
        }
    }

    if ((opt->n_count_fa == 0 && opt->n_count_fq == 0) || (opt->n_target_fa == 0 && opt->n_target_fq == 0)) {
        fputs("error: no input file\n\n", stderr);
        return -1;
    }

	return 0;
}

void correct(option_correct_t *opt)
{
    uint64_t i;
    corrector_t ct;

    corrector_init(&ct, opt->m);

    ct.n_target_file = 0;
    for (i = 0; i < opt->n_target_fa; ++i) {
        strcpy(ct.fname[i], opt->target_fa_name[i]);
        strcat(ct.fname[i], opt->o);
        corrector_read_fasta(&ct, opt->target_fa_fp[i], true);
    }
    for (i = 0; i < opt->n_target_fq; ++i) {
        strcpy(ct.fname[opt->n_target_fa + i], opt->target_fq_name[i]);
        strcat(ct.fname[opt->n_target_fa + i], opt->o);
        corrector_read_fastq(&ct, opt->target_fq_fp[i], true);
    }
    for (i = 0; i < opt->n_count_fa; ++i) {
        corrector_read_fasta(&ct, opt->count_fa_fp[i], false);
    }
    for (i = 0; i < opt->n_count_fq; ++i) {
        corrector_read_fastq(&ct, opt->count_fq_fp[i], false);
    }


    for (i = 0; i < opt->n_k; ++i) {
        corrector_set_th(&ct, opt->c);
        corrector_set_upper_th(&ct, opt->d);
        corrector_set_max_edit(&ct, opt->e);

		if (!(opt->no_ungap)) {
			corrector_make_table(&ct, opt->k[i]);
			fprintf(stderr, "ungap correction...\neach dot below indicates %ld reads in a thread (thread ID, 0).\n", DISPLAY_NUM_READ_UNIT);
			corrector_ungap_correct(&ct);
			fprintf(stderr, "ungap correction: K = %lu, THRESHOLD = %lu, CORRECTED_READS = %lu, CORRECTED_BASES = %lu\n", opt->k[i], ct.th, ct.n_cor_read, ct.n_cor_base);
			corrector_destroy_table(&ct);
		}

		if (!(opt->no_gapped)) {
			corrector_make_table(&ct, opt->k[i]);
			fprintf(stderr, "gapped correction...\neach dot below indicates %ld reads in a thread (thread ID, 0).\n", DISPLAY_NUM_READ_UNIT);
			corrector_gapped_correct(&ct);
			fprintf(stderr, "gapped correction: K = %lu, THRESHOLD = %lu, CORRECTED = %lu, CORRECTED_BASES = %lu\n", opt->k[i], ct.th, ct.n_cor_read, ct.n_cor_base);
			corrector_destroy_table(&ct);
		}
    }

    corrector_show_seq(&ct);
    corrector_destroy(&ct);

    fputs("correction completed!\n", stderr);
}


void correct_mt(option_correct_t *opt)
{
    uint64_t i;
    corrector_threads_t ctt;

    corrector_threads_init(&ctt, opt->t, opt->m);

    for (i = 0; i < opt->n_target_fa; ++i) {
        strcpy(ctt.fname[i], opt->target_fa_name[i]);
        strcat(ctt.fname[i], opt->o);
        corrector_threads_read_fasta(&ctt, opt->target_fa_fp[i], true);
    }
    for (i = 0; i < opt->n_target_fq; ++i) {
        strcpy(ctt.fname[opt->n_target_fa + i], opt->target_fq_name[i]);
        strcat(ctt.fname[opt->n_target_fa + i], opt->o);
        corrector_threads_read_fastq(&ctt, opt->target_fq_fp[i], true);
    }
    for (i = 0; i < opt->n_count_fa; ++i) {
        corrector_threads_read_fasta(&ctt, opt->count_fa_fp[i], false);
    }
    for (i = 0; i < opt->n_count_fq; ++i) {
        corrector_threads_read_fastq(&ctt, opt->count_fq_fp[i], false);
    }

    for (i = 0; i < opt->n_k; ++i) {
        corrector_threads_set_th(&ctt, opt->c);
        corrector_threads_set_upper_th(&ctt, opt->d);
        corrector_threads_set_max_edit(&ctt, opt->e);

		if (!(opt->no_ungap)) {
			corrector_threads_make_table(&ctt, opt->k[i]);
			fprintf(stderr, "ungap correction...\neach dot below indicates %ld reads in a thread (thread ID, 0).\n", DISPLAY_NUM_READ_UNIT);
			corrector_threads_ungap_correct(&ctt);
			fprintf(stderr, "ungap correction: K = %lu, THRESHOLD = %lu, CORRECTED = %lu, CORRECTED_BASES = %lu\n", opt->k[i], ctt.th, ctt.n_cor_read, ctt.n_cor_base);
			corrector_threads_destroy_table(&ctt);
		}

		if (!(opt->no_gapped)) {
			corrector_threads_make_table(&ctt, opt->k[i]);
			fprintf(stderr, "gapped correction...\neach dot below indicates %ld reads in a thread (thread ID, 0).\n", DISPLAY_NUM_READ_UNIT);
			corrector_threads_gapped_correct(&ctt);
			fprintf(stderr, "gapped correction: K = %lu, THRESHOLD = %lu, CORRECTED = %lu, CORRECTED_BASES = %lu\n", opt->k[i], ctt.th, ctt.n_cor_read, ctt.n_cor_base);
			corrector_threads_destroy_table(&ctt);
		}
    }

    corrector_threads_show_seq(&ctt);
    corrector_threads_destroy(&ctt);

    fputs("correction completed!\n", stderr);
}

void corrector_init(corrector_t *ct, uint64_t mem_lim)
{
    memset(ct, 0, sizeof(corrector_t));
    ct->mem_lim = mem_lim;
}

void corrector_destroy_table(corrector_t *ct)
{
    if (ct->val_table != NULL)
        my_free(ct->val_table);
    ct->val_table = NULL;
    if (ct->key_table != NULL)
        my_free(ct->key_table);
    ct->key_table = NULL;
}

void corrector_destroy(corrector_t *ct)
{
    corrector_destroy_table(ct);
    if (ct->count_seq_file != NULL)
        fclose(ct->count_seq_file);
    if (ct->target_seq_file != NULL)
        fclose(ct->target_seq_file);
    memset(ct, 0, sizeof(corrector_t));
}

void corrector_read_fasta(corrector_t *ct, FILE *fp, bool is_target)
{
    int8_t c;
    seq_t seq = {0};

	FILE **out_fpp = &(ct->count_seq_file);
	if (is_target)
		out_fpp = &(ct->target_seq_file);

    if (*out_fpp == NULL && (*out_fpp = tmpfile_open()) == NULL) {
        fputs("error: cannot create tmpfile_open!\n", stderr);
        my_exit(1);
    }
    fseek(*out_fpp, 0, SEEK_END);
    while ((c = getc(fp)) != '>') {
        if (c == EOF) {
            fputs("error: invalid file format!\n", stderr);
            my_exit(1);
        }
    }
    while (c != EOF) {
        while ((c = getc(fp)) != '\n' && c != EOF)
            ;
        seq.len = 0;
        seq.n_unknown = 0;
        while ((c = getc(fp)) != '>' && c != EOF) {
            if (c == '\n')
                continue;
            if (seq.len >= MAX_READ_LEN) {
                fputs("error: too long sequence!\n", stderr);
                my_exit(1);
            }
            if (ascii2num(c) == 4) {
                seq.unknown_pos[seq.n_unknown] = seq.len;
                ++seq.n_unknown;
            }
            else
                seq.base[seq.len] = ascii2num(c);
            ++seq.len;
        }
        seq_write(&seq, *out_fpp);
		if (is_target)
			++ct->n_target_read[ct->n_target_file];
    }

	if (is_target)
		++ct->n_target_file;
}

void corrector_read_fastq(corrector_t *ct, FILE *fp, bool is_target)
{
    int64_t i;
    int8_t c;
    seq_t seq = {0};

	FILE **out_fpp = &(ct->count_seq_file);
	if (is_target)
		out_fpp = &(ct->target_seq_file);

    if (*out_fpp == NULL && (*out_fpp = tmpfile_open()) == NULL) {
        fputs("error: cannot create tmpfile_open!\n", stderr);
        my_exit(1);
    }
    fseek(*out_fpp, 0, SEEK_END);
    while ((c = getc(fp)) != '@') {
        if (c == EOF) {
            fputs("error: invalid file format!\n", stderr);
            my_exit(1);
        }
    }
    while (c != EOF) {
        while ((c = getc(fp)) != '\n' && c != EOF)
            ;
        seq.len = 0;
        seq.n_unknown = 0;
        while ((c = getc(fp)) != '+') {
            if (c == '\n')
                continue;
            if (seq.len >= MAX_READ_LEN) {
                fputs("error: too long sequence!\n", stderr);
                my_exit(1);
            }
            seq.base[seq.len] = ascii2num(c);
            if (seq.base[seq.len] == 4) {
                seq.unknown_pos[seq.n_unknown] = seq.len;
                ++seq.n_unknown;
                seq.base[i] = 0;
            }
            ++seq.len;
        }

        while ((c = getc(fp)) != '\n' && c != EOF)
            ;

        c = fastq_skip_qual(seq.len, fp);

        seq_write(&seq, *out_fpp);
		if (is_target)
			++ct->n_target_read[ct->n_target_file];
    }

	if (is_target)
		++ct->n_target_file;
}

void corrector_make_table(corrector_t *ct, uint64_t k_len)
{
    uint64_t i;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    uint64_t n_used;
    uint64_t max;
    uint64_t *distr;
    uint16_t occ;
    kmer_t kmer;
    seq_t seq;
    FILE *stored_fp;
    FILE *unstored_fp;
    FILE *tmp_fp;

    fprintf(stderr, "K=%lu\nmaking hash table...\n", k_len);

    if (ct->val_table != NULL)
        corrector_destroy_table(ct);
    ct->k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    ct->k_len = k_len;
    ct->table_size = 1ull;
    for (i = 0; i < k_len*2; ++i) {
        ct->table_size *= 2ull;
        if (ct->table_size*(sizeof(uint16_t)+sizeof(uint64_t)) > ct->mem_lim * (1.0-MEM_MARGIN))
            break;
    }
    ct->index_len = i;
    ct->table_size /= 2ull;

    ct->val_table = (uint16_t *)my_calloc(ct->table_size, sizeof(uint16_t));
    check_alloc(ct->val_table);
    ct->key_table = (uint64_t *)my_malloc(ct->table_size * sizeof(uint64_t));
    check_alloc(ct->key_table);

    stored_fp = tmpfile_open();
    check_tmpfile(stored_fp);
    unstored_fp = tmpfile_open();
    check_tmpfile(unstored_fp);

    rewind(ct->count_seq_file);
    n_used = 0;

    while (seq_read(&seq, ct->count_seq_file)) {
        if (seq.len < k_len)
            continue;
        seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
        seq.n_unknown = 0;
        kmer.fwd = kmer.rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer.fwd = kmer.fwd << 2 | (uint64_t)seq.base[i];
            kmer.rev = (kmer.rev | (uint64_t)(3^seq.base[k_len-i-2])) << 2;
        }
        for (i = 0; i < seq.len - k_len + 1; ++i) {
            kmer.fwd = ((kmer.fwd << 2) & ct->k_mask) | (uint64_t)seq.base[i+k_len-1];
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^seq.base[i+k_len-1]) << (2*(k_len-1)));
            if (seq.unknown_pos[seq.n_unknown] < i + k_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                continue;
            }

            key = min_u64(kmer.fwd, kmer.rev);
            index = hash64(key, ct->index_len) & (ct->table_size - 1);
            if (!(ct->val_table[index])) {
                if (n_used <= ct->table_size*MAX_LOAD_FACTOR) {
                    ct->val_table[index] = 1;
                    ct->key_table[index] = key;
                    ++n_used;
                }
                else
                    fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
            }
            else if (ct->key_table[index] == key) {
                if (ct->val_table[index] < UINT16_MAX-1)
                    ++ct->val_table[index];
            }
            else {
                step = rehash64(key, ct->index_len);
                while (1) {
                    index = (index+step) & (ct->table_size-1);
                    if (!(ct->val_table[index])) {
                        if (n_used <= ct->table_size*MAX_LOAD_FACTOR) {
                            ct->val_table[index] = 1;
                            ct->key_table[index] = key;
                            ++n_used;
                        }
                        else
                            fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
                        break;
                    }
                    else if (ct->key_table[index] == key) {
                        if (ct->val_table[index] < UINT16_MAX-1)
                            ++ct->val_table[index];
                        break;
                    }
                }
            }
        }
    }

    distr = (uint64_t *)my_calloc(UINT16_MAX, sizeof(uint64_t));
    check_alloc(distr);

    max = 0;
    while (1) {
        for (i = 0; i < ct->table_size; ++i) {
            if (ct->val_table[i]) {
                occ = ct->val_table[i];
                fwrite(&ct->key_table[i], sizeof(uint64_t), 1, stored_fp);
                fwrite(&occ, sizeof(uint16_t), 1, stored_fp);
                ++distr[occ];
                ct->val_table[i] = 0;
                if (occ > max) 
                    max = occ;
            }
        }

        if (!(ftell(unstored_fp) / sizeof(uint64_t)))
            break;

        rewind(unstored_fp);
        n_used = 0;
        tmp_fp = tmpfile_open();
        check_tmpfile(tmp_fp);
        while (fread(&key, sizeof(uint64_t), 1, unstored_fp)) {
            index = hash64(key, ct->index_len) & (ct->table_size - 1);
            if (!(ct->val_table[index])) {
                if (n_used <= ct->table_size*MAX_LOAD_FACTOR) {
                    ct->val_table[index] = 1;
                    ct->key_table[index] = key;
                    ++n_used;
                }
                else
                    fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
            }
            else if (ct->key_table[index] == key) {
                if (ct->val_table[index] < UINT16_MAX-1)
                    ++ct->val_table[index];
            }
            else {
                step = rehash64(key, ct->index_len);
                while (1) {
                    index = (index+step) & (ct->table_size-1);
                    if (!(ct->val_table[index])) {
                        if (n_used <= ct->table_size*MAX_LOAD_FACTOR) {
                            ct->val_table[index] = 1;
                            ct->key_table[index] = key;
                            ++n_used;
                        }
                        else
                            fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
                        break;
                    }
                    else if (ct->key_table[index] == key) {
                        if (ct->val_table[index] < UINT16_MAX-1)
                            ++ct->val_table[index];
                        break;
                    }
                }
            }
        }
        fclose(unstored_fp);
        unstored_fp = tmp_fp;
    }
    fclose(unstored_fp);

    if (ct->th < 0)
		ct->th = get_left_minimal_smooth(distr, max, SMOOTHING_WINDOW);

    n_used = 0;
    for (i = ct->th; i <= max; ++i)
        n_used += distr[i];

    fprintf(stderr, "KMER OCCURRENCE THRESHOLD = %lu\n", ct->th);
    fputs("KMER OCCURRENCE DISTRIBUTION:\nBEGIN\n", stderr);
    for (i = 1; i <= max; ++i)
        fprintf(stderr, "%lu\t%lu\n", i, distr[i]);
    fputs("END\n", stderr);

    my_free(distr);

    while ((double)ct->table_size * MAX_LOAD_FACTOR < n_used) {
        ct->table_size *= 2;
        ++ct->index_len;
    }
    check_mem_usage(ct->table_size*(sizeof(uint64_t)+sizeof(uint16_t)), ct->mem_lim);
    fprintf(stderr, "BUCKET_SIZE=%lu, ", ct->table_size); 
    fprintf(stderr, "COV_CUTOFF=%lu, ", ct->th); 
    fprintf(stderr, "NUM_OCCUPIED_SLOT=%lu, ", n_used);  
    fprintf(stderr, "LOAD_FACTOR=%f, ", (double)n_used / (double)ct->table_size);
    fprintf(stderr, "MEM_USAGE=%luB\n", ct->table_size * (sizeof(uint16_t)+sizeof(uint64_t))); 

    ct->val_table = my_realloc(ct->val_table, ct->table_size*sizeof(uint16_t));
    check_alloc(ct->val_table);
    memset(ct->val_table, 0, ct->table_size*sizeof(uint16_t));
    ct->key_table = my_realloc(ct->key_table, ct->table_size*sizeof(uint64_t));
    check_alloc(ct->key_table);

    rewind(stored_fp);
    while (fread(&key, sizeof(uint64_t), 1, stored_fp)) {
        fread(&occ, sizeof(uint16_t), 1, stored_fp);
        if (occ < ct->th)
            continue;
        index = hash64(key, ct->index_len) & (ct->table_size - 1);
        if (!(ct->val_table[index])) {
            ct->val_table[index] = occ;
            ct->key_table[index] = key;
        }
        else {
            step = rehash64(key, ct->index_len);
            index = (index+step) & (ct->table_size-1);
            while (ct->val_table[index])
                index = (index+step) & (ct->table_size-1);
            ct->val_table[index] = occ;
            ct->key_table[index] = key;
        }
    }

    fclose(stored_fp);
}

void corrector_ungap_correct(corrector_t *ct)
{
    int64_t i;
    int64_t j;
    int64_t k_len;
    int64_t occ;
    int64_t max_occ;
    int64_t st;
    int64_t ed;
    int64_t maxd;
    int64_t min_score;
    int64_t score;
    int64_t left_score;
    uint8_t best_base;
	int64_t next_display_n = DISPLAY_NUM_READ_UNIT;
	int64_t n_read = 0;
    seq_t seq;
    seq_t raw_seq;
    seq_t new_seq;
    kmer_t kmer[MAX_READ_LEN];
    kmer_t best_kmer;
    FILE *new_file;

    if (ct->th <= 1)
        return;
    k_len = ct->k_len;
    ct->n_cor_read = 0;
    ct->n_cor_base = 0;
    new_file = tmpfile_open();
    check_tmpfile(new_file);

    rewind(ct->target_seq_file);
    while (seq_read(&raw_seq, ct->target_seq_file)) {
		if (ct->id == 0) {
			++n_read;
			if (n_read >= next_display_n) {
				fputs(".", stderr);
				fflush(stderr);
				next_display_n += DISPLAY_NUM_READ_UNIT;
			}
		}

        if (raw_seq.len < k_len) {
            seq_write(&raw_seq, new_file);
            continue;
        }
        seq = raw_seq;
        seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
        seq.n_unknown = 0;
        kmer[0].fwd = kmer[0].rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer[0].fwd = kmer[0].fwd << 2 | (uint64_t)seq.base[i];
            kmer[0].rev = (kmer[0].rev | (uint64_t)(3^seq.base[k_len-i-2])) << 2;
        }
        j = st = ed = 0;
        for (i = 0; i < seq.len - k_len + 1; ++i) {
            kmer[i].fwd = ((kmer[i].fwd << 2) & ct->k_mask) | (uint64_t)seq.base[i+k_len-1];
            kmer[i].rev = (kmer[i].rev >> 2) | ((uint64_t)(3^seq.base[i+k_len-1]) << (2*(k_len-1)));
            kmer[i+1] = kmer[i];
            if (seq.unknown_pos[seq.n_unknown] < i + k_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                occ = 0;
            }
            else
                occ = corrector_occ(ct, min_u64(kmer[i].fwd, kmer[i].rev));

            if (!occ) {
                if (i - j > ed - st) {
                    st = j;
                    ed = i;
                }
                j = i + 1;
            }
        }
        if (i - j > ed - st) {
            st = j;
            ed = i;
        }
        ed += k_len-1;
        if (ed - st == seq.len || ed < k_len) {
            seq_write(&raw_seq, new_file);
            continue;
        }
        for (i = 0; i < seq.n_unknown; ++i)
            seq.base[seq.unknown_pos[i]] = 4;
        maxd = (int64_t)(seq.len * ct->max_edit);
/*LEFT*/
/************************************************************************************************/
        score = 0;
        min_score = 0;
        for (i = st; i > 0; --i) {
            min_score = maxd + 1;
			bool valid_flag = false;
            for (j = 0; j < 4; ++j) {
                kmer[i-1].fwd = (kmer[i].fwd >> 2) | ((uint64_t)j << (2*(k_len-1)));
                kmer[i-1].rev = ((kmer[i].rev << 2) & ct->k_mask) | (uint64_t)(3^j);
                occ = corrector_occ(ct, min_u64(kmer[i-1].fwd, kmer[i-1].rev));
                if (!occ || occ > ct->upper_th)
                    continue;
				valid_flag = true;
                if (seq.base[i-1] == j) {
                    if (score < min_score || (score == min_score && occ > max_occ)) {
                        min_score = score;
                        max_occ = occ;
                        best_base = j;
                        best_kmer = kmer[i-1];
                    }
                    break;
                }
                else {
                    if (score+1 < min_score || (score+1 == min_score && occ > max_occ)) {
                        min_score = score+1;
                        max_occ = occ;
                        best_base = j;
                        best_kmer = kmer[i-1];
                    }
                }
            }
			if (!valid_flag) {
				min_score = score;
				best_base = seq.base[i-1];
				best_kmer = kmer[i-1];
			}

            score = min_score;
            kmer[i-1] = best_kmer;
            new_seq.base[i-1] = best_base;
        }
		left_score = min_score;
/************************************************************************************************/
        for (i = st; i < ed; ++i)
            new_seq.base[i] = seq.base[i];
/*RIGHT*/
/************************************************************************************************/
        maxd = min_64(seq.len-ed, (int64_t)(seq.len*ct->max_edit) - min_score);
        for (i = ed; i < seq.len; ++i) {
            min_score = maxd + 1;
			bool valid_flag = false;
            for (j = 0; j < 4; ++j) {
                kmer[i-k_len+1].fwd = ((kmer[i-k_len].fwd << 2) & ct->k_mask) | (uint64_t)j;
                kmer[i-k_len+1].rev = (kmer[i-k_len].rev >> 2) | ((uint64_t)(3^j) << (2*(k_len-1)));
                occ = corrector_occ(ct, min_u64(kmer[i-k_len+1].fwd, kmer[i-k_len+1].rev));
                if (!occ || occ > ct->upper_th) 
                    continue;
				valid_flag = true;
                if (seq.base[i] == j) {
                    if (score < min_score || (score == min_score && occ > max_occ)) {
                        min_score = score;
                        max_occ = occ;
                        best_base = j;
                        best_kmer = kmer[i-k_len+1];
                    }
                    break;
                }
                else {
                    if (score+1 < min_score || (score+1 == min_score && occ > max_occ)) {
                        min_score = score+1;
                        max_occ = occ;
                        best_base = j;
                        best_kmer = kmer[i-k_len+1];
                    }
                }
            }
            if (!valid_flag) {
				min_score = score;
				best_base = seq.base[i];
				best_kmer = kmer[i-k_len+1];
			}

            score = min_score;
            kmer[i-k_len+1] = best_kmer;
            new_seq.base[i] = best_base;
        }
/************************************************************************************************/
        new_seq.len = seq.len;
        new_seq.n_unknown = 0;
        seq_write(&new_seq, new_file);
        ++ct->n_cor_read;
		ct->n_cor_base += left_score + min_score;
    }
    fclose(ct->target_seq_file);
    ct->target_seq_file = new_file;

	if (ct->id == 0) {
		fputs("\n", stderr);
		fflush(stderr);
	}
}

void corrector_gapped_correct(corrector_t *ct)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t k_len;
    int64_t occ;
    int64_t max_occ;
    int64_t st;
    int64_t ed;
    int64_t maxd;
    int64_t min_score;
    int64_t min_dist;
    int64_t left_dist;
    int64_t min_pos;
    int64_t best_track;
    int64_t n_col; 
    int16_t *score_mat;
    uint8_t *track_mat;
	int64_t next_display_n = DISPLAY_NUM_READ_UNIT;
	int64_t n_read = 0;
    seq_t seq;
    seq_t raw_seq;
    seq_t new_seq;
    kmer_t *kmer_mat;
    kmer_t kmer[MAX_READ_LEN];
    kmer_t tmp_kmer;
    kmer_t best_kmer;
    FILE *new_file;

    if (ct->th <= 1)
        return;
    k_len = ct->k_len;
    ct->n_cor_read = 0;
    ct->n_cor_base = 0;
    new_file = tmpfile_open();
    check_tmpfile(new_file);

    n_col = 2 * (int64_t)(MAX_READ_LEN * ct->max_edit + 0.5) + 1; 

    check_mem_usage(n_col*(MAX_READ_LEN+1)*(sizeof(int16_t)+sizeof(uint8_t)+sizeof(kmer_t)), ct->mem_lim);
    score_mat = (int16_t *)my_malloc((n_col*(MAX_READ_LEN+1))*sizeof(int16_t));
    check_alloc(score_mat);
    track_mat = (uint8_t *)my_malloc((n_col*(MAX_READ_LEN+1))*sizeof(uint8_t));
    check_alloc(track_mat);
    kmer_mat = (kmer_t *)my_malloc((n_col*(MAX_READ_LEN+1))*sizeof(kmer_t));
    check_alloc(kmer_mat);

    rewind(ct->target_seq_file);
    while (seq_read(&raw_seq, ct->target_seq_file)) {
		if (ct->id == 0) {
			++n_read;
			if (n_read >= next_display_n) {
				fputs(".", stderr);
				fflush(stderr);
				next_display_n += DISPLAY_NUM_READ_UNIT;
			}
		}

        if (raw_seq.len < k_len) {
            seq_write(&raw_seq, new_file);
            continue;
        }
        seq = raw_seq;
        seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
        seq.n_unknown = 0;
        kmer[0].fwd = kmer[0].rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer[0].fwd = kmer[0].fwd << 2 | (uint64_t)seq.base[i];
            kmer[0].rev = (kmer[0].rev | (uint64_t)(3^seq.base[k_len-i-2])) << 2;
        }
        j = st = ed = 0;
        for (i = 0; i < seq.len - k_len + 1; ++i) {
            kmer[i].fwd = ((kmer[i].fwd << 2) & ct->k_mask) | (uint64_t)seq.base[i+k_len-1];
            kmer[i].rev = (kmer[i].rev >> 2) | ((uint64_t)(3^seq.base[i+k_len-1]) << (2*(k_len-1)));
            kmer[i+1] = kmer[i];
            if (seq.unknown_pos[seq.n_unknown] < i + k_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                occ = 0;
            }
            else
                occ = corrector_occ(ct, min_u64(kmer[i].fwd, kmer[i].rev));

            if (!occ) {
                if (i - j > ed - st) {
                    st = j;
                    ed = i;
                }
                j = i + 1;
            }
        }
        if (i - j > ed - st) {
            st = j;
            ed = i;
        }
        ed += k_len-1;
        if (ed - st == seq.len || ed < k_len) {
            seq_write(&raw_seq, new_file);
            continue;
        }
        for (i = 0; i < seq.n_unknown; ++i)
            seq.base[seq.unknown_pos[i]] = 4;
/*LEFT*/
/************************************************************************************************/
        maxd = min_64(st, (int64_t)(seq.len*ct->max_edit+0.5));
        for (i = 0; i < maxd; ++i)
            score_mat[i] = maxd+1;
        score_mat[maxd] = 0;
        track_mat[maxd] = 9;
        kmer_mat[maxd] = kmer[st];
        min_dist = 0;
        min_pos = maxd;
        for (i = maxd; i < 2*maxd; ++i) {
            if (score_mat[i] > maxd) {
                score_mat[i+1] = maxd + 1;
                continue;
            }
            max_occ = 0;
            for (j = 0; j < 4; ++j) {
                tmp_kmer.fwd = (kmer_mat[i].fwd >> 2) | ((uint64_t)j << (2*(k_len-1)));
                tmp_kmer.rev = ((kmer_mat[i].rev << 2) & ct->k_mask) | (uint64_t)(3^j);
                occ = corrector_occ(ct, min_u64(tmp_kmer.fwd, tmp_kmer.rev));
                if (occ && occ > max_occ) {
                    max_occ = occ;
                    best_track = 4 + j;
                    best_kmer = tmp_kmer;
                }
            }
            if (max_occ > 0) {
                score_mat[i+1] = score_mat[i] + 1;
                track_mat[i+1] = best_track;
                kmer_mat[i+1] = best_kmer;
            }
            else
                score_mat[i+1] = maxd + 1;
        }
        for (i = 0; i < st; ++i) {
            min_dist = maxd + 1;
            for (j = 0; j <= 2*maxd; ++j) {
                min_score = maxd + 1;
                if (score_mat[i*n_col + j] <= maxd) {
                    for (k = 0; k < 4; ++k) {
                        tmp_kmer.fwd = (kmer_mat[i*n_col + j].fwd >> 2) | ((uint64_t)k << (2*(k_len-1)));
                        tmp_kmer.rev = ((kmer_mat[i*n_col + j].rev << 2) & ct->k_mask) | (uint64_t)(3^k);
                        occ = corrector_occ(ct, min_u64(tmp_kmer.fwd, tmp_kmer.rev));
                        if (!occ) 
                            continue;
                        if (seq.base[st-i-1] == k) {
                            if (score_mat[i*n_col + j] < min_score || (score_mat[i*n_col + j] == min_score && occ > max_occ)) {
                                min_score = score_mat[i*n_col + j];
                                max_occ = occ;
                                best_track = k;
                                best_kmer = tmp_kmer;
                            }
                            break;
                        }
                        else {
                            if (score_mat[i*n_col + j]+1 < min_score || (score_mat[i*n_col + j]+1 == min_score && occ > max_occ)) {
                                min_score = score_mat[i*n_col + j]+1;
                                max_occ = occ;
                                best_track = k;
                                best_kmer = tmp_kmer;
                            }
                        }
                    }
                }
                if (j > 0 && score_mat[(i+1)*n_col + (j-1)] < maxd) {
                    for (k = 0; k < 4; ++k) {
                        tmp_kmer.fwd = (kmer_mat[(i+1)*n_col + (j-1)].fwd >> 2) | ((uint64_t)k << (2*(k_len-1)));
                        tmp_kmer.rev = ((kmer_mat[(i+1)*n_col + (j-1)].rev << 2) & ct->k_mask) | (uint64_t)(3^k);
                        occ = corrector_occ(ct, min_u64(tmp_kmer.fwd, tmp_kmer.rev));
                        if (!occ) 
                            continue;
                        if (score_mat[(i+1)*n_col + (j-1)]+1 < min_score || (score_mat[(i+1)*n_col + (j-1)]+1 == min_score && occ > max_occ)) {
                            min_score = score_mat[(i+1)*n_col + (j-1)]+1;
                            max_occ = occ;
                            best_track = 4 + k;
                            best_kmer = tmp_kmer;
                        }
                    }
                }
                if (j < 2*maxd && score_mat[i*n_col + (j+1)] < maxd) {
                    occ = corrector_occ(ct, min_u64(kmer_mat[i*n_col+(j+1)].fwd, kmer_mat[i*n_col+(j+1)].rev));
                    if (!occ) 
                        continue;
                    if (score_mat[i*n_col + (j+1)]+1 < min_score || (score_mat[i*n_col + (j+1)]+1 == min_score && occ > max_occ)) {
                        min_score = score_mat[i*n_col + (j+1)]+1;
                        max_occ = occ;
                        best_track = 8;
                        best_kmer = kmer_mat[i*n_col + (j+1)];
                    }
                }
                if (min_score <= maxd) {
                    score_mat[(i+1)*n_col + j] = min_score;
                    track_mat[(i+1)*n_col + j] = best_track;
                    kmer_mat[(i+1)*n_col + j] = best_kmer;
                    if (min_score < min_dist) {
                        min_dist = min_score;
                        min_pos = j;
                    }
                }
                else
                    score_mat[(i+1)*n_col + j] = maxd + 1;
            }
            if (min_dist > maxd)
                break;
        }

        if (i != st || seq.len == MAX_READ_LEN) {
            seq_write(&raw_seq, new_file);
            continue;
        }
		left_dist = min_dist;

        new_seq.len = 0;
        j = min_pos;
        while (track_mat[i*n_col + j] != 9) {
            if (track_mat[i*n_col + j] < 4) {
                new_seq.base[new_seq.len] = track_mat[i*n_col + j];    
                --i;
                ++new_seq.len;
            }
            else if (track_mat[i*n_col + j] < 8) {
                new_seq.base[new_seq.len] = track_mat[i*n_col + j] - 4;
                --j;
                ++new_seq.len;
            }
            else if (track_mat[i*n_col + j] == 8) {
                --i;
                ++j;
            }
        }
/************************************************************************************************/
        for (i = st; i < ed; ++i) {
            new_seq.base[new_seq.len] = seq.base[i];
            ++new_seq.len;
        }
/*RIGHT*/
/************************************************************************************************/
        maxd = min_64(seq.len - ed,  (int64_t)(seq.len * ct->max_edit+0.5) - min_dist);
        for (i = 0; i < maxd; ++i)
            score_mat[i] = maxd+1;
        score_mat[maxd] = 0;
        track_mat[maxd] = 9;
        kmer_mat[maxd] = kmer[ed-k_len];
        min_dist = 0;
        min_pos = maxd;
        for (i = maxd; i < 2*maxd; ++i) {
            if (score_mat[i] > maxd) {
                score_mat[i+1] = maxd + 1;
                continue;
            }
            max_occ = 0;
            for (j = 0; j < 4; ++j) {
                tmp_kmer.fwd = ((kmer_mat[i].fwd << 2) & ct->k_mask) | (uint64_t)j;
                tmp_kmer.rev = (kmer_mat[i].rev >> 2) | ((uint64_t)(3^j) << (2*(k_len-1)));
                occ = corrector_occ(ct, min_u64(tmp_kmer.fwd, tmp_kmer.rev));
                if (occ && occ > max_occ) {
                    max_occ = occ;
                    best_track = 4 + j;
                    best_kmer = tmp_kmer;
                }
            }
            if (max_occ > 0) {
                score_mat[i+1] = score_mat[i] + 1;
                track_mat[i+1] = best_track;
                kmer_mat[i+1] = best_kmer;
            }
            else
                score_mat[i+1] = maxd + 1;
        }
        for (i = 0; i < seq.len-ed; ++i) {
            min_dist = maxd + 1;
            for (j = 0; j <= 2*maxd; ++j) {
                min_score = maxd + 1;
                if (score_mat[i*n_col + j] <= maxd) {
                    for (k = 0; k < 4; ++k) {
                        tmp_kmer.fwd = ((kmer_mat[i*n_col + j].fwd << 2) & ct->k_mask) | (uint64_t)k;
                        tmp_kmer.rev = (kmer_mat[i*n_col + j].rev >> 2) | ((uint64_t)(3^k) << (2*(k_len-1)));
                        occ = corrector_occ(ct, min_u64(tmp_kmer.fwd, tmp_kmer.rev));
                        if (!occ)
                            continue;
                        if (seq.base[ed+i] == k) {
                            if (score_mat[i*n_col + j] < min_score || (score_mat[i*n_col + j] == min_score && occ > max_occ)) {
                                min_score = score_mat[i*n_col + j];
                                max_occ = occ;
                                best_track = k;
                                best_kmer = tmp_kmer;
                            }
                            break;
                        }
                        else {
                            if (score_mat[i*n_col + j]+1 < min_score || (score_mat[i*n_col + j]+1 == min_score && occ > max_occ)) {
                                min_score = score_mat[i*n_col + j]+1;
                                max_occ = occ;
                                best_track = k;
                                best_kmer = tmp_kmer;
                            }
                        }
                    }
                }
                if (j > 0 && score_mat[(i+1)*n_col + (j-1)] < maxd) {
                    for (k = 0; k < 4; ++k) {
                        tmp_kmer.fwd = ((kmer_mat[(i+1)*n_col + (j-1)].fwd << 2) & ct->k_mask) | (uint64_t)k;
                        tmp_kmer.rev = (kmer_mat[(i+1)*n_col + (j-1)].rev >> 2) | ((uint64_t)(3^k) << (2*(k_len-1)));
                        occ = corrector_occ(ct, min_u64(tmp_kmer.fwd, tmp_kmer.rev));
                        if (!occ) 
                            continue;
                        if (score_mat[(i+1)*n_col + (j-1)]+1 < min_score || (score_mat[(i+1)*n_col + (j-1)]+1 == min_score && occ > max_occ)) {
                            min_score = score_mat[(i+1)*n_col + (j-1)]+1;
                            max_occ = occ;
                            best_track = 4 + k;
                            best_kmer = tmp_kmer;
                        }
                    }
                }
                if (j < 2*maxd && score_mat[i*n_col + (j+1)] < maxd) {
                    occ = corrector_occ(ct, min_u64(kmer_mat[i*n_col+(j+1)].fwd, kmer_mat[i*n_col+(j+1)].rev));
                    if (!occ) 
                        continue;
                    if (score_mat[i*n_col + (j+1)]+1 < min_score || (score_mat[i*n_col + (j+1)]+1 == min_score && occ > max_occ)) {
                        min_score = score_mat[i*n_col + (j+1)]+1;
                        max_occ = occ;
                        best_track = 8;
                        best_kmer = kmer_mat[i*n_col + (j+1)];
                    }
                }
                if (min_score <= maxd) {
                    score_mat[(i+1)*n_col + j] = min_score;
                    track_mat[(i+1)*n_col + j] = best_track;
                    kmer_mat[(i+1)*n_col + j] = best_kmer;
                    if (min_score < min_dist) {
                        min_dist = min_score;
                        min_pos = j;
                    }
                }
                else
                    score_mat[(i+1)*n_col + j] = maxd + 1;
            }
            if (min_dist > maxd)
                break;
        }

        if (i != seq.len-ed) {
            seq_write(&raw_seq, new_file);
            continue;
        }
        new_seq.len += seq.len-ed+min_pos-maxd;
        j = min_pos;
        while (track_mat[i*n_col + j] != 9) {
            if (track_mat[i*n_col + j] < 4) {
                --new_seq.len;
                new_seq.base[new_seq.len] = track_mat[i*n_col + j];    
                --i;
            }
            else if (track_mat[i*n_col + j] < 8) {
                --new_seq.len;
                new_seq.base[new_seq.len] = track_mat[i*n_col + j] - 4;
                --j;
            }
            else if (track_mat[i*n_col + j] == 8) {
                --i;
                ++j;
            }
        }
/************************************************************************************************/
        new_seq.len += seq.len-ed+min_pos-maxd;
        new_seq.n_unknown = 0;
        seq_write(&new_seq, new_file);
        ++ct->n_cor_read;
		ct->n_cor_base += left_dist + min_dist;
    }
    fclose(ct->target_seq_file);
    ct->target_seq_file = new_file;
    my_free(score_mat);
    my_free(track_mat);
    my_free(kmer_mat);

	if (ct->id == 0) {
		fputs("\n", stderr);
		fflush(stderr);
	}
}

void corrector_show_seq(corrector_t *ct)
{
    uint64_t i;
    uint64_t j;
    uint64_t k;
    uint8_t base_name[] = {'A', 'C', 'G', 'T', 'N'};
    seq_t seq;
    FILE *out[MAX_FILE_NUM * 2];
    
    for (i = 0; i < ct->n_target_file; ++i) {
        if((out[i] = fopen(ct->fname[i], "w")) == NULL) {
            fprintf(stderr, "error: cannot open \"%s\"\n", ct->fname[i]);
            my_exit(1);
        }
    }

    rewind(ct->target_seq_file);
    for (i = 0; i < ct->n_target_file; ++i) {
        for (j = 0; j < ct->n_target_read[i]; ++j) {
            seq_read(&seq, ct->target_seq_file);
            for (k = 0; k < seq.n_unknown; ++k)
                seq.base[seq.unknown_pos[k]] = 4;
            fprintf(out[i], ">s%lu\n", j+1);
            for (k = 0; k < seq.len; ++k)
                putc(base_name[seq.base[k]], out[i]);
            putc('\n', out[i]);
        }
    }

    for (i = 0; i < ct->n_target_file; ++i)
        fclose(out[i]);
}

double corrector_get_ave_cov(corrector_t *ct)
{
    uint64_t i;
    uint64_t num;
    uint64_t total;

    total = num = 0;
    for (i = 0; i < ct->table_size; ++i) {
        if (!ct->val_table[i])
            continue;
        ++num;
        total += ct->val_table[i];
    }

    return (double)total / num;
}

void corrector_set_id(corrector_t *ct, int64_t id)
{
    ct->id = id;
}

void corrector_set_th(corrector_t *ct, int64_t th)
{
    ct->th = th;
}

void corrector_set_upper_th(corrector_t *ct, int64_t upper_th)
{
    ct->upper_th = upper_th;
}

void corrector_set_max_edit(corrector_t *ct, double max_edit)
{
    ct->max_edit = max_edit;
}

void corrector_threads_init(corrector_threads_t *ctt, uint64_t n_thread, uint64_t mem_lim)
{
    uint64_t i;

    memset(ctt, 0, sizeof(corrector_threads_t));
    ctt->n_thread = n_thread;
    ctt->mem_lim = mem_lim;
    for (i = 0; i < n_thread; ++i) {
        corrector_init(&ctt->ct[i], mem_lim);
        ctt->ct[i].id = i;
    }
}
        
void corrector_threads_destroy_table(corrector_threads_t *ctt)
{
    if (ctt->val_table != NULL)
        my_free(ctt->val_table);
    ctt->val_table = NULL;
    if (ctt->key_table != NULL)
        my_free(ctt->key_table);
    ctt->key_table = NULL;
}

void corrector_threads_destroy(corrector_threads_t *ctt)
{
    uint64_t i;

    corrector_threads_destroy_table(ctt);
    for (i = 0; i < ctt->n_thread; ++i) {
        if (ctt->ct[i].count_seq_file != NULL)
            fclose(ctt->ct[i].count_seq_file);
        if (ctt->ct[i].target_seq_file != NULL)
            fclose(ctt->ct[i].target_seq_file);
	}
    memset(ctt, 0, sizeof(corrector_t));
}

void corrector_threads_read_fasta(corrector_threads_t *ctt, FILE *fp, bool is_target)
{
    uint64_t i;
    int8_t c;
    seq_t seq = {0};

	FILE **out_fpp[MAX_THREAD];
    for (i = 0; i < ctt->n_thread; ++i) {
		if (is_target)
			out_fpp[i] = &(ctt->ct[i].target_seq_file);
		else
			out_fpp[i] = &(ctt->ct[i].count_seq_file);

        if (*(out_fpp[i]) == NULL && (*(out_fpp[i]) = tmpfile_open()) == NULL) {
            fputs("error: cannot create tmpfile_open!\n", stderr);
            my_exit(1);
        }
        fseek(*(out_fpp[i]), 0, SEEK_END);
    }
    while ((c = getc(fp)) != '>') {
        if (c == EOF) {
            fputs("error: invalid file format!\n", stderr);
            my_exit(1);
        }
    }
    while (c != EOF) {
        while ((c = getc(fp)) != '\n' && c != EOF)
            ;
        seq.len = 0;
        seq.n_unknown = 0;
        while ((c = getc(fp)) != '>' && c != EOF) {
            if (c == '\n')
                continue;
            if (seq.len >= MAX_READ_LEN) {
                fputs("error: too long sequence!\n", stderr);
                my_exit(1);
            }
            if (ascii2num(c) == 4) {
                seq.unknown_pos[seq.n_unknown] = seq.len;
                ++seq.n_unknown;
            }
            else
                seq.base[seq.len] = ascii2num(c);
            ++seq.len;
        }
        seq_write(&seq, *(out_fpp[ctt->fid]));
        ctt->fid = (ctt->fid+1) % ctt->n_thread;
		if (is_target)
			++ctt->n_target_read[ctt->n_target_file];
    }
	if (is_target)
		++ctt->n_target_file;
}

void corrector_threads_read_fastq(corrector_threads_t *ctt, FILE *fp, bool is_target)
{
    uint16_t i;
    int8_t c;
    uint8_t *u8_p;
    uint8_t base[MAX_READ_LEN];
    uint16_t len;
    uint16_t n_unknown;
    uint16_t unknown_pos[MAX_READ_LEN+1];

	FILE **out_fpp[MAX_THREAD];
    for (i = 0; i < ctt->n_thread; ++i) {
		if (is_target)
			out_fpp[i] = &(ctt->ct[i].target_seq_file);
		else
			out_fpp[i] = &(ctt->ct[i].count_seq_file);

        if (*(out_fpp[i]) == NULL && (*(out_fpp[i]) = tmpfile_open()) == NULL) {
            fputs("error: cannot create tmpfile_open!\n", stderr);
            my_exit(1);
        }
        fseek(*(out_fpp[i]), 0, SEEK_END);
    }
    while ((c = getc(fp)) != '@') {
        if (c == EOF) {
            fputs("error: invalid file format!\n", stderr);
            my_exit(1);
        }
    }
    while (c != EOF) {
        while ((c = getc(fp)) != '\n' && c != EOF)
            ;
        len = 0;
        n_unknown = 0;
        while ((c = getc(fp)) != '+') {
            if (c == '\n')
                continue;
            if (len >= MAX_READ_LEN) {
                fputs("error: too long sequence!\n", stderr);
                my_exit(1);
            }
            base[len] = ascii2num(c);
            if (base[len] == 4) {
                unknown_pos[n_unknown] = len;
                ++n_unknown;
                base[len] = 0;
            }
            ++len;
        }

        while ((c = getc(fp)) != '\n' && c != EOF)
            ;

        c = fastq_skip_qual(len, fp);
        u8_p = base;

        fwrite(&n_unknown, sizeof(uint16_t), 1, *(out_fpp[ctt->fid])); 
        fwrite(unknown_pos, sizeof(uint16_t), n_unknown, *(out_fpp[ctt->fid])); 
        fwrite(&len, sizeof(uint16_t), 1, *(out_fpp[ctt->fid])); 
        fwrite(u8_p, sizeof(uint8_t), len, *(out_fpp[ctt->fid])); 

        ctt->fid = (ctt->fid+1) % ctt->n_thread;
		if (is_target)
			++ctt->n_target_read[ctt->n_target_file];
    }

	if (is_target)
		++ctt->n_target_file;
}

void corrector_threads_make_table(corrector_threads_t *ctt, uint64_t k_len)
{
    uint64_t i;
    uint64_t j;
    uint64_t max;
    uint64_t n_used;
    uint64_t *distr;
    uint16_t occ;

    fprintf(stderr, "K=%lu\nmaking hash table...\n", k_len);

    if (ctt->val_table != NULL)
        corrector_threads_destroy_table(ctt);
    if (k_len > 32 || k_len < 0) {
        fprintf(stderr, "error: function \'corrector_threads_make_table\' cannot accept k_len=%lu\n", k_len);
        my_exit(1);
    }
    ctt->k_len = k_len;
    ctt->table_size = 1ull;
    for (i = 0; i < k_len*2; ++i) {
        ctt->table_size *= 2ull;
        if (ctt->table_size*(sizeof(uint16_t)+sizeof(uint64_t))+(ctt->table_size/GRANULE+1)*sizeof(pthread_mutex_t) > ctt->mem_lim * (1.0-MEM_MARGIN))
            break;
    }
    ctt->index_len = i;
    ctt->table_size /= 2ull;

    ctt->val_table = (uint16_t *)my_calloc(ctt->table_size, sizeof(uint16_t));
    check_alloc(ctt->val_table);
    ctt->key_table = (uint64_t *)my_malloc(ctt->table_size * sizeof(uint64_t));
    check_alloc(ctt->key_table);

    g_ctt_table_mutex = (pthread_mutex_t *)my_malloc((ctt->table_size/GRANULE+1) * sizeof(pthread_mutex_t));
    check_alloc(g_ctt_table_mutex);
    for (i = 0; i < ctt->table_size/GRANULE+1; ++i)
        pthread_mutex_init(&g_ctt_table_mutex[i], NULL);


    for (i = 0; i < ctt->n_thread; ++i) {
        ctt->ct[i].val_table = ctt->val_table;
        ctt->ct[i].key_table = ctt->key_table;
        ctt->ct[i].table_size = ctt->table_size;
        ctt->ct[i].k_len = k_len;
        ctt->ct[i].k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
        ctt->ct[i].index_len = ctt->index_len;
        g_ctt_unstored_fp[i] = tmpfile_open();
        check_tmpfile(g_ctt_unstored_fp[i]);
        pthread_create(&ctt->thread[i], NULL, (void *)corrector_threads_count, &ctt->ct[i]);
    }
    for (i = 0; i < ctt->n_thread; ++i) {
        pthread_join(ctt->thread[i], NULL); 
        g_ctt_stored_fp[i] = tmpfile_open();
        check_tmpfile(g_ctt_stored_fp[i]);
    }

    distr = (uint64_t *)my_calloc(UINT16_MAX-1, sizeof(uint64_t));
    check_alloc(distr);

    max = 0;
    while (1) {
        j = 0;
        for (i = 0; i < ctt->table_size; ++i) {
            if(ctt->val_table[i]) {
                occ = ctt->val_table[i];
                fwrite(&ctt->key_table[i], sizeof(uint64_t), 1, g_ctt_stored_fp[j]);
                fwrite(&occ, sizeof(uint16_t), 1, g_ctt_stored_fp[j]);
                ++distr[occ];
                ctt->val_table[i] = 0;
                if (occ > max) 
                    max = occ;
                j = (j+1) % ctt->n_thread;
            }
        }

        for (i = 0; i < ctt->n_thread; ++i) {
            if(ftell(g_ctt_unstored_fp[i]) / sizeof(uint64_t))
                break;
        }
        if (i == ctt->n_thread)
            break;

        for (i = 0; i < ctt->n_thread; ++i)
            pthread_create(&ctt->thread[i], NULL, (void *)corrector_threads_recount, &ctt->ct[i]);
        for (i = 0; i < ctt->n_thread; ++i)
            pthread_join(ctt->thread[i], NULL); 
    }

	if (ctt->th < 0)
		ctt->th = get_left_minimal_smooth(distr, max, SMOOTHING_WINDOW);

    n_used = 0;
    for (i = ctt->th; i <= max; ++i)
        n_used += distr[i];


    fprintf(stderr, "KMER OCCURRENCE THRESHOLD = %lu\n", ctt->th);
    fputs("KMER OCCURRENCE DISTRIBUTION:\nBEGIN\n", stderr);
    for (i = 1; i <= max; ++i)
        fprintf(stderr, "%lu\t%lu\n", i, distr[i]);
    fputs("END\n", stderr);

    my_free(distr);

    while ((double)ctt->table_size * MAX_LOAD_FACTOR < n_used) {
        ctt->table_size *= 2;
        ++ctt->index_len;
    }
    check_mem_usage(ctt->table_size*(sizeof(uint64_t)+sizeof(uint16_t))+(ctt->table_size/GRANULE+1)*sizeof(pthread_mutex_t), ctt->mem_lim);
    fprintf(stderr, "BUCKET_SIZE=%lu, ", ctt->table_size); 
    fprintf(stderr, "COV_CUTOFF=%lu, ", ctt->th); 
    fprintf(stderr, "NUM_OCCUPIED_SLOT=%lu, ", n_used);  
    fprintf(stderr, "LOAD_FACTOR=%f, ", (double)n_used / (double)ctt->table_size);
    fprintf(stderr, "MEM_USAGE=%luB\n", ctt->table_size * (sizeof(uint16_t)+sizeof(uint64_t))+(ctt->table_size/GRANULE+1)*sizeof(pthread_mutex_t)); 

    if (ctt->table_size*(sizeof(uint16_t)+sizeof(uint64_t))+(ctt->table_size/GRANULE+1)*sizeof(pthread_mutex_t) > ctt->mem_lim) {
        ctt->val_table = (uint16_t *)my_realloc(ctt->val_table, ctt->table_size*sizeof(uint16_t));
        check_alloc(ctt->val_table);
        memset(ctt->val_table, 0, ctt->table_size*sizeof(uint16_t));
        ctt->key_table = (uint64_t *)my_realloc(ctt->key_table, ctt->table_size*sizeof(uint64_t));
        check_alloc(ctt->key_table);
        g_ctt_table_mutex = (pthread_mutex_t *)my_realloc(g_ctt_table_mutex, (ctt->table_size/GRANULE+1)*sizeof(pthread_mutex_t));
        check_alloc(g_ctt_table_mutex);
    }

    for (i = 0; i < ctt->n_thread; ++i) {
        ctt->ct[i].val_table = ctt->val_table;
        ctt->ct[i].key_table = ctt->key_table;
        ctt->ct[i].table_size = ctt->table_size;
        ctt->ct[i].index_len = ctt->index_len;
        ctt->ct[i].th = ctt->th;
        fclose(g_ctt_unstored_fp[i]);
        pthread_create(&ctt->thread[i], NULL, (void *)corrector_threads_finish_table, &ctt->ct[i]);
    }
    for (i = 0; i < ctt->n_thread; ++i) {
        pthread_join(ctt->thread[i], NULL); 
        fclose(g_ctt_stored_fp[i]);
    }
    
    my_free(g_ctt_table_mutex);
}

void corrector_threads_count(corrector_t *ct)
{
    uint64_t i;
    uint64_t key;
    uint64_t index;
    uint64_t k_len;
    kmer_t kmer;
    seq_t seq;

    k_len = ct->k_len;
    rewind(ct->count_seq_file);
    while (seq_read(&seq, ct->count_seq_file)) {
        if (seq.len < k_len)
            continue;
        seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
        seq.n_unknown = 0;
        kmer.fwd = kmer.rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer.fwd = kmer.fwd << 2 | (uint64_t)seq.base[i];
            kmer.rev = (kmer.rev | (uint64_t)(3^seq.base[k_len-i-2])) << 2;
        }
        for (i = 0; i < seq.len - k_len + 1; ++i) {
            kmer.fwd = ((kmer.fwd << 2) & ct->k_mask) | (uint64_t)seq.base[i+k_len-1];
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^seq.base[i+k_len-1]) << (2*(k_len-1)));
            if (seq.unknown_pos[seq.n_unknown] < i + k_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                continue;
            }
            key = kmer.fwd <= kmer.rev ? kmer.fwd : kmer.rev;
            index = hash64(key, ct->index_len) & (ct->table_size - 1);
            pthread_mutex_lock(&g_ctt_table_mutex[index/GRANULE]);
            if (!(ct->val_table[index])) {
                ct->val_table[index] = 1;
                ct->key_table[index] = key;
            }
            else if (ct->key_table[index] == key) {
                if (ct->val_table[index] < UINT16_MAX-1)
                    ++ct->val_table[index];
            }
            else
                fwrite(&key, sizeof(uint64_t), 1, g_ctt_unstored_fp[ct->id]);
            pthread_mutex_unlock(&g_ctt_table_mutex[index/GRANULE]);
        }
    }
}

void corrector_threads_recount(corrector_t *ct)
{
    uint64_t key;
    uint64_t index;
    FILE *tmp_fp;

    tmp_fp = tmpfile_open();
    check_tmpfile(tmp_fp);

    rewind(g_ctt_unstored_fp[ct->id]);
    while (fread(&key, sizeof(uint64_t), 1, g_ctt_unstored_fp[ct->id])) {
        index = hash64(key, ct->index_len) & (ct->table_size - 1);
        pthread_mutex_lock(&g_ctt_table_mutex[index/GRANULE]);
        if (!(ct->val_table[index])) {
            ct->val_table[index] = 1;
            ct->key_table[index] = key;
        }
        else if (ct->key_table[index] == key) {
            if (ct->val_table[index] < UINT16_MAX-1)
                ++ct->val_table[index];
        }
        else
            fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
        pthread_mutex_unlock(&g_ctt_table_mutex[index/GRANULE]);
    }
    fclose(g_ctt_unstored_fp[ct->id]);
    g_ctt_unstored_fp[ct->id] = tmp_fp;
}

void corrector_threads_finish_table(corrector_t *ct)
{
    uint64_t key;
    uint64_t index;
    uint64_t step;
    uint16_t occ;

    rewind(g_ctt_stored_fp[ct->id]);
    while (fread(&key, sizeof(uint64_t), 1, g_ctt_stored_fp[ct->id])) {
        fread(&occ, sizeof(uint16_t), 1, g_ctt_stored_fp[ct->id]);
        if (occ < ct->th)
            continue;
        index = hash64(key, ct->index_len) & (ct->table_size - 1);
        pthread_mutex_lock(&g_ctt_table_mutex[index/GRANULE]);
        if (!(ct->val_table[index])) {
            ct->val_table[index] = occ;
            ct->key_table[index] = key;
            pthread_mutex_unlock(&g_ctt_table_mutex[index/GRANULE]);
        }
        else {
            pthread_mutex_unlock(&g_ctt_table_mutex[index/GRANULE]);
            step = rehash64(key, ct->index_len);
            index = (index+step) & (ct->table_size-1);
            while (1) {
                pthread_mutex_lock(&g_ctt_table_mutex[index/GRANULE]);
                if (!ct->val_table[index]) {
                    ct->val_table[index] = occ;
                    ct->key_table[index] = key;
                    pthread_mutex_unlock(&g_ctt_table_mutex[index/GRANULE]);
                    break;
                }
                pthread_mutex_unlock(&g_ctt_table_mutex[index/GRANULE]);
                index = (index+step) & (ct->table_size-1);
            }
        }
    }
}

void corrector_threads_ungap_correct(corrector_threads_t *ctt)
{
    uint64_t i;

    for (i = 0; i < ctt->n_thread; ++i) {
        corrector_set_max_edit(&ctt->ct[i], ctt->max_edit);
        corrector_set_upper_th(&ctt->ct[i], ctt->upper_th);
        corrector_set_id(&ctt->ct[i], i);
        pthread_create(&ctt->thread[i], NULL, (void *)corrector_ungap_correct, &ctt->ct[i]);
    }
    for (i = 0; i < ctt->n_thread; ++i) {
        pthread_join(ctt->thread[i], NULL);
	}
    ctt->n_cor_read = 0;
    ctt->n_cor_base = 0;
    for (i = 0; i < ctt->n_thread; ++i) {
        ctt->n_cor_read += ctt->ct[i].n_cor_read;
        ctt->n_cor_base += ctt->ct[i].n_cor_base;
    }
}

void corrector_threads_gapped_correct(corrector_threads_t *ctt)
{
    uint64_t i;

    for (i = 0; i < ctt->n_thread; ++i) {
        corrector_set_max_edit(&ctt->ct[i], ctt->max_edit);
        corrector_set_upper_th(&ctt->ct[i], ctt->upper_th);
        corrector_set_id(&ctt->ct[i], i);
        pthread_create(&ctt->thread[i], NULL, (void *)corrector_gapped_correct, &ctt->ct[i]);
    }
    for (i = 0; i < ctt->n_thread; ++i) {
        pthread_join(ctt->thread[i], NULL);
	}
    ctt->n_cor_read = 0;
    ctt->n_cor_base = 0;
    for (i = 0; i < ctt->n_thread; ++i) {
        ctt->n_cor_read += ctt->ct[i].n_cor_read;
        ctt->n_cor_base += ctt->ct[i].n_cor_base;
    }
}

void corrector_threads_show_seq(corrector_threads_t *ctt)
{
    int64_t i;
    uint64_t j;
    uint64_t k;
    uint64_t n_seq;
    int8_t base_name[] = {'A', 'C', 'G', 'T', 'N'};
    seq_t seq;
    FILE *out;
    
    k = 0;
    for (i = 0; i < ctt->n_thread; ++i)
        rewind(ctt->ct[i].target_seq_file);
    
    n_seq = 0;
    for (i = 0; i < ctt->n_target_file; ++i) {
        if((out = fopen(ctt->fname[i], "w")) == NULL)
            fprintf(stderr, "error: cannot open \"%s\"\n", ctt->fname[i]);

        for (j = 0; j < ctt->n_target_read[i]; ++j) {
            seq_read(&seq, ctt->ct[n_seq % ctt->n_thread].target_seq_file);
            for (k = 0; k < seq.n_unknown; ++k)
                seq.base[seq.unknown_pos[k]] = 4;
            fprintf(out, ">s%lu\n", j+1);
            for (k = 0; k < seq.len; ++k)
                putc(base_name[seq.base[k]], out);
            putc('\n', out);
            ++n_seq;
        }

        fclose(out);
    }
}

void corrector_threads_set_th(corrector_threads_t *ctt, int64_t th)
{
    ctt->th = th;
}

void corrector_threads_set_upper_th(corrector_threads_t *ctt, int64_t upper_th)
{
    ctt->upper_th = upper_th;
}

void corrector_threads_set_max_edit(corrector_threads_t *ctt, double max_edit)
{
    ctt->max_edit = max_edit;
}
