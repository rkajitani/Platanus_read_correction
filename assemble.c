#include "assemble.h"

void print_assemble_usage(void)
{
    option_assemble_t def;

    option_assemble_init(&def);

    fprintf(stderr, "\nUsage platanus assemble [Options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR               : prefix of output files (default %s, length <= %u)\n", def.o, MAX_FILE_LEN);
    fprintf(stderr, "    -f FILE1 [FILE2 ...] : reads file (fasta or fastq, number <= %u)\n", MAX_FILE_NUM);
    fprintf(stderr, "    -k INT               : initial k-mer size (default %lu)\n", def.k);
    fprintf(stderr, "    -s INT               : step size of k-mer extension (>= 1, default %lu)\n", def.s);
	fprintf(stderr, "    -n INT               : initial k-mer coverage cutoff (default %lu, 0 means auto)\n", def.n);
    fprintf(stderr, "    -c INT               : minimun k-mer coverage (default %lu)\n", def.c);
    fprintf(stderr, "    -a FLOAT             : k-mer extension safety level (default %.1f)\n", def.a);
    fprintf(stderr, "    -u FLOAT             : maximum difference for bubble crush (identity, default %.1f)\n", def.u);
    fprintf(stderr, "    -d FLOAT             : maximum difference for branch cutting (coverage ratio, default %.1f)\n", def.d);
    fprintf(stderr, "    -t INT               : number of threads (<= %u, default %lu)\n", MAX_THREAD, def.t);
    fprintf(stderr, "    -m INT               : memory limit(GB, >= 1, default %llu)\n", def.m / GIBIBYTE);

    fprintf(stderr, "Outputs:\n");
    fprintf(stderr, "    PREFIX_contig.fa\n");
    fprintf(stderr, "    PREFIX_contigBubble.fa\n");
    fprintf(stderr, "    PREFIX_kmerFrq.tsv\n");

    option_assemble_destroy(&def);
}

void option_assemble_init(option_assemble_t *opt)
{

    strcpy(opt->o, "out");
    opt->k = 32;
    opt->s = 10;
	opt->n = 0;
    opt->c = 1;
    opt->a = 10.0;
    opt->u = 0.1;
    opt->d = 0.5;
    opt->t = 1;
    opt->m = 4 * GIBIBYTE;
    opt->n_f = 0;
}

void option_assemble_destroy(const option_assemble_t *opt)
{
    uint64_t i;

    for (i = 0; i < opt->n_f; ++i) 
        fclose(opt->f[i].fp);
}

uint8_t option_assemble_parse(option_assemble_t *opt, int argc, char **argv)
{
    uint64_t i;

    for (i = 2; i < argc;) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help") || !strcmp(argv[i], "--help"))
            return -1;
        else if (!strcmp(argv[i], "-o"))
            i = option_string(argc, argv, i, MAX_FILE_LEN, opt->o);
        else if (!strcmp(argv[i], "-f"))
            i = option_multi_file(argc, argv, i, opt->f, &(opt->n_f));
        else if (!strcmp(argv[i], "-k"))
            i = option_int(argc, argv, i, 1, INT64_MAX, &(opt->k));
        else if (!strcmp(argv[i], "-s"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->s));
		else if (!strcmp(argv[i], "-n"))
			i = option_int(argc, argv, i, 0, UINT16_MAX, &(opt->n));
        else if (!strcmp(argv[i], "-c"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->c));
        else if (!strcmp(argv[i], "-a"))
            i = option_float(argc, argv, i, -DBL_MAX, DBL_MAX, &(opt->a));
        else if (!strcmp(argv[i], "-u"))
            i = option_float(argc, argv, i, 0.0, 1.0, &(opt->u));
        else if (!strcmp(argv[i], "-d"))
            i = option_float(argc, argv, i, 0.0, 1.0, &(opt->d));
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

    if (opt->n_f == 0) {
        fputs("error: no input file\n", stderr);
        return -1;
    }

    return 0;
}

void assemble(const option_assemble_t *opt)
{
    char out_name[MAX_FILE_LEN + MAX_EXTENSION_NAME + 1];
    uint64_t i;
    uint64_t k = opt->k;
    uint64_t pre_k;
    uint64_t max_k;
    uint64_t cov_cut;
    uint64_t len_cut;
    const uint64_t len_step = opt->s;
    const uint64_t min_cov = opt->c;
    double ave_len = 0.0;
    double ave_cov = 0.0;
    const double min_log_p_join = log(1.0 - pow(10.0, -(opt->a)));
    const double bubble_th = opt->u;
    const double branch_th = opt->d;
    FILE *read_fp;
    FILE *contig_fp;
    FILE *sorted_key_fp;
    mcounter_t mc;
    mgraph_t mg;
    lcounter_t lc;
    lgraph_t lg;


    mcounter_init(&mc, opt->m);
    lcounter_init(&lc, opt->m);
    mgraph_init(&mg, branch_th, bubble_th);
    lgraph_init(&lg, branch_th, bubble_th);

    read_fp = tmpfile_open();
    check_tmpfile(read_fp);

    for (i = 0; i < opt->n_f; ++i) {
		if (opt->f[i].fmt == FASTA)
			read_fasta(opt->f[i].fp, read_fp);
		else
			read_fastq(opt->f[i].fp, read_fp);
	}

    if (k <= 32) {
        fprintf(stderr, "K = %lu, saving kmers from reads...\n", k);
        mcounter_save_kmer_read(&mc, k, read_fp);
		if (opt->n == 0)
			cov_cut = get_left_minimal_smooth(mc.occ_distr, mc.max_occ, SMOOTHING_WINDOW);
		else
			cov_cut = opt->n;
        ave_len = get_ave(mc.len_distr, 0, MAX_READ_LEN);
        ave_cov = get_ave(mc.occ_distr, cov_cut, mc.max_occ);
        ave_cov = ave_cov * ave_len / (ave_len - k + 1.0);
        max_k = max_kmer_len(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step);

        fprintf(stderr, "AVE_READ_LEN=%f\n", ave_len);
        show_kmer_extension(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step, cov_cut);

        strcpy(out_name, opt->o);
        sprintf(out_name, "%s_%lumerFrq.tsv", opt->o, k);
        show_distribution(mc.occ_distr, 1, mc.max_occ, out_name);

        sorted_key_fp = sorted_key_from_kmer_file(mc.kmer_fp, cov_cut);
        mcounter_load_kmer(&mc, cov_cut);
        mgraph_save(&mg, &mc, sorted_key_fp);
        fclose(sorted_key_fp);

        mcounter_destroy_table(&mc);
        mgraph_load(&mg);
        mgraph_cut_branch_iterative(&mg);
    }
    else {
        fprintf(stderr, "K = %lu, saving kmers from reads...\n", k);
        lcounter_save_kmer_read(&lc, k, read_fp);
		if (opt->n == 0)
			cov_cut = get_left_minimal_smooth(lc.occ_distr, lc.max_occ, SMOOTHING_WINDOW);
		else
			cov_cut = opt->n;
        ave_len = get_ave(lc.len_distr, 0, MAX_READ_LEN);
        ave_cov = get_ave(lc.occ_distr, cov_cut, lc.max_occ);
        ave_cov = ave_cov * ave_len / (ave_len - k + 1.0);
        max_k = max_kmer_len(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step);

        fprintf(stderr, "AVE_READ_LEN=%f\n", ave_len);
        show_kmer_extension(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step, cov_cut);

        strcpy(out_name, opt->o);
        sprintf(out_name, "%s_%lumerFrq.tsv", opt->o, k);
        show_distribution(lc.occ_distr, 1, lc.max_occ, out_name);

        sorted_key_fp = sorted_lkey_from_lmer_file(lc.kmer_fp, k, cov_cut);
        lcounter_load_kmer(&lc, cov_cut);
        lgraph_save(&lg, &lc, sorted_key_fp);
        fclose(sorted_key_fp);
        lcounter_destroy_table(&lc);
        lgraph_load(&lg);
        lgraph_cut_branch_iterative(&lg);
    }

    pre_k = k;
    if (k < max_k && k + len_step > max_k)
        k = max_k;
    else
        k += len_step;

    while (k <= max_k) {
        if (k <= 32) {
            mgraph_save_edge_kmer(&mg, &mc, k);
            fputs("extracting reads...\n", stderr);
            mcounter_load_kmer(&mc, 1);
            mcounter_extract_read(&mc, &read_fp);
            mcounter_destroy_table(&mc);
            contig_fp = mgraph_save_contig(&mg, k);
            mgraph_destroy_table(&mg);
            mcounter_load_contig(&mc, k, (ave_len - k + 1.0) / (ave_len - pre_k + 1.0), contig_fp);
            mcounter_save_additional_kmer_read(&mc, k, read_fp);
            fclose(contig_fp);
            cov_cut = decrease_cov_cut(cov_cut, ave_cov, ave_len, min_log_p_join, k, pre_k);
            if (cov_cut < min_cov)
                cov_cut = min_cov;
            fprintf(stderr, "COVERAGE_CUTOFF = %lu\n", cov_cut);
            sorted_key_fp = sorted_key_from_kmer_file(mc.kmer_fp, cov_cut);
            mcounter_load_kmer(&mc, cov_cut);
            mgraph_save(&mg, &mc, sorted_key_fp);
            fclose(sorted_key_fp);
            mcounter_destroy_table(&mc);
            mgraph_load(&mg);
            mgraph_cut_branch_iterative(&mg);
        }
        else {
            if (pre_k <= 32) {
                mgraph_save_edge_kmer(&mg, &mc, k);
                fputs("extracting reads...\n", stderr);
                mcounter_load_kmer(&mc, 1);
                mcounter_extract_read(&mc, &read_fp);
                mcounter_destroy_table(&mc);
                contig_fp = mgraph_save_contig(&mg, k);
                mgraph_destroy_table(&mg);
            }
            else {
                lgraph_save_edge_kmer(&lg, &lc, k);
                fputs("extracting reads...\n", stderr);
                lcounter_load_kmer(&lc, 1);
                lcounter_extract_read(&lc, &read_fp);
                lcounter_destroy_table(&lc);
                contig_fp = lgraph_save_contig(&lg, k);
                lgraph_destroy_table(&lg);
            }
            lcounter_load_contig(&lc, k, (ave_len - k + 1.0) / (ave_len - pre_k + 1.0), contig_fp);
            lcounter_save_additional_kmer_read(&lc, k, read_fp);
            fclose(contig_fp);
            cov_cut = decrease_cov_cut(cov_cut, ave_cov, ave_len, min_log_p_join, k, pre_k);
            if (cov_cut < min_cov)
                cov_cut = min_cov;
            fprintf(stderr, "COVERAGE_CUTOFF = %lu\n", cov_cut);
            sorted_key_fp = sorted_lkey_from_lmer_file(lc.kmer_fp, k, cov_cut);
            lcounter_load_kmer(&lc, cov_cut);
            lgraph_save(&lg, &lc, sorted_key_fp);
            fclose(sorted_key_fp);
            lcounter_destroy_table(&lc);
            lgraph_load(&lg);
            lgraph_cut_branch_iterative(&lg);
        }

        pre_k = k;
        if (k < max_k && k + len_step > max_k)
            k = max_k;
        else
            k += len_step;
    }

    k = pre_k;

    fclose(read_fp);

    if (k <= 32) {
        mgraph_straight_stats(&mg, &len_cut, &cov_cut, &ave_cov);
        fprintf(stderr, "LENGTH_CUTOFF = %lu\nCOVERAGE_CUTOFF = %lu\nAVERAGE_COVERAGE = %f\n", len_cut, cov_cut, ave_cov);
        mgraph_delete_erroneous_straight_node(&mg, len_cut, cov_cut);
    }
    else {
        lgraph_straight_stats(&lg, &len_cut, &cov_cut, &ave_cov);
        fprintf(stderr, "LENGTH_CUTOFF = %lu\nCOVERAGE_CUTOFF = %lu\nAVERAGE_COVERAGE = %f\n", len_cut, cov_cut, ave_cov);
        lgraph_delete_erroneous_straight_node(&lg, len_cut, cov_cut);
    }

    strcpy(out_name, opt->o);
    strcat(out_name, "_contigBubble.fa");
    if (k <= 32)
        contig_fp = mgraph_crush_bubble_iterative(&mg, ave_cov);
    else
        contig_fp = lgraph_crush_bubble_iterative(&lg, ave_cov);
    print_contig(contig_fp, out_name, ave_len / (ave_len - k + 1.0));

    strcpy(out_name, opt->o);
    strcat(out_name, "_contig.fa");
    if (k <= 32)
        mgraph_show_contig(&mg, ave_len / (ave_len - k + 1.0), out_name);
    else
        lgraph_show_contig(&lg, ave_len / (ave_len - k + 1.0), out_name);

    mgraph_destroy(&mg);
    lgraph_destroy(&lg);

    fputs("assemble completed!\n", stderr);
}

void assemble_mt(const option_assemble_t *opt)
{
    char out_name[MAX_FILE_LEN + MAX_EXTENSION_NAME + 1];
    uint64_t i;
    int64_t j;
    uint64_t k = opt->k;
    uint64_t pre_k;
    uint64_t max_k;
    const uint64_t n_thread = opt->t;
    const uint64_t len_step = opt->s;
    const uint64_t min_cov = opt->c;
    uint64_t cov_cut;
    uint64_t len_cut;
    double ave_len = 0.0;
    double ave_cov = 0.0;
    const double min_log_p_join = log(1.0 - pow(10.0, -(opt->a)));
    const double bubble_th = opt->u;
    const double branch_th = opt->d;
    FILE *read_fp[MAX_THREAD];
    FILE *contig_fp;
    FILE *sorted_key_fp;
    mcounter_t mc;
    mgraph_t mg;
    lcounter_t lc;
    lgraph_t lg;


    mcounter_init(&mc, opt->m);
    lcounter_init(&lc, opt->m);
    mgraph_init(&mg, branch_th, bubble_th);
    lgraph_init(&lg, branch_th, bubble_th);

    omp_set_num_threads(n_thread);

    for (i = 0; i < n_thread; ++i) {
        read_fp[i] = tmpfile_open();
        check_tmpfile(read_fp[i]);
    }

    for (i = 0; i < opt->n_f; ++i) {
		if (opt->f[i].fmt == FASTA)
			read_fasta_mt(opt->f[i].fp, read_fp, n_thread);
		else
			read_fastq_mt(opt->f[i].fp, read_fp, n_thread);
	}
    
    if (k <= 32) {
        fprintf(stderr, "K = %lu, saving kmers from reads...\n", k);
        mcounter_save_kmer_read_mt(&mc, k, read_fp, n_thread);
		if (opt->n == 0)
			cov_cut = get_left_minimal_smooth(mc.occ_distr, mc.max_occ, SMOOTHING_WINDOW);
		else
			cov_cut = opt->n;
        ave_len = get_ave(mc.len_distr, 0, MAX_READ_LEN);
        ave_cov = get_ave(mc.occ_distr, cov_cut, mc.max_occ);
        ave_cov = ave_cov * ave_len / (ave_len - k + 1.0);
        max_k = max_kmer_len(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step);

        fprintf(stderr, "AVE_READ_LEN=%f\n", ave_len);
        show_kmer_extension(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step, cov_cut);

        strcpy(out_name, opt->o);
        sprintf(out_name, "%s_%lumerFrq.tsv", opt->o, k);
        show_distribution(mc.occ_distr, 1, mc.max_occ, out_name);

        sorted_key_fp = sorted_key_from_kmer_file(mc.kmer_fp, cov_cut);
        mcounter_load_kmer(&mc, cov_cut);
        mcounter_save_table(&mc);
        mgraph_save(&mg, &mc, sorted_key_fp);
        fclose(sorted_key_fp);
        mcounter_destroy_table(&mc);
        mgraph_load(&mg);
        mgraph_cut_branch_iterative(&mg);
    }
    else {
        fprintf(stderr, "K = %lu, saving kmers from reads...\n", k);
        lcounter_save_kmer_read_mt(&lc, k, read_fp, n_thread);
		if (opt->n == 0)
			cov_cut = get_left_minimal_smooth(lc.occ_distr, lc.max_occ, SMOOTHING_WINDOW);
		else
			cov_cut = opt->n;
        ave_len = get_ave(lc.len_distr, 0, MAX_READ_LEN);
        ave_cov = get_ave(lc.occ_distr, cov_cut, lc.max_occ);
        ave_cov = ave_cov * ave_len / (ave_len - k + 1.0);
        max_k = max_kmer_len(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step);

        fprintf(stderr, "AVE_READ_LEN=%f\n", ave_len);
        show_kmer_extension(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step, cov_cut);

        strcpy(out_name, opt->o);
        sprintf(out_name, "%s_%lumerFrq.tsv", opt->o, k);
        show_distribution(lc.occ_distr, 1, lc.max_occ, out_name);

        sorted_key_fp = sorted_lkey_from_lmer_file(lc.kmer_fp, k, cov_cut);
        lcounter_load_kmer(&lc, cov_cut);
        lgraph_save(&lg, &lc, sorted_key_fp);
        fclose(sorted_key_fp);
        lcounter_destroy_table(&lc);
        lgraph_load(&lg);
        lgraph_cut_branch_iterative(&lg);
    }

    pre_k = k;
    if (k < max_k && k + len_step > max_k)
        k = max_k;
    else
        k += len_step;

    while (k <= max_k) {
        if (k <= 32) {
            mgraph_save_edge_kmer(&mg, &mc, k);
            fputs("extracting reads...\n", stderr);
            mcounter_load_kmer(&mc, 1);
#            pragma omp parallel for schedule(static, 1)
            for (j = 0; j < n_thread; ++j) {
                mcounter_extract_read(&mc, &read_fp[j]);
            }
            mcounter_destroy_table(&mc);
            contig_fp = mgraph_save_contig(&mg, k);
            mgraph_destroy_table(&mg);
            mcounter_load_contig(&mc, k, (ave_len - k + 1.0) / (ave_len - pre_k + 1.0), contig_fp);
            mcounter_save_additional_kmer_read_mt(&mc, k, read_fp, n_thread);
            fclose(contig_fp);
            cov_cut = decrease_cov_cut(cov_cut, ave_cov, ave_len, min_log_p_join, k, pre_k);
            if (cov_cut < min_cov)
                cov_cut = min_cov;
            fprintf(stderr, "COVERAGE_CUTOFF = %lu\n", cov_cut);
            sorted_key_fp = sorted_key_from_kmer_file(mc.kmer_fp, cov_cut);
            mcounter_load_kmer(&mc, cov_cut);
            mgraph_save(&mg, &mc, sorted_key_fp);
            fclose(sorted_key_fp);
            mcounter_destroy_table(&mc);
            mgraph_load(&mg);
            mgraph_cut_branch_iterative(&mg);
        }
        else {
            if (pre_k <= 32) {
                mgraph_save_edge_kmer(&mg, &mc, k);
                fputs("extracting reads...\n", stderr);
                mcounter_load_kmer(&mc, 1);
#                pragma omp parallel for schedule(static, 1)
                for (j = 0; j < n_thread; ++j) {
                    mcounter_extract_read(&mc, &read_fp[j]);
                }
                mcounter_destroy_table(&mc);
                contig_fp = mgraph_save_contig(&mg, k);
                mgraph_destroy_table(&mg);
            }
            else {
                lgraph_save_edge_kmer(&lg, &lc, k);
                fputs("extracting reads...\n", stderr);
                lcounter_load_kmer(&lc, 1);
#                pragma omp parallel for schedule(static, 1)
                for (j = 0; j < n_thread; ++j) {
                    lcounter_extract_read(&lc, &read_fp[j]);
                }
                lcounter_destroy_table(&lc);
                contig_fp = lgraph_save_contig(&lg, k);
                lgraph_destroy_table(&lg);
            }
            lcounter_load_contig(&lc, k, (ave_len - k + 1.0) / (ave_len - pre_k + 1.0), contig_fp);
            lcounter_save_additional_kmer_read_mt(&lc, k, read_fp, n_thread);
            fclose(contig_fp);
            cov_cut = decrease_cov_cut(cov_cut, ave_cov, ave_len, min_log_p_join, k, pre_k);
            if (cov_cut < min_cov)
                cov_cut = min_cov;
            fprintf(stderr, "COVERAGE_CUTOFF = %lu\n", cov_cut);
            sorted_key_fp = sorted_lkey_from_lmer_file(lc.kmer_fp, k, cov_cut);
            lcounter_load_kmer(&lc, cov_cut);
            lgraph_save(&lg, &lc, sorted_key_fp);
            fclose(sorted_key_fp);
            lcounter_destroy_table(&lc);
            lgraph_load(&lg);
            lgraph_cut_branch_iterative(&lg);
        }

        pre_k = k;
        if (k < max_k && k + len_step > max_k)
            k = max_k;
        else
            k += len_step;
    }

    k = pre_k;

    for (i = 0; i < n_thread; ++i)
        fclose(read_fp[i]);

    if (k <= 32) {
        mgraph_straight_stats(&mg, &len_cut, &cov_cut, &ave_cov);
        fprintf(stderr, "LENGTH_CUTOFF = %lu\nCOVERAGE_CUTOFF = %lu\nAVERAGE_COVERAGE = %f\n", len_cut, cov_cut, ave_cov);
        mgraph_delete_erroneous_straight_node(&mg, len_cut, cov_cut);
    }
    else {
        lgraph_straight_stats(&lg, &len_cut, &cov_cut, &ave_cov);
        fprintf(stderr, "LENGTH_CUTOFF = %lu\nCOVERAGE_CUTOFF = %lu\nAVERAGE_COVERAGE = %f\n", len_cut, cov_cut, ave_cov);
        lgraph_delete_erroneous_straight_node(&lg, len_cut, cov_cut);
    }

    strcpy(out_name, opt->o);
    strcat(out_name, "_contigBubble.fa");
    if (k <= 32)
        contig_fp = mgraph_crush_bubble_iterative(&mg, ave_cov);
    else
        contig_fp = lgraph_crush_bubble_iterative(&lg, ave_cov);
    print_contig(contig_fp, out_name, ave_len / (ave_len - k + 1.0));

    strcpy(out_name, opt->o);
    strcat(out_name, "_contig.fa");
    if (k <= 32)
        mgraph_show_contig(&mg, ave_len / (ave_len - k + 1.0), out_name);
    else
        lgraph_show_contig(&lg, ave_len / (ave_len - k + 1.0), out_name);

    mgraph_destroy(&mg);
    lgraph_destroy(&lg);

    fputs("assemble completed!\n", stderr);
}

void show_distribution(const uint64_t *distr, const uint64_t min, const uint64_t max, const char *out_name)
{
    uint64_t i;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    for (i = min; i <= max; ++i)
        fprintf(out, "%lu\t%lu\n", i, distr[i]);

    fclose(out);
}

void show_kmer_extension(const double min_log_p_join, const double ave_cov, const double ave_len, const uint64_t min_cov, uint64_t k, uint64_t const len_step, uint64_t cov_cut)
{
    const uint64_t max_k = max_kmer_len(min_log_p_join, ave_cov, ave_len, min_cov, k, len_step);
    uint64_t pre_k;


    fprintf(stderr, "\nKMER_EXTENSION:\n");
    fprintf(stderr, "K=%lu", k);
    fprintf(stderr, ", KMER_COVERAGE=%.2f (>= %lu)", ave_cov * (ave_len - k + 1.0) / ave_len, cov_cut);
    fprintf(stderr, ", COVERAGE_CUTOFF=%lu\n", cov_cut);

    pre_k = k;
    if (k < max_k && k + len_step > max_k)
        k = max_k;
    else
        k += len_step;

    while (k <= max_k) {
        cov_cut = decrease_cov_cut(cov_cut, ave_cov, ave_len, min_log_p_join, k, pre_k);
        if (cov_cut < min_cov)
            cov_cut = min_cov;

        fprintf(stderr, "K=%lu", k);
        fprintf(stderr, ", KMER_COVERAGE=%.2f", ave_cov * (ave_len - k + 1.0) / ave_len);
        fprintf(stderr, ", COVERAGE_CUTOFF=%lu", cov_cut);
        fprintf(stderr, ", PROB_SPLIT=10e%.6f\n", log10(1.0 - exp(log_prob_join(cov_cut, ave_cov, ave_len, k, pre_k))));

        pre_k = k;
        if (k < max_k && k + len_step > max_k)
            k = max_k;
        else
            k += len_step;
    }
    putc('\n', stderr);
}

double decrease_coverage(const double cov, const double ave_len, const uint64_t large_k, const uint64_t small_k)
{
    return cov * (ave_len - large_k + 1.0) / (ave_len - small_k + 1.0);
}

double log_prob_join(const uint64_t cov_cut, const double ave_cov, const double ave_len, const uint64_t large_k, const uint64_t small_k)
{
    int64_t i;
    int64_t j;
    double p;
    double s;
    double tmp_ave_cov = ave_cov  * (ave_len - large_k + 1.0) / ave_len;

    s = 0.0;
    for (i = 0; i < cov_cut; ++i) {
        p = 0.0;
        for (j = 1; j <= i; ++j)
            p += log(tmp_ave_cov) - log((double)j);
        s += exp(p);
    }
    s = exp(-tmp_ave_cov + log(s));

    return ((large_k - small_k) + 1.0)*(-s);
}

uint64_t decrease_cov_cut(uint64_t cov_cut, double ave_cov, double ave_len, double min_log_p_join, uint64_t large_k, uint64_t small_k)
{
    uint64_t i;

    if (cov_cut <= 1)
        return 1;

    for (i = cov_cut; i > 1; --i) {
        if (log_prob_join(i, ave_cov, ave_len, large_k, small_k) > min_log_p_join)
            break;
    }
    return i;
}

uint64_t max_kmer_len(const double min_log_p_join, const double ave_cov, const double ave_len, const uint64_t min_cov, uint64_t k_len, const uint64_t len_step)
{
    uint64_t i;

    for (k_len += len_step; k_len < (uint64_t)(ave_len + 0.5); k_len += len_step) {
        if (log_prob_join(min_cov, ave_cov, ave_len, k_len, k_len - len_step) < min_log_p_join)
            break;
    }
    k_len -= len_step;

    for (i = 1; i < len_step; ++i) {
        if (log_prob_join(min_cov, ave_cov, ave_len, k_len + i, k_len) < min_log_p_join)
            break;
    }
    --i;

    return k_len + i;
}
