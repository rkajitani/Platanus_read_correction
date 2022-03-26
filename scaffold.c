#include "scaffold.h"

void print_scaffold_usage(void)
{
    option_scaffold_t def = {};

    option_scaffold_init(&def);

    fprintf(stderr, "\nUsage: platanus scaffold [Options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR                             : prefix of output file (default %s, length <= %u)\n", def.o, MAX_FILE_LEN);
    fprintf(stderr, "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n");
    fprintf(stderr, "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n");
    fprintf(stderr, "    -ip{INT} PAIR1 [PAIR2 ...]         : INT = lib_id, PAIR1 = inward_pair_file (reads in 1 file, fasta or fastq)\n");
    fprintf(stderr, "    -op{INT} PAIR1 [PAIR2 ...]         : INT = lib_id, PAIR1 = outward_pair_file (reads in 1 file, fasta or fastq)\n");
    fprintf(stderr, "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : INT = lib_id, (FWD1, REV1) = inward_pair_files (reads in 2 files, fasta or fastq)\n");
    fprintf(stderr, "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : INT = lib_id, (FWD1, REV1) = outward_pair_files (reads in 2 files, fasta or fastq)\n");
    fprintf(stderr, "    -n{INT1} INT2                      : INT1 = lib_id, INT2 = minimum_insert_size\n");
    fprintf(stderr, "    -a{INT1} INT2                      : INT1 = lib_id, INT2 = average_insert_size\n");
    fprintf(stderr, "    -d{INT1} INT2                      : INT1 = lib_id, INT2 = standard_deviation_of_insert_size\n");
    fprintf(stderr, "    -s INT                             : mapping seed length (default %ld)\n", def.s);
    fprintf(stderr, "    -v INT                             : minimum overlap length (default %ld)\n", def.v);
    fprintf(stderr, "    -l INT                             : minimum number of link (default %ld)\n", def.l);
    fprintf(stderr, "    -u FLOAT                           : maximum difference for bubble crush (identity, default %.1f)\n", def.u);
    fprintf(stderr, "    -t INT                             : number of threads (<= %ld, default 1)\n", def.t);
    fprintf(stderr, "    -m INT                             : memory limit (GB, >= 1, default %llu)\n", def.m / GIBIBYTE);

    fprintf(stderr, "Outputs:\n");
    fprintf(stderr, "    PREFIX_scaffold.fa\n");
    fprintf(stderr, "    PREFIX_scaffoldBubble.fa\n");

    option_scaffold_destroy(&def);
}

void option_scaffold_init(option_scaffold_t *opt)
{
	int i;

    strcpy(opt->o, "out");
    opt->s = 32;
    opt->v = 32;
    opt->l = 3;
    opt->u = 0.1;
    opt->m = 4 * GIBIBYTE;
    opt->t = 1;
    opt->n_c = 0;
    opt->n_pair = 0;
    opt->n_b = 0;

	for (i = 0; i < MAX_FILE_NUM; ++i) {
		opt->n[i] = 0;
		opt->a[i] = -1;
		opt->d[i] = -1;
	}
}

void option_scaffold_destroy(const option_scaffold_t *opt)
{
    uint64_t i;

    for (i = 0; i < opt->n_c; ++i)
        fclose(opt->c[i]);
    for (i = 0; i < opt->n_b; ++i)
        fclose(opt->b[i]);
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

uint8_t option_scaffold_parse(option_scaffold_t *opt, int argc, char **argv)
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
        else if (!strcmp(argv[i], "-s"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->s));
        else if (!strcmp(argv[i], "-v"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->v));
        else if (!strcmp(argv[i], "-l"))
            i = option_int(argc, argv, i, 0, INT64_MAX, &(opt->l));
        else if (!strcmp(argv[i], "-u"))
            i = option_float(argc, argv, i, 0.0, 1.0, &(opt->u));
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
        else if (strstr(argv[i], "-n") == argv[i])
            i = option_indexed_int(argc, argv, i, opt->n);
        else if (strstr(argv[i], "-a") == argv[i])
            i = option_indexed_int(argc, argv, i, opt->a);
        else if (strstr(argv[i], "-d") == argv[i])
            i = option_indexed_int(argc, argv, i, opt->d);
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

void scaffold_mt(option_scaffold_t *opt)
{
    char out_name[MAX_FILE_LEN + MAX_EXTENSION_NAME + 1];
    const int64_t n_thread = opt->t;
    int64_t i;
    int64_t j;
	int64_t pre_n_node;
    const int64_t seed_len = opt->s;
    int64_t key_len = min_64(seed_len, SCA_HASH_OVERLAP);
    const int64_t min_overlap = opt->v;
    int64_t link_th;
    const int64_t min_link = opt->l;
    int64_t n_lib;
    int64_t n_pair[MAX_FILE_NUM] = {0};
    const double bubble_th = opt->u;
    double ave_cov;
    contig_t con;
    contig_t bub;
    hetero_mapper_t mp;
    scaffolder_t sc;
    seqlib_t *lib[MAX_FILE_NUM];
    seqlib_input_t *pair[MAX_FILE_NUM];
    FILE *bubble_fp;

    if (key_len > 32)
        key_len = 32;

    bubble_fp = tmpfile_open();

    check_tmpfile(bubble_fp);
    contig_init(&con);
    contig_init(&bub);

    hetero_mapper_init(&mp, seed_len, key_len, opt->m);
    scaffolder_init(&sc, seed_len, min_overlap, opt->m);

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
    for (i = -2; i < n_lib; ++i) {
        if (i == -2) {
            for (j = 0; j < opt->n_c; ++j)
                contig_read_fasta_cov(&con, opt->c[j]);
            contig_load_seq(&con);
            contig_load_cov(&con);
            ave_cov = contig_ave_cov(&con, contig_median_len(&con));
            mapper_set_contig(&(mp.cmp), &con);
            mapper_make_table(&(mp.cmp));
        }
        else if (i == -1 && opt->n_b > 0) {
            for (j = 0; j < opt->n_b; ++j)
                contig_read_fasta(&bub, opt->b[j]);
            contig_load_seq(&bub);
            mapper_set_contig(&(mp.bmp), &bub);
            mapper_make_table(&(mp.bmp));
        }
        else if (i >= 0) {
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

    fprintf(stderr, "CONTIG_AVERAGE_COVERAGE = %f\n", ave_cov);
    
    hetero_mapper_merge_mt(&mp, n_thread);

    for (i = 0; i < n_lib; ++i) {
        fprintf(stderr, "[LIBRARY %ld]\n", i + 1);
        lib[i]->ave_len = (double)(lib[i]->total_len) / (2 * lib[i]->n_pair) + 0.5;
        hetero_mapper_map_pair_mt(&mp, lib[i], opt->n[i], n_thread);
        seqlib_estimate_ins(lib[i]);
        fprintf(stderr, "Insert size estimated: AVE = %ld, SD = %ld\n", lib[i]->ave_ins, lib[i]->sd_ins);

		if (opt->a[i] >= 0 || opt->d[i] >= 0) {
			if (opt->a[i] >= 0) {
				lib[i]->ave_ins = opt->a[i];
				fprintf(stderr, "Average insert size spesified: AVE = %ld\n", lib[i]->ave_ins);
			}
			if (opt->d[i] >= 0)
				lib[i]->sd_ins = opt->d[i];
			else
				lib[i]->sd_ins = (int64_t)((double)lib[i]->ave_ins / 10.0 + 0.5);
			fprintf(stderr, "SD insert size spesified: SD = %ld\n", lib[i]->sd_ins);
		}
    }

    qsort(lib, n_lib, sizeof(seqlib_t *), cmp_seqlib_ptr);
    scaffolder_set_seqlib(&sc, lib[0]);
    scaffolder_save_overlap(&sc, &(mp.cmp), SCA_HASH_OVERLAP, lib[0]->sd_ins * MIN_TOL_FACTOR);
    scaffolder_init_sca(&sc, &con, mp.cmp.seed_len, ave_cov);
    hetero_mapper_destroy(&mp);
    contig_destroy(&bub);
    scaffolder_load_overlap(&sc);

	pre_n_node = sc.n_node;
	while (1) {
		for (i = 0; i < n_lib; ++i) {
			fprintf(stderr, "[LIBRARY %ld]\nAVE_INS = %ld, SD_INS = %ld\n", i + 1, lib[i]->ave_ins, lib[i]->sd_ins);
			scaffolder_set_seqlib(&sc, lib[i]);
			for (j = MIN_TOL_FACTOR; j <= MAX_TOL_FACTOR; ++j) {
				scaffolder_set_tol(&sc, j * (lib[i])->sd_ins);
				fprintf(stderr, "TOLERENCE = %ld\n", sc.tol);
				link_th = scaffolder_estimate_link(&sc);
				if (link_th < min_link)
					link_th = min_link;
				scaffolder_link(&sc, link_th);
				scaffolder_make_graph(&sc);
				scaffolder_crush_hetero_bubble(&sc, &bubble_fp);
				scaffolder_crush_bubble_iterative(&sc, bubble_th, &bubble_fp);
				scaffolder_delete_erroneous_edge_iterative(&sc);
				if (i > 0) {
					scaffolder_delete_repeat_edge(&sc);
					scaffolder_split(&sc);
					scaffolder_link(&sc, link_th);
					scaffolder_make_graph(&sc);
					scaffolder_crush_hetero_bubble(&sc, &bubble_fp);
					scaffolder_crush_bubble_iterative(&sc, bubble_th, &bubble_fp);
					scaffolder_delete_erroneous_edge_iterative(&sc);
				}
				scaffolder_detect_repeat(&sc);
				scaffolder_scaffold(&sc);

				scaffolder_link(&sc, link_th);
				scaffolder_make_graph(&sc);
				scaffolder_crush_hetero_bubble(&sc, &bubble_fp);
				scaffolder_crush_bubble_iterative(&sc, bubble_th, &bubble_fp);
				scaffolder_detect_repeat(&sc);
				scaffolder_scaffold(&sc);
			}
		}

		for (i = 0; i < n_lib; ++i) {
			scaffolder_set_seqlib(&sc, lib[i]);
			link_th = scaffolder_estimate_link(&sc);
			if (link_th <= min_link)
				continue;
			fprintf(stderr, "[LIBRARY %ld]\nAVE_INS = %ld, SD_INS = %ld\n", i + 1, lib[i]->ave_ins, lib[i]->sd_ins);
			for (j = MIN_TOL_FACTOR; j <= MAX_TOL_FACTOR; ++j) {
				scaffolder_set_tol(&sc, j * (lib[i])->sd_ins);
				fprintf(stderr, "TOLERENCE = %ld\n", sc.tol);
				scaffolder_link(&sc, min_link);
				scaffolder_make_graph(&sc);
				scaffolder_crush_hetero_bubble(&sc, &bubble_fp);
				scaffolder_crush_bubble_iterative(&sc, bubble_th, &bubble_fp);
				scaffolder_delete_erroneous_edge_iterative(&sc);
				scaffolder_delete_repeat_edge(&sc);
				scaffolder_detect_repeat(&sc);
				scaffolder_scaffold(&sc);

				scaffolder_link(&sc, min_link);
				scaffolder_make_graph(&sc);
				scaffolder_crush_hetero_bubble(&sc, &bubble_fp);
				scaffolder_crush_bubble_iterative(&sc, bubble_th, &bubble_fp);
				scaffolder_detect_repeat(&sc);
				scaffolder_scaffold(&sc);
			}
		}

		if (pre_n_node <= sc.n_node)
			break;
		
		pre_n_node = sc.n_node;
	}

    strcpy(out_name, opt->o);
    strcat(out_name, "_scaffold.fa");
    scaffolder_cut_and_print_seq(&sc, MIN_SCA_LEN, out_name);

    strcpy(out_name, opt->o);
    strcat(out_name, "_scaffoldBubble.fa");
    print_scaffold_bubble(bubble_fp, out_name);
    fclose(bubble_fp);

    scaffolder_destroy(&sc);
    contig_destroy(&con);

    for (i = 0; i < n_lib; ++i) {
        seqlib_destroy(lib[i]);
        my_free(lib[i]);
    }

    fputs("scaffold completed!\n", stderr);
}

int64_t contig_median_len(const contig_t *con)
{
    int64_t i;
    int64_t median;
	varr64_t buf;

	varr64_init(&buf);
    for (i = 0; i < con->n_seq; ++i) {
		varr64_push(&buf, con->seq[i].len);
	}
    qsort(buf.ele, buf.size, sizeof(uint64_t), cmp_int64);
	median = buf.ele[buf.size / 2];
	varr64_destroy(&buf);
	return median;
}

int64_t contig_n50_len(const contig_t *con)
{
    int64_t i;
    int64_t total = 0;
    int64_t n = 0;
	varr64_t buf;

	varr64_init(&buf);
    for (i = 0; i < con->n_seq; ++i) {
		total += con->seq[i].len;
		varr64_push(&buf, con->seq[i].len);
	}
    qsort(buf.ele, buf.size, sizeof(uint64_t), cmp_int64);
    for (i = 0; i < con->n_seq; ++i) {
		n += buf.ele[i];
		if (n > total / 2)
			break;
	}
	n = buf.ele[i];
	varr64_destroy(&buf);
	return n;
}

double contig_ave_cov(const contig_t *con, int64_t min_len)
{
    int64_t i;
    int64_t sum;
    int64_t num;

    sum = num = 0;
    for (i = 0; i < con->n_seq; ++i) {
		if (con->seq[i].len >= min_len) {
			sum += con->cov[i] * con->seq[i].len;
			num += con->seq[i].len;
		}
    }
    return (double)sum / num;

/*
    int64_t i;
    int64_t q1;
    int64_t q3;
    int64_t sum;
    int64_t num;
	varr64_t buf;

	varr64_init(&buf);
    for (i = 0; i < con->n_seq; ++i)
		varr64_push(&buf, con->cov[i]);
    qsort(buf.ele, buf.size, sizeof(uint64_t), cmp_int64);
	q1 = buf.ele[buf.size / 4];
	q3 = buf.ele[buf.size * 3 / 4];
	varr64_destroy(&buf);

    sum = num = 0;
    for (i = 0; i < con->n_seq; ++i) {
		if (con->cov[i] >= q1 && con->cov[i] <= q3) {
			sum += con->cov[i] * con->seq[i].len;
			num += con->seq[i].len;
		}
    }
    return (double)sum / num;
*/
}

void scaffolder_init(scaffolder_t *sc, const int64_t s, const int64_t o, const int64_t mem_lim)
{
    memset(sc, 0, sizeof(scaffolder_t));
    sc->seed_len = s;
    sc->min_overlap = o;
    sc->mem_lim = mem_lim;
}

void scaffolder_destroy_graph(scaffolder_t *sc)
{
    if (sc->node != NULL)
        my_free(sc->node);
    sc->mem_usage -=  sc->n_node * sizeof(sc_node_t);
    sc->node = NULL;

    if (sc->edge_pool != NULL )
        my_free(sc->edge_pool);
    sc->mem_usage -=  sc->edge_pool_size * sizeof(sc_edge_t);
    sc->edge_pool = NULL;

    if (sc->bkdn_pool != NULL )
        my_free(sc->bkdn_pool);
    sc->mem_usage -=  sc->bkdn_pool_size * sizeof(int64_t);
    sc->bkdn_pool = NULL;

    if (sc->con_pool != NULL )
        my_free(sc->con_pool);
    sc->mem_usage -=  sc->con_pool_size * sizeof(scaffold_part_t);
    sc->con_pool = NULL;
}

void scaffolder_destroy(scaffolder_t *sc)
{
    scaffolder_destroy_graph(sc);

    if (sc->con_inc != NULL)
        my_free(sc->con_inc);
    if (sc->ol_table != NULL)
        my_free(sc->ol_table);

    if (sc->con_fp != NULL)
        fclose(sc->con_fp);
    if (sc->link_fp != NULL)
        fclose(sc->link_fp);
    if (sc->ol_fp != NULL)
        fclose(sc->ol_fp);

    memset(sc, 0, sizeof(scaffolder_t));
}

void scaffolder_save_overlap(scaffolder_t *sc, const mapper_t *mp, const int64_t hash_overlap, const int64_t len_cutoff)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t buf_size;
    position_t *buf;
    bool is_l_unknown;
    bool is_r_unknown;
    kmer_t lkmer;
    kmer_t rkmer;
    overlap_rec_t ol;
    overlap_rec_t pre_ol;

    fprintf(stderr, "saving overlaps... (LEN_CUTOFF=%ld)\n", len_cutoff);

    sc->hash_overlap = hash_overlap;

    buf = (position_t *)my_malloc(mp->max_occ * 2 * sizeof(position_t));

    if (sc->ol_fp != NULL)
        fclose(sc->ol_fp);
    sc->ol_fp = tmpfile_open();
    check_tmpfile(sc->ol_fp);

    sc->n_con = mp->n_seq;
    sc->seq_pool_size = mp->seq_pool_size;
    sc->seq_pool = mp->seq_pool;
    sc->con = mp->seq;
    sc->seed_len = mp->seed_len;

    for (i = 0; i < sc->n_con; ++i) {
        if (sc->con[i].len < len_cutoff)
            continue;
        lkmer.fwd = lkmer.rev = rkmer.fwd = rkmer.rev = 0;
        is_l_unknown = is_r_unknown = false;
        for (j = 0;    j < mp->key_len; ++j) {
            if (sc->con[i].seq[j] == 4)
                is_l_unknown = true;
            if (sc->con[i].seq[sc->con[i].len - j - 1] == 4)
                is_r_unknown = true;
            lkmer.fwd = (lkmer.fwd << 2) | (uint64_t)sc->con[i].seq[j];
            lkmer.rev = (lkmer.rev << 2) | (uint64_t)(3 ^ sc->con[i].seq[mp->key_len - j - 1]);
            rkmer.fwd = (rkmer.fwd << 2) | (uint64_t)sc->con[i].seq[sc->con[i].len - mp->key_len + j];
            rkmer.rev = (rkmer.rev << 2) | (uint64_t)(3 ^ sc->con[i].seq[sc->con[i].len - j - 1]);
        }

        if (!is_l_unknown) {
            buf_size = mapper_map_seed(mp, &lkmer, buf);
            pre_ol.id1 = 0;
            for (j = 0; j < buf_size; ++j) {
                if (abs(buf[j].id) - 1 <= i || (buf[j].id - 1 == i && buf[j].ofst == 0))
                    continue;

                if (buf[j].id > 0) {
                    ol.len = sc->con[buf[j].id - 1].len - buf[j].ofst;
                    if (ol.len < sc->hash_overlap || ol.len > sc->con[i].len || sc->con[buf[j].id - 1].len < len_cutoff)
                        continue;
                    for (k = mp->key_len; k < ol.len; ++k) {
                        if (sc->con[i].seq[k] != sc->con[buf[j].id - 1].seq[buf[j].ofst + k])
                            break;
                    }
                }
                else {
                    ol.len = buf[j].ofst + mp->key_len;
                    if (ol.len < sc->hash_overlap || ol.len > sc->con[i].len || sc->con[-(buf[j].id) - 1].len < len_cutoff)
                        continue;
                    for (k = mp->key_len; k < ol.len; ++k) {
                        if (sc->con[i].seq[k] != (3 ^ sc->con[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1]))
                            break;
                    }
                }
                if (k < ol.len)
                    continue;

                ol.id1 = -(i + 1);
                ol.id2 = -(buf[j].id);

                if (ol.id1 == pre_ol.id1 && ol.id2 == pre_ol.id2) {
                    if (ol.len > pre_ol.len) 
                        pre_ol.len = ol.len;
                }
                else {
                    if (pre_ol.id1 != 0)
                        fwrite(&pre_ol, sizeof(overlap_rec_t), 1, sc->ol_fp);
                    pre_ol = ol;
                }
            }
            if (pre_ol.id1 != 0)
                fwrite(&pre_ol, sizeof(overlap_rec_t), 1, sc->ol_fp);
        }

        if (!is_r_unknown) {
            buf_size = mapper_map_seed(mp, &rkmer, buf);
            pre_ol.id1 = 0;
            for (j = 0; j < buf_size; ++j) {
                if (abs(buf[j].id) - 1 <= i || (buf[j].id - 1 == i && buf[j].ofst == sc->con[i].len - mp->key_len))
                    continue;

                if (buf[j].id > 0) {
                    ol.len = buf[j].ofst + mp->key_len;
                    if (ol.len < sc->hash_overlap || ol.len > sc->con[i].len || sc->con[buf[j].id - 1].len < len_cutoff)
                        continue;
                    for (k = mp->key_len; k < ol.len; ++k) {
                        if (sc->con[i].seq[sc->con[i].len - k - 1] != sc->con[buf[j].id - 1].seq[buf[j].ofst + mp->key_len -  k - 1])
                            break;
                    }
                }
                else {
                    ol.len = sc->con[-(buf[j].id) - 1].len - buf[j].ofst;
                    if (ol.len < sc->hash_overlap || ol.len > sc->con[i].len || sc->con[-(buf[j].id) - 1].len < len_cutoff)
                        continue;
                    for (k = mp->key_len; k < ol.len; ++k) {
                        if (sc->con[i].seq[sc->con[i].len - k - 1] != (3 ^ sc->con[-(buf[j].id) - 1].seq[buf[j].ofst + k]))
                            break;
                    }
                }
                if (k < ol.len)
                    continue;

                ol.id1 = i + 1;
                ol.id2 = buf[j].id;

                if (ol.id1 == pre_ol.id1 && ol.id2 == pre_ol.id2) {
                    if (ol.len > pre_ol.len) 
                        pre_ol.len = ol.len;
                }
                else {
                    if (pre_ol.id1 != 0)
                        fwrite(&pre_ol, sizeof(overlap_rec_t), 1, sc->ol_fp);
                    pre_ol = ol;
                }
            }
            if (pre_ol.id1 != 0)
                fwrite(&pre_ol, sizeof(overlap_rec_t), 1, sc->ol_fp);
        }
    }

    my_free(buf);
}

void scaffolder_load_overlap(scaffolder_t *sc)
{
    int64_t i;
    int64_t n_ol;
    overlap_rec_t ol;

    fputs("loading overlaps...\n", stderr);

    n_ol = 0;
    rewind(sc->ol_fp);
    while (fread(&ol, sizeof(overlap_rec_t), 1, sc->ol_fp))
        ++n_ol;


    sc->ol_table_size = 1;
    sc->ol_index_len = 0;
    for (i = 0; i < 64; ++i) {
        sc->ol_table_size *= 2;
        ++sc->ol_index_len;
        if (sc->ol_table_size * MAX_LOAD_FACTOR > n_ol)
            break;
    }
    sc->mem_usage += sc->ol_table_size*sizeof(overlap_rec_t);
    check_mem_usage(sc->mem_usage, sc->mem_lim);
    fprintf(stderr, "MEM_USAGE=%ld\n", sc->mem_usage);

    sc->ol_table = my_calloc(sc->ol_table_size, sizeof(overlap_rec_t));

    rewind(sc->ol_fp);
    while (fread(&ol, sizeof(overlap_rec_t), 1, sc->ol_fp))
        scaffolder_insert_overlap(sc, ol);
}

void scaffolder_insert_overlap(scaffolder_t *sc, overlap_rec_t ol)
{
    uint64_t key;
    int64_t index;
    int64_t step;

    key = (int64_t)ol.id1 * (int64_t)ol.id2;
    index = hash64(key, sc->ol_index_len) & (sc->ol_table_size-1);
    if (overlap_rec_is_emp(sc->ol_table[index]))
        sc->ol_table[index] = ol;
    else {
        step = rehash64(key, sc->ol_index_len);
        index = (index+step) & (sc->ol_table_size-1);
        while (!(overlap_rec_is_emp(sc->ol_table[index]))) {
            index = (index+step) & (sc->ol_table_size-1);
        }
        sc->ol_table[index] = ol;
    }
}

int64_t scaffolder_get_overlap(scaffolder_t *sc, int64_t id1, int64_t id2)
{
    uint64_t key;
    int64_t index;
    int64_t step;

    if (abs64(id1) > abs64(id2)) {
        index = id1;
        id1 = -id2;
        id2 = -index;
    }

    key = id1 * id2;
    index = hash64(key, sc->ol_index_len) & (sc->ol_table_size-1);
    if (!(overlap_rec_is_emp(sc->ol_table[index])) && sc->ol_table[index].id1 == id1 && sc->ol_table[index].id2 == id2)
        return sc->ol_table[index].len;
    else {
        step = rehash64(key, sc->ol_index_len);
        index = (index+step) & (sc->ol_table_size-1);
        while (!(overlap_rec_is_emp(sc->ol_table[index]))) {
            if (sc->ol_table[index].id1 == id1 && sc->ol_table[index].id2 == id2)
                return sc->ol_table[index].len;
            index = (index+step) & (sc->ol_table_size-1);
        }
        return scaffolder_get_short_overlap(sc, id1, id2);
    }
}

int64_t scaffolder_get_short_overlap(scaffolder_t *sc, int64_t id1, int64_t id2)
{
    int64_t i;
    int64_t j;

    if (id1 > 0) {
        if (id2 > 0) {
            for (i = sc->hash_overlap - 1; i >= sc->min_overlap; --i) {
                for (j = 0; j < i; ++j) {
                    if (sc->con[id1 - 1].seq[sc->con[id1 - 1].len - j - 1] != sc->con[id2 - 1].seq[i - j - 1])
                        break;
                }
                if (j == i)
                    return i;
            }
        }
        else {
            for (i = sc->hash_overlap - 1; i >= sc->min_overlap; --i) {
                for (j = 0; j < i; ++j) {
                    if (sc->con[id1 - 1].seq[sc->con[id1 - 1].len - j - 1] != (3^(sc->con[-id2 - 1].seq[sc->con[-id2 - 1].len - i + j])))
                        break;
                }
                if (j == i)
                    return i;
            }
        }
    }
    else {
        if (id2 > 0) {
            for (i = sc->hash_overlap - 1; i >= sc->min_overlap; --i) {
                for (j = 0; j < i; ++j) {
                    if ((3^(sc->con[-id1 - 1].seq[j])) != sc->con[id2 - 1].seq[i - j - 1])
                        break;
                }
                if (j == i)
                    return i;
            }
        }
        else {
            for (i = sc->hash_overlap - 1; i >= sc->min_overlap; --i) {
                for (j = 0; j < i; ++j) {
                    if (sc->con[-id1 - 1].seq[j] != sc->con[-id2 - 1].seq[sc->con[-id2 - 1].len - i + j])
                        break;
                }
                if (j == i)
                    return i;
            }
        }
    }

    return 0;
}
                
int64_t scaffolder_get_sca_overlap(scaffolder_t *sc, int64_t id1, int64_t id2)
{
    if (id1 > 0)
        id1 = sc->node[id1 - 1].con[sc->node[id1 - 1].n_con - 1].id;
    else
        id1 = -(sc->node[-id1 - 1].con[0].id);

    if (id2 > 0)
        id2 = sc->node[id2 - 1].con[0].id;
    else
        id2 = -(sc->node[-id2 - 1].con[sc->node[-id2 - 1].n_con - 1].id);
    
    return scaffolder_get_overlap(sc, id1, id2);
}

int64_t scaffolder_estimate_link(scaffolder_t *sc)
{
    return expected_link(sc->lib, (double)sc->seq_pool_size / sc->n_con, (double)sc->seq_pool_size / sc->n_con, (double)(sc->lib->ave_ins) - 2.0 * sc->lib->ave_len);
}

void scaffolder_init_sca(scaffolder_t *sc, contig_t *con, const int64_t seed_len, const double ave_cov)
{
    int64_t i;

    sc->seed_len = seed_len;
    sc->n_con = con->n_seq;
    sc->con = con->seq;
    sc->cov = con->cov;
    sc->seq_pool_size = con->seq_pool_size;
    sc->seq_pool = con->seq_pool;
    sc->mem_usage += con->mem_usage;
	sc->ave_cov = ave_cov;

    sc->n_node = sc->con_pool_size = sc->n_con;
    sc->mem_usage += sc->n_node * (sizeof(sc_node_t) + sizeof(scaffold_part_t));
    check_mem_usage(sc->mem_usage, sc->mem_lim);
    sc->node = (sc_node_t *)my_calloc(sc->n_node, sizeof(sc_node_t));
    sc->con_pool = (scaffold_part_t *)my_calloc(sc->con_pool_size, sizeof(scaffold_part_t));

    for (i = 0; i < sc->n_node; ++i) {
        sc->node[i].n_con = 1;
        sc->node[i].len = sc->con_pool[i].ed = sc->con[i].len;
        sc->node[i].con = &(sc->con_pool[i]);
        sc->con_pool[i].id = i + 1;
    }

    sc->mem_usage += sc->n_con * sizeof(position_t);
    check_mem_usage(sc->mem_usage, sc->mem_lim);
    sc->con_inc = (position_t *)my_malloc(sc->n_con * sizeof(position_t));

    for (i = 0; i < sc->n_con; ++i) {
        sc->con_inc[i].id = i + 1;
        sc->con_inc[i].ofst = 0;
    }
}

void scaffolder_set_seqlib(scaffolder_t *sc, seqlib_t *lib)
{
    sc->lib = lib;
}

void scaffolder_set_tol(scaffolder_t *sc, int64_t tol)
{
    sc->tol = tol;
}

void scaffolder_link(scaffolder_t *sc, int64_t min_link)
{
    int64_t i;
    int64_t j;
    int64_t len_cutoff;
    int64_t num;
    int64_t sum;
    int64_t total_link;
    FILE *tmp_fp;
    link_t link;
    link_t *link_pool;
    position_t fres;
    position_t rres;
    varr64_t bkdn1;
    varr64_t bkdn2;

    fprintf(stderr, "linking scaffolds (MIN_LINK = %ld)\n", min_link);

    sc->min_link = min_link;

    if (sc->tol > sc->seed_len)
        len_cutoff = sc->tol * 2;
    else
        len_cutoff = sc->seed_len * 2;

    tmp_fp = tmpfile_open();
    check_tmpfile(tmp_fp);

    total_link = 0;
    rewind(sc->lib->mapped_fp);
    while (fread(&fres, sizeof(position_t), 1, sc->lib->mapped_fp)) {
        fread(&rres, sizeof(position_t), 1, sc->lib->mapped_fp);

        i = abs64(fres.id) - 1;
        if (sc->con_inc[i].id == 0)
            continue;

        fres.id = fres.id > 0 ? sc->con_inc[i].id : -(sc->con_inc[i].id);
        fres.ofst = sc->con_inc[i].id > 0 ? fres.ofst : sc->con[i].len - fres.ofst - 1;
        fres.ofst += sc->node[abs64(fres.id) - 1].con[sc->con_inc[i].ofst].st;

        j = abs64(rres.id) - 1;
        if (sc->con_inc[j].id == 0)
            continue;
        rres.id = rres.id > 0 ? sc->con_inc[j].id : -(sc->con_inc[j].id);
        rres.ofst = sc->con_inc[j].id > 0 ? rres.ofst : sc->con[j].len - rres.ofst - 1;
        rres.ofst += sc->node[abs64(rres.id) - 1].con[sc->con_inc[j].ofst].st;

        if (abs64(fres.id) < abs64(rres.id)) {
            link.id1 = fres.id;
            link.ofst1 = sc->con_inc[i].ofst;
            link.id2 = -(rres.id);
            link.ofst2 = sc->con_inc[j].ofst;
        }
        else {
            link.id1 = rres.id;
            link.ofst1 = sc->con_inc[j].ofst;
            link.id2 = -(fres.id);
            link.ofst2 = sc->con_inc[i].ofst;
        }

        link.gap = sc->lib->ave_ins;
        if (fres.id > 0) {
            if (sc->node[fres.id - 1].len < len_cutoff)
                continue;
            link.gap -= sc->node[fres.id-1].len - fres.ofst;
            if (rres.id > 0) {
                if (sc->node[rres.id - 1].len < len_cutoff || fres.id == rres.id)
                    continue;
                link.gap -= sc->node[rres.id-1].len - rres.ofst;
            }
            else {
                if (sc->node[-(rres.id) - 1].len < len_cutoff || fres.id == -(rres.id))
                    continue;
                link.gap -= rres.ofst + 1;
            }
        }
        else {
            if (sc->node[-(fres.id) - 1].len < len_cutoff)
                continue;
            link.gap -= fres.ofst + 1;
            if (rres.id > 0) {
                if (sc->node[rres.id - 1].len < len_cutoff || fres.id == -(rres.id))
                    continue;
                link.gap -= sc->node[rres.id-1].len - rres.ofst;
            }
            else {
                if (sc->node[-(rres.id) - 1].len < len_cutoff || fres.id == rres.id)
                    continue;
                link.gap -= rres.ofst + 1;
            }
        }

        if (-link.gap > sc->tol + scaffolder_get_sca_overlap(sc, link.id1, link.id2))
            continue;

        fwrite(&link, sizeof(link_t), 1, tmp_fp);
        ++total_link;
    }

    link_pool = (link_t *)my_malloc((total_link + 1) * sizeof(link_t));

    rewind(tmp_fp);
    for (i = 0; i < total_link; ++i)
        fread(&(link_pool[i]), sizeof(link_t), 1, tmp_fp);
    fclose(tmp_fp);

    fputs("sorting links...\n", stderr);
    qsort(link_pool, total_link, sizeof(link_t), cmp_link);

    if (sc->link_fp != NULL)
        fclose(sc->link_fp);
    sc->link_fp = tmpfile_open();
    check_tmpfile(sc->link_fp);

    varr64_init(&bkdn1);
    varr64_init(&bkdn2);

    link_pool[total_link].id1 = 0;
    i = 0;
    while (i < total_link) {
        varr64_resize(&bkdn1, sc->node[abs64(link_pool[i].id1) - 1].n_con);
        memset(bkdn1.ele, 0, bkdn1.size * sizeof(uint64_t));
        varr64_resize(&bkdn2, sc->node[abs64(link_pool[i].id2) - 1].n_con);
        memset(bkdn2.ele, 0, bkdn2.size * sizeof(uint64_t));
        num = sum = 0;
        do {
            ++num;
            sum += link_pool[i].gap;
            ++(bkdn1.ele[link_pool[i].ofst1]);
            ++(bkdn2.ele[link_pool[i].ofst2]);
            ++i;
        } while (link_pool[i-1].id1 == link_pool[i].id1 && link_pool[i-1].id2 == link_pool[i].id2);

        if (num < min_link)
            continue;

        link = link_pool[i - 1];
        link.gap = (double)sum / num + 0.5;

        fwrite(&num, sizeof(int64_t), 1, sc->link_fp);
        fwrite(&link, sizeof(link_t), 1, sc->link_fp);
        fwrite(bkdn1.ele, sizeof(int64_t), bkdn1.size, sc->link_fp);
        fwrite(bkdn2.ele, sizeof(int64_t), bkdn2.size, sc->link_fp);

        ++(sc->node[abs64(link.id1) - 1].n_edge);
        ++(sc->node[abs64(link.id2) - 1].n_edge);
    }

    varr64_destroy(&bkdn1);
    varr64_destroy(&bkdn2);

    my_free(link_pool);
}

void scaffolder_make_graph(scaffolder_t *sc)
{
    int64_t i;
    int64_t j;
    int64_t m;
    int64_t n;
    int64_t n_link;
    link_t link;

    fprintf(stderr, "constructing scaffold graph\n");

    sc->edge_pool_size = sc->bkdn_pool_size = 0;
    for (i = 0; i < sc->n_node; ++i) {
        sc->edge_pool_size += sc->node[i].n_edge;
        sc->bkdn_pool_size += sc->node[i].n_edge * sc->node[i].n_con;
    }

    sc->mem_usage +=  sc->edge_pool_size * sizeof(sc_edge_t) + sc->bkdn_pool_size * sizeof(int64_t);
    check_mem_usage(sc->mem_usage, sc->mem_lim);
    sc->edge_pool = my_malloc(sc->edge_pool_size * sizeof(sc_edge_t));
    sc->bkdn_pool = my_calloc(sc->bkdn_pool_size, sizeof(int64_t));

    sc->edge_pool_size = sc->bkdn_pool_size = 0;
    for (i = 0; i < sc->n_node; ++i) {
        if (sc->node[i].n_edge == 0) {
            sc->node[i].edge = NULL;
            continue;
        }

        sc->node[i].edge = &(sc->edge_pool[sc->edge_pool_size]);
        sc->edge_pool_size += sc->node[i].n_edge;
        for (j = 0; j < sc->node[i].n_edge; ++j) {
            sc->node[i].edge[j].bkdn = &(sc->bkdn_pool[sc->bkdn_pool_size]);
            sc->bkdn_pool_size += sc->node[i].n_con;
        }

        sc->node[i].n_edge = 0;
    }

    rewind(sc->link_fp);
    while (fread(&n_link, sizeof(int64_t), 1, sc->link_fp)) {
        fread(&link, sizeof(link_t), 1, sc->link_fp);
        i = abs64(link.id1) - 1;
        j = abs64(link.id2) - 1;

        i = abs64(link.id1) - 1;
        j = abs64(link.id2) - 1;
        m = sc->node[i].n_edge;
        n = sc->node[j].n_edge;

        sc->node[i].edge[m].n_link = n_link;
        sc->node[j].edge[n].n_link = n_link;

        fread(sc->node[i].edge[m].bkdn, sizeof(int64_t), sc->node[i].n_con, sc->link_fp);
        fread(sc->node[j].edge[n].bkdn, sizeof(int64_t), sc->node[j].n_con, sc->link_fp);

        sc->node[i].edge[m].len = link.gap;
        sc->node[j].edge[n].len = link.gap;
        sc->node[i].edge[m].dir = (int8_t)(link.id1 / (i + 1));
        sc->node[j].edge[n].dir = (int8_t)(-(link.id2) / (j + 1));
        if (link.id1 * link.id2 > 0) {
            sc->node[i].edge[m].ed = j + 1;
            sc->node[j].edge[n].ed = i + 1;
        }
        else {
            sc->node[i].edge[m].ed = -(j + 1);
            sc->node[j].edge[n].ed = -(i + 1);
        }
        ++(sc->node[i].n_edge);
        ++(sc->node[j].n_edge);
    }
}

int64_t scaffolder_delete_erroneous_edge(scaffolder_t *sc)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t n_delete;
    double e_rate;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
    varr64_t ids;

    varr64_init(&ids);
    n_delete = 0;

    for (i = 0; i < sc->n_node; ++i) {
        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {

                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if (e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if (e2->len  + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

				e_rate = (double)(e1->n_link) / expected_link(sc->lib, sc->node[i].len, n1->len, e1->len);
				e_rate /= (double)(e2->n_link) / expected_link(sc->lib, sc->node[i].len, n2->len, e2->len);
				if (e_rate <= SC_EDGE_EXPECTED_RATE_TH) {
					varr64_push(&ids, i + 1);
					varr64_push(&ids, e1->ed);
					++n_delete;
				}
				else {
					e_rate = (double)(e2->n_link) / expected_link(sc->lib, sc->node[i].len, n2->len, e2->len);
					e_rate /= (double)(e1->n_link) / expected_link(sc->lib, sc->node[i].len, n1->len, e1->len);
					if (e_rate <= SC_EDGE_EXPECTED_RATE_TH) {
						varr64_push(&ids, i + 1);
						varr64_push(&ids, e2->ed);
						++n_delete;
					}
				}
            }
        }
    }

    scaffolder_delete_edges(sc, &ids);
    varr64_destroy(&ids);

    return n_delete;
}

void scaffolder_crush_bubble_iterative(scaffolder_t *sc, double bubble_th, FILE **bubble_fpp)
{
    uint64_t n_crush;
    uint64_t total_crush;

    fprintf(stderr, "removing bubbles... (MAX_BUBBLE_IDENTITY = %f)\n", bubble_th);
    total_crush = 0;
    do {
        n_crush = scaffolder_crush_bubble(sc, bubble_th, bubble_fpp);
        total_crush += n_crush;
        fprintf(stderr, "NUM_CRUSH=%lu\n", n_crush);
    } while (n_crush > 0);

    fprintf(stderr, "TOTAL_NUM_CRUSH=%lu\n", total_crush);
}

void scaffolder_delete_erroneous_edge_iterative(scaffolder_t *sc)
{
    uint64_t n_delete;
    uint64_t total_delete;

    fputs("removing erroneous edges...\n", stderr);
    fprintf(stderr, "EXPECTED_RATE_THRESHOLD=%f\n", SC_EDGE_EXPECTED_RATE_TH);
    total_delete = 0;
    do {
        n_delete = scaffolder_delete_erroneous_edge(sc);
        total_delete += n_delete;
        fprintf(stderr, "NUM_DELETE=%lu\n", n_delete);
    } while (n_delete > 0);

    fprintf(stderr, "TOTAL_NUM_DELETE=%lu\n", total_delete);
}

void scaffolder_delete_repeat_edge(scaffolder_t *sc)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t m;
    int64_t n;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
    varr64_t ids;

    varr64_init(&ids);

    fputs("deleting edges from repeat contigs...\n", stderr);

    for (i = 0; i < sc->n_node; ++i) {
        if (sc->node[i].n_con == 1)
            continue;
        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {
                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if (e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if (e2->len + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

                for (m = 0; m < sc->node[i].n_con; ++m) {
					if (e1->bkdn[m] < sc->min_link || e2->bkdn[m] < sc->min_link)
                        continue;

                    for (n = 0; n < sc->node[i].n_edge; ++n) {
                        sc->node[i].edge[n].n_link -= sc->node[i].edge[n].bkdn[m];
                        sc->node[i].edge[n].bkdn[m] = 0;
                    }
                    sc->con_inc[abs64(sc->node[i].con[m].id) - 1].id = 0;
                }
            }
        }
    }

    for (i = 0; i < sc->n_node; ++i) {
        for (j = 0; j < sc->node[i].n_edge; ++j) {
            if (sc->node[i].edge[j].n_link < sc->min_link) {
                varr64_push(&ids, i + 1);
                varr64_push(&ids, sc->node[i].edge[j].ed);
            }
        }
    }

    scaffolder_delete_edges(sc, &ids);

    varr64_destroy(&ids);
}

void scaffolder_detect_repeat(scaffolder_t *sc)
{
    int64_t i;
    int64_t j;
    int64_t k;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
	const double cov_th = sc->ave_cov * 2.0;
//	const double len_th = (int64_t)(scaffolder_ave_len(sc) * 2.0 + 0.5);

    for (i = 0; i < sc->n_node; ++i) {
//		if (sc->node[i].n_con == 1 && sc->node[i].len < len_th && scaffolder_node_cov(sc, &(sc->node[i])) > cov_th) {
		if (sc->node[i].n_con == 1 && scaffolder_node_cov(sc, &(sc->node[i])) > cov_th) {
			sc->node[i].state |= SC_REP;
			continue;
		}

        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {
                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if (e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if (e2->len + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

                sc->node[i].state |= SC_REP;
                j = sc->node[i].n_edge;
                break;
            }
        }
    }
}

void scaffolder_delete_edges(scaffolder_t *sc, varr64_t *ids)
{
    int64_t i;
    int64_t j;
    int64_t id1;
    int64_t id2;
    sc_node_t *n;

    for (i = 0; i < ids->size; i += 2) {
        id1 = ids->ele[i];
        id2 = ids->ele[i + 1];
        n = &(sc->node[id1 - 1]);
        for (j = 0; j < n->n_edge; ++j) {
            if (n->edge[j].ed == id2) {
                n->edge[j] = n->edge[n->n_edge - 1];
                --(n->n_edge);
                break;
            }
        }

        if (id2 > 0) {
            n = &(sc->node[id2 - 1]);
            for (j = 0; j < n->n_edge; ++j) {
                if (n->edge[j].ed == id1) {
                    n->edge[j] = n->edge[n->n_edge - 1];
                    --(n->n_edge);
                    break;
                }
            }
        }
        else {
            n = &(sc->node[-id2 - 1]);
            for (j = 0; j < n->n_edge; ++j) {
                if (n->edge[j].ed == -id1) {
                    n->edge[j] = n->edge[n->n_edge - 1];
                    --(n->n_edge);
                    break;
                }
            }
        }
    }
}

void scaffolder_split(scaffolder_t *sc)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t m;
    int64_t st;
    int64_t new_n_node;
    int64_t new_con_pool_size;
    bool is_split;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
    varr8_t breakpoint;
    FILE *sca_fp;

    new_n_node = new_con_pool_size = 0;
    varr8_init(&breakpoint);
    sca_fp = tmpfile_open();
    check_tmpfile(sca_fp);

    fputs("splitting erroneous scaffolds...\n", stderr);

    for (i = 0; i < sc->n_node; ++i) {
        if (sc->node[i].n_con == 1) {
            j = 1;
            fwrite(&j, sizeof(int64_t), 1, sca_fp);
            fwrite(sc->node[i].con, sizeof(scaffold_part_t), 1, sca_fp);
            ++new_n_node;
            ++new_con_pool_size;
            continue;
        }
        is_split = false;
        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {
                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if (e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if (e2->len + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

                varr8_resize(&breakpoint, sc->node[i].n_con + 1);
                memset(breakpoint.ele, 0, sc->node[i].n_con + 1);
                m = 0;
                while (m < sc->node[i].n_con) {
                    while (m < sc->node[i].n_con && e1->bkdn[m] < sc->min_link && e1->bkdn[m] < sc->min_link)
                        ++m;
                    breakpoint.ele[m] = 1;
                    while (m < sc->node[i].n_con && e1->bkdn[m] >= sc->min_link)
                        ++m;
                    breakpoint.ele[m] = 1;
                    while (m < sc->node[i].n_con && e2->bkdn[m] >= sc->min_link)
                        ++m;
                    breakpoint.ele[m] = 1;
                }
                is_split = true;
            }
        }

        if (is_split) {
            j = 0;
            while (j < sc->node[i].n_con) {
                st = sc->node[i].con[j].st;
                k = j;
                while (breakpoint.ele[j + 1] == 0) {
                    sc->node[i].con[j].st -= st;
                    sc->node[i].con[j].ed -= st;
                    ++j;
                }
                sc->node[i].con[j].st -= st;
                sc->node[i].con[j].ed -= st;
                ++j;
                m = j - k;

                fwrite(&m, sizeof(int64_t), 1, sca_fp);
                fwrite(&(sc->node[i].con[k]), sizeof(scaffold_part_t), m, sca_fp);
                ++new_n_node;
                new_con_pool_size += m;
            }
        }
        else {
            fwrite(&(sc->node[i].n_con), sizeof(int64_t), 1, sca_fp);
            fwrite(sc->node[i].con, sizeof(scaffold_part_t), sc->node[i].n_con ,sca_fp);
            ++new_n_node;
            new_con_pool_size += sc->node[i].n_con;
        }

    }

    varr8_destroy(&breakpoint);

    scaffolder_remake(sc, new_n_node, new_con_pool_size, sca_fp);
    fclose(sca_fp);
}

void scaffolder_scaffold(scaffolder_t *sc)
{
    FILE *sca_fp;
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t tmp;
    int64_t n_con;
    int64_t new_n_node;
    int64_t new_con_pool_size;
    int64_t min_dist;
    int64_t min_n_link;
    int64_t min_i;
    int64_t min_st;
    int64_t ol_len;
    sc_node_t *new_node;
    scaffold_part_t tmp_con;
    sc_layout_t inc;
    sc_layout_t cand;

    fputs("scaffolding...\n", stderr);

    sc_layout_init(&inc);
    sc_layout_init(&cand);

    new_n_node = new_con_pool_size = 0;
    sca_fp = tmpfile_open();
    check_tmpfile(sca_fp);

    for (i = 0; i < sc->n_node; ++i) {
        if (sc->node[i].state & (SC_INC | SC_REP | SC_DEL))
            continue;

        new_node = &(sc->node[i]);
        sc_layout_resize(&inc, 1);
        sc_layout_resize(&cand, 0);
        inc.part[0].st =  inc.part[0].dist = 0;
        inc.part[0].ed = new_node->len;
        inc.part[0].id = i + 1;
        n_con = new_node->n_con;
        new_node->state |= SC_INC;

        for (k = 0; k < new_node->n_edge; ++k) {
            tmp = abs64(new_node->edge[k].ed) - 1;
            if ((sc->node[tmp].state & SC_INC) && !(sc->node[tmp].state & SC_REP))
                continue;

            sc_layout_resize(&cand, cand.size + 1);
            if (new_node->edge[k].dir > 0) {
                cand.part[cand.size - 1].st = inc.part[0].ed + new_node->edge[k].len;
                cand.part[cand.size - 1].ed = cand.part[cand.size - 1].st + sc->node[tmp].len;
            }
            else {
                cand.part[cand.size - 1].ed = -(new_node->edge[k].len);
                cand.part[cand.size - 1].st = cand.part[cand.size - 1].ed - sc->node[tmp].len;
            }
            cand.part[cand.size - 1].id = new_node->edge[k].ed;
            cand.part[cand.size - 1].dist = 1;
            cand.part[cand.size - 1].n_link = new_node->edge[k].n_link;
        }

        while (cand.size > 0) {
            min_dist = cand.part[0].dist;

            min_n_link = cand.part[0].st;
            min_i = 0;
            for (j = 1; j < cand.size; ++j) {
                if (cand.part[j].dist < min_dist) {
                    min_dist = cand.part[j].dist;
                    min_n_link = cand.part[j].st;
                    min_i = j;
                }
                else if (cand.part[j].dist == min_dist && cand.part[j].st < min_n_link) {
                    min_n_link = cand.part[j].st;
                    min_i = j;
                }
            }

            tmp = abs64(cand.part[min_i].id) - 1;
            if ((sc->node[tmp].state & SC_INC) && !((sc->node[tmp].state & SC_INC) & SC_REP)) {
                sc_layout_delete(&cand, min_i);
                continue;
            }

            for (j = 0; j < inc.size; ++j) {
                if (cand.part[min_i].ed <= inc.part[j].st)
                    continue;
                if (cand.part[min_i].st >= inc.part[j].ed)
                    continue;
                if (abs64(inc.part[j].ed - cand.part[min_i].st) <= sc->tol + scaffolder_get_sca_overlap(sc, inc.part[j].id, cand.part[min_i].id))
                    continue;
                if (abs64(cand.part[min_i].ed - inc.part[j].st) <= sc->tol + scaffolder_get_sca_overlap(sc, cand.part[min_i].id, inc.part[j].id))
                    continue;
                break;
            }

            if (j == inc.size) {
                sc_layout_resize(&inc, inc.size + 1);
                inc.part[inc.size - 1] = cand.part[min_i];

                new_node = &(sc->node[abs64(inc.part[inc.size - 1].id) - 1]);
                if (~(new_node->state) & SC_REP) {
                    for (k = 0; k < new_node->n_edge; ++k) {
                        tmp = abs64(new_node->edge[k].ed) - 1;
                        if ((sc->node[tmp].state & SC_INC) && !(sc->node[tmp].state & SC_REP))
                            continue;

                        sc_layout_resize(&cand, cand.size + 1);
                        if (inc.part[inc.size - 1].id * new_node->edge[k].dir > 0) {
                            cand.part[cand.size - 1].st = inc.part[inc.size - 1].ed + new_node->edge[k].len;
                            cand.part[cand.size - 1].ed = cand.part[cand.size - 1].st + sc->node[tmp].len;
                        }
                        else {
                            cand.part[cand.size - 1].ed = inc.part[inc.size - 1].st - new_node->edge[k].len;
                            cand.part[cand.size - 1].st = cand.part[cand.size - 1].ed - sc->node[tmp].len;
                        }
                        cand.part[cand.size - 1].id = inc.part[inc.size - 1].id > 0 ? new_node->edge[k].ed : -(new_node->edge[k].ed);
                        cand.part[cand.size - 1].dist = inc.part[inc.size - 1].dist + 1;
                        cand.part[cand.size - 1].n_link = new_node->edge[k].n_link;
                    }
                }

                n_con += new_node->n_con;
                if (!(new_node->state & SC_REP))
                    new_node->state |= SC_INC;
            }

            sc_layout_delete(&cand, min_i);
        }

        qsort(inc.part, inc.size, sizeof(sc_layout_part_t), cmp_sc_layout_part);

        for (j = 0; sc->node[abs64(inc.part[j].id) - 1].state & SC_REP; ++j)
            n_con -= sc->node[abs64(inc.part[j].id) - 1].n_con;
        for (; sc->node[abs64(inc.part[inc.size - 1].id) - 1].state & SC_REP; --inc.size)
            n_con -= sc->node[abs64(inc.part[inc.size - 1].id) - 1].n_con;

        fwrite(&n_con, sizeof(int64_t), 1, sca_fp);
        min_st = inc.part[j].st;
        for (; j < inc.size; ++j) {
            tmp = abs64(inc.part[j].id) - 1;

            sc->node[tmp].state |= SC_INC;

            inc.part[j].st -= min_st;
            inc.part[j].ed -= min_st;

            if (inc.part[j].st != 0) {
                ol_len = scaffolder_get_sca_overlap(sc, inc.part[j-1].id, inc.part[j].id);
                if (ol_len > 0 && ol_len + inc.part[j].st - inc.part[j-1].ed <= sc->tol) {
                    ol_len = inc.part[j-1].ed - inc.part[j].st - ol_len;
                    for (k = j; k < inc.size; ++k) {
                        inc.part[k].ed += ol_len;
                        inc.part[k].st += ol_len;
                    }
                }
                else if (inc.part[j].st < inc.part[j-1].ed) {
                    ol_len = inc.part[j-1].ed - inc.part[j].st + 1;
                    for (k = j; k < inc.size; ++k) {
                        inc.part[k].ed += ol_len;
                        inc.part[k].st += ol_len;
                    }
                }
            }

            if (inc.part[j].id > 0) {
                for (k = 0; k < sc->node[tmp].n_con; ++k) {
                    tmp_con.st = inc.part[j].st + sc->node[tmp].con[k].st;
                    tmp_con.ed = inc.part[j].st + sc->node[tmp].con[k].ed;
                    tmp_con.id = sc->node[tmp].con[k].id;
                    fwrite(&tmp_con, sizeof(scaffold_part_t), 1, sca_fp);
                }
            }
            else {
                for (k = sc->node[tmp].n_con - 1; k >= 0; --k) {
                    tmp_con.st = inc.part[j].st + sc->node[tmp].len - sc->node[tmp].con[k].ed;
                    tmp_con.ed = inc.part[j].st + sc->node[tmp].len - sc->node[tmp].con[k].st;
                    tmp_con.id = -(sc->node[tmp].con[k].id);
                    fwrite(&tmp_con, sizeof(scaffold_part_t), 1, sca_fp);
                }
            }
        }

        new_con_pool_size += n_con;
        ++new_n_node;
    }

    sc_layout_destroy(&inc);
    sc_layout_destroy(&cand);

    for (i = 0; i < sc->n_node; ++i) {
        if (!(sc->node[i].state & SC_REP) || (sc->node[i].state & (SC_INC | SC_DEL)))
            continue;
		fwrite(&(sc->node[i].n_con), sizeof(int64_t), 1, sca_fp);
		for (j = 0; j < sc->node[i].n_con; ++j)
			fwrite(&(sc->node[i].con[j]), sizeof(scaffold_part_t), 1, sca_fp);
		new_con_pool_size += sc->node[i].n_con;
		++new_n_node;
    }

    scaffolder_remake(sc, new_n_node, new_con_pool_size, sca_fp);
    fclose(sca_fp);
}

void scaffolder_remake(scaffolder_t *sc, int64_t new_n_node, int64_t new_con_pool_size, FILE *sca_fp)
{
    int64_t i;
    int64_t j;
    int64_t tmp;
	int64_t *n_occ;

	n_occ = (int64_t *)my_calloc(sc->n_con, sizeof(int64_t));

    scaffolder_destroy_graph(sc);

    sc->con_pool_size = new_con_pool_size;
    sc->n_node = new_n_node;
    sc->mem_usage += sc->n_node * sizeof(sc_node_t) + sc->con_pool_size * sizeof(scaffold_part_t);
    check_mem_usage(sc->mem_usage, sc->mem_lim);
    sc->node = (sc_node_t *)my_realloc(sc->node, sc->n_node * sizeof(sc_node_t));
    memset(sc->node, 0, sc->n_node * sizeof(sc_node_t));
    sc->con_pool = (scaffold_part_t *)my_malloc(sc->con_pool_size * sizeof(scaffold_part_t));

    sc->con_pool_size = 0;
    rewind(sca_fp);
    for (i = 0; i < sc->n_node; ++i) {
        fread(&(sc->node[i].n_con), sizeof(int64_t), 1, sca_fp);
        sc->node[i].con = &(sc->con_pool[sc->con_pool_size]);
        fread(sc->node[i].con, sizeof(scaffold_part_t), sc->node[i].n_con, sca_fp);
        sc->node[i].len = sc->node[i].con[sc->node[i].n_con - 1].ed;
        sc->con_pool_size += sc->node[i].n_con;

        for (j = 0; j < sc->node[i].n_con; ++j) {
            tmp = abs64(sc->node[i].con[j].id) - 1;
			++n_occ[tmp];
			sc->con_inc[tmp].ofst = j;
			if (sc->node[i].con[j].id > 0)
				sc->con_inc[tmp].id = i + 1;
			else
				sc->con_inc[tmp].id = -(i + 1);
        }
    }

	for (i = 0; i < sc->n_con; ++i) {
		if (n_occ[i] != 1)
			sc->con_inc[i].id = 0;
	}

	my_free(n_occ);
}

int cmp_link(const void *a, const void *b)
{
    if (((link_t *)a)->id1 != ((link_t *)b)->id1)
        return ((link_t *)a)->id1 - ((link_t *)b)->id1;
    else if (((link_t *)a)->id2 != ((link_t *)b)->id2)
        return ((link_t *)a)->id2 - ((link_t *)b)->id2;
    else
        return ((link_t *)a)->gap - ((link_t *)b)->gap;
}

int cmp_sc_layout_part(const void *a, const void *b)
{
    return ((sc_layout_part_t *)a)->st - ((sc_layout_part_t *)b)->ed - (((sc_layout_part_t *)b)->st - ((sc_layout_part_t *)a)->ed);
}

double expected_link(seqlib_t *lib, double l1, double l2, double g)
{
    double r, c, ave_ins, sd_ins, n_link;

    ave_ins = (double)(lib->ave_ins);
    sd_ins = (double)(lib->sd_ins);
    r = (double)lib->ave_len;
    c = lib->ave_cov;

    n_link = 0.0;

    n_link += (l1 + g - ave_ins + l2) * erf((l1 + g - ave_ins + l2) / (M_SQRT2 * sd_ins));
    n_link += (M_SQRT2 * sd_ins / sqrt(M_PI)) * exp(-pow(((l1 + g - ave_ins + l2) / (M_SQRT2 * sd_ins)), 2.0));
    n_link -= (r + g - ave_ins + l2) * erf((r + g - ave_ins + l2) / (M_SQRT2 * sd_ins));
    n_link -= (M_SQRT2 * sd_ins / sqrt(M_PI)) * exp(-pow(((r + g - ave_ins + l2) / (M_SQRT2 * sd_ins)), 2.0));

    n_link -= (l1 + g - ave_ins + r) * erf((l1 + g - ave_ins + r) / (M_SQRT2 * sd_ins));
    n_link -= (M_SQRT2 * sd_ins / sqrt(M_PI)) * exp(-pow(((l1 + g - ave_ins + r) / (M_SQRT2 * sd_ins)), 2.0));
    n_link += (r + g - ave_ins + r) * erf((r + g - ave_ins + r) / (M_SQRT2 * sd_ins));
    n_link += (M_SQRT2 * sd_ins / sqrt(M_PI)) * exp(-pow(((r + g - ave_ins + r) / (M_SQRT2 * sd_ins)), 2.0));

    n_link *= c / (4.0 * r);

    if (n_link < 1.0)
        n_link = 1.0;

    return n_link;
}

void scaffolder_print_layout(scaffolder_t *sc, char *out_name)
{
    int64_t i;
    int64_t j;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    for (i = 0; i < sc->n_node; ++i) {
        fprintf(out, "scaffold%ld_len%ld\n{\n", i + 1, sc->node[i].con[sc->node[i].n_con - 1].ed);
        for (j = 0; j < sc->node[i].n_con; ++j) {
            if (sc->node[i].con[j].id > 0)
                fprintf(out, "contig%ld\t+\t%ld\t%ld\n", sc->node[i].con[j].id, sc->node[i].con[j].st, sc->node[i].con[j].ed);
            else    
                fprintf(out, "contig%ld\t-\t%ld\t%ld\n", -(sc->node[i].con[j].id), sc->node[i].con[j].st, sc->node[i].con[j].ed);
        }
        fputs("}\n\n", out);
    }

    fclose(out);
}

void scaffolder_print_seq(scaffolder_t *sc, char *out_name)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t l;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    for (i = 0; i < sc->n_node; ++i) {

        fprintf(out, ">scaffold%ld_len%ld\n", i + 1, sc->node[i].con[sc->node[i].n_con - 1].ed);
        l = 0;
        for (j = 0; j < sc->node[i].n_con; ++j) {

            if (j > 0 && sc->node[i].con[j-1].ed > sc->node[i].con[j].st)
                k = sc->node[i].con[j-1].ed - sc->node[i].con[j].st;
            else
                k = 0;

            if (sc->node[i].con[j].id > 0) {
                for (; k < sc->con[sc->node[i].con[j].id - 1].len; ++k) {
                    putc(num2ascii(sc->con[sc->node[i].con[j].id - 1].seq[k]), out);
                    ++l;
                    if (l % LINE_LENGTH == 0)
                        putc('\n', out);
                }
            }
            else {
                for (k = sc->con[-(sc->node[i].con[j].id) - 1].len - 1 - k; k >= 0; --k) {
                    putc(num2ascii(3 ^ sc->con[-(sc->node[i].con[j].id) - 1].seq[k]), out);
                    ++l;
                    if (l % LINE_LENGTH == 0)
                        putc('\n', out);
                }
            }
            if (j == sc->node[i].n_con - 1)
                continue;

            for (k = sc->node[i].con[j].ed; k < sc->node[i].con[j+1].st; ++k) {
                putc('N', out);
                ++l;
                if (l % LINE_LENGTH == 0)
                    putc('\n', out);
            }
        }
        if (l % LINE_LENGTH != 0)
            putc('\n', out);
    }

    fclose(out);
}

void scaffolder_cut_and_print_seq(scaffolder_t *sc, int64_t min_len, char *out_name)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t l;
    int64_t len;
    int64_t max_n_con;
    int64_t min_overlap;
    int64_t n_sca;
    FILE *out;
    int64_t *max_l_ol;
    int64_t *max_r_ol;
    int64_t *l_cut;
    int64_t *r_cut;
    int64_t *gap;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    max_l_ol = (int64_t *)my_calloc(sc->n_con, sizeof(int64_t));
    max_r_ol = (int64_t *)my_calloc(sc->n_con, sizeof(int64_t));

    min_overlap = sc->min_overlap;
    sc->min_overlap = 1;

    for (i = 0; i < sc->ol_table_size; ++i) {
        if (overlap_rec_is_emp(sc->ol_table[i]))
            continue;
        if (sc->ol_table[i].id1 > 0) {
            if (sc->ol_table[i].len > max_r_ol[sc->ol_table[i].id1 - 1])
                max_r_ol[sc->ol_table[i].id1 - 1] = sc->ol_table[i].len;
        }
        else if (sc->ol_table[i].len > max_l_ol[-(sc->ol_table[i].id1) - 1])
            max_l_ol[-(sc->ol_table[i].id1) - 1] = sc->ol_table[i].len;

        if (sc->ol_table[i].id2 > 0) {
            if (sc->ol_table[i].len > max_l_ol[sc->ol_table[i].id2 - 1])
                max_l_ol[sc->ol_table[i].id2 - 1] = sc->ol_table[i].len;
        }
        else if (sc->ol_table[i].len > max_r_ol[-(sc->ol_table[i].id2) - 1])
            max_r_ol[-(sc->ol_table[i].id2) - 1] = sc->ol_table[i].len;
    }

    max_n_con = 0;
    for (i = 0; i < sc->n_node; ++i) {
        if (sc->node[i].n_con > max_n_con)
            max_n_con = sc->node[i].n_con;
    }
    l_cut = (int64_t *)my_malloc(max_n_con * sizeof(int64_t));
    r_cut = (int64_t *)my_malloc(max_n_con * sizeof(int64_t));
    gap = (int64_t *)my_malloc(max_n_con * sizeof(int64_t));

    n_sca = 0;
    for (i = 0; i < sc->n_node; ++i) {
        for (j = 0; j < sc->node[i].n_con; ++j) {
            if (sc->con_inc[abs64(sc->node[i].con[j].id) - 1].id != 0)
                break;
        }
        if (j == sc->node[i].n_con)
            continue;

        len = 0;

        if (sc->node[i].con[0].id > 0)
            l_cut[0] = max_l_ol[sc->node[i].con[0].id - 1] / 2;
        else 
            l_cut[0] = max_r_ol[-(sc->node[i].con[0].id) - 1] / 2;

        for (j = 0; j < sc->node[i].n_con - 1; ++j) {
            if (sc->node[i].con[j].ed > sc->node[i].con[j + 1].st) {
                r_cut[j] = 0;
                l_cut[j + 1] = sc->node[i].con[j].ed - sc->node[i].con[j + 1].st;
                gap[j] = 0;
            }
            else {
                r_cut[j] = l_cut[j + 1] = scaffolder_get_overlap(sc, sc->node[i].con[j].id, sc->node[i].con[j + 1].id) / 2;
                gap[j] = sc->node[i].con[j + 1].st - sc->node[i].con[j].ed + 2 * r_cut[j];
            }
            len += (sc->node[i].con[j].ed - sc->node[i].con[j].st) - l_cut[j] - r_cut[j] + gap[j];
        }

        if (sc->node[i].con[j].id > 0)
            r_cut[j] = max_r_ol[sc->node[i].con[j].id - 1] / 2;
        else 
            r_cut[j] = max_l_ol[-(sc->node[i].con[j].id) - 1] / 2;
        len += (sc->node[i].con[j].ed - sc->node[i].con[j].st) - l_cut[j] - r_cut[j];
        gap[j] = 0;

        if (len < min_len)
            continue;

        ++n_sca;

		if (sc->node[i].n_con > 1)
			fprintf(out, ">scaffold%ld_len%ld_cov%ld\n", n_sca, len, (int64_t)(scaffolder_node_cov(sc, &(sc->node[i])) + 0.5));
		else
			fprintf(out, ">scaffold%ld_len%ld_cov%ld_single\n", n_sca, len, (int64_t)(scaffolder_node_cov(sc, &(sc->node[i])) + 0.5));

        l = 0;
        for (j = 0; j < sc->node[i].n_con; ++j) {
            if (sc->node[i].con[j].id > 0) {
                for (k = l_cut[j]; k < sc->con[sc->node[i].con[j].id - 1].len - r_cut[j]; ++k) {
                    putc(num2ascii(sc->con[sc->node[i].con[j].id - 1].seq[k]), out);
                    ++l;
                    if (l % LINE_LENGTH == 0)
                        putc('\n', out);
                }
            }
            else {
                for (k = sc->con[-(sc->node[i].con[j].id) - 1].len - l_cut[j] - 1; k >= (int64_t)r_cut[j]; --k) {
                    if (sc->con[-(sc->node[i].con[j].id) - 1].seq[k] != 4)
                        putc(num2ascii(3 ^ sc->con[-(sc->node[i].con[j].id) - 1].seq[k]), out);
                    else
                        putc('N', out);
                    ++l;
                    if (l % LINE_LENGTH == 0)
                        putc('\n', out);
                }
            }

            for (k = 0; k < gap[j]; ++k) {
                putc('N', out);
                ++l;
                if (l % LINE_LENGTH == 0)
                    putc('\n', out);
            }
        }
        if (l % LINE_LENGTH != 0)
            putc('\n', out);
    }

    fclose(out);
    my_free(l_cut);
    my_free(r_cut);
    my_free(gap);
    my_free(max_l_ol);
    my_free(max_r_ol);

    sc->min_overlap = min_overlap;
}

void sc_layout_init(sc_layout_t *lo)
{
    memset(lo, 0, sizeof(sc_layout_t));
}

void sc_layout_destroy(sc_layout_t *lo)
{
    if (lo->part != NULL)
        my_free(lo->part);
    memset(lo, 0, sizeof(sc_layout_t));
}

void sc_layout_resize(sc_layout_t *lo, int64_t size)
{
    if (size > lo->capa) {
        lo->capa = size;
        lo->part = (sc_layout_part_t *)my_realloc(lo->part, lo->capa * sizeof(sc_layout_part_t)); 
    }
    lo->size = size;
}

void sc_layout_delete(sc_layout_t *lo, int64_t i)
{
    --(lo->size);
    lo->part[i] = lo->part[lo->size];
}

int64_t scaffolder_crush_bubble(scaffolder_t *sc, double bubble_th, FILE **bubble_fpp)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t m;
    int64_t n;
    int64_t redge;
    int64_t ledge;
    int64_t n_crush;
    double cov1;
    double cov2;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
    sc_layout_t layout1;
    sc_layout_t layout2;
    sc_layout_t *lp;
    sc_layout_t work;
    varr64_t ids;
    varr8_t sca_seq1;
    varr8_t sca_seq2;
    varr64_t al_mat;
    const double cov_th = sc->ave_cov * 2.0;

    varr64_init(&ids);
    sc_layout_init(&layout1);
    sc_layout_init(&layout2);
    sc_layout_init(&work);
    varr8_init(&sca_seq1);
    varr8_init(&sca_seq2);
    varr64_init(&al_mat);
    n_crush = 0;

    scaffolder_detect_repeat(sc);

    for (i = 0; i < sc->n_node; ++i) {
        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {
                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if ((n1->state & SC_DEL) || e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if ((n2->state & SC_DEL) || e2->len + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

                scaffolder_layout_nodes(sc, n1, &layout1, &work);
                scaffolder_layout_nodes(sc, n2, &layout2, &work);

                redge = ledge =  min_64(layout1.size, layout2.size);
                for (m = 0; m < ledge; ++m) {
                    if (layout1.part[m].id != layout2.part[m].id)
                        break;
                }
                if (m == 0 || m == ledge)
                    continue;
                ledge = m - 1;

                for (m = 1; m <= redge; ++m) {
                    if (layout1.part[layout1.size - m].id != layout2.part[layout2.size - m].id)
                        break;
                }
                if (m == 1)
                    continue;
                redge = m - 1;

                cov1 = scaffolder_layout_ave_cov(sc, &(layout1.part[ledge + 1]), layout1.size - redge - ledge - 1);
                cov2 = scaffolder_layout_ave_cov(sc, &(layout2.part[ledge + 1]), layout2.size - redge - ledge - 1);

                lp = cov1 < cov2 ? &layout1 : &layout2;
                if (redge + ledge + 1 >= lp->size || cov1 + cov2 > cov_th)
                    continue;


                scaffolder_layout2seq(sc, &(layout1.part[ledge + 1]), layout1.size - redge - ledge - 1, &sca_seq1);
                scaffolder_layout2seq(sc, &(layout2.part[ledge + 1]), layout2.size - redge - ledge - 1, &sca_seq2);
                if (abs64((int64_t)sca_seq1.size - (int64_t)sca_seq2.size) > max_64(sca_seq1.size, sca_seq2.size) * bubble_th)
                    continue;
                if (align_scaffold(&sca_seq1, &sca_seq2, &al_mat) > max_64(sca_seq1.size, sca_seq2.size) * bubble_th)
                    continue;

                for (m = ledge + 1; m < lp->size - redge; ++m) {
                    n1 = &(sc->node[abs64(lp->part[m].id) - 1]);
                    for (n = 0; n < n1->n_edge; ++n) {
                        varr64_push(&ids, (int64_t)(n1 - sc->node) + 1);
                        varr64_push(&ids, n1->edge[n].ed);
                    }
                    for (n = 0; n < n1->n_con; ++n)
                        sc->con_inc[abs64(n1->con[n].id) - 1].id = 0;
                    n1->state |= SC_DEL;
                }

                scaffolder_layout2seq(sc, &(lp->part[ledge + 1]), lp->size - redge - ledge, &sca_seq1);
                fwrite(&(sca_seq1.size), sizeof(uint64_t), 1, *bubble_fpp);
                fwrite(sca_seq1.ele, sizeof(uint8_t), sca_seq1.size, *bubble_fpp);
                fwrite(cov1 < cov2 ? &cov1 : &cov2, sizeof(double), 1, *bubble_fpp);

                ++n_crush;
            }
        }
    }

    sc_layout_destroy(&layout1);
    sc_layout_destroy(&layout2);
    sc_layout_destroy(&work);
    varr8_destroy(&sca_seq1);
    varr8_destroy(&sca_seq2);
    varr64_destroy(&al_mat);

    scaffolder_delete_edges(sc, &ids);
    varr64_destroy(&ids);

    for (i = 0; i < sc->n_node; ++i)
        sc->node[i].state &= ~SC_REP;

    return n_crush;
}

void scaffolder_layout_nodes(scaffolder_t *sc, sc_node_t *new_node, sc_layout_t *ret, sc_layout_t *work)
{
    int64_t j;
    int64_t k;
    int64_t tmp;
    int64_t min_dist;
    int64_t min_n_link;
    int64_t min_i;

    sc_layout_resize(ret, 1);
    sc_layout_resize(work, 0);
    ret->part[0].st =  ret->part[0].dist = 0;
    ret->part[0].ed = new_node->len;
    ret->part[0].id = (int64_t)(new_node - sc->node) + 1;
    new_node->state |= SC_INC;

    for (k = 0; k < new_node->n_edge; ++k) {
        tmp = abs64(new_node->edge[k].ed) - 1;
        if ((sc->node[tmp].state & SC_INC) && !(sc->node[tmp].state & SC_REP))
            continue;

        sc_layout_resize(work, work->size + 1);
        if (new_node->edge[k].dir > 0) {
            work->part[work->size - 1].st = ret->part[0].ed + new_node->edge[k].len;
            work->part[work->size - 1].ed = work->part[work->size - 1].st + sc->node[tmp].len;
        }
        else {
            work->part[work->size - 1].ed = -(new_node->edge[k].len);
            work->part[work->size - 1].st = work->part[work->size - 1].ed - sc->node[tmp].len;
        }
        work->part[work->size - 1].id = new_node->edge[k].ed;
        work->part[work->size - 1].dist = 1;
        work->part[work->size - 1].n_link = new_node->edge[k].n_link;
    }

    while (work->size > 0) {
        min_dist = work->part[0].dist;

        min_n_link = work->part[0].st;
        min_i = 0;
        for (j = 1; j < work->size; ++j) {
            if (work->part[j].dist < min_dist) {
                min_dist = work->part[j].dist;
                min_n_link = work->part[j].st;
                min_i = j;
            }
            else if (work->part[j].dist == min_dist && work->part[j].st < min_n_link) {
                min_n_link = work->part[j].st;
                min_i = j;
            }
        }

        tmp = abs64(work->part[min_i].id) - 1;
        if ((sc->node[tmp].state & SC_INC) && !((sc->node[tmp].state & SC_INC) & SC_REP)) {
            sc_layout_delete(work, min_i);
            continue;
        }

        for (j = 0; j < ret->size; ++j) {
            if (work->part[min_i].ed <= ret->part[j].st)
                continue;
            if (work->part[min_i].st >= ret->part[j].ed)
                continue;
            if (abs64(ret->part[j].ed - work->part[min_i].st) <= sc->tol + scaffolder_get_sca_overlap(sc, ret->part[j].id, work->part[min_i].id))
                continue;
            if (abs64(work->part[min_i].ed - ret->part[j].st) <= sc->tol + scaffolder_get_sca_overlap(sc, work->part[min_i].id, ret->part[j].id))
                continue;
            break;
        }

        if (j == ret->size) {
            sc_layout_resize(ret, ret->size + 1);
            ret->part[ret->size - 1] = work->part[min_i];

            new_node = &(sc->node[abs64(ret->part[ret->size - 1].id) - 1]);
            if (~(new_node->state) & SC_REP) {
                for (k = 0; k < new_node->n_edge; ++k) {
                    tmp = abs64(new_node->edge[k].ed) - 1;
                    if ((sc->node[tmp].state & SC_INC) && !(sc->node[tmp].state & SC_REP))
                        continue;

                    sc_layout_resize(work, work->size + 1);
                    if (ret->part[ret->size - 1].id * new_node->edge[k].dir > 0) {
                        work->part[work->size - 1].st = ret->part[ret->size - 1].ed + new_node->edge[k].len;
                        work->part[work->size - 1].ed = work->part[work->size - 1].st + sc->node[tmp].len;
                    }
                    else {
                        work->part[work->size - 1].ed = ret->part[ret->size - 1].st - new_node->edge[k].len;
                        work->part[work->size - 1].st = work->part[work->size - 1].ed - sc->node[tmp].len;
                    }
                    work->part[work->size - 1].id = ret->part[ret->size - 1].id > 0 ? new_node->edge[k].ed : -(new_node->edge[k].ed);
                    work->part[work->size - 1].dist = ret->part[ret->size - 1].dist + 1;
                    work->part[work->size - 1].n_link = new_node->edge[k].n_link;
                }
            }

            if (!(new_node->state & SC_REP))
                new_node->state |= SC_INC;
        }

        sc_layout_delete(work, min_i);
    }

    qsort(ret->part, ret->size, sizeof(sc_layout_part_t), cmp_sc_layout_part);

    for (j = 1; j < ret->size; ++j) {
        sc->node[abs64(ret->part[j].id) - 1].state &= ~SC_INC;
        ret->part[j].st -= ret->part[0].st;
        ret->part[j].ed -= ret->part[0].st;
        if (ret->part[j].st != 0) {
            tmp = scaffolder_get_sca_overlap(sc, ret->part[j-1].id, ret->part[j].id);
            if (tmp + ret->part[j].st - ret->part[j-1].ed <= sc->tol) {
                tmp = ret->part[j-1].ed - ret->part[j].st - tmp;
                for (k = j; k < ret->size; ++k) {
                    ret->part[k].ed += tmp;
                    ret->part[k].st += tmp;
                }
            }
            else if (ret->part[j].st < ret->part[j-1].ed) {
                tmp = ret->part[j-1].ed - ret->part[j].st + 1;
                for (k = j; k < ret->size; ++k) {
                    ret->part[k].ed += tmp;
                    ret->part[k].st += tmp;
                }
            }
        }
    }
    sc->node[abs64(ret->part[0].id) - 1].state &= ~SC_INC;
    ret->part[0].ed -= ret->part[0].st;
    ret->part[0].st = 0;
}

void scaffolder_layout2seq(scaffolder_t *sc, sc_layout_part_t *lo, int64_t lo_size, varr8_t *ret)
{
    int64_t i;
    int64_t j;
    int64_t k;
    sc_node_t *node;
    lseq_t *con;

    varr8_clear(ret);

    for (i = 0; i < lo_size; ++i) {
        if (lo[i].id > 0) {
            node = &(sc->node[lo[i].id - 1]);
            for (j = 0; j < node->n_con; ++j) {

                if (i == 0)
                    k = 0;
                else if (j == 0)
                    k = lo[i - 1].ed - lo[i].st;
                else
                    k = node->con[j - 1].ed - node->con[j].st;

                for (; k < 0; ++k)
                    varr8_push(ret, 4);

                if (node->con[j].id > 0) {
                    con = &(sc->con[node->con[j].id - 1]);
                    for (; k < con->len; ++k)
                        varr8_push(ret, con->seq[k]);
                }
                else {
                    con = &(sc->con[-(node->con[j].id) - 1]);
                    for (; k < con->len; ++k) {
                        if (con->seq[con->len - k - 1] >= 4)
                            varr8_push(ret, 4);
                        else 
                            varr8_push(ret, 3^(con->seq[con->len - k - 1]));
                    }
                }
            }
        }
        else {
            node = &(sc->node[-(lo[i].id) - 1]);
            for (j = node->n_con - 1; j >= 0; --j) {

                if (i == 0)
                    k = 0;
                else if (j == node->n_con - 1)
                    k = lo[i - 1].ed - lo[i].st;
                else
                    k = node->con[j].ed - node->con[j + 1].st;

                for (; k < 0; ++k) {
                    varr8_push(ret, 4);
                }

                if (node->con[j].id > 0) {
                    con = &(sc->con[node->con[j].id - 1]);
                    for (; k < con->len; ++k) {
                        if (con->seq[con->len - k - 1] >= 4)
                            varr8_push(ret, 4);
                        else 
                            varr8_push(ret, 3^(con->seq[con->len - k - 1]));
                    }
                }
                else {
                    con = &(sc->con[-(node->con[j].id) - 1]);
                    for (; k < con->len; ++k)
                        varr8_push(ret, con->seq[k]);
                }
            }
        }
    }
}

int64_t align_scaffold(varr8_t *sca1, varr8_t *sca2, varr64_t *work)
{
    int64_t m;
    int64_t n;
    int64_t wsize;

    wsize = sca2->size + 1;
    varr64_resize(work, 2 * wsize);
    for (n = 0; n < sca2->size + 1; ++n)
        work->ele[n] = n;
    for (m = 0; m < sca1->size; ++m) {
        work->ele[(m&1)*wsize] = m;
        work->ele[(~m&1)*wsize] = m+1;
        for (n = 0; n < sca2->size; ++n) {
            if (sca1->ele[m] == sca2->ele[n]) {
                work->ele[(~m&1)*wsize + n+1] = work->ele[(m&1)*wsize + n];
                continue;
            }
            work->ele[(~m&1)*wsize + n+1] = work->ele[(m&1)*wsize + n] + 1;
            if (work->ele[(~m&1)*wsize + n+1] > work->ele[(m&1)*wsize + n+1] + 1)
                work->ele[(~m&1)*wsize + n+1] = work->ele[(m&1)*wsize + n+1] + 1;
            if (work->ele[(~m&1)*wsize + n+1] > work->ele[(~m&1)*wsize + n] + 1)
                work->ele[(~m&1)*wsize + n+1] = work->ele[(~m&1)*wsize + n] + 1;
        }
    }
    return (int64_t)(work->ele[(m&1)*wsize + n]);
}

double scaffolder_layout_ave_cov(scaffolder_t *sc, sc_layout_part_t *lo, int64_t lo_size)
{
    int64_t i;
    int64_t j;
    int64_t num;
    int64_t sum;
    sc_node_t *node;

    sum = num = 0;
    for (i = 0; i < lo_size; ++i) {
        node = &(sc->node[abs64(lo[i].id) - 1]);
        for (j = 0; j < node->n_con; ++j) {
            if (node->con[j].id > 0) {
                num += sc->con[node->con[j].id - 1].len;
                sum += sc->cov[node->con[j].id - 1] * sc->con[node->con[j].id - 1].len;
            }
            else {
                num += sc->con[-(node->con[j].id) - 1].len;
                sum += sc->cov[-(node->con[j].id) - 1] * sc->con[-(node->con[j].id) - 1].len;
            }
        }
    }

    if (num != 0)
        return (double)sum / num;
    else
        return 0.0;
}

double scaffolder_node_cov(scaffolder_t *sc, sc_node_t *node)
{
    int64_t j;
    int64_t num;
    int64_t sum;

    sum = num = 0;
	for (j = 0; j < node->n_con; ++j) {
		if (node->con[j].id > 0) {
			num += sc->con[node->con[j].id - 1].len;
			sum += sc->cov[node->con[j].id - 1] * sc->con[node->con[j].id - 1].len;
		}
		else {
			num += sc->con[-(node->con[j].id) - 1].len;
			sum += sc->cov[-(node->con[j].id) - 1] * sc->con[-(node->con[j].id) - 1].len;
		}
	}

    if (num != 0)
        return (double)sum / num;
    else
        return 0.0;
}

void print_scaffold_bubble(FILE *bubble_fp, char *out_name)
{
    int64_t i;
    int64_t n_seq;
    uint64_t len;
    double cov;
    varr8_t seq;
    FILE *out;

    varr8_init(&seq);

    out = fopen(out_name, "w");
    check_file_open(out, out_name);

    n_seq = 0;
    rewind(bubble_fp);
    while (fread(&len, sizeof(uint64_t), 1, bubble_fp)) {
        varr8_resize(&seq, len);
        fread(seq.ele, sizeof(uint8_t), len, bubble_fp);
        fread(&cov, sizeof(double), 1, bubble_fp);

        ++n_seq;
        fprintf(out, ">seq%lu_len%lu_cov%u\n", n_seq, len, (uint16_t)(cov + 0.5));
        for (i = 0; i < len; ++i) {
            putc(num2ascii(seq.ele[i]), out);
            if ((i + 1) % LINE_LENGTH == 0)
                putc('\n', out);
        }
        if (i % LINE_LENGTH != 0)
            putc('\n', out);
    }

    varr8_destroy(&seq);
    fclose(out);
} 

int64_t scaffolder_delete_hetero_edge(scaffolder_t *sc)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t m;
    int64_t n_delete;
	double cov1;
	double cov2;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
    varr64_t ids;
	const double homo_cov_th = (int64_t)(sc->ave_cov * 1.5 + 0.5);
	const double hetero_cov_th = (int64_t)(sc->ave_cov * 0.75 + 0.5);

    varr64_init(&ids);
    n_delete = 0;

    for (i = 0; i < sc->n_node; ++i) {
        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {

                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if (e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if (e2->len  + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

				if (scaffolder_node_cov(sc, &(sc->node[i])) > homo_cov_th)
					continue;

				cov1 = scaffolder_node_cov(sc, n1);
				cov2 = scaffolder_node_cov(sc, n2);
				if (n1->len > n2->len) {
					n1 = n2;
					cov1 = cov2;
				}

				if (cov1 > hetero_cov_th || cov2 > hetero_cov_th)
					continue;

				++n_delete;
				n1->state |= SC_DEL;
				for (m = 0; m < n1->n_edge; ++m) {
					varr64_push(&ids, (int64_t)(n1 - sc->node) + 1);
					varr64_push(&ids, n1->edge[m].ed);
				}
				for (m = 0; m < n1->n_con; ++m)
					sc->con_inc[abs64(n1->con[m].id) - 1].id = 0;
            }
        }
    }

    scaffolder_delete_edges(sc, &ids);
    varr64_destroy(&ids);

    fprintf(stderr, "NUM_DELETE=%lu (COVERAGE_THRESHOLD)\n", n_delete);

    return n_delete;
}

int64_t scaffolder_crush_hetero_bubble(scaffolder_t *sc, FILE **bubble_fpp)
{
    int64_t i;
    int64_t j;
    int64_t k;
    int64_t m;
    int64_t n;
    int64_t redge;
    int64_t ledge;
    int64_t n_crush;
    double cov1;
    double cov2;
    sc_edge_t *e1;
    sc_edge_t *e2;
    sc_node_t *n1;
    sc_node_t *n2;
    sc_layout_t layout1;
    sc_layout_t layout2;
    sc_layout_t *lp;
    sc_layout_t work;
    varr64_t ids;
    varr8_t sca_seq;
	const double homo_cov_th = (int64_t)(sc->ave_cov * 1.5 + 0.5);
	const double hetero_cov_th = (int64_t)(sc->ave_cov * 0.75 + 0.5);

    varr64_init(&ids);
    sc_layout_init(&layout1);
    sc_layout_init(&layout2);
    sc_layout_init(&work);
    varr8_init(&sca_seq);
    n_crush = 0;

    scaffolder_detect_repeat(sc);

    for (i = 0; i < sc->n_node; ++i) {
        for (j = 0; j < sc->node[i].n_edge - 1; ++j) {
            for (k = j + 1; k < sc->node[i].n_edge; ++k) {
                e1 = &(sc->node[i].edge[j]);
                e2 = &(sc->node[i].edge[k]);

                if (e1->dir * e2->dir < 0)
                    continue;
                n1 = &(sc->node[abs64(e1->ed) - 1]);
                if ((n1->state & SC_DEL) || e1->len + n1->len <= e2->len)
                    continue;
                n2 = &(sc->node[abs64(e2->ed) - 1]);
                if ((n2->state & SC_DEL) || e2->len + n2->len <= e1->len)
                    continue;

                if (e1->dir > 0) {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed))

                        continue;
                }
                else {
                    if (abs64(e1->len + n1->len - e2->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e2->ed, e1->ed) ||
                        abs64(e2->len + n2->len - e1->len) <= sc->tol + scaffolder_get_sca_overlap(sc, e1->ed, e2->ed))

                        continue;
                }

                scaffolder_layout_nodes(sc, n1, &layout1, &work);
                scaffolder_layout_nodes(sc, n2, &layout2, &work);

                redge = ledge =  min_64(layout1.size, layout2.size);

                for (m = 0; m < ledge; ++m) {
                    if (layout1.part[m].id != layout2.part[m].id)
                        break;
                }
                if (m == 0 || m == ledge)
                    continue;
                ledge = m - 1;

				if (scaffolder_node_cov(sc, &(sc->node[abs64(layout1.part[ledge].id) - 1])) >= homo_cov_th)
					continue;

                for (m = 1; m <= redge; ++m) {
                    if (layout1.part[layout1.size - m].id != layout2.part[layout2.size - m].id)
                        break;
                }
                if (m == 1)
                    continue;
                redge = m - 1;

				if (scaffolder_node_cov(sc, &(sc->node[abs64(layout1.part[layout1.size - redge].id) - 1])) > homo_cov_th)
					continue;

                cov1 = scaffolder_layout_ave_cov(sc, &(layout1.part[ledge + 1]), layout1.size - redge - ledge - 1);
                cov2 = scaffolder_layout_ave_cov(sc, &(layout2.part[ledge + 1]), layout2.size - redge - ledge - 1);

                lp = cov1 < cov2 ? &layout1 : &layout2;

                if (redge + ledge + 1 >= lp->size || cov1 > hetero_cov_th || cov2 > hetero_cov_th)
                    continue;

                for (m = ledge + 1; m < lp->size - redge; ++m) {
                    n1 = &(sc->node[abs64(lp->part[m].id) - 1]);
                    for (n = 0; n < n1->n_edge; ++n) {
                        varr64_push(&ids, (int64_t)(n1 - sc->node) + 1);
                        varr64_push(&ids, n1->edge[n].ed);
                    }
                    for (n = 0; n < n1->n_con; ++n)
                        sc->con_inc[abs64(n1->con[n].id) - 1].id = 0;
                    n1->state |= SC_DEL;
                }

                scaffolder_layout2seq(sc, &(lp->part[ledge + 1]), lp->size - redge - ledge, &sca_seq);
                fwrite(&(sca_seq.size), sizeof(uint64_t), 1, *bubble_fpp);
                fwrite(sca_seq.ele, sizeof(uint8_t), sca_seq.size, *bubble_fpp);
                fwrite(cov1 < cov2 ? &cov1 : &cov2, sizeof(double), 1, *bubble_fpp);

                ++n_crush;
            }
        }
    }

    sc_layout_destroy(&layout1);
    sc_layout_destroy(&layout2);
    sc_layout_destroy(&work);
    varr8_destroy(&sca_seq);

    scaffolder_delete_edges(sc, &ids);
    varr64_destroy(&ids);

    for (i = 0; i < sc->n_node; ++i)
        sc->node[i].state &= ~SC_REP;

    fprintf(stderr, "NUM_CRUSH=%lu (COVERAGE_THRESHOLD)\n", n_crush);

    return n_crush;
}

double scaffolder_ave_len(scaffolder_t *sc)
{
	int64_t i;
	int64_t sum = 0;

	for (i = 0; i < sc->n_node; ++i)
		sum += sc->node[i].len;
	return (double)sum / sc->n_node;
}

int64_t scaffolder_get_similer_overlap(const scaffolder_t *sc, int64_t id1, int64_t id2)
{
    int64_t j;
    int64_t k;
	int64_t overlap = 0;
    int64_t max_ol;
    int64_t tol_mis;
    int64_t n_mis;
	double l_rate = 1.0;
	double max_mis_rate = 0.05;

    if (id1 > 0) {
        if (id2 > 0) {
			max_ol = min_64(sc->con[id1 - 1].len, sc->con[id2 - 1].len);
			for (j = max_ol; j >= sc->min_overlap; --j) {
				tol_mis = j * max_mis_rate + 0.5;
				n_mis = 0;
				for (k = 0; k < j; ++k) {
                    if (sc->con[id1 - 1].seq[sc->con[id1 - 1].len - k - 1] != sc->con[id2 - 1].seq[j - k - 1])
						continue;
					++n_mis;
					if (n_mis > tol_mis)
						break;
				}
				if (n_mis <= tol_mis && (double)n_mis / k < l_rate) {
					l_rate = (double)n_mis / k;
					overlap = j;
				}
			}
        }
        else {
			max_ol = min_64(sc->con[id1 - 1].len, sc->con[-id2 - 1].len);
			for (j = max_ol; j >= sc->min_overlap; --j) {
				tol_mis = j * max_mis_rate + 0.5;
				n_mis = 0;
				for (k = 0; k < j; ++k) {
                    if (sc->con[id1 - 1].seq[sc->con[id1 - 1].len - k - 1] != (3^(sc->con[-id2 - 1].seq[sc->con[-id2 - 1].len - j + k])))
						continue;
					++n_mis;
					if (n_mis > tol_mis)
						break;
				}
				if (n_mis <= tol_mis && (double)n_mis / k < l_rate) {
					l_rate = (double)n_mis / k;
					overlap = j;
				}
			}
        }
    }
    else {
        if (id2 > 0) {
			max_ol = min_64(sc->con[-id1 - 1].len, sc->con[id2 - 1].len);
			for (j = max_ol; j >= sc->min_overlap; --j) {
				tol_mis = j * max_mis_rate + 0.5;
				n_mis = 0;
				for (k = 0; k < j; ++k) {
                    if ((3^(sc->con[-id1 - 1].seq[k])) != sc->con[id2 - 1].seq[j - k - 1])
						continue;
					++n_mis;
					if (n_mis > tol_mis)
						break;
				}
				if (n_mis <= tol_mis && (double)n_mis / k < l_rate) {
					l_rate = (double)n_mis / k;
					overlap = j;
				}
			}
        }
        else {
			max_ol = min_64(sc->con[-id1 - 1].len, sc->con[-id2 - 1].len);
			for (j = max_ol; j >= sc->min_overlap; --j) {
				tol_mis = j * max_mis_rate + 0.5;
				n_mis = 0;
				for (k = 0; k < j; ++k) {
					if (sc->con[-id1 - 1].seq[k] != sc->con[-id2 - 1].seq[sc->con[-id2 - 1].len - j + k])
						continue;
					++n_mis;
					if (n_mis > tol_mis)
						break;
				}
				if (n_mis <= tol_mis && (double)n_mis / k < l_rate) {
					l_rate = (double)n_mis / k;
					overlap = j;
				}
            }
        }
    }

	return overlap;
}

int64_t scaffolder_get_similer_sca_overlap(scaffolder_t *sc, int64_t id1, int64_t id2)
{
    if (id1 > 0)
        id1 = sc->node[id1 - 1].con[sc->node[id1 - 1].n_con - 1].id;
    else
        id1 = -(sc->node[-id1 - 1].con[0].id);

    if (id2 > 0)
        id2 = sc->node[id2 - 1].con[0].id;
    else
        id2 = -(sc->node[-id2 - 1].con[sc->node[-id2 - 1].n_con - 1].id);
    
    return scaffolder_get_similer_overlap(sc, id1, id2);
}
