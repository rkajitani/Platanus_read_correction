#include "lcounter.h"

void lcounter_init(lcounter_t *lc, const uint64_t mem_lim)
{
    memset(lc, 0, sizeof(lcounter_t));
    lc->mem_lim = mem_lim;
}

void lcounter_destroy_table(lcounter_t *lc)
{
    if (lc->val_table != NULL)
        my_free(lc->val_table);
    lc->val_table = NULL;
    if (lc->key_table != NULL)
        my_free(lc->key_table);
    lc->key_table = NULL;
    lc->table_size = 0;
}

void lcounter_destroy(lcounter_t *lc)
{
    lcounter_destroy_table(lc);
    if (lc->occ_distr != NULL)
        my_free(lc->occ_distr);
    if (lc->kmer_fp != NULL)
        fclose(lc->kmer_fp);
    memset(lc, 0, sizeof(lcounter_t));
}

void lcounter_save_kmer_read(lcounter_t *lc, const uint64_t l_len, FILE *read_fp)
{
    uint64_t i;
    uint64_t bl;
    uint64_t *key;
    uint64_t index;
    uint64_t step;
    uint64_t n_used;
    uint64_t max;
    uint64_t fwd[(MAX_READ_LEN + 31) / 32];
    uint64_t rev[(MAX_READ_LEN + 31) / 32];
    uint16_t occ;
    seq_t seq;
    FILE *unstored_fp;
    FILE *tmp_fp;

    if (lc->val_table != NULL)
        lcounter_destroy_table(lc);
    lc->l_len = l_len;
    lc->bl = bl = (l_len + 31) / 32;

    lc->table_size = 1ull;
    for (i = 0; i < l_len*2; ++i) {
        lc->table_size *= 2ull;
        if (lc->table_size*(sizeof(uint16_t) + bl*sizeof(uint64_t)) > lc->mem_lim / 2)
            break;
    }
    lc->index_len = i;
    lc->table_size /= 2ull;

    lc->val_table = (uint16_t *)my_calloc(lc->table_size, sizeof(uint16_t));
    lc->key_table = (uint64_t *)my_malloc(lc->table_size * bl * sizeof(uint64_t));

    if (lc->kmer_fp != NULL)
        fclose(lc->kmer_fp);
    lc->kmer_fp = tmpfile_open();
    check_tmpfile(lc->kmer_fp);
    unstored_fp = tmpfile_open();
    check_tmpfile(unstored_fp);

    rewind(read_fp);
    n_used = 0;
    while (seq_read(&seq, read_fp)) {
        ++lc->len_distr[seq.len];
        if (seq.len < l_len)
            continue;
        seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
        seq.n_unknown = 0;
        for (i = 0; i < l_len - 1; ++i) {
            lmer_set(fwd, i+1, seq.base[i]);
            lmer_set(rev, l_len-i-2, 3^seq.base[i]);
        }
        for (i = 0; i < seq.len - l_len + 1; ++i) {
            lmer_lshift(fwd, l_len);
            lmer_set(fwd, l_len-1, seq.base[i+l_len-1]);
            lmer_rshift(rev, l_len);
            lmer_set(rev, 0, 3^seq.base[i+l_len-1]);
            if (seq.unknown_pos[seq.n_unknown] < i + l_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                continue;
            }
            key = lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev;
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
            if (!(lc->val_table[index])) {
                if (n_used <= lc->table_size*MAX_LOAD_FACTOR) {
                    lc->val_table[index] = 1;
                    lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                    ++n_used;
                }
                else
                    fwrite(key, sizeof(uint64_t), bl, unstored_fp);
            }
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                if (lc->val_table[index] < UINT16_MAX-1)
                    ++lc->val_table[index];
            }
            else {
                step = lmer_rehash(key, l_len, lc->index_len);
                while (1) {
                    index = (index+step) & (lc->table_size-1);
                    if (!(lc->val_table[index])) {
                        if (n_used <= lc->table_size*MAX_LOAD_FACTOR) {
                            lc->val_table[index] = 1;
                            lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                            ++n_used;
                        }
                        else
                            fwrite(key, sizeof(uint64_t), bl, unstored_fp);
                        break;
                    }
                    else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                        if (lc->val_table[index] < UINT16_MAX-1)
                            ++lc->val_table[index];
                        break;
                    }
                }
            }
        }
    }

    if (lc->occ_distr != NULL)
        my_free(lc->occ_distr);
    lc->occ_distr = (uint64_t *)my_calloc(UINT16_MAX, sizeof(uint64_t));

    key = fwd;
    max = 0;
    while (1) {
        for (i = 0; i < lc->table_size; ++i) {
            if(lc->val_table[i]) {
                occ = lc->val_table[i];
                fwrite(&(lc->key_table[i * bl]), sizeof(uint64_t), bl, lc->kmer_fp);
                fwrite(&occ, sizeof(uint16_t), 1, lc->kmer_fp);
                ++lc->occ_distr[occ];
                lc->val_table[i] = 0;
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
        while (fread(key, sizeof(uint64_t), bl, unstored_fp)) {
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size-1);
            if (!(lc->val_table[index])) {
                if (n_used <= lc->table_size*MAX_LOAD_FACTOR) {
                    lc->val_table[index] = 1;
                    lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                    ++n_used;
                }
                else
                    fwrite(key, sizeof(uint64_t), bl, tmp_fp);
            }
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                if (lc->val_table[index] < UINT16_MAX-1)
                    ++lc->val_table[index];
            }
            else {
                step = lmer_rehash(key, l_len, lc->index_len);
                while (1) {
                    index = (index+step) & (lc->table_size-1);
                    if (!(lc->val_table[index])) {
                        if (n_used <= lc->table_size*MAX_LOAD_FACTOR) {
                            lc->val_table[index] = 1;
                            lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                            ++n_used;
                        }
                        else
                            fwrite(key, sizeof(uint64_t), bl, tmp_fp);
                        break;
                    }
                    else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                        if (lc->val_table[index] < UINT16_MAX-1)
                            ++lc->val_table[index];
                        break;
                    }
                }
            }
        }
        fclose(unstored_fp);
        unstored_fp = tmp_fp;
    }

    fclose(unstored_fp);

    for (i = UINT16_MAX - 1; i > 0; --i) {
        if (lc->occ_distr[i] > 0)
            break;
    }
    lc->max_occ = i;
}

void lcounter_load_kmer(lcounter_t *lc, const uint64_t min_occ)
{
    uint64_t n_key;
    uint64_t bl;
    uint64_t l_len;
    uint64_t index;
    uint64_t step;
    uint16_t occ;
    lmer_t key;
    FILE *new_kmer_fp;

    fprintf(stderr, "loading kmers...\n");

    bl = lc->bl;
    l_len = lc->l_len; 

    n_key = 0;
    rewind(lc->kmer_fp);
    while(fread(key, sizeof(uint64_t), bl, lc->kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, lc->kmer_fp);
        if (occ >= min_occ)
            ++n_key;
    }

    lc->table_size = 2;
    lc->index_len = 1;
    while ((double)lc->table_size * MAX_LOAD_FACTOR < n_key) {
        lc->table_size *= 2;
        ++lc->index_len;
    }
    while (2 * lc->table_size*(sizeof(uint16_t) + bl*sizeof(uint64_t)) < lc->mem_lim / 2) {
        lc->table_size *= 2;
        ++lc->index_len;
    }

    check_mem_usage(lc->table_size*(sizeof(uint64_t)+sizeof(uint16_t)), lc->mem_lim);
/*
    fprintf(stderr, "BUCKET_SIZE=%lu, ", lc->table_size); 
    fprintf(stderr, "MIN_OCCURRENCE=%lu, ", min_occ);
    fprintf(stderr, "NUM_OCCUPIED_SLOT=%lu, ", n_key / bl);  
    fprintf(stderr, "LOAD_FACTOR=%f, ", (double)n_key / (double)lc->table_size);
    fprintf(stderr, "MEM_USAGE=%luB\n", lc->table_size * (sizeof(uint16_t) + bl*sizeof(uint64_t))); 
*/

    lc->key_table = my_realloc(lc->key_table, lc->table_size * bl * sizeof(uint64_t));
    lc->val_table = my_realloc(lc->val_table, lc->table_size*sizeof(uint16_t));
    memset(lc->val_table, 0, lc->table_size*sizeof(uint16_t));
    new_kmer_fp = tmpfile_open();
    check_tmpfile(new_kmer_fp);

    rewind(lc->kmer_fp);
    while (fread(key, sizeof(uint64_t), bl, lc->kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, lc->kmer_fp);
        if (occ < min_occ) {
            fwrite(key, sizeof(uint64_t), bl, new_kmer_fp);
            fwrite(&occ, sizeof(uint16_t), 1, new_kmer_fp);
            continue;
        }

        index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
        if (!(lc->val_table[index])) {
            lc->val_table[index] = occ;
            lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
        }
        else {
            step = lmer_rehash(key, l_len, lc->index_len);
            index = (index+step) & (lc->table_size-1);
            while (lc->val_table[index])
                index = (index+step) & (lc->table_size-1);
            lc->val_table[index] = occ;
            lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
        }
    }

    fclose(lc->kmer_fp);
    lc->kmer_fp = new_kmer_fp;
}

void lcounter_load_contig(lcounter_t *lc, const uint64_t l_len, const double occ_ratio, FILE *contig_fp)
{
    uint64_t i;
    uint64_t bl;
    uint64_t *key;
    uint64_t index;
    uint64_t step;
    uint64_t total_key;
    uint64_t clen;
    lmer_t fwd;
    lmer_t rev;
    uint16_t occ;
    varr8_t cseq;

    fprintf(stderr, "K = %lu, loading kmers from contigs...\n", l_len);

    lc->l_len = l_len;
    lc->bl = bl = (l_len + 31) / 32;

    total_key = 0;
    rewind(contig_fp);
    while (fread(&clen, sizeof(uint64_t), 1, contig_fp)) {
        fseek(contig_fp, sizeof(uint16_t) + ((clen+3)/4)*sizeof(uint8_t), SEEK_CUR);
        if (clen >= l_len)
            total_key += clen - l_len + 1;
    }

    if (lc->val_table != NULL)
        lcounter_destroy_table(lc);
    lc->table_size = 1;
    for (i = 0; i < l_len*2; ++i) {
        lc->table_size *= 2ull;
        if (lc->table_size*(sizeof(uint16_t) + bl*sizeof(uint64_t)) > lc->mem_lim / 2)
            break;
    }
    lc->index_len = i;
    lc->table_size /= 2ull;
    while ((double)lc->table_size * MAX_LOAD_FACTOR < total_key) {
        lc->table_size *= 2;
        ++lc->index_len;
    }

    lc->val_table = (uint16_t *)my_calloc(lc->table_size, sizeof(uint16_t));
    lc->key_table = (uint64_t *)my_malloc(lc->table_size * bl * sizeof(uint64_t));

    rewind(contig_fp);
    varr8_init(&cseq);
    while (fread(&clen, sizeof(uint64_t), 1, contig_fp)) {
        fread(&occ, sizeof(uint16_t), 1, contig_fp);
        occ = occ * occ_ratio + 0.5;
        varr8_resize(&cseq, (clen+3)/4);
        fread(cseq.ele, sizeof(uint8_t), (clen+3)/4, contig_fp);
        if (clen < l_len)
            continue;
        for (i = 0; i < l_len - 1; ++i) {
            lmer_set(fwd, i+1, bseq_get(cseq.ele, i));
            lmer_set(rev, l_len-i-2, 3^bseq_get(cseq.ele, i));
        }
        for (i = 0; i < clen - l_len + 1; ++i) {
            lmer_lshift(fwd, l_len);
            lmer_set(fwd, l_len-1, bseq_get(cseq.ele, i+l_len-1));
            lmer_rshift(rev, l_len);
            lmer_set(rev, 0, 3^bseq_get(cseq.ele, i+l_len-1));
            key = lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev;
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
            if (!lc->val_table[index]) {
                lc->val_table[index] = occ;
                lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
            }
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                if (lc->val_table[index] < occ)
                    lc->val_table[index] = occ;
            }
            else {
                step = lmer_rehash(key, l_len, lc->index_len);
                index = (index + step) & (lc->table_size-1);
                while (lc->val_table[index] && lmer_cmp(&(lc->key_table[index * bl]), key, l_len) != 0)
                    index = (index+step) & (lc->table_size-1);
                if (lc->val_table[index] < occ)
                    lc->val_table[index] = occ;
                lmer_cpy(&(lc->key_table[index * bl]), key, l_len);
            }
        }
    }
    varr8_destroy(&cseq);
}

void lcounter_save_additional_kmer_read(lcounter_t *lc, const uint64_t l_len, FILE *read_fp)
{
    uint64_t i;
    uint64_t bl;
    uint64_t *key;
    uint64_t index;
    uint64_t step;
    uint64_t n_used;
    lmer_t fwd;
    lmer_t rev;
    seq_t seq = {0};
    FILE *unstored_fp;
    FILE *tmp_fp;

    fprintf(stderr, "K = %lu, saving additional kmers from reads...\n", l_len);

    bl = lc->bl;

    if (lc->kmer_fp != NULL)
        fclose(lc->kmer_fp);
    lc->kmer_fp = tmpfile_open();
    check_tmpfile(lc->kmer_fp);
    unstored_fp = tmpfile_open();
    check_tmpfile(unstored_fp);

    rewind(read_fp);
    while (seq_read(&seq, read_fp)) {
        if (seq.len < l_len)
            continue;
        seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
        seq.n_unknown = 0;
        for (i = 0; i < l_len - 1; ++i) {
            lmer_set(fwd, i+1, seq.base[i]);
            lmer_set(rev, l_len-i-2, 3^seq.base[i]);
        }
        for (i = 0; i < seq.len - l_len + 1; ++i) {
            lmer_lshift(fwd, l_len);
            lmer_set(fwd, l_len-1, seq.base[i+l_len-1]);
            lmer_rshift(rev, l_len);
            lmer_set(rev, 0, 3^seq.base[i+l_len-1]);
            if (seq.unknown_pos[seq.n_unknown] < i + l_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                continue;
            }
            key = lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev;
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
            if (!(lc->val_table[index]))
                fwrite(key, sizeof(uint64_t), bl, unstored_fp);
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) != 0) {
                step = lmer_rehash(key, l_len, lc->index_len);
                while (1) {
                    index = (index+step) & (lc->table_size-1);
                    if (!(lc->val_table[index])) {
                        fwrite(key, sizeof(uint64_t), bl, unstored_fp);
                        break;
                    }
                    else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                        break;
                    }
                }
            }
        }
    }

    for (i = 0; i < lc->table_size; ++i) {
        if (lc->val_table[i]) {
            fwrite(&(lc->key_table[i * bl]), sizeof(uint64_t), bl, lc->kmer_fp);
            fwrite(&(lc->val_table[i]), sizeof(uint16_t), 1, lc->kmer_fp);
            lc->val_table[i] = 0;
        }
    }

    key = fwd;

    while (ftell(unstored_fp) / sizeof(uint64_t)) {
        rewind(unstored_fp);
        n_used = 0;
        tmp_fp = tmpfile_open();
        check_tmpfile(tmp_fp);
        while (fread(key, sizeof(uint64_t), bl, unstored_fp)) {
            index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size-1);
            if (!(lc->val_table[index])) {
                if (n_used <= lc->table_size*MAX_LOAD_FACTOR) {
                    lc->val_table[index] = 1;
                    lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                    ++n_used;
                }
                else
                    fwrite(key, sizeof(uint64_t), bl, tmp_fp);
            }
            else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                if (lc->val_table[index] < UINT16_MAX-1)
                    ++lc->val_table[index];
            }
            else {
                step = lmer_rehash(key, l_len, lc->index_len);
                while (1) {
                    index = (index+step) & (lc->table_size-1);
                    if (!(lc->val_table[index])) {
                        if (n_used <= lc->table_size*MAX_LOAD_FACTOR) {
                            lc->val_table[index] = 1;
                            lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                            ++n_used;
                        }
                        else
                            fwrite(key, sizeof(uint64_t), bl, tmp_fp);
                        break;
                    }
                    else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                        if (lc->val_table[index] < UINT16_MAX-1)
                            ++lc->val_table[index];
                        break;
                    }
                }
            }
        }
        for (i = 0; i < lc->table_size; ++i) {
            if (!(lc->val_table[i]))
                continue;
            fwrite(&(lc->key_table[i * bl]), sizeof(uint64_t), bl, lc->kmer_fp);
            fwrite(&(lc->val_table[i]), sizeof(uint16_t), 1, lc->kmer_fp);
            lc->val_table[i] = 0;
        }
        fclose(unstored_fp);
        unstored_fp = tmp_fp;
    }

    fclose(unstored_fp);
}

void lcounter_save_kmer_read_mt(lcounter_t *lc, const uint64_t l_len, FILE **read_fp, const uint64_t n_thread)
{
    int64_t i;
    int64_t j;
    uint16_t occ;
    FILE *unstored_fp[MAX_THREAD];
    FILE *tmp_fp[MAX_THREAD];
    lcounter_t lcs[MAX_THREAD];
    pthread_t thread[MAX_THREAD];
    pthread_mutex_t *mutex;
    arg_lcounter_lock_and_count_t arg[MAX_THREAD];

    omp_set_num_threads(n_thread);

#    pragma omp parallel for schedule(static, 1)
    for (i = 0; i < n_thread; ++i) {
        lcounter_init(&(lcs[i]), lc->mem_lim / n_thread);
        lcounter_save_kmer_read(&(lcs[i]), l_len, read_fp[i]);
        lcounter_destroy_table(&(lcs[i]));
    }

    if (lc->val_table != NULL)
        lcounter_destroy_table(lc);
    lc->l_len = l_len;
    lc->bl = (l_len + 31) / 32;

    lc->table_size = 1ull;
    for (i = 0; i < l_len*2; ++i) {
        lc->table_size *= 2ull;
        if (lc->table_size*(sizeof(uint16_t) + lc->bl*sizeof(uint64_t)) + (lc->table_size/GRANULE+1) * sizeof(pthread_mutex_t)> lc->mem_lim / 2)
            break;
    }
    lc->index_len = i;
    lc->table_size /= 2ull;

    lc->val_table = (uint16_t *)my_calloc(lc->table_size, sizeof(uint16_t));
    lc->key_table = (uint64_t *)my_malloc(lc->table_size * lc->bl * sizeof(uint64_t));

    if (lc->kmer_fp != NULL)
        fclose(lc->kmer_fp);
    lc->kmer_fp = tmpfile_open();
    check_tmpfile(lc->kmer_fp);

    mutex = (pthread_mutex_t *)my_malloc((lc->table_size / GRANULE + 1) * sizeof(pthread_mutex_t));
    for (i = 0; i < lc->table_size/GRANULE+1; ++i)
        pthread_mutex_init(&mutex[i], NULL);

    for (i = 0; i < n_thread; ++i) {
        arg[i].lc = lc;
        arg[i].mutex = mutex;
        arg[i].in_kmer_fp = lcs[i].kmer_fp;
        unstored_fp[i] = tmpfile_open();
        check_tmpfile(unstored_fp[i]);
        arg[i].unstored_fp = unstored_fp[i];
        pthread_create(&thread[i], NULL, (void *)lcounter_lock_and_count1, &arg[i]);
    }
    for (i = 0; i < n_thread; ++i) {
        pthread_join(thread[i], NULL); 
        for (j = 1; j <= MAX_READ_LEN; ++j)
            lc->len_distr[j] += lcs[i].len_distr[j];
        lcounter_destroy(&(lcs[i]));
    }

    if (lc->occ_distr != NULL)
        my_free(lc->occ_distr);
    lc->occ_distr = (uint64_t *)my_calloc(UINT16_MAX, sizeof(uint64_t));

    while (1) {
        for (i = 0; i < lc->table_size; ++i) {
            if (lc->val_table[i]) {
                occ = lc->val_table[i];
                fwrite(&lc->key_table[i * lc->bl], sizeof(uint64_t), lc->bl, lc->kmer_fp);
                fwrite(&occ, sizeof(uint16_t), 1, lc->kmer_fp);
                ++lc->occ_distr[occ];
                lc->val_table[i] = 0;
            }
        }

        for (i = 0; i < n_thread; ++i) {
            if(ftell(unstored_fp[i]) / sizeof(uint64_t))
                break;
        }
        if (i == n_thread)
            break;

        for (i = 0; i < n_thread; ++i) {
            tmp_fp[i] = tmpfile_open();
            check_tmpfile(tmp_fp[i]);
            arg[i].in_kmer_fp = unstored_fp[i];
            arg[i].unstored_fp = tmp_fp[i];
            pthread_create(&thread[i], NULL, (void *)lcounter_lock_and_count1, &arg[i]);
        }
        for (i = 0; i < n_thread; ++i) {
            pthread_join(thread[i], NULL); 
            fclose(unstored_fp[i]);
            unstored_fp[i] = tmp_fp[i];
        }
    }

    for (i = 0; i < n_thread; ++i)
        fclose(unstored_fp[i]);

    my_free(mutex);

    for (i = UINT16_MAX - 1; i > 0; --i) {
        if (lc->occ_distr[i] > 0)
            break;
    }
    lc->max_occ = i;
}

void lcounter_save_additional_kmer_read_mt(lcounter_t *lc, const uint64_t l_len, FILE **read_fp, const uint64_t n_thread)
{
    int64_t i;
    uint64_t j;
    uint64_t bl;
    uint64_t *key;
    uint64_t index;
    uint64_t step;
    lmer_t fwd;
    lmer_t rev;
    seq_t seq;
    FILE *unstored_fp[MAX_THREAD];
    FILE *tmp_fp[MAX_THREAD];
    pthread_t thread[MAX_THREAD];
    pthread_mutex_t *mutex;
    arg_lcounter_lock_and_count_t arg[MAX_THREAD];

    fprintf(stderr, "K = %lu, saving additional kmers from reads...\n", l_len);

    bl = lc->bl;
    if (lc->kmer_fp != NULL)
        fclose(lc->kmer_fp);
    lc->kmer_fp = tmpfile_open();
    check_tmpfile(lc->kmer_fp);

    omp_set_num_threads(n_thread);

#    pragma omp parallel for private(j, index, step, key, fwd, rev, seq) schedule(static, 1)
    for (i = 0; i < n_thread; ++i) {
        unstored_fp[i] = tmpfile_open();
        check_tmpfile(unstored_fp[i]);
        rewind(read_fp[i]);
        while (seq_read(&seq, read_fp[i])) {
            if (seq.len < l_len)
                continue;
            seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
            seq.n_unknown = 0;
            for (j = 0; j < l_len - 1; ++j) {
                lmer_set(fwd, j+1, seq.base[j]);
                lmer_set(rev, l_len-j-2, 3^seq.base[j]);
            }
            for (j = 0; j < seq.len - l_len + 1; ++j) {
                lmer_lshift(fwd, l_len);
                lmer_set(fwd, l_len-1, seq.base[j+l_len-1]);
                lmer_rshift(rev, l_len);
                lmer_set(rev, 0, 3^seq.base[j+l_len-1]);
                if (seq.unknown_pos[seq.n_unknown] < j + l_len) {
                    if (seq.unknown_pos[seq.n_unknown] <= j)
                        ++seq.n_unknown;
                    continue;
                }
                key = lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev;
                index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
                if (!(lc->val_table[index]))
                    fwrite(key, sizeof(uint64_t), bl, unstored_fp[i]);
                else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) != 0) {
                    step = lmer_rehash(key, l_len, lc->index_len);
                    while (1) {
                        index = (index+step) & (lc->table_size-1);
                        if (!(lc->val_table[index])) {
                            fwrite(key, sizeof(uint64_t), bl, unstored_fp[i]);
                            break;
                        }
                        else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                            break;
                        }
                    }
                }
            }
        }
    }

    mutex = (pthread_mutex_t *)my_malloc((lc->table_size / GRANULE + 1) * sizeof(pthread_mutex_t));
    for (i = 0; i < lc->table_size/GRANULE+1; ++i)
        pthread_mutex_init(&mutex[i], NULL);

    for (i = 0; i < lc->table_size; ++i) {
        if (lc->val_table[i]) {
            fwrite(&(lc->key_table[i * bl]), sizeof(uint64_t), bl, lc->kmer_fp);
            fwrite(&(lc->val_table[i]), sizeof(uint16_t), 1, lc->kmer_fp);
            lc->val_table[i] = 0;
        }
    }

    while (1) {
        for (i = 0; i < n_thread; ++i) {
            if(ftell(unstored_fp[i]) / sizeof(uint64_t))
                break;
        }
        if (i == n_thread)
            break;

        for (i = 0; i < n_thread; ++i) {
            tmp_fp[i] = tmpfile_open();
            check_tmpfile(tmp_fp[i]);
            arg[i].lc = lc;
            arg[i].mutex = mutex;
            arg[i].in_kmer_fp = unstored_fp[i];
            arg[i].unstored_fp = tmp_fp[i];
            pthread_create(&thread[i], NULL, (void *)lcounter_lock_and_count2, &arg[i]);
        }
        for (i = 0; i < n_thread; ++i) {
            pthread_join(thread[i], NULL); 
            fclose(unstored_fp[i]);
            unstored_fp[i] = tmp_fp[i];
        }

        for (i = 0; i < lc->table_size; ++i) {
            if (lc->val_table[i]) {
                fwrite(&lc->key_table[i * bl], sizeof(uint64_t), bl, lc->kmer_fp);
                fwrite(&(lc->val_table[i]), sizeof(uint16_t), 1, lc->kmer_fp);
                lc->val_table[i] = 0;
            }
        }
    }

    for (i = 0; i < n_thread; ++i)
        fclose(unstored_fp[i]);

    my_free(mutex);
}

void lcounter_lock_and_count1(const arg_lcounter_lock_and_count_t *arg)
{
    uint64_t bl;
    uint64_t index;
    uint16_t occ;
    FILE *in_kmer_fp;
    FILE *unstored_fp;
    lmer_t key;
    lcounter_t *lc;
    pthread_mutex_t *mutex;

    lc = arg->lc;
    in_kmer_fp = arg->in_kmer_fp;
    unstored_fp = arg->unstored_fp;
    mutex = arg->mutex;
    bl = (lc->l_len + 31) / 32;

    rewind(in_kmer_fp);
    while (fread(&key, sizeof(uint64_t), bl, in_kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, in_kmer_fp);
        index = lmer_hash(key, lc->l_len, lc->index_len) & (lc->table_size - 1);

        pthread_mutex_lock(&mutex[index/GRANULE]);
        if (!(lc->val_table[index])) {
            lc->val_table[index] = occ;
            lmer_cpy(&(lc->key_table[index * bl]), key, lc->l_len);;
        }
        else if (lmer_cmp(&(lc->key_table[index * bl]), key, lc->l_len) == 0) {
            if ((uint64_t)occ + lc->val_table[index] < UINT16_MAX - 1)
                lc->val_table[index] += occ;
            else
                lc->val_table[index] = UINT16_MAX - 1;
        }
        else {
            fwrite(key, sizeof(uint64_t), bl, unstored_fp);
            fwrite(&occ, sizeof(uint16_t), 1, unstored_fp);
        }
        pthread_mutex_unlock(&mutex[index/GRANULE]);
    }
}

void lcounter_lock_and_count2(const arg_lcounter_lock_and_count_t *arg)
{
    uint64_t bl;
    uint64_t index;
    FILE *in_kmer_fp;
    FILE *unstored_fp;
    lmer_t key;
    lcounter_t *lc;
    pthread_mutex_t *mutex;

    lc = arg->lc;
    in_kmer_fp = arg->in_kmer_fp;
    unstored_fp = arg->unstored_fp;
    mutex = arg->mutex;
    bl = (lc->l_len + 31) / 32;

    rewind(in_kmer_fp);
    while (fread(&key, sizeof(uint64_t), bl, in_kmer_fp)) {
        index = lmer_hash(key, lc->l_len, lc->index_len) & (lc->table_size - 1);

        pthread_mutex_lock(&mutex[index/GRANULE]);
        if (!(lc->val_table[index])) {
            lc->val_table[index] = 1;
            lmer_cpy(&(lc->key_table[index * bl]), key, lc->l_len);;
        }
        else if (lmer_cmp(&(lc->key_table[index * bl]), key, lc->l_len) == 0) {
            if (lc->val_table[index] < UINT16_MAX - 1)
                ++lc->val_table[index];
        }
        else {
            fwrite(key, sizeof(uint64_t), bl, unstored_fp);
        }
        pthread_mutex_unlock(&mutex[index/GRANULE]);
    }
}

uint64_t lcounter_sorted_lmer_list(const lcounter_t *lc, uint64_t **buf)
{
    uint64_t i;
    uint64_t buf_size;

    buf_size = 0;
    for (i = 0; i < lc->table_size; ++i) {
        if (!lc->val_table[i])
            continue;
        ++buf_size;
    }
    *buf = (uint64_t *)my_malloc(buf_size * lc->bl * sizeof(uint64_t));

    buf_size = 0;
    for (i = 0; i < lc->table_size; ++i) {
        if (!lc->val_table[i])
            continue;
        lmer_cpy(&((*buf)[buf_size * lc->bl]), &(lc->key_table[i * lc->bl]), lc->l_len);    
        ++buf_size;
    }
    g_l_len = lc->l_len;
    qsort(*buf, buf_size, sizeof(uint64_t) * lc->bl, cmp_lmer_uniq);

    return buf_size;
}

void lcounter_extract_read(const lcounter_t *lc, FILE **read_fpp)
{
    uint64_t i;
    uint64_t l_len;
    lmer_t fwd;
    lmer_t rev;
    uint64_t *key;
    seq_t seq;
    seq_t tmp_seq;
    FILE *new_read_fp;

    l_len = lc->l_len;
    
    new_read_fp = tmpfile_open();
    check_tmpfile(new_read_fp);
    rewind(*read_fpp);
    while (seq_read(&seq, *read_fpp)) {
        if (seq.len < l_len)
            continue;
        tmp_seq = seq;
        tmp_seq.unknown_pos[tmp_seq.n_unknown] = MAX_READ_LEN+1;
        tmp_seq.n_unknown = 0;
        for (i = 0; i < l_len - 1; ++i) {
            lmer_set(fwd, i+1, tmp_seq.base[i]);
            lmer_set(rev, l_len-i-2, 3^tmp_seq.base[i]);
        }
        for (i = 0; i < tmp_seq.len - l_len + 1; ++i) {
            lmer_lshift(fwd, l_len);
            lmer_set(fwd, l_len-1, tmp_seq.base[i+l_len-1]);
            lmer_rshift(rev, l_len);
            lmer_set(rev, 0, 3^tmp_seq.base[i+l_len-1]);
            if (tmp_seq.unknown_pos[tmp_seq.n_unknown] < i + l_len) {
                if (tmp_seq.unknown_pos[tmp_seq.n_unknown] <= i)
                    ++tmp_seq.n_unknown;
                continue;
            }

            key = lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev;
            if (lcounter_occ(lc, key) > 0) {
                seq_write(&seq, new_read_fp);
                break;
            }
        }
    }

    fclose(*read_fpp);
    *read_fpp = new_read_fp;
}

double lcounter_estimate_error_rate(lcounter_t *lc, const uint64_t l_len, FILE *read_fp)
{
    double err_rate; 
    double ave_len;
    double cov_cut;
    lcounter_t tmp_lc;
    lcounter_init(&tmp_lc, lc->mem_lim);

    fprintf(stderr, "K = %lu, saving kmers from reads...\n", l_len);
    lcounter_save_kmer_read(lc, l_len, read_fp);
    
    fprintf(stderr, "K = %lu, saving kmers from reads...\n", l_len - 1);
    lcounter_save_kmer_read(&tmp_lc, l_len - 1, read_fp);
    cov_cut = (double)get_left_minimal(tmp_lc.occ_distr, tmp_lc.max_occ);
    lcounter_load_kmer(&tmp_lc, (uint64_t)(cov_cut + 0.5));

    ave_len = get_ave(lc->len_distr, 0, MAX_READ_LEN);
    err_rate = lcounter_calculate_error_rate(lc, &tmp_lc, ave_len);
    lcounter_destroy_table(lc);
    lcounter_destroy(&tmp_lc);
    return err_rate;
}

double lcounter_estimate_error_rate_mt(lcounter_t *lc, const uint64_t l_len, FILE **read_fp, const uint64_t n_thread)
{
    double err_rate; 
    double ave_len;
    double cov_cut;
    lcounter_t tmp_lc;

    lcounter_init(&tmp_lc, lc->mem_lim);

    fprintf(stderr, "K = %lu, saving kmers from reads...\n", l_len);
    lcounter_save_kmer_read_mt(lc, l_len, read_fp, n_thread);
    
    fprintf(stderr, "K = %lu, saving kmers from reads...\n", l_len - 1);
    lcounter_save_kmer_read_mt(&tmp_lc, l_len - 1, read_fp, n_thread);
    cov_cut = (double)get_left_minimal(tmp_lc.occ_distr, tmp_lc.max_occ);
    lcounter_load_kmer(&tmp_lc, (uint64_t)(cov_cut + 0.5));

    ave_len = get_ave(lc->len_distr, 0, MAX_READ_LEN);
    err_rate = lcounter_calculate_error_rate(lc, &tmp_lc, ave_len);
    lcounter_destroy_table(lc);
    lcounter_destroy(&tmp_lc);
    return err_rate;
}

double lcounter_calculate_error_rate(const lcounter_t *lc1, const lcounter_t *lc2, const double ave_len)
{
    uint64_t i;
    uint64_t cov_cut;
    uint64_t base = 0;
    uint64_t tmp;
    uint64_t occ2;
    uint16_t occ1;
    double sum1;
    double sum2;
    lmer_t key;
    lmer_t fwd;
    lmer_t rev;

    cov_cut = get_left_minimal(lc1->occ_distr, lc1->max_occ);

    sum1 = sum2 = 0.0;
    rewind(lc1->kmer_fp);
    while (fread(key, sizeof(uint64_t), lc1->bl, lc1->kmer_fp)) {
        fread(&occ1, sizeof(uint16_t), 1, lc1->kmer_fp);
        if (occ1 < cov_cut || occ1 == UINT16_MAX - 1)
            continue;

        for (i = 0; i < lc2->l_len - 1; ++i) {
            base = lmer_get(key, i);
            lmer_set(fwd, i, base);
            lmer_set(rev, lc2->l_len - i - 1, 3^base);
        }

        occ2 = lcounter_occ(lc2, lmer_cmp(fwd, rev, lc2->l_len) < 0 ? fwd : rev);
        if (occ2 == 0 || occ2 == UINT16_MAX - 1)
            continue;
        for (i = (base + 1) % 4; i != base; i = (i + 1) % 4) {
            lmer_set(fwd, 0, i);
            lmer_set(rev, lc2->l_len - 1, 3^i);
            if (lcounter_occ(lc2, lmer_cmp(fwd, rev, lc2->l_len) < 0 ? fwd : rev) > 0)
                break;
        }
        if (i != base)
            continue;

        base = lmer_get(key, lc1->l_len - 1);
        lmer_lshift(fwd, lc2->l_len);
        lmer_set(fwd, lc2->l_len - 1, base);
        lmer_rshift(rev, lc2->l_len);
        lmer_set(rev, 0, 3^base);
        tmp = lcounter_occ(lc2, lmer_cmp(fwd, rev, lc2->l_len) < 0 ? fwd : rev);
        if (tmp == 0 || tmp == UINT16_MAX - 1)
            continue;
        occ2 += tmp;
        for (i = (base + 1) % 4; i != base; i = (i + 1) % 4) {
            lmer_set(fwd, lc2->l_len - 1, i);
            lmer_set(rev, 0, 3^i);
            if (lcounter_occ(lc2, lmer_cmp(fwd, rev, lc2->l_len) < 0 ? fwd : rev) > 0)
                break;
        }
        if (i != base)
            continue;

        sum1 += (double)occ1;
        sum2 += (double)occ2 / 2.0;
    }
    
//fprintf(stderr, "sum1=%f, sum2=%f\n", sum1, sum2);
    return (1.0 - (sum1 / sum2 * (ave_len - lc2->l_len + 1) / (ave_len - lc1->l_len + 1)));
}

uint64_t lcounter_get_coverage_cutoff(const lcounter_t *lc, const double err_rate)
{
    uint64_t i;
    uint64_t ave_len;
    uint64_t total_len;
    uint64_t n_read;
    uint64_t total_kmer;
    uint64_t cum;
    double e_err_kmer;

    total_len = n_read = 0;
    for (i = 1; i <= MAX_READ_LEN; ++i) {
        n_read += lc->len_distr[i];
        total_len += i * lc->len_distr[i];
    }
    ave_len = (uint64_t)((double)total_len / n_read + 0.5);

    total_kmer = 0;
    for (i = 1; i <= lc->max_occ; ++i)
        total_kmer += i * lc->occ_distr[i];
    
    e_err_kmer = 0.0;
    for (i = 0; i < ave_len; ++i) {
        if (i >= lc->l_len - 1 && i + lc->l_len <= ave_len)
            e_err_kmer += (double)(lc->l_len);
        else
            e_err_kmer += (double)min_u64(i + 1, ave_len - i);
    }
    e_err_kmer /= (double)ave_len;
    e_err_kmer *= (double)total_len;
    e_err_kmer *= err_rate;

fprintf(stderr, "EXPECTED_ERRONEOUS_KMER = %f\n", e_err_kmer);

    cum = 0;
    for (i = 1; i <= lc->max_occ; ++i) {
        if (e_err_kmer - cum < 0.5 * (i * lc->occ_distr[i]))
            break;
        cum += i * lc->occ_distr[i];
    }

    return i;
}

FILE *sorted_lkey_from_lmer_file(FILE *lmer_fp, const uint64_t l_len, const uint64_t min_occ)
{
    uint16_t occ;
    uint64_t bl;
    lmer_t key;
    varr64_t buf;
    FILE *sorted_key_fp;

    bl = (l_len + 31) / 32;

    varr64_init(&buf);
    rewind(lmer_fp);
    while(fread(key, sizeof(uint64_t), bl, lmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, lmer_fp);
        if (occ < min_occ)
            continue;
        varr64_resize(&buf, buf.size + bl); 
        lmer_cpy(&(buf.ele[buf.size - bl]), key, l_len);
    }

    g_l_len = l_len;
    qsort(buf.ele, buf.size / bl, sizeof(uint64_t) * bl, cmp_lmer_uniq);
    sorted_key_fp = tmpfile_open();
    fwrite(buf.ele, sizeof(uint64_t), buf.size, sorted_key_fp);
    varr64_destroy(&buf);

    return sorted_key_fp;
}

FILE *sorted_lkey_from_contig_file(FILE *contig_fp, const uint64_t l_len)
{
    uint64_t i;
    uint64_t clen;
    uint64_t bl;
    uint64_t *key;
    lmer_t fwd;
    lmer_t rev;
    varr8_t cseq;
    varr64_t buf;
    FILE *sorted_key_fp;

    bl = (l_len + 31) / 32;

    varr64_init(&buf);
    rewind(contig_fp);
    varr8_init(&cseq);
    while (fread(&clen, sizeof(uint64_t), 1, contig_fp)) {
        fseek(contig_fp, sizeof(uint16_t), SEEK_CUR);
        varr8_resize(&cseq, (clen+3)/4);
        fread(cseq.ele, sizeof(uint8_t), (clen+3)/4, contig_fp);
        if (clen < l_len)
            continue;
        for (i = 0; i < l_len - 1; ++i) {
            lmer_set(fwd, i+1, bseq_get(cseq.ele, i));
            lmer_set(rev, l_len-i-2, 3^bseq_get(cseq.ele, i));
        }
        for (i = 0; i < clen - l_len + 1; ++i) {
            lmer_lshift(fwd, l_len);
            lmer_set(fwd, l_len-1, bseq_get(cseq.ele, i+l_len-1));
            lmer_rshift(rev, l_len);
            lmer_set(rev, 0, 3^bseq_get(cseq.ele, i+l_len-1));
            key = lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev;
            varr64_resize(&buf, buf.size + bl); 
            lmer_cpy(&(buf.ele[buf.size - bl]), key, l_len);
        }
    }
    varr8_destroy(&cseq);

    g_l_len = l_len;
    qsort(buf.ele, buf.size / bl, sizeof(uint64_t) * bl, cmp_lmer_uniq);
    sorted_key_fp = tmpfile_open();
    fwrite(buf.ele, sizeof(uint64_t), buf.size, sorted_key_fp);
    varr64_destroy(&buf);

    return sorted_key_fp;
}

void lcounter_count_quick(lcounter_t *lc, const uint8_t *seq, const uint64_t len, const uint64_t l_len)
{
    uint64_t i;
    uint64_t n_key;
    uint64_t bl;
    uint64_t st;
    uint64_t index;
    uint64_t step;
    lmer_t key;
    bool is_init;

    lc->l_len = l_len;
    bl = lc->bl = (l_len + 31) / 32;

    if (l_len >= len + 1) {
        memset(lc->val_table, 0, lc->table_size*sizeof(uint16_t));
        return;
    }

    st = n_key = 0;
    is_init = true;
    while (st < len - l_len + 1) {
        if (is_init) {
            for (i = 0; i < l_len - 1; ++i) {
                if (seq[st + i] == 4)
                    break;
            }
            if (i == l_len - 1)
                is_init = false;
            else {
                st += i + 1;
                continue;
            }
        }
        if (seq[st + l_len - 1] == 4) {
            st += l_len;
            is_init = true;
            continue;
        }

        ++st;
        ++n_key;
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

    st = 0;
    is_init = true;
    while (st < len - l_len + 1) {
        if (is_init) {
            for (i = 0; i < l_len - 1; ++i) {
                if (seq[st + i] == 4)
                    break;
                lmer_set(key, i + 1, seq[st + i]);
            }
            if (i == l_len - 1)
                is_init = false;
            else {
                st += i + 1;
                continue;
            }
        }
        if (seq[st + l_len - 1] == 4) {
            st += l_len;
            is_init = true;
            continue;
        }
        lmer_lshift(key, l_len);
        lmer_set(key, l_len - 1, seq[st + l_len - 1]);

        index = lmer_hash(key, l_len, lc->index_len) & (lc->table_size - 1);
        if (!(lc->val_table[index])) {
            lc->val_table[index] = 1;
            lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
        }
        else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
            if (lc->val_table[index] < UINT16_MAX-1)
                ++lc->val_table[index];
        }
        else {
            step = lmer_rehash(key, l_len, lc->index_len);
            while (1) {
                index = (index+step) & (lc->table_size-1);
                if (!(lc->val_table[index])) {
                    lc->val_table[index] = 1;
                    lmer_cpy(&(lc->key_table[index * bl]), key, l_len);;
                    break;
                }
                else if (lmer_cmp(&(lc->key_table[index * bl]), key, l_len) == 0) {
                    if (lc->val_table[index] < UINT16_MAX-1)
                        ++lc->val_table[index];
                    break;
                }
            }
        }

        ++st;
    }
}
