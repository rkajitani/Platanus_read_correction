#include "mcounter.h"

void mcounter_init(mcounter_t *mc, const uint64_t mem_lim)
{
    memset(mc, 0, sizeof(mcounter_t));
    mc->mem_lim = mem_lim;
}

void mcounter_destroy_table(mcounter_t *mc)
{
    if (mc->val_table != NULL)
        my_free(mc->val_table);
    mc->val_table = NULL;
    if (mc->key_table != NULL)
        my_free(mc->key_table);
    mc->key_table = NULL;
    mc->table_size = 0;
}

void mcounter_destroy(mcounter_t *mc)
{
    mcounter_destroy_table(mc);
    if (mc->occ_distr != NULL)
        my_free(mc->occ_distr);
    if(mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    if(mc->table_fp != NULL)
        fclose(mc->table_fp);
    memset(mc, 0, sizeof(mcounter_t));
}

void mcounter_save_kmer_read(mcounter_t *mc, const uint64_t k_len, FILE *read_fp)
{
    uint64_t i;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    uint64_t n_used;
    uint16_t occ;
    kmer_t kmer;
    seq_t seq;
    FILE *unstored_fp;
    FILE *tmp_fp;

    mc->k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    mc->k_len = k_len;

    if (mc->val_table != NULL)
        mcounter_destroy_table(mc);

    mc->table_size = 1ull;
    for (i = 0; i < k_len*2; ++i) {
        mc->table_size *= 2ull;
        if (mc->table_size*(sizeof(uint16_t)+sizeof(uint64_t)) > mc->mem_lim / 2)
            break;
    }
    mc->index_len = i;
    mc->table_size /= 2ull;

    mc->val_table = (uint16_t *)my_calloc(mc->table_size, sizeof(uint16_t));
    mc->key_table = (uint64_t *)my_malloc(mc->table_size * sizeof(uint64_t));

    if (mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    mc->kmer_fp = tmpfile_open();
    check_tmpfile(mc->kmer_fp);
    unstored_fp = tmpfile_open();
    check_tmpfile(unstored_fp);

    rewind(read_fp);
    n_used = 0;

    while (seq_read(&seq, read_fp)) {
        ++mc->len_distr[seq.len];
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
            kmer.fwd = ((kmer.fwd << 2) & mc->k_mask) | (uint64_t)seq.base[i+k_len-1];
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^seq.base[i+k_len-1]) << (2*(k_len-1)));
            if (seq.unknown_pos[seq.n_unknown] < i + k_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                continue;
            }

            key = min_u64(kmer.fwd, kmer.rev);
            index = hash64(key, mc->index_len) & (mc->table_size - 1);
            if (!(mc->val_table[index])) {
                if (n_used <= mc->table_size*MAX_LOAD_FACTOR) {
                    mc->val_table[index] = 1;
                    mc->key_table[index] = key;
                    ++n_used;
                }
                else
                    fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
            }
            else if (mc->key_table[index] == key) {
                if (mc->val_table[index] < UINT16_MAX-1)
                    ++mc->val_table[index];
            }
            else {
                step = rehash64(key, mc->index_len);
                while (1) {
                    index = (index+step) & (mc->table_size-1);
                    if (!(mc->val_table[index])) {
                        if (n_used <= mc->table_size*MAX_LOAD_FACTOR) {
                            mc->val_table[index] = 1;
                            mc->key_table[index] = key;
                            ++n_used;
                        }
                        else
                            fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
                        break;
                    }
                    else if (mc->key_table[index] == key) {
                        if (mc->val_table[index] < UINT16_MAX-1)
                            ++mc->val_table[index];
                        break;
                    }
                }
            }
        }
    }

    if (mc->occ_distr != NULL)
        my_free(mc->occ_distr);
    mc->occ_distr = (uint64_t *)my_calloc(UINT16_MAX, sizeof(uint64_t));

    while (1) {
        for (i = 0; i < mc->table_size; ++i) {
            if (mc->val_table[i]) {
                occ = mc->val_table[i];
                fwrite(&mc->key_table[i], sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&occ, sizeof(uint16_t), 1, mc->kmer_fp);
                ++mc->occ_distr[occ];
                mc->val_table[i] = 0;
            }
        }

        if (!(ftell(unstored_fp) / sizeof(uint64_t)))
            break;

        rewind(unstored_fp);
        n_used = 0;
        tmp_fp = tmpfile_open();
        check_tmpfile(tmp_fp);
        while (fread(&key, sizeof(uint64_t), 1, unstored_fp)) {
            index = hash64(key, mc->index_len) & (mc->table_size - 1);
            if (!(mc->val_table[index])) {
                if (n_used <= mc->table_size*MAX_LOAD_FACTOR) {
                    mc->val_table[index] = 1;
                    mc->key_table[index] = key;
                    ++n_used;
                }
                else
                    fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
            }
            else if (mc->key_table[index] == key) {
                if (mc->val_table[index] < UINT16_MAX-1)
                    ++mc->val_table[index];
            }
            else {
                step = rehash64(key, mc->index_len);
                while (1) {
                    index = (index+step) & (mc->table_size-1);
                    if (!(mc->val_table[index])) {
                        if (n_used <= mc->table_size*MAX_LOAD_FACTOR) {
                            mc->val_table[index] = 1;
                            mc->key_table[index] = key;
                            ++n_used;
                        }
                        else
                            fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
                        break;
                    }
                    else if (mc->key_table[index] == key) {
                        if (mc->val_table[index] < UINT16_MAX-1)
                            ++mc->val_table[index];
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
        if (mc->occ_distr[i] > 0)
            break;
    }
    mc->max_occ = i;
}

void mcounter_load_kmer(mcounter_t *mc, const uint64_t min_occ)
{
    uint64_t key;
    uint64_t n_key;
    uint64_t index;
    uint64_t step;
    uint16_t occ;
    FILE *new_kmer_fp;

    fprintf(stderr, "loading kmers...\n");

    n_key = 0;
    rewind(mc->kmer_fp);
    while(fread(&key, sizeof(uint64_t), 1, mc->kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, mc->kmer_fp);
        if (occ >= min_occ)
            ++n_key;
    }

    mc->table_size = 2;
    mc->index_len = 1;
    while ((double)mc->table_size * MAX_LOAD_FACTOR < n_key) {
        mc->table_size *= 2;
        ++mc->index_len;
    }
    while (2 * mc->table_size * (sizeof(uint16_t)+sizeof(uint64_t)) < mc->mem_lim / 2) {
        mc->table_size *= 2;
        ++mc->index_len;
    }

    check_mem_usage(mc->table_size*(sizeof(uint64_t)+sizeof(uint16_t)), mc->mem_lim);
/*
    fprintf(stderr, "BUCKET_SIZE=%lu, ", mc->table_size);
    fprintf(stderr, "MIN_OCCURRENCE=%lu, ", min_occ);
    fprintf(stderr, "NUM_OCCUPIED_SLOT=%lu, ", n_key);
    fprintf(stderr, "LOAD_FACTOR=%f, ", (double)(n_key) / (double)mc->table_size);
    fprintf(stderr, "MEM_USAGE=%luB\n", mc->table_size * (sizeof(uint16_t)+sizeof(uint64_t))); 
*/

    mc->key_table = my_realloc(mc->key_table, mc->table_size * sizeof(uint64_t));
    mc->val_table = my_realloc(mc->val_table, mc->table_size * sizeof(uint16_t));
    memset(mc->val_table, 0, mc->table_size * sizeof(uint16_t));
    new_kmer_fp = tmpfile_open();
    check_tmpfile(new_kmer_fp);

    rewind(mc->kmer_fp);
    while (fread(&key, sizeof(uint64_t), 1, mc->kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, mc->kmer_fp);
        if (occ < min_occ) {
            fwrite(&key, sizeof(uint64_t), 1, new_kmer_fp);
            fwrite(&occ, sizeof(uint16_t), 1, new_kmer_fp);
            continue;
        }

        index = hash64(key, mc->index_len) & (mc->table_size - 1);
        if (!(mc->val_table[index])) {
            mc->val_table[index] = occ;
            mc->key_table[index] = key;
        }
        else {
            step = rehash64(key, mc->index_len);
            index = (index+step) & (mc->table_size-1);
            while (mc->val_table[index])
                index = (index+step) & (mc->table_size-1);
            mc->val_table[index] = occ;
            mc->key_table[index] = key;
        }
    }

    fclose(mc->kmer_fp);
    mc->kmer_fp = new_kmer_fp;
}

void mcounter_load_contig(mcounter_t *mc, const uint64_t k_len, const double occ_ratio, FILE *contig_fp)
{
    uint64_t i;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    uint64_t total_key;
    uint64_t clen;
    uint16_t occ;
    kmer_t kmer;
    varr8_t cseq;

    fprintf(stderr, "K = %lu, loading kmers from contigs...\n", k_len);

    mc->k_len = k_len;
    mc->k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;

    if (mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    mc->kmer_fp = tmpfile_open();
    check_tmpfile(mc->kmer_fp);

    rewind(contig_fp);
    total_key = 0;
    while (fread(&clen, sizeof(uint64_t), 1, contig_fp)) {
        fseek(contig_fp, sizeof(uint16_t) + ((clen+3)/4)*sizeof(uint8_t), SEEK_CUR);
        if (clen >= k_len)
            total_key += clen - k_len + 1;
    }

    if (mc->val_table != NULL)
        mcounter_destroy_table(mc);

    mc->table_size = 1;
    for (i = 0; i < k_len*2; ++i) {
        mc->table_size *= 2ull;
        if (mc->table_size*(sizeof(uint16_t)+sizeof(uint64_t)) > mc->mem_lim / 2)
            break;
    }
    mc->index_len = i;
    mc->table_size /= 2ull;
    while ((double)mc->table_size * MAX_LOAD_FACTOR < total_key) {
        mc->table_size *= 2;
        ++mc->index_len;
    }

    mc->val_table = (uint16_t *)my_calloc(mc->table_size, sizeof(uint16_t));
    mc->key_table = (uint64_t *)my_malloc(mc->table_size * sizeof(uint64_t));

    rewind(contig_fp);
    varr8_init(&cseq);
    while (fread(&clen, sizeof(uint64_t), 1, contig_fp)) {
        fread(&occ, sizeof(uint16_t), 1, contig_fp);
        occ = occ * occ_ratio + 0.5;
        varr8_resize(&cseq, (clen+3)/4);
        fread(cseq.ele, sizeof(uint8_t), (clen+3)/4, contig_fp);
        if (clen < k_len)
            continue;
        kmer.fwd = kmer.rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer.fwd = kmer.fwd << 2 | (uint64_t)bseq_get(cseq.ele, i);
            kmer.rev = (kmer.rev | (uint64_t)(3^bseq_get(cseq.ele, k_len-i-2))) << 2;
        }
        for (i = 0; i < clen - k_len + 1; ++i) {
            kmer.fwd = ((kmer.fwd << 2) & mc->k_mask) | (uint64_t)bseq_get(cseq.ele, i+k_len-1);
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^bseq_get(cseq.ele, i+k_len-1)) << (2*(k_len-1)));
            key = min_u64(kmer.fwd, kmer.rev);
            index = hash64(key, mc->index_len) & (mc->table_size - 1);
            if (!mc->val_table[index]) {
                mc->val_table[index] = occ;
                mc->key_table[index] = key;
            }
            else if (mc->key_table[index] == key) {
                if (mc->val_table[index] < occ)
                    mc->val_table[index] = occ;
            }
            else {
                step = rehash64(key, mc->index_len);
                index = (index + step) & (mc->table_size-1);
                while (mc->val_table[index] && mc->key_table[index] != key)
                    index = (index+step) & (mc->table_size-1);
                if (mc->val_table[index] < occ)
                    mc->val_table[index] = occ;
                mc->key_table[index] = key;
            }
        }
    }
    varr8_destroy(&cseq);
}

void mcounter_save_additional_kmer_read(mcounter_t *mc, const uint64_t k_len, FILE *read_fp)
{
    uint64_t i;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    uint64_t n_used;
    kmer_t kmer;
    seq_t seq;
    FILE *unstored_fp;
    FILE *tmp_fp;

    fprintf(stderr, "K = %lu, saving additional kmers from reads...\n", k_len);

    if (mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    mc->kmer_fp = tmpfile_open();
    check_tmpfile(mc->kmer_fp);
    unstored_fp = tmpfile_open();
    check_tmpfile(unstored_fp);

    rewind(read_fp);
    while (seq_read(&seq, read_fp)) {
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
            kmer.fwd = ((kmer.fwd << 2) & mc->k_mask) | (uint64_t)seq.base[i+k_len-1];
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^seq.base[i+k_len-1]) << (2*(k_len-1)));
            if (seq.unknown_pos[seq.n_unknown] < i + k_len) {
                if (seq.unknown_pos[seq.n_unknown] <= i)
                    ++seq.n_unknown;
                continue;
            }

            key = min_u64(kmer.fwd, kmer.rev);
            index = hash64(key, mc->index_len) & (mc->table_size - 1);
            if (!(mc->val_table[index])) {
                fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
            }
            else if (mc->key_table[index] != key) {
                step = rehash64(key, mc->index_len);
                while (1) {
                    index = (index+step) & (mc->table_size-1);
                    if (!(mc->val_table[index])) {
                        fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
                        break;
                    }
                    else if (mc->key_table[index] == key) {
                        break;
                    }
                }
            }
        }
    }

    for (i = 0; i < mc->table_size; ++i) {
        if (mc->val_table[i]) {
            fwrite(&mc->key_table[i], sizeof(uint64_t), 1, mc->kmer_fp);
            fwrite(&(mc->val_table[i]), sizeof(uint16_t), 1, mc->kmer_fp);
            mc->val_table[i] = 0;
        }
    }

    while (ftell(unstored_fp) / sizeof(uint64_t)) {

        rewind(unstored_fp);
        n_used = 0;
        tmp_fp = tmpfile_open();
        check_tmpfile(tmp_fp);
        while (fread(&key, sizeof(uint64_t), 1, unstored_fp)) {
            index = hash64(key, mc->index_len) & (mc->table_size - 1);
            if (!(mc->val_table[index])) {
                if (n_used <= mc->table_size*MAX_LOAD_FACTOR) {
                    mc->val_table[index] = 1;
                    mc->key_table[index] = key;
                    ++n_used;
                }
                else
                    fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
            }
            else if (mc->key_table[index] == key) {
                if (mc->val_table[index] < UINT16_MAX-1)
                    ++mc->val_table[index];
            }
            else {
                step = rehash64(key, mc->index_len);
                while (1) {
                    index = (index+step) & (mc->table_size-1);
                    if (!(mc->val_table[index])) {
                        if (n_used <= mc->table_size*MAX_LOAD_FACTOR) {
                            mc->val_table[index] = 1;
                            mc->key_table[index] = key;
                            ++n_used;
                        }
                        else
                            fwrite(&key, sizeof(uint64_t), 1, tmp_fp);
                        break;
                    }
                    else if (mc->key_table[index] == key) {
                        if (mc->val_table[index] < UINT16_MAX-1)
                            ++mc->val_table[index];
                        break;
                    }
                }
            }
        }
        for (i = 0; i < mc->table_size; ++i) {
            if (mc->val_table[i]) {
                fwrite(&mc->key_table[i], sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&(mc->val_table[i]), sizeof(uint16_t), 1, mc->kmer_fp);
                mc->val_table[i] = 0;
            }
        }
        fclose(unstored_fp);
        unstored_fp = tmp_fp;
    }

    fclose(unstored_fp);
}

void mcounter_save_kmer_read_mt(mcounter_t *mc, const uint64_t k_len, FILE **read_fp, const uint64_t n_thread)
{
    int64_t i;
    int64_t j;
    uint16_t occ;
    FILE *unstored_fp[MAX_THREAD];
    FILE *tmp_fp[MAX_THREAD];
    mcounter_t mcs[MAX_THREAD];
    pthread_t thread[MAX_THREAD];
    pthread_mutex_t *mutex;
    arg_mcounter_lock_and_count_t arg[MAX_THREAD];

    omp_set_num_threads(n_thread);

#    pragma omp parallel for schedule(static, 1)
    for (i = 0; i < n_thread; ++i) {
        mcounter_init(&(mcs[i]), mc->mem_lim / n_thread);
        mcounter_save_kmer_read(&(mcs[i]), k_len, read_fp[i]);
        mcounter_destroy_table(&(mcs[i]));
    }

    mc->k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    mc->k_len = k_len;

    if (mc->val_table != NULL)
        mcounter_destroy_table(mc);

    mc->table_size = 1ull;
    for (i = 0; i < k_len*2; ++i) {
        mc->table_size *= 2ull;
        if (mc->table_size*(sizeof(uint16_t)+sizeof(uint64_t)) + (mc->table_size/GRANULE+1) * sizeof(pthread_mutex_t) > mc->mem_lim / 2)
            break;
    }
    mc->index_len = i;
    mc->table_size /= 2ull;

    mc->val_table = (uint16_t *)my_calloc(mc->table_size, sizeof(uint16_t));
    mc->key_table = (uint64_t *)my_malloc(mc->table_size * sizeof(uint64_t));

    if (mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    mc->kmer_fp = tmpfile_open();
    check_tmpfile(mc->kmer_fp);

    mutex = (pthread_mutex_t *)my_malloc((mc->table_size / GRANULE + 1) * sizeof(pthread_mutex_t));
    for (i = 0; i < mc->table_size/GRANULE+1; ++i)
        pthread_mutex_init(&mutex[i], NULL);

    for (i = 0; i < n_thread; ++i) {
        arg[i].mc = mc;
        arg[i].mutex = mutex;
        arg[i].in_kmer_fp = mcs[i].kmer_fp;
        unstored_fp[i] = tmpfile_open();
        check_tmpfile(unstored_fp[i]);
        arg[i].unstored_fp = unstored_fp[i];
        pthread_create(&thread[i], NULL, (void *)mcounter_lock_and_count1, &arg[i]);
    }
    for (i = 0; i < n_thread; ++i) {
        pthread_join(thread[i], NULL); 
        for (j = 1; j <= MAX_READ_LEN; ++j)
            mc->len_distr[j] += mcs[i].len_distr[j];
        mcounter_destroy(&(mcs[i]));
    }

    if (mc->occ_distr != NULL)
        my_free(mc->occ_distr);
    mc->occ_distr = (uint64_t *)my_calloc(UINT16_MAX, sizeof(uint64_t));

    while (1) {
        for (i = 0; i < mc->table_size; ++i) {
            if (mc->val_table[i]) {
                occ = mc->val_table[i];
                fwrite(&mc->key_table[i], sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&occ, sizeof(uint16_t), 1, mc->kmer_fp);
                ++mc->occ_distr[occ];
                mc->val_table[i] = 0;
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
            pthread_create(&thread[i], NULL, (void *)mcounter_lock_and_count1, &arg[i]);
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
        if (mc->occ_distr[i] > 0)
            break;
    }
    mc->max_occ = i;
}

void mcounter_save_additional_kmer_read_mt(mcounter_t *mc, const uint64_t k_len, FILE **read_fp, const uint64_t n_thread)
{
    int64_t i;
    uint64_t j;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    kmer_t kmer;
    seq_t seq;
    FILE *unstored_fp[MAX_THREAD];
    FILE *tmp_fp[MAX_THREAD];
    pthread_t thread[MAX_THREAD];
    pthread_mutex_t *mutex;
    arg_mcounter_lock_and_count_t arg[MAX_THREAD];

    fprintf(stderr, "K = %lu, saving additional kmers from reads...\n", k_len);

    if (mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    mc->kmer_fp = tmpfile_open();
    check_tmpfile(mc->kmer_fp);

    omp_set_num_threads(n_thread);

#    pragma omp parallel for private(j, index, step, key, kmer, seq) schedule(static, 1)
    for (i = 0; i < n_thread; ++i) {
        unstored_fp[i] = tmpfile_open();
        check_tmpfile(unstored_fp[i]);
        rewind(read_fp[i]);
        while (seq_read(&seq, read_fp[i])) {
            if (seq.len < k_len)
                continue;
            seq.unknown_pos[seq.n_unknown] = MAX_READ_LEN+1;
            seq.n_unknown = 0;
            kmer.fwd = kmer.rev = 0ull;
            for (j = 0; j < k_len - 1; ++j) {
                kmer.fwd = kmer.fwd << 2 | (uint64_t)seq.base[j];
                kmer.rev = (kmer.rev | (uint64_t)(3^seq.base[k_len-j-2])) << 2;
            }
            for (j = 0; j < seq.len - k_len + 1; ++j) {
                kmer.fwd = ((kmer.fwd << 2) & mc->k_mask) | (uint64_t)seq.base[j+k_len-1];
                kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^seq.base[j+k_len-1]) << (2*(k_len-1)));
                if (seq.unknown_pos[seq.n_unknown] < j + k_len) {
                    if (seq.unknown_pos[seq.n_unknown] <= j)
                        ++seq.n_unknown;
                    continue;
                }

                key = min_u64(kmer.fwd, kmer.rev);
                index = hash64(key, mc->index_len) & (mc->table_size - 1);
                if (!(mc->val_table[index])) {
                    fwrite(&key, sizeof(uint64_t), 1, unstored_fp[i]);
                }
                else if (mc->key_table[index] != key) {
                    step = rehash64(key, mc->index_len);
                    while (1) {
                        index = (index+step) & (mc->table_size-1);
                        if (!(mc->val_table[index])) {
                            fwrite(&key, sizeof(uint64_t), 1, unstored_fp[i]);
                            break;
                        }
                        else if (mc->key_table[index] == key) {
                            break;
                        }
                    }
                }
            }
        }
    }

    mutex = (pthread_mutex_t *)my_malloc((mc->table_size / GRANULE + 1) * sizeof(pthread_mutex_t));
    for (i = 0; i < mc->table_size/GRANULE+1; ++i)
        pthread_mutex_init(&mutex[i], NULL);

    for (i = 0; i < mc->table_size; ++i) {
        if (mc->val_table[i]) {
            fwrite(&mc->key_table[i], sizeof(uint64_t), 1, mc->kmer_fp);
            fwrite(&(mc->val_table[i]), sizeof(uint16_t), 1, mc->kmer_fp);
            mc->val_table[i] = 0;
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
            arg[i].mc = mc;
            arg[i].mutex = mutex;
            arg[i].in_kmer_fp = unstored_fp[i];
            arg[i].unstored_fp = tmp_fp[i];
            pthread_create(&thread[i], NULL, (void *)mcounter_lock_and_count2, &arg[i]);
        }
        for (i = 0; i < n_thread; ++i) {
            pthread_join(thread[i], NULL); 
            fclose(unstored_fp[i]);
            unstored_fp[i] = tmp_fp[i];
        }

        for (i = 0; i < mc->table_size; ++i) {
            if (mc->val_table[i]) {
                fwrite(&mc->key_table[i], sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&(mc->val_table[i]), sizeof(uint16_t), 1, mc->kmer_fp);
                mc->val_table[i] = 0;
            }
        }
    }

    for (i = 0; i < n_thread; ++i)
        fclose(unstored_fp[i]);

    my_free(mutex);
}

void mcounter_lock_and_count1(const arg_mcounter_lock_and_count_t *arg)
{
    uint64_t key;
    uint64_t index;
    uint16_t occ;
    FILE *in_kmer_fp;
    FILE *unstored_fp;
    mcounter_t *mc;
    pthread_mutex_t *mutex;

    mc = arg->mc;
    in_kmer_fp = arg->in_kmer_fp;
    unstored_fp = arg->unstored_fp;
    mutex = arg->mutex;

    rewind(in_kmer_fp);
    while (fread(&key, sizeof(uint64_t), 1, in_kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, in_kmer_fp);
        index = hash64(key, mc->index_len) & (mc->table_size - 1);

        pthread_mutex_lock(&mutex[index/GRANULE]);
        if (!(mc->val_table[index])) {
            mc->val_table[index] = occ;
            mc->key_table[index] = key;
        }
        else if (mc->key_table[index] == key) {
            if ((uint64_t)occ + mc->val_table[index] < UINT16_MAX - 1)
                mc->val_table[index] += occ;
            else
                mc->val_table[index] = UINT16_MAX - 1;
        }
        else {
            fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
            fwrite(&occ, sizeof(uint16_t), 1, unstored_fp);
        }
        pthread_mutex_unlock(&mutex[index/GRANULE]);
    }
}

void mcounter_lock_and_count2(const arg_mcounter_lock_and_count_t *arg)
{
    uint64_t key;
    uint64_t index;
    FILE *in_kmer_fp;
    FILE *unstored_fp;
    mcounter_t *mc;
    pthread_mutex_t *mutex;

    mc = arg->mc;
    in_kmer_fp = arg->in_kmer_fp;
    unstored_fp = arg->unstored_fp;
    mutex = arg->mutex;

    rewind(in_kmer_fp);
    while (fread(&key, sizeof(uint64_t), 1, in_kmer_fp)) {
        index = hash64(key, mc->index_len) & (mc->table_size - 1);

        pthread_mutex_lock(&mutex[index/GRANULE]);
        if (!(mc->val_table[index])) {
            mc->val_table[index] = 1;
            mc->key_table[index] = key;
        }
        else if (mc->key_table[index] == key) {
            if (mc->val_table[index] < UINT16_MAX - 1)
                ++mc->val_table[index];
        }
        else {
            fwrite(&key, sizeof(uint64_t), 1, unstored_fp);
        }
        pthread_mutex_unlock(&mutex[index/GRANULE]);
    }
}

uint64_t mcounter_sorted_kmer_list(const mcounter_t *mc, uint64_t ** buf)
{
    uint64_t i;
    uint64_t buf_size;

    buf_size = 0;
    for (i = 0; i < mc->table_size; ++i) {
        if (!mc->val_table[i])
            continue;
        ++buf_size;
    }
    *buf = (uint64_t *)my_malloc(buf_size * sizeof(uint64_t));
    buf_size = 0;
    for (i = 0; i < mc->table_size; ++i) {
        if (!mc->val_table[i])
            continue;
        (*buf)[buf_size] = mc->key_table[i];
        ++buf_size;
    }
    qsort(*buf, buf_size, sizeof(uint64_t), cmp_kmer_uniq);

    return buf_size;
}

void mcounter_extract_read(const mcounter_t *mc, FILE **read_fpp)
{
    uint64_t i;
    uint64_t k_len;
    uint64_t k_mask;
    uint64_t key;
    kmer_t kmer;
    seq_t seq;
    seq_t tmp_seq;
    FILE *new_read_fp;

    k_len = mc->k_len;
    k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    
    new_read_fp = tmpfile_open();
    check_tmpfile(new_read_fp);
    rewind(*read_fpp);

    while (seq_read(&seq, *read_fpp)) {
        if (seq.len < k_len)
            continue;
        tmp_seq = seq;
        tmp_seq.unknown_pos[tmp_seq.n_unknown] = MAX_READ_LEN+1;
        tmp_seq.n_unknown = 0;
        kmer.fwd = kmer.rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer.fwd = kmer.fwd << 2 | (uint64_t)tmp_seq.base[i];
            kmer.rev = (kmer.rev | (uint64_t)(3^tmp_seq.base[k_len-i-2])) << 2;
        }
        for (i = 0; i < tmp_seq.len - k_len + 1; ++i) {
            kmer.fwd = ((kmer.fwd << 2) & k_mask) | (uint64_t)tmp_seq.base[i+k_len-1];
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^tmp_seq.base[i+k_len-1]) << (2*(k_len-1)));
            if (tmp_seq.unknown_pos[tmp_seq.n_unknown] < i + k_len) {
                if (tmp_seq.unknown_pos[tmp_seq.n_unknown] <= i)
                    ++tmp_seq.n_unknown;
                continue;
            }

            key = min_u64(kmer.fwd, kmer.rev);
            if (mcounter_occ(mc, key) > 0) {
                seq_write(&seq, new_read_fp);
                break;
            }
        }
    }

    fclose(*read_fpp);
    *read_fpp = new_read_fp;
}

double mcounter_estimate_error_rate(mcounter_t *mc, const uint64_t k_len, FILE *read_fp)
{
    double err_rate; 
    double ave_len;
    double cov_cut;
    mcounter_t tmp_mc;
    mcounter_init(&tmp_mc, mc->mem_lim);

    fprintf(stderr, "K = %lu, saving kmers from reads...\n", k_len);
    mcounter_save_kmer_read(mc, k_len, read_fp);
    
    fprintf(stderr, "K = %lu, saving kmers from reads...\n", k_len - 1);
    mcounter_save_kmer_read(&tmp_mc, k_len - 1, read_fp);
    cov_cut = (double)get_left_minimal(tmp_mc.occ_distr, tmp_mc.max_occ);
    mcounter_load_kmer(&tmp_mc, (uint64_t)(cov_cut + 0.5));

    ave_len = get_ave(mc->len_distr, 0, MAX_READ_LEN);
    err_rate = mcounter_calculate_error_rate(mc, &tmp_mc, ave_len);
    mcounter_destroy_table(mc);
    mcounter_destroy(&tmp_mc);
    return err_rate;
}

double mcounter_estimate_error_rate_mt(mcounter_t *mc, const uint64_t k_len, FILE **read_fp, const uint64_t n_thread)
{
    double err_rate; 
    double ave_len;
    double cov_cut;
    mcounter_t tmp_mc;

    mcounter_init(&tmp_mc, mc->mem_lim);

    fprintf(stderr, "K = %lu, saving kmers from reads...\n", k_len);
    mcounter_save_kmer_read_mt(mc, k_len, read_fp, n_thread);
    
    fprintf(stderr, "K = %lu, saving kmers from reads...\n", k_len - 1);
    mcounter_save_kmer_read_mt(&tmp_mc, k_len - 1, read_fp, n_thread);
    cov_cut = (double)get_left_minimal(tmp_mc.occ_distr, tmp_mc.max_occ);
    mcounter_load_kmer(&tmp_mc, (uint64_t)(cov_cut + 0.5));

    ave_len = get_ave(mc->len_distr, 0, MAX_READ_LEN);
    err_rate = mcounter_calculate_error_rate(mc, &tmp_mc, ave_len);
    mcounter_destroy_table(mc);
    mcounter_destroy(&tmp_mc);
    return err_rate;
}

double mcounter_calculate_error_rate(const mcounter_t *mc1, const mcounter_t *mc2, const double ave_len)
{
    uint64_t i;
    uint64_t k_len;
    uint64_t mask;
    uint64_t key;
    uint64_t cov_cut;
    uint64_t base = 0;
    uint64_t tmp;
    uint64_t occ2;
    uint16_t occ1;
    double sum1;
    double sum2;
    kmer_t lkmer;
    kmer_t rkmer;

    k_len = mc1->k_len;
    mask = ~(~0ull << (2*(k_len-1)));
    cov_cut = get_left_minimal(mc1->occ_distr, mc1->max_occ);

    sum1 = sum2 = 0.0;
    rewind(mc1->kmer_fp);
    while (fread(&key, sizeof(uint64_t), 1, mc1->kmer_fp)) {
        fread(&occ1, sizeof(uint16_t), 1, mc1->kmer_fp);
        if (occ1 < cov_cut || occ1 == UINT16_MAX - 1)
            continue;

        lkmer.fwd = key >> 2;
        lkmer.rev = 0;
        for (i = 0; i < k_len - 1; ++i) {
            base = ((key >> (2 * (i + 1))) & 0x3);
            lkmer.rev |= (3^base) << (2 * (k_len - i - 2));
        }

        occ2 = mcounter_occ(mc2, min_u64(lkmer.fwd, lkmer.rev));
        if (occ2 == 0 || occ2 == UINT16_MAX - 1)
            continue;
        for (i = (base + 1) % 4; i != base; i = (i + 1) % 4) {
            lkmer.fwd = (lkmer.fwd & (mask >> 2)) | (i << (2 * (k_len - 2)));
            lkmer.rev = (lkmer.rev & ~0x3) | (3^i);
            if (mcounter_occ(mc2, min_u64(lkmer.fwd, lkmer.rev)) > 0)
                break;
        }
        if (i != base)
            continue;

        base = (key & 0x3);
        rkmer.fwd = (key & mask);
        rkmer.rev = (lkmer.rev >> 2) | (3^base << (2 * (k_len - 2)));
        tmp = mcounter_occ(mc2, min_u64(rkmer.fwd, rkmer.rev));
        if (tmp == 0 || tmp == UINT16_MAX - 1)
            continue;
        occ2 += tmp;
        for (i = (base + 1) % 4; i != base; i = (i + 1) % 4) {
            rkmer.fwd = (rkmer.fwd & ~0x3) | i;
            rkmer.rev = (rkmer.rev & (mask >> 2)) | ((3^i) << (2 * (k_len - 2)));
            if (mcounter_occ(mc2, min_u64(rkmer.fwd, rkmer.rev)) > 0)
                break;
        }
        if (i != base)
            continue;

        sum1 += (double)occ1;
        sum2 += (double)occ2 / 2.0;
    }
    
    return (1.0 - (sum1 / sum2 * (ave_len - mc2->k_len + 1) / (ave_len - mc1->k_len + 1)));
}

uint64_t mcounter_get_coverage_cutoff(const mcounter_t *mc, const double err_rate)
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
        n_read += mc->len_distr[i];
        total_len += i * mc->len_distr[i];
    }
    ave_len = (uint64_t)((double)total_len / n_read + 0.5);

    total_kmer = 0;
    for (i = 1; i <= mc->max_occ; ++i)
        total_kmer += i * mc->occ_distr[i];
    
    e_err_kmer = 0.0;
    for (i = 0; i < ave_len; ++i) {
        if (i >= mc->k_len - 1 && i + mc->k_len <= ave_len)
            e_err_kmer += (double)(mc->k_len);
        else
            e_err_kmer += (double)min_u64(i + 1, ave_len - i);
    }
    e_err_kmer /= (double)ave_len;
    e_err_kmer *= (double)total_len;
    e_err_kmer *= err_rate;

fprintf(stderr, "EXPECTED_ERRONEOUS_KMER = %f\n", e_err_kmer);

    cum = 0;
    for (i = 1; i <= mc->max_occ; ++i) {
        if (e_err_kmer - cum < 0.5 * (i * mc->occ_distr[i]))
            break;
        cum += i * mc->occ_distr[i];
    }

    return i;
}

FILE *sorted_key_from_kmer_file(FILE *kmer_fp, const uint64_t min_occ)
{
    uint16_t occ;
    uint64_t key;
    varr64_t buf;
    FILE *sorted_key_fp;

    varr64_init(&buf);
    rewind(kmer_fp);
    while(fread(&key, sizeof(uint64_t), 1, kmer_fp)) {
        fread(&occ, sizeof(uint16_t), 1, kmer_fp);
        if (occ < min_occ)
            continue;
        varr64_push(&buf, key);
    }

    qsort(buf.ele, buf.size, sizeof(uint64_t), cmp_kmer_uniq);
    sorted_key_fp = tmpfile_open();
    fwrite(buf.ele, sizeof(uint64_t), buf.size, sorted_key_fp);
    varr64_destroy(&buf);

    return sorted_key_fp;
}

FILE *sorted_key_from_contig_file(FILE *contig_fp, const uint64_t k_len)
{
    uint64_t i;
    uint64_t clen;
    uint64_t mask;
    kmer_t kmer;
    varr8_t cseq;
    varr64_t buf;
    FILE *sorted_key_fp;

    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;

    varr64_init(&buf);
    varr8_init(&cseq);
    rewind(contig_fp);
    while (fread(&clen, sizeof(uint64_t), 1, contig_fp)) {
        fseek(contig_fp, sizeof(uint16_t), SEEK_CUR);
        varr8_resize(&cseq, (clen+3)/4);
        fread(cseq.ele, sizeof(uint8_t), (clen+3)/4, contig_fp);
        if (clen < k_len)
            continue;
        kmer.fwd = kmer.rev = 0ull;
        for (i = 0; i < k_len - 1; ++i) {
            kmer.fwd = kmer.fwd << 2 | (uint64_t)bseq_get(cseq.ele, i);
            kmer.rev = (kmer.rev | (uint64_t)(3^bseq_get(cseq.ele, k_len-i-2))) << 2;
        }
        for (i = 0; i < clen - k_len + 1; ++i) {
            kmer.fwd = ((kmer.fwd << 2) & mask) | (uint64_t)bseq_get(cseq.ele, i+k_len-1);
            kmer.rev = (kmer.rev >> 2) | ((uint64_t)(3^bseq_get(cseq.ele, i+k_len-1)) << (2*(k_len-1)));
            varr64_push(&buf, min_u64(kmer.fwd, kmer.rev));
        }
    }
    varr8_destroy(&cseq);

    qsort(buf.ele, buf.size, sizeof(uint64_t), cmp_kmer_uniq);
    sorted_key_fp = tmpfile_open();
    fwrite(buf.ele, sizeof(uint64_t), buf.size, sorted_key_fp);
    varr64_destroy(&buf);

    return sorted_key_fp;
}

void mcounter_count_quick(mcounter_t *mc, const uint8_t *seq, const uint64_t len, const uint64_t k_len)
{
    uint64_t i;
    uint64_t n_key;
    uint64_t mask;
    uint64_t st;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    bool is_init;

    mc->k_len = k_len;
    mask = mc->k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;

    if (k_len >= len + 1) {
        memset(mc->val_table, 0, mc->table_size * sizeof(uint16_t));
        return;
    }

    st = n_key = 0;
    is_init = true;
    while (st < len - k_len + 1) {
        if (is_init) {
            for (i = 0; i < k_len - 1; ++i) {
                if (seq[st + i] == 4)
                    break;
            }
            if (i == k_len - 1)
                is_init = false;
            else {
                st += i + 1;
                continue;
            }
        }
        if (seq[st + k_len - 1] == 4) {
            st += k_len;
            is_init = true;
            continue;
        }

        ++st;
        ++n_key;
    }

    if ((double)mc->table_size * MAX_LOAD_FACTOR <= n_key) {
        mc->table_size = 2;
        mc->index_len = 1;
        while ((double)mc->table_size * MAX_LOAD_FACTOR <= n_key) {
            mc->table_size *= 2;
            ++mc->index_len;
        }
        mc->key_table = my_realloc(mc->key_table, mc->table_size * sizeof(uint64_t));
        mc->val_table = my_realloc(mc->val_table, mc->table_size * sizeof(uint16_t));
    }

    memset(mc->val_table, 0, mc->table_size * sizeof(uint16_t));

    st = key = 0;
    is_init = true;
    while (st < len - k_len + 1) {
        if (is_init) {
            key = 0;
            for (i = 0; i < k_len - 1; ++i) {
                if (seq[st + i] == 4)
                    break;
                key = (key << 2) | (uint64_t)seq[st + i];
            }
            if (i == k_len - 1)
                is_init = false;
            else {
                st += i + 1;
                continue;
            }
        }
        if (seq[st + k_len - 1] == 4) {
            st += k_len;
            is_init = true;
            continue;
        }
        key = ((key << 2) & mask) | (uint64_t)seq[st + k_len - 1];

        index = hash64(key, mc->index_len) & (mc->table_size - 1);
        if (!(mc->val_table[index])) {
            mc->val_table[index] = 1;
            mc->key_table[index] = key;
        }
        else if (mc->key_table[index] == key) {
            if (mc->val_table[index] < UINT16_MAX-1)
                ++mc->val_table[index];
        }
        else {
            step = rehash64(key, mc->index_len);
            index = (index + step) & (mc->table_size-1);
            while (mc->val_table[index] && mc->key_table[index] != key)
                index = (index+step) & (mc->table_size-1);
            if (mc->val_table[index] < UINT16_MAX-1)
                ++mc->val_table[index];
            mc->key_table[index] = key;
        }

        ++st;
    }
}

void mcounter_save_table(mcounter_t *mc)
{
	if (mc->table_fp != NULL)
		fclose(mc->table_fp);
	mc->table_fp = tmpfile_open();
	check_tmpfile(mc->table_fp);

	fwrite(&(mc->k_len), sizeof(uint64_t), 1, mc->table_fp);
	fwrite(&(mc->k_mask), sizeof(uint64_t), 1, mc->table_fp);
	fwrite(&(mc->index_len), sizeof(uint64_t), 1, mc->table_fp);
	fwrite(&(mc->table_size), sizeof(uint64_t), 1, mc->table_fp);
	fwrite(mc->key_table, sizeof(uint64_t), mc->table_size, mc->table_fp);
	fwrite(mc->val_table, sizeof(uint16_t), mc->table_size, mc->table_fp);
}
