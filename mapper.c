#include "mapper.h"

void mapper_init(mapper_t *mp, const int64_t seed_len, const int64_t key_len, const int64_t mem_lim)
{
    memset(mp, 0, sizeof(mapper_t));
    mp->seed_len = seed_len;
    mp->key_len = key_len;
    if (mp->key_len > 32)
        mp->key_len = 32;
    mp->mask = mp->key_len < 32 ? ~(~0ull << (2*mp->key_len)) : ~0ull;
    mp->mem_lim = mem_lim;
}

void mapper_destroy_table(mapper_t *mp)
{
    if (mp->table != NULL) {
        my_free(mp->table);
        my_free(mp->pos_pool);
        mp->mem_usage -= mp->table_size * sizeof(map_rec_t) + mp->pos_pool_size * sizeof(position_t);
    }
    mp->table_size = mp->pos_pool_size = 0;;
    mp->table = NULL;
    mp->pos_pool = NULL;
}

void mapper_destroy(mapper_t *mp)
{
    mapper_destroy_table(mp);
    memset(mp, 0, sizeof(mapper_t));
}

void mapper_set_contig(mapper_t *mp, const contig_t *con)
{
    mp->mem_usage += con->mem_usage;
    mp->n_seq = con->n_seq;
    mp->seq_pool_size = con->seq_pool_size;
    mp->seq = con->seq;
}

void mapper_make_table(mapper_t *mp)
{
    int64_t i;
    int64_t j;
    int64_t len;
    int64_t st;
    int64_t n_key;
    uint64_t key;
    uint64_t index;
    uint64_t step;
    int8_t *seq;
    bool is_init;
    FILE *key_occ_fp;
    FILE *key_pos_fp;
    FILE *uniq_fp;
    position_t pos;

    fprintf(stderr, "K=%ld, making hash table...\n", mp->key_len);

    mp->table_size = 1;
    mp->index_len = 0;
    for (i = 0; i < 2 * mp->key_len; ++i) {
        mp->table_size *= 2;
        ++mp->index_len;
        if (mp->table_size * MAX_LOAD_FACTOR > mp->seq_pool_size - mp->key_len + 1)
            break;
    }

    mp->table = my_calloc(mp->table_size, sizeof(map_rec_t));
    key_pos_fp = tmpfile_open();
    check_tmpfile(key_pos_fp);

    for (i = 0; i < mp->n_seq; ++i) {
        len = mp->seq[i].len;
        if (len < mp->key_len)
            continue;
        seq = mp->seq[i].seq;
        st = 0;
        is_init = true;
        while (st < len - mp->key_len + 1) {
            if (is_init) {
                key = 0;
                for (j = 0; j < mp->key_len - 1; ++j) {
                    if (seq[st + j] == 4)
                        break;
                    key = (key << 2) | (uint64_t)seq[st + j];
                }
                if (j == mp->key_len - 1)
                    is_init = false;
                else {
                    st += j + 1;
                    continue;
                }
            }
            if (seq[st + mp->key_len - 1] == 4) {
                st += mp->key_len;
                is_init = true;
                continue;
            }
            key = ((key << 2) & mp->mask) | (uint64_t)seq[st + mp->key_len - 1];

            pos.id = i + 1;
            pos.ofst = st;
            fwrite(&key, sizeof(uint64_t), 1, key_pos_fp);
            fwrite(&pos, sizeof(position_t), 1, key_pos_fp);

            index = hash64(key, mp->index_len) & (mp->table_size-1);
            if (!(map_rec_is_emp(mp->table[index])) && mp->table[index].key != key) {
                step = rehash64(key, mp->index_len);
                index = (index+step) & (mp->table_size-1);
                while (!(map_rec_is_emp(mp->table[index])) && mp->table[index].key != key) {
                    index = (index+step) & (mp->table_size-1);
                }
            }
            mp->table[index].key = key;
            ++(mp->table[index].val);

            ++st;
        }
    }

    key_occ_fp = tmpfile_open();
    check_tmpfile(key_occ_fp);
    mp->pos_pool_size = mp->max_occ = n_key = 0;
    for (i = 0; i < mp->table_size; ++i) {
        if (map_rec_is_emp(mp->table[i]))
            continue;
        if (mp->table[i].val > mp->max_occ)
            mp->max_occ = mp->table[i].val;
        if (mp->table[i].val > 1) {
            mp->pos_pool_size += mp->table[i].val; 
            fwrite(&(mp->table[i].key), sizeof(uint64_t), 1, key_occ_fp);
            fwrite(&(mp->table[i].val), sizeof(int64_t), 1, key_occ_fp);
        }
        ++n_key;
    }

    mp->table_size = 1;
    mp->index_len = 0;
    for (i = 0; i < 2 * mp->key_len; ++i) {
        mp->table_size *= 2;
        ++mp->index_len;
        if (mp->table_size * MAX_LOAD_FACTOR > n_key)
            break;
    }

    mp->mem_usage += mp->table_size*sizeof(map_rec_t) + mp->pos_pool_size*sizeof(position_t);
    check_mem_usage(mp->mem_usage, mp->mem_lim);
    fprintf(stderr, "MEM_USAGE=%ld\n", mp->mem_usage);

    mp->table = (map_rec_t *)my_realloc(mp->table, mp->table_size * sizeof(map_rec_t));
    memset(mp->table, 0, mp->table_size * sizeof(map_rec_t)); 
    mp->pos_pool = (position_t *)my_calloc(mp->pos_pool_size + 1, sizeof(position_t));

    rewind(key_occ_fp);
    while (fread(&key, sizeof(uint64_t), 1, key_occ_fp)) {
        index = hash64(key, mp->index_len) & (mp->table_size-1);
        if (!(map_rec_is_emp(mp->table[index]))) {
            step = rehash64(key, mp->index_len);
            index = (index+step) & (mp->table_size-1);
            while (!(map_rec_is_emp(mp->table[index])))
                index = (index+step) & (mp->table_size-1);
        }
        mp->table[index].key = key;
        fread(&(mp->table[index].val), sizeof(int64_t), 1, key_occ_fp);
    }
    fclose(key_occ_fp);

    mp->pos_pool_size = 1;
    for (i = 0; i < mp->table_size; ++i) {
        if (map_rec_is_emp(mp->table[i]))
            continue;
        mp->pos_pool_size += mp->table[i].val;
        mp->table[i].val = -(mp->pos_pool_size);
    }

    uniq_fp = tmpfile_open();
    check_tmpfile(uniq_fp);
    rewind(key_pos_fp);
    while (fread(&key, sizeof(uint64_t), 1, key_pos_fp)) {
        fread(&pos, sizeof(position_t), 1, key_pos_fp);
        index = hash64(key, mp->index_len) & (mp->table_size-1);
        if (!(map_rec_is_emp(mp->table[index])) && mp->table[index].key != key) {
            step = rehash64(key, mp->index_len);
            index = (index+step) & (mp->table_size-1);
            while (!(map_rec_is_emp(mp->table[index])) && mp->table[index].key != key)
                index = (index+step) & (mp->table_size-1);
        }

        if (map_rec_is_emp(mp->table[index])) {
            fwrite(&key, sizeof(uint64_t), 1, uniq_fp);
            fwrite(&pos, sizeof(position_t), 1, uniq_fp);
            continue;
        }

        if (mp->table[index].val < 0) {
            mp->table[index].val *= -1; 
            pos.id *= -1;
        }

        --(mp->table[index].val);
        mp->pos_pool[mp->table[index].val] = pos;
    }
    fclose(key_pos_fp);

    rewind(uniq_fp);
    while (fread(&key, sizeof(uint64_t), 1, uniq_fp)) {
        index = hash64(key, mp->index_len) & (mp->table_size-1);
        if (!(map_rec_is_emp(mp->table[index])) && mp->table[index].key != key) {
            step = rehash64(key, mp->index_len);
            index = (index+step) & (mp->table_size-1);
            while (!(map_rec_is_emp(mp->table[index])) && mp->table[index].key != key)
                index = (index+step) & (mp->table_size-1);
        }
        mp->table[index].key = key;
        fread(&(mp->table[index].pos), sizeof(position_t), 1, uniq_fp);
        mp->table[index].pos.id *= -1;
    }
    fclose(uniq_fp);
}

int64_t mapper_map_seed(const mapper_t *mp, const kmer_t *kmer, position_t *buf)
{
    uint64_t index;
    uint64_t step;
    int64_t n_mapped;

    n_mapped = 0;

    index = hash64(kmer->fwd, mp->index_len) & (mp->table_size-1);
    if (!map_rec_is_emp(mp->table[index])) {
        if (kmer->fwd != mp->table[index].key) {
            step = rehash64(kmer->fwd, mp->index_len);
            index = (index + step) & (mp->table_size-1);
            while (!map_rec_is_emp(mp->table[index])) {
                if (kmer->fwd == mp->table[index].key) {
                    if (mp->table[index].pos.id < 0)
                        buf[n_mapped] = mp->table[index].pos;
                    else {
                        for (index = mp->table[index].val; mp->pos_pool[index].id > 0; ++index) {
                            buf[n_mapped] = mp->pos_pool[index];
                            ++n_mapped;
                        }
                        buf[n_mapped] = mp->pos_pool[index];
                    }
                    buf[n_mapped].id *= -1;
                    ++n_mapped;
                    break;
                }
                index = (index+step) & (mp->table_size-1);
            }
        }
        else {
            if (mp->table[index].pos.id < 0)
                buf[n_mapped] = mp->table[index].pos;
            else {
                for (index = mp->table[index].val; mp->pos_pool[index].id > 0; ++index) {
                    buf[n_mapped] = mp->pos_pool[index];
                    ++n_mapped;
                }
                buf[n_mapped] = mp->pos_pool[index];
            }
            buf[n_mapped].id *= -1;
            ++n_mapped;
        }
    }

    index = hash64(kmer->rev, mp->index_len) & (mp->table_size-1);
    if (!map_rec_is_emp(mp->table[index])) {
        if (kmer->rev != mp->table[index].key) {
            step = rehash64(kmer->rev, mp->index_len);
            index = (index + step) & (mp->table_size-1);
            while (!map_rec_is_emp(mp->table[index])) {
                if (kmer->rev == mp->table[index].key) {
                    if (mp->table[index].pos.id < 0)
                        buf[n_mapped] = mp->table[index].pos;
                    else {
                        for (index = mp->table[index].val; mp->pos_pool[index].id > 0; ++index) {
                            buf[n_mapped] = mp->pos_pool[index];
                            buf[n_mapped].id *= -1;
                            ++n_mapped;
                        }
                        buf[n_mapped] = mp->pos_pool[index];
                    }
                    ++n_mapped;
                    break;
                }
                index = (index+step) & (mp->table_size-1);
            }
        }
        else {
            if (mp->table[index].pos.id < 0)
                buf[n_mapped] = mp->table[index].pos;
            else {
                for (index = mp->table[index].val; mp->pos_pool[index].id > 0; ++index) {
                    buf[n_mapped] = mp->pos_pool[index];
                    buf[n_mapped].id *= -1;
                    ++n_mapped;
                }
                buf[n_mapped] = mp->pos_pool[index];
            }
            ++n_mapped;
        }
    }

    return n_mapped;
}

position_t mapper_map_read(const mapper_t *mp, const seq_t *read, position_t *buf)
{
    int64_t i;
    int64_t j;
    int64_t k = 0;
    int64_t buf_size;
    int64_t max;
    int64_t sec_max;
    int64_t n_res;
    position_t res[MAX_READ_LEN];
    kmer_t kmer;

    n_res = 0;

    for (i = read->len - mp->seed_len; i > -(mp->seed_len); i -= mp->seed_len) {
        if (i < 0)
            i = 0;
        kmer.fwd = kmer.rev = 0;
        for (j = 0; j < mp->key_len; ++j) {
            if (read->base[i + j] == 4)
                break;
            kmer.fwd = (kmer.fwd << 2) | (uint64_t)read->base[i + j];
            kmer.rev = (kmer.rev << 2) | (uint64_t)(3^read->base[i + mp->key_len - j - 1]);
        }
        if (j != mp->key_len)
            continue;

        buf_size = mapper_map_seed(mp, &kmer, buf);

        res[n_res].id = 0;
        for (j = 0; j < buf_size; ++j) {
            if (buf[j].id > 0) {
                if (buf[j].ofst > mp->seq[buf[j].id - 1].len - mp->seed_len)
                    continue;
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k])
                        break;
                }
                buf[j].ofst -= i;
            }
            else {
                if (buf[j].ofst < mp->seed_len - mp->key_len)
                    continue;
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k]))
                        break;
                }
                buf[j].ofst += i + mp->key_len - 1;
            }

            if (k == mp->seed_len) {
                if (res[n_res].id == 0)
                    res[n_res] = buf[j];
                else
                    break;
            }
        }
        if (j == buf_size && res[n_res].id != 0)
            ++n_res;
    }
    if (n_res == 0) {
        res[0].id = 0;
        return res[0];
    }

    qsort(res, n_res, sizeof(position_t), cmp_position);

    i = max = sec_max = res[n_res].id = 0;
    while (res[i].id != 0) {
        j = i;
        while (res[i].id == res[i+1].id)
            ++i;
        ++i;
        if (i - j >= max) {
            sec_max = max;
            max = i - j;
            k = j;
        }
    }

    if (max == sec_max) {
        res[0].id = 0;
        return res[0];
    }

    j = 0;
    for (i = 0; i < max; ++i)
        j += res[k + i].ofst;
    res[k].ofst = (double)j / max + 0.5;
    return res[k];
}

void mapper_map_pair_mt(const mapper_t *mp, seqlib_t *lib, const int64_t min_ins, const int64_t n_thread)
{
    int c;
    int64_t i;
    int64_t sum;
    int64_t total_pair;
    int64_t n_map_same;
    int64_t n_map_diff;
    int64_t ins_len;
    position_t fpos;
    position_t rpos;
    seq_t fwd;
    seq_t rev;
    position_t *buf;

    fputs("mapping reads...\n", stderr);
    fprintf(stderr, "MIN_INS_LEN = %ld\n", min_ins);

    omp_set_num_threads(n_thread);

    sum = total_pair = n_map_same = n_map_diff = 0;

#    pragma omp parallel for  schedule(static, 1) \
        private(buf, fwd, rev, fpos, rpos, ins_len) \
        reduction(+: sum, total_pair, n_map_same, n_map_diff)
    for (i = 0; i < n_thread; ++i) {
        buf = (position_t *)my_malloc(mp->max_occ * 2 * sizeof(position_t));

        if (lib[i].ins_len_fp != NULL)
            fclose(lib[i].ins_len_fp);
        lib[i].ins_len_fp = tmpfile_open();
        check_tmpfile(lib[i].ins_len_fp);

        if (lib[i].mapped_fp != NULL)
            fclose(lib[i].mapped_fp);
        lib[i].mapped_fp = tmpfile_open();
        check_tmpfile(lib[i].mapped_fp);

        rewind(lib[i].pair_fp);
        while (seq_read(&fwd, lib[i].pair_fp)) {
            seq_read(&rev, lib[i].pair_fp);
            ++total_pair;

            if (fwd.len < mp->seed_len || rev.len < mp->seed_len)
                continue;

            fpos = mapper_map_read(mp, &fwd, buf);
            rpos = mapper_map_read(mp, &rev, buf);

            if (fpos.id == 0 || rpos.id == 0)
                continue;

            if (fpos.id == -rpos.id) {
                if (fpos.id > 0 && fpos.ofst < rpos.ofst)
                    ins_len = (int64_t)rpos.ofst - fpos.ofst;
                else if (rpos.id > 0 && rpos.ofst < fpos.ofst)
                    ins_len = (int64_t)fpos.ofst - rpos.ofst;
                else
                    continue;

                if (ins_len < min_ins || ins_len < fwd.len || ins_len < rev.len)
                    continue;

                fwrite(&ins_len, sizeof(int64_t), 1, lib[i].ins_len_fp);
                ++n_map_same;
            }
            else if (fpos.id != rpos.id) {
                fwrite(&fpos, sizeof(position_t), 1, lib[i].mapped_fp);
                fwrite(&rpos, sizeof(position_t), 1, lib[i].mapped_fp);
                ++n_map_diff;
            }

            sum += fwd.len + rev.len;
        }
        my_free(buf);
    }

    if (n_map_same == 0)
        fprintf(stderr, "WARNING: no read mapped in the same contig!\n");
    if (n_map_diff == 0)
        fprintf(stderr, "WARNING: no read links contigs!\n");

    lib->ave_cov = (double)sum / mp->seq_pool_size;

#    pragma omp parallel sections private(i, c)
    {
#        pragma omp section
        {
            for (i = 1; i < n_thread; ++i) {
                rewind(lib[i].ins_len_fp);
                while (fread(&c, sizeof(int8_t), 1, lib[i].ins_len_fp))
                    putc(c, lib->ins_len_fp);
                fclose(lib[i].ins_len_fp);
            }
        }

#        pragma omp section
        {
            for (i = 1; i < n_thread; ++i) {
                rewind(lib[i].mapped_fp);
                while(fread(&c, sizeof(int8_t), 1, lib[i].mapped_fp))
                    putc(c, lib->mapped_fp);
                fclose(lib[i].mapped_fp);
            }
        }
    }

    fprintf(stderr, "TOTAL_PAIR = %ld\n", total_pair);
    fprintf(stderr, "MAPPED_PAIR = %ld (%f)\n", n_map_same + n_map_diff, (double)(n_map_same + n_map_diff) / total_pair);
    fprintf(stderr, "MAPPED_IN_DIFFERENT_CONTIGS = %ld (%f)\n", n_map_diff, (double)n_map_diff / total_pair);
    fprintf(stderr, "MAPPED_IN_SAME_CONTIG = %ld (%f)\n", n_map_same, (double)n_map_same / total_pair);
    fprintf(stderr, "AVERAGE_COVERAGE = %f\n", lib->ave_cov);
}

void mapper_map_small_gap_mt(const mapper_t *mp, seqlib_t *lib, FILE **gap_seq_fp, const int64_t n_thread)
{
    int8_t c;
    int8_t gap_seq[MAX_READ_LEN];
    int32_t gap_seq_len;
    int64_t i;
    int64_t j;
    int64_t t;
    int64_t st;
    int64_t ed;
    int64_t lbuf_size;
    int64_t rbuf_size;
    position_t lres;
    position_t rres;
    position_t *lbuf;
    position_t *rbuf;
    kmer_t lkmer;
    kmer_t rkmer;
    seq_t read;
    FILE *tmp_fp[MAX_THREAD];

    fputs("mapping reads that cover small gaps...\n", stderr);

    omp_set_num_threads(n_thread);

#    pragma omp parallel for  schedule(static, 1) \
        private(i, j, st, ed, gap_seq, gap_seq_len, lbuf, lbuf_size, rbuf, rbuf_size, lkmer, rkmer, read, lres, rres)
    for (t = 0; t < n_thread; ++t) {
rres.ofst = 0;
ed = 0;
        lbuf = (position_t *)my_malloc(mp->max_occ * 2 * sizeof(position_t));
        rbuf = (position_t *)my_malloc(mp->max_occ * 2 * sizeof(position_t));

        tmp_fp[t] = tmpfile_open();
        check_tmpfile(tmp_fp[t]);

        rewind(lib[t].pair_fp);
        while (seq_read(&read, lib[t].pair_fp)) {
            if (read.len < 2 * mp->seed_len)
                continue;

            lkmer.fwd = lkmer.rev = rkmer.fwd = rkmer.rev = 0;
            for (i = 0; i < mp->key_len; ++i) {
                if (read.base[i] == 4 || read.base[read.len - mp->seed_len + i] == 4)
                    break;
                lkmer.fwd = (lkmer.fwd << 2) | (uint64_t)read.base[i];
                lkmer.rev = (lkmer.rev << 2) | (uint64_t)(3^read.base[mp->key_len - i - 1]);
                rkmer.fwd = (rkmer.fwd << 2) | (uint64_t)read.base[read.len - mp->seed_len + i];
                rkmer.rev = (rkmer.rev << 2) | (uint64_t)(3^read.base[read.len - mp->seed_len + mp->key_len - i - 1]);
            }
            if (i != mp->key_len)
                continue;

            lbuf_size = mapper_map_seed(mp, &lkmer, lbuf);
            lres.id = 0;
            for (i = 0; i < lbuf_size; ++i) {
                if (lbuf[i].id > 0) {
                    if (lbuf[i].ofst > mp->seq[lbuf[i].id - 1].len - mp->seed_len)
                        continue;
                    for (j = mp->key_len; j < mp->seed_len; ++j) {
                        if (mp->seq[lbuf[i].id - 1].seq[lbuf[i].ofst + j] != read.base[j])
                            break;
                    }
                }
                else {
                    if (lbuf[i].ofst < mp->seed_len - mp->key_len)
                        continue;
                    for (j = mp->key_len; j < mp->seed_len; ++j) {
                        if (mp->seq[-(lbuf[i].id) - 1].seq[lbuf[i].ofst + mp->key_len - j - 1] != (3^read.base[j]))
                            break;
                    }
                    lbuf[i].ofst -= mp->seed_len - mp->key_len;
                }
                if (j != mp->seed_len)
                    lbuf[i].id = 0;;
            }

            rbuf_size = mapper_map_seed(mp, &rkmer, rbuf);
            rres.id = 0;
            for (i = 0; i < rbuf_size; ++i) {
                if (rbuf[i].id > 0) {
                    if (rbuf[i].ofst > mp->seq[rbuf[i].id - 1].len - mp->seed_len)
                        continue;
                    for (j = mp->key_len; j < mp->seed_len; ++j) {
                        if (mp->seq[rbuf[i].id - 1].seq[rbuf[i].ofst + j] != read.base[read.len - mp->seed_len + j])
                            break;
                    }
                }
                else {
                    if (rbuf[i].ofst < mp->seed_len - mp->key_len)
                        continue;
                    for (j = mp->key_len; j < mp->seed_len; ++j) {
                        if (mp->seq[-(rbuf[i].id) - 1].seq[rbuf[i].ofst + mp->key_len - j - 1] != (3^read.base[read.len - mp->seed_len + j]))
                            break;
                    }
                    rbuf[i].ofst -= mp->seed_len - mp->key_len;
                }
                if (j != mp->seed_len)
                    rbuf[i].id = 0;;
            }

            lres.id = 0;
            for (i = 0; i < lbuf_size; ++i) {
                for (j = 0; j < rbuf_size; ++j) {
                    if (lbuf[i].id == 0 || rbuf[j].id == 0 || 
                        lbuf[i].id != rbuf[j].id || 
                        (lbuf[i].id > 0 && (lbuf[i].ofst >= rbuf[j].ofst || rbuf[j].ofst - lbuf[i].ofst + mp->seed_len > 2*read.len )) ||
                        (lbuf[i].id < 0 && (rbuf[j].ofst >= lbuf[i].ofst || lbuf[i].ofst - rbuf[j].ofst + mp->seed_len > 2*read.len ))) {

                        continue;
                    }
                    if (lres.id != 0)
                        break;
                    lres = lbuf[i];
                    rres = rbuf[j];
                }
                if (j != rbuf_size)
                    break;
            }
            if (i != lbuf_size || lres.id == 0)
                continue;

            if (lres.id > 0) {
                st = 0;
                for (i = lres.ofst + mp->seed_len; i < rres.ofst; ++i) {
                    if (mp->seq[lres.id - 1].seq[i] == 4) {
                        st = i - lres.ofst;
                        break;
                    }
                }
                if (st == 0)
                    continue;
                for (i = rres.ofst - 1; i >= lres.ofst + mp->seed_len; --i) {
                    if (mp->seq[lres.id - 1].seq[i] == 4) {
                        ed = read.len - (rres.ofst + mp->seed_len - 1 - i);
                        break;
                    }
                }
                for (i = lres.ofst + st; i < rres.ofst - (read.len - ed - mp->seed_len); ++i) {
                    if (mp->seq[lres.id - 1].seq[i] != 4)
                        break;
                }
                if (i != rres.ofst - (read.len - ed - mp->seed_len))
                    continue;
                gap_seq_len = 0;
                for (i = st; i < ed; ++i) {
                    gap_seq[gap_seq_len] = read.base[i];
                    ++gap_seq_len;
                }
                lres.ofst = lres.ofst + st;
            }
            else {
                st = 0;
                for (i = rres.ofst + mp->seed_len; i < lres.ofst; ++i) {
                    if (mp->seq[-(lres.id) - 1].seq[i] == 4) {
                        st = i - rres.ofst;
                        break;
                    }
                }
                if (st == 0)
                    continue;
                for (i = lres.ofst - 1; i >= rres.ofst + mp->seed_len; --i) {
                    if (mp->seq[-(lres.id) - 1].seq[i] == 4) {
                        ed = read.len - (lres.ofst + mp->seed_len - 1 - i);
                        break;
                    }
                }
                for (i = rres.ofst + st; i < lres.ofst - (read.len - ed - mp->seed_len); ++i) {
                    if (mp->seq[-(lres.id) - 1].seq[i] != 4)
                        break;
                }
                if (i != lres.ofst - (read.len - ed - mp->seed_len))
                    continue;
                gap_seq_len = 0;
                for (i = read.len - st - 1; i >= read.len - ed; --i) {
                    gap_seq[gap_seq_len] = 3^(read.base[i]);
                    ++gap_seq_len;
                }
                lres.ofst = rres.ofst + st;
                lres.id *= -1;
            }

            fwrite(&lres, sizeof(position_t), 1, tmp_fp[t]);
            fwrite(&gap_seq_len, sizeof(int32_t), 1, tmp_fp[t]);
            fwrite(gap_seq, sizeof(int8_t), gap_seq_len, tmp_fp[t]);
        }

        my_free(lbuf);
        my_free(rbuf);
    }

    fseek(*gap_seq_fp, 0, SEEK_END);
    for (t = 0; t < n_thread; ++t) {
        rewind(tmp_fp[t]);
        while(fread(&c, sizeof(int8_t), 1, tmp_fp[t]))
            putc(c, *gap_seq_fp);
        fclose(tmp_fp[t]);
    }
}

void mapper_map_pair_nolink_mt(const mapper_t *mp, seqlib_t *lib, const int64_t n_thread)
{
    int c;
    int64_t i;
    int64_t j;
    int64_t total_pair;
    int64_t n_map_same;
    int64_t ins_len;
    position_t fpos;
    position_t rpos;
    seq_t fwd;
    seq_t rev;
    position_t *buf;

    fputs("mapping reads...\n", stderr);

    omp_set_num_threads(n_thread);

    total_pair = n_map_same = 0;

#    pragma omp parallel for  schedule(static, 1) \
        private(j, buf, fwd, rev, fpos, rpos, ins_len) reduction(+: total_pair, n_map_same)
    for (i = 0; i < n_thread; ++i) {
        buf = (position_t *)my_malloc(mp->max_occ * 2 * sizeof(position_t));

        if (lib[i].ins_len_fp != NULL)
            fclose(lib[i].ins_len_fp);
        lib[i].ins_len_fp = tmpfile_open();
        check_tmpfile(lib[i].ins_len_fp);

        if (lib[i].mapped_fp != NULL)
            fclose(lib[i].mapped_fp);
        lib[i].mapped_fp = tmpfile_open();
        check_tmpfile(lib[i].mapped_fp);

        rewind(lib[i].pair_fp);
        while (seq_read(&fwd, lib[i].pair_fp)) {
            seq_read(&rev, lib[i].pair_fp);
            ++total_pair;

            if (fwd.len < mp->seed_len || rev.len < mp->seed_len)
                continue;

            fpos = mapper_map_read(mp, &fwd, buf);
            rpos = mapper_map_read(mp, &rev, buf);
            
            if (fpos.id == 0 && rpos.id == 0)
                continue;

            fwrite(&fpos, sizeof(position_t), 1, lib[i].mapped_fp);
            for (j = 0; j < fwd.n_unknown; ++j)
                fwd.base[fwd.unknown_pos[j]] = 4;
            fwrite(&(fwd.len), sizeof(uint16_t), 1, lib[i].mapped_fp);
            fwrite(fwd.base, sizeof(uint8_t), fwd.len, lib[i].mapped_fp);

            fwrite(&rpos, sizeof(position_t), 1, lib[i].mapped_fp);
            for (j = 0; j < rev.n_unknown; ++j)
                rev.base[rev.unknown_pos[j]] = 4;
            fwrite(&(rev.len), sizeof(uint16_t), 1, lib[i].mapped_fp);
            fwrite(rev.base, sizeof(uint8_t), rev.len, lib[i].mapped_fp);

            if (fpos.id == -rpos.id) {
                if (fpos.id > 0 && fpos.ofst < rpos.ofst) {
                    ins_len = (int64_t)rpos.ofst - fpos.ofst;
                    if (ins_len - fwd.len - rev.len > 0 && 
                        memchr(&(mp->seq[fpos.id - 1].seq[fpos.ofst + fwd.len]), 4, ins_len - fwd.len - rev.len) != NULL)
                        continue;
                }
                else if (rpos.id > 0 && rpos.ofst < fpos.ofst) {
                    ins_len = (int64_t)fpos.ofst - rpos.ofst;
                    if (ins_len - fwd.len - rev.len > 0 && 
                        memchr(&(mp->seq[rpos.id - 1].seq[rpos.ofst + rev.len]), 4, ins_len - fwd.len - rev.len) != NULL)
                        continue;
                }
                else
                    continue;

                if (ins_len < fwd.len || ins_len < rev.len)
					continue;

                fwrite(&ins_len, sizeof(int64_t), 1, lib[i].ins_len_fp);
                ++n_map_same;
            }
        }
        my_free(buf);
    }

    if (n_map_same == 0) {
        fprintf(stderr, "error: no read mapped in the same contig!\n");
        my_exit(1);
    }

    for (i = 1; i < n_thread; ++i) {
        rewind(lib[i].ins_len_fp);
        while (fread(&c, sizeof(int8_t), 1, lib[i].ins_len_fp))
            putc(c, lib->ins_len_fp);
        fclose(lib[i].ins_len_fp);
    }

    fprintf(stderr, "TOTAL_PAIR = %ld\n", total_pair);
    fprintf(stderr, "MAPPED_IN_SAME_CONTIG = %ld (%f)\n", n_map_same, (double)n_map_same / total_pair);
}

void hetero_mapper_init(hetero_mapper_t *mp, const int64_t seed_len, const int64_t key_len, const int64_t mem_lim)
{
    memset(mp, 0, sizeof(hetero_mapper_t));
    mp->seed_len = mp->cmp.seed_len = mp->bmp.seed_len = seed_len;
    mp->key_len = mp->cmp.key_len = mp->bmp.key_len = key_len;
    if (mp->cmp.key_len > 32)
        mp->key_len = mp->cmp.key_len = mp->bmp.key_len = 32;
    mp->cmp.mask = mp->bmp.mask = mp->cmp.key_len < 32 ? ~(~0ull << (2*mp->cmp.key_len)) : ~0ull;
    mp->cmp.mem_lim = mp->bmp.mem_lim = mem_lim;
}

void hetero_mapper_destroy(hetero_mapper_t *mp)
{
    mapper_destroy(&(mp->cmp));
    mapper_destroy(&(mp->bmp));
    if (mp->bubble_map != NULL)
        my_free(mp->bubble_map);
    memset(mp, 0, sizeof(hetero_mapper_t));
}

void hetero_mapper_merge_mt(hetero_mapper_t *mp, const int64_t n_thread)
{
    int64_t i;
    int64_t j;
    int64_t t;
    int64_t max_l_len;
    int64_t max_r_len;
    int64_t buf_size;
    position_t *buf;
    position_t lres = {0};
    position_t rres = {0};
    kmer_t lkmer;
    kmer_t rkmer;
    lseq_t *bub;

    if (mp->bmp.n_seq == 0)
        return;

    fputs("merging contigs and bubbles...\n", stderr);

    if (mp->bubble_map != NULL)
        my_free(mp->bubble_map);
    mp->bubble_map = (position_t *)my_calloc(mp->bmp.n_seq, sizeof(position_t));
    
    omp_set_num_threads(n_thread);

#    pragma omp parallel for  schedule(static, 1) \
        private(i, j, buf, buf_size, bub, lkmer, rkmer, max_l_len, max_r_len) \
        firstprivate(lres, rres)
    for (t = 0; t < n_thread; ++t) {
        buf = (position_t *)my_malloc(mp->cmp.max_occ * 2 * sizeof(position_t));

        for (bub = &(mp->bmp.seq[t]); bub - mp->bmp.seq < mp->bmp.n_seq; bub += n_thread) {
            if (bub->len < 2 * mp->key_len)
                continue;

            lkmer.fwd = lkmer.rev = rkmer.fwd = rkmer.rev = 0;
            for (i = 0; i < mp->key_len; ++i) {
                if (bub->seq[i] == 4 || bub->seq[bub->len - mp->seed_len + i] == 4)
                    break;
                lkmer.fwd = (lkmer.fwd << 2) | (uint64_t)bub->seq[i];
                lkmer.rev = (lkmer.rev << 2) | (uint64_t)(3^bub->seq[mp->key_len - i - 1]);
                rkmer.fwd = (rkmer.fwd << 2) | (uint64_t)bub->seq[bub->len - mp->key_len + i];
                rkmer.rev = (rkmer.rev << 2) | (uint64_t)(3^bub->seq[bub->len - i - 1]);
            }
            if (i != mp->key_len)
                continue;

            buf_size = mapper_map_seed(&(mp->cmp), &lkmer, buf);
            max_l_len = 0;
            for (i = 0; i < buf_size; ++i) {
                if (buf[i].id > 0) {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + j >= mp->cmp.seq[buf[i].id - 1].len ||
                            mp->cmp.seq[buf[i].id - 1].seq[buf[i].ofst + j] != bub->seq[j])
                            break;
                    }
                }
                else {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + mp->key_len - j <= 0  || 
                            mp->cmp.seq[-(buf[i].id) - 1].seq[buf[i].ofst + mp->key_len - j - 1] != (3^bub->seq[j]))
                            break;
                    }
                    buf[i].ofst += mp->key_len - 1;
                }
                if (j > max_l_len) {
                    max_l_len = j;
                    lres =  buf[i];
                }
            }
            if (max_l_len < mp->seed_len)
                continue;

            buf_size = mapper_map_seed(&(mp->cmp), &rkmer, buf);
            max_r_len = 0;
            for (i = 0; i < buf_size; ++i) {
                if (buf[i].id > 0) {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + mp->key_len - j <= 0 ||
                            mp->cmp.seq[buf[i].id - 1].seq[buf[i].ofst + mp->key_len - j - 1] != bub->seq[bub->len - j - 1])
                            break;
                    }
                    buf[i].ofst += mp->key_len - 1;
                }
                else {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + j >= mp->cmp.seq[-(buf[i].id) - 1].len || 
                            mp->cmp.seq[-(buf[i].id) - 1].seq[buf[i].ofst + j] != (3^bub->seq[bub->len - j - 1]))
                            break;
                    }
                }
                if (j > max_r_len) {
                    max_r_len = j;
                    rres =  buf[i];
                }
            }
            if (max_r_len < mp->seed_len || rres.id != lres.id)
                continue;

            mp->bubble_map[bub - mp->bmp.seq].id = lres.id;
            mp->bubble_map[bub - mp->bmp.seq].ofst = (lres.ofst + rres.ofst) / 2;
        }
        my_free(buf);
    }
}

position_t hetero_mapper_map_read(const hetero_mapper_t *mp, const seq_t *read, position_t *buf)
{
    int64_t i;
    int64_t j;
    int64_t k = 0;
    int64_t buf_size;
    int64_t max;
    int64_t sec_max;
    int64_t n_res;
    position_t res[MAX_READ_LEN];
    kmer_t kmer;

    n_res = 0;

    for (i = read->len - mp->seed_len; i > -(mp->seed_len); i -= mp->seed_len) {
        if (i < 0)
            i = 0;
        kmer.fwd = kmer.rev = 0;
        for (j = 0; j < mp->key_len; ++j) {
            if (read->base[i + j] == 4)
                break;
            kmer.fwd = (kmer.fwd << 2) | (uint64_t)read->base[i + j];
            kmer.rev = (kmer.rev << 2) | (uint64_t)(3^read->base[i + mp->key_len - j - 1]);
        }
        if (j != mp->key_len)
            continue;

        buf_size = mapper_map_seed(&(mp->cmp), &kmer, buf);
        res[n_res].id = 0;
        for (j = 0; j < buf_size; ++j) {
            if (buf[j].id > 0) {
                if (buf[j].ofst > mp->cmp.seq[buf[j].id - 1].len - mp->seed_len)
                    continue;
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->cmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k])
                        break;
                }
                buf[j].ofst -= i;
            }
            else {
                if (buf[j].ofst < mp->seed_len - mp->key_len)
                    continue;
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->cmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k]))
                        break;
                }
                buf[j].ofst += i + mp->key_len - 1;
            }

            if (k == mp->seed_len) {
                if (res[n_res].id == 0)
                    res[n_res] = buf[j];
                else
                    break;
            }
        }
        if (j == buf_size) {
            if (res[n_res].id != 0) {
                ++n_res;
                continue;
            }
        }
        else
            continue;

        buf_size = mapper_map_seed(&(mp->bmp), &kmer, buf);
        res[n_res].id = 0;
        for (j = 0; j < buf_size; ++j) {
            if (buf[j].id > 0) {
                if (buf[j].ofst > mp->bmp.seq[buf[j].id - 1].len - mp->seed_len)
                    continue;
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->bmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k])
                        break;
                }
                buf[j].ofst -= i;
                if (mp->bubble_map[buf[j].id - 1].id > 0) {
                    buf[j].ofst = mp->bubble_map[buf[j].id - 1].ofst + buf[j].ofst - mp->bmp.seq[buf[j].id - 1].len / 2;
                    buf[j].id = mp->bubble_map[buf[j].id - 1].id;
                }
                else {
                    buf[j].ofst = mp->bubble_map[buf[j].id - 1].ofst - buf[j].ofst + mp->bmp.seq[buf[j].id - 1].len / 2;
                    buf[j].id = mp->bubble_map[buf[j].id - 1].id;
                }
            }
            else {
                if (buf[j].ofst < mp->seed_len - mp->key_len)
                    continue;
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->bmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k]))
                        break;
                }
                buf[j].ofst += i + mp->key_len - 1;
                if (mp->bubble_map[-(buf[j].id) - 1].id > 0) {
                    buf[j].ofst = mp->bubble_map[-(buf[j].id) - 1].ofst + buf[j].ofst - mp->bmp.seq[-(buf[j].id) - 1].len / 2;
                    buf[j].id = -(mp->bubble_map[-(buf[j].id) - 1].id);
                }
                else {
                    buf[j].ofst = mp->bubble_map[-(buf[j].id) - 1].ofst - buf[j].ofst + mp->bmp.seq[-(buf[j].id) - 1].len / 2;
                    buf[j].id = -(mp->bubble_map[-(buf[j].id) - 1].id);
                }
            }

            if (k == mp->seed_len) {
                if (res[n_res].id == 0)
                    res[n_res] = buf[j];
                else
                    break;
            }
        }
        if (j == buf_size && res[n_res].id != 0)
            ++n_res;
    }

    if (n_res == 0) {
        res[0].id = 0;
        return res[0];
    }

    qsort(res, n_res, sizeof(position_t), cmp_position);

    i = max = sec_max = res[n_res].id = 0;
    while (res[i].id != 0) {
        j = i;
        while (res[i].id == res[i+1].id)
            ++i;
        ++i;
        if (i - j >= max) {
            sec_max = max;
            max = i - j;
            k = j;
        }
    }

    if (max == sec_max) {
        res[0].id = 0;
        return res[0];
    }

    j = 0;
    for (i = 0; i < max; ++i)
        j += res[k + i].ofst;
    res[k].ofst = (double)j / max + 0.5;

    return res[k];
}

void hetero_mapper_map_pair_mt(const hetero_mapper_t *mp, seqlib_t *lib, const int64_t min_ins, const int64_t n_thread)
{
    int c;
    int64_t i;
    int64_t sum;
    int64_t total_pair;
    int64_t n_map_same;
    int64_t n_map_diff;
    int64_t ins_len;
    position_t fpos;
    position_t rpos;
    seq_t fwd;
    seq_t rev;
    position_t *buf;

    if (mp->bmp.n_seq == 0) {
        mapper_map_pair_mt(&(mp->cmp), lib, min_ins, n_thread);
        return;
    }

    fputs("mapping reads...\n", stderr);
    fprintf(stderr, "MIN_INS_LEN = %ld\n", min_ins);

    omp_set_num_threads(n_thread);

    sum = total_pair = n_map_same = n_map_diff = 0;

#    pragma omp parallel for  schedule(static, 1) \
        private(buf, fwd, rev, fpos, rpos, ins_len) \
        reduction(+: sum, total_pair, n_map_same, n_map_diff)
    for (i = 0; i < n_thread; ++i) {
        if (mp->cmp.max_occ > mp->bmp.max_occ)
            buf = (position_t *)my_malloc(mp->cmp.max_occ * 2 * sizeof(position_t));
        else
            buf = (position_t *)my_malloc(mp->bmp.max_occ * 2 * sizeof(position_t));

        if (lib[i].ins_len_fp != NULL)
            fclose(lib[i].ins_len_fp);
        lib[i].ins_len_fp = tmpfile_open();
        check_tmpfile(lib[i].ins_len_fp);

        if (lib[i].mapped_fp != NULL)
            fclose(lib[i].mapped_fp);
        lib[i].mapped_fp = tmpfile_open();
        check_tmpfile(lib[i].mapped_fp);

        rewind(lib[i].pair_fp);
        while (seq_read(&fwd, lib[i].pair_fp)) {
            seq_read(&rev, lib[i].pair_fp);
            ++total_pair;

            if (fwd.len < mp->seed_len || rev.len < mp->seed_len)
                continue;

            fpos = hetero_mapper_map_read(mp, &fwd, buf);
            rpos = hetero_mapper_map_read(mp, &rev, buf);

            if (fpos.id == 0 || rpos.id == 0)
                continue;

            if (fpos.id == -rpos.id) {
                if (fpos.id > 0 && fpos.ofst < rpos.ofst)
                    ins_len = (int64_t)rpos.ofst - fpos.ofst;
                else if (rpos.id > 0 && rpos.ofst < fpos.ofst)
                    ins_len = (int64_t)fpos.ofst - rpos.ofst;
                else
                    continue;

                if (ins_len < min_ins || ins_len < fwd.len || ins_len < rev.len)
                    continue;

                fwrite(&ins_len, sizeof(int64_t), 1, lib[i].ins_len_fp);
                ++n_map_same;
            }
            else if (fpos.id != rpos.id) {
                fwrite(&fpos, sizeof(position_t), 1, lib[i].mapped_fp);
                fwrite(&rpos, sizeof(position_t), 1, lib[i].mapped_fp);
                ++n_map_diff;
            }

            sum += fwd.len + rev.len;
        }
        my_free(buf);
    }

    if (n_map_same == 0)
        fprintf(stderr, "WARNING: no read mapped in the same contig!\n");
    if (n_map_diff == 0)
        fprintf(stderr, "WARNING: no read links contigs!\n");

    lib->ave_cov = (double)sum / mp->cmp.seq_pool_size;

#    pragma omp parallel sections private(i, c)
    {
#        pragma omp section
        {
            for (i = 1; i < n_thread; ++i) {
                rewind(lib[i].ins_len_fp);
                while (fread(&c, sizeof(int8_t), 1, lib[i].ins_len_fp))
                    putc(c, lib->ins_len_fp);
                fclose(lib[i].ins_len_fp);
            }
        }

#        pragma omp section
        {
            for (i = 1; i < n_thread; ++i) {
                rewind(lib[i].mapped_fp);
                while(fread(&c, sizeof(int8_t), 1, lib[i].mapped_fp))
                    putc(c, lib->mapped_fp);
                fclose(lib[i].mapped_fp);
            }
        }
    }

    fprintf(stderr, "TOTAL_PAIR = %ld\n", total_pair);
    fprintf(stderr, "MAPPED_PAIR = %ld (%f)\n", n_map_same + n_map_diff, (double)(n_map_same + n_map_diff) / total_pair);
    fprintf(stderr, "MAPPED_IN_DIFFERENT_CONTIGS = %ld (%f)\n", n_map_diff, (double)n_map_diff / total_pair);
    fprintf(stderr, "MAPPED_IN_SAME_CONTIG = %ld (%f)\n", n_map_same, (double)n_map_same / total_pair);
    fprintf(stderr, "AVERAGE_COVERAGE = %f\n", lib->ave_cov);
}

void mapper_map_bubble(const mapper_t *mp, const contig_t *bubble, const int64_t min_match_len, FILE **pos_fp, const int64_t n_thread)
{
    int64_t i;
    int64_t j;
    int64_t max_l_len;
    int64_t max_r_len;
    int64_t t;
    int64_t buf_size;
    bool is_mapped;
    position_t *buf;
    position_t lres = {0};
    position_t rres = {0};
    kmer_t lkmer;
    kmer_t rkmer;
    lseq_t *bub;
    varr8_t tmp_seq;

    fputs("mapping bubbles...\n", stderr);

    omp_set_num_threads(n_thread);

#    pragma omp parallel for  schedule(static, 1) \
        private(i, j, is_mapped, buf, buf_size, bub, lkmer, rkmer, max_l_len, max_r_len, tmp_seq) \
        firstprivate(lres, rres)
    for (t = 0; t < n_thread; ++t) {

        pos_fp[t] = tmpfile_open();
        check_tmpfile(pos_fp[t]);

        varr8_init(&tmp_seq);

        buf = (position_t *)my_malloc(mp->max_occ * 2 * sizeof(position_t));

        for (bub = &(bubble->seq[t]); bub - bubble->seq < bubble->n_seq; bub += n_thread) {
            is_mapped = false;
            if (bub->len < 2 * mp->key_len) {
                fwrite(&is_mapped, sizeof(bool), 1, pos_fp[t]);
                continue;
            }

            lkmer.fwd = lkmer.rev = rkmer.fwd = rkmer.rev = 0;
            for (i = 0; i < mp->key_len; ++i) {
                if (bub->seq[i] == 4 || bub->seq[bub->len - mp->seed_len + i] == 4)
                    break;
                lkmer.fwd = (lkmer.fwd << 2) | (uint64_t)bub->seq[i];
                lkmer.rev = (lkmer.rev << 2) | (uint64_t)(3^bub->seq[mp->key_len - i - 1]);
                rkmer.fwd = (rkmer.fwd << 2) | (uint64_t)bub->seq[bub->len - mp->key_len + i];
                rkmer.rev = (rkmer.rev << 2) | (uint64_t)(3^bub->seq[bub->len - i - 1]);
            }
            if (i != mp->key_len) {
                fwrite(&is_mapped, sizeof(bool), 1, pos_fp[t]);
                continue;
            }

            buf_size = mapper_map_seed(mp, &lkmer, buf);
            max_l_len = 0;
            for (i = 0; i < buf_size; ++i) {
                if (buf[i].id > 0) {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + j >= mp->seq[buf[i].id - 1].len ||
                            mp->seq[buf[i].id - 1].seq[buf[i].ofst + j] != bub->seq[j] ||
                            bub->seq[j] == 4)
                            break;
                    }
                    buf[i].ofst += j;
                }
                else {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + mp->key_len - j <= 0  || 
                            mp->seq[-(buf[i].id) - 1].seq[buf[i].ofst + mp->key_len - j - 1] != (3^bub->seq[j]))
                            break;
                    }
                    buf[i].ofst += mp->key_len - 1 - j;
                }
                if (j >= max_l_len) {
                    lres = buf[i];
                    if (j > max_l_len)
                        max_l_len = j;
                    else
                        lres.id = 0;
                }
            }
            if (max_l_len < min_match_len || lres.id == 0) {
                fwrite(&is_mapped, sizeof(bool), 1, pos_fp[t]);
                continue;
            }

            buf_size = mapper_map_seed(mp, &rkmer, buf);
            max_r_len = 0;
            for (i = 0; i < buf_size; ++i) {
                if (buf[i].id > 0) {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + mp->key_len - j <= 0 ||
                            mp->seq[buf[i].id - 1].seq[buf[i].ofst + mp->key_len - j - 1] != bub->seq[bub->len - j - 1] ||
                            bub->seq[bub->len - j - 1] == 4)
                            break;
                    }
                    buf[i].ofst += mp->key_len - 1 - j;
                }
                else {
                    for (j = mp->key_len; j < bub->len; ++j) {
                        if (buf[i].ofst + j >= mp->seq[-(buf[i].id) - 1].len || 
                            mp->seq[-(buf[i].id) - 1].seq[buf[i].ofst + j] != (3^bub->seq[bub->len - j - 1]))
                            break;
                    }
                    buf[i].ofst += j;
                }
                if (j >= max_r_len) {
                    rres = buf[i];
                    if (j > max_r_len)
                        max_r_len = j;
                    else
                        rres.id = 0;
                }
            }
            if (max_r_len < min_match_len || rres.id == 0 || rres.id != lres.id) {
                fwrite(&is_mapped, sizeof(bool), 1, pos_fp[t]);
                continue;
            }

            if (lres.id > 0) {
                if (bub->len < max_r_len + max_l_len || lres.ofst > rres.ofst + 1) {
                    j = max_64(max_r_len + max_l_len - bub->len, lres.ofst - rres.ofst - 1);
                    max_l_len -= j;
                    lres.ofst -= j;
                }
                varr8_resize(&tmp_seq, bub->len - max_r_len - max_l_len);
                for (j = 0; j < tmp_seq.size; ++j)
                    tmp_seq.ele[j] = bub->seq[max_l_len + j];
            }
            else {
                j = lres.ofst;
                lres.ofst = rres.ofst;
                rres.ofst = j;
                lres.id *= -1;

                if (bub->len < max_r_len + max_l_len || lres.ofst > rres.ofst + 1) {
                    j = max_64(max_r_len + max_l_len - bub->len, lres.ofst - rres.ofst - 1);
                    max_l_len -= j;
                    lres.ofst -= j;
                }
                varr8_resize(&tmp_seq, bub->len - max_r_len - max_l_len);
                for (j = 0; j < tmp_seq.size; ++j) {
                    if (bub->seq[max_l_len + j] != 4)
                        tmp_seq.ele[tmp_seq.size - j - 1] = (3^bub->seq[max_l_len + j]);
                    else
                        tmp_seq.ele[tmp_seq.size - j - 1] = 4;
                }
            }

            is_mapped = true;
            fwrite(&is_mapped, sizeof(bool), 1, pos_fp[t]);
            fwrite(&(tmp_seq.size), sizeof(uint64_t), 1, pos_fp[t]);
            fwrite(tmp_seq.ele, sizeof(int8_t), tmp_seq.size, pos_fp[t]);

            varr8_resize(&tmp_seq, rres.ofst - lres.ofst + 1);
            for (j = 0; j < tmp_seq.size; ++j)
                tmp_seq.ele[j] = mp->seq[lres.id - 1].seq[lres.ofst + j];
            fwrite(&(tmp_seq.size), sizeof(uint64_t), 1, pos_fp[t]);
            fwrite(tmp_seq.ele, sizeof(int8_t), tmp_seq.size, pos_fp[t]);
            fwrite(&lres, sizeof(position_t), 1, pos_fp[t]);
        }

        varr8_destroy(&tmp_seq);
        my_free(buf);
    }
}

int8_t hetero_mapper_map_read_ungap(const hetero_mapper_t *mp, const seq_t *read, const double min_idt, position_t *buf)
{
    int64_t i;
    int64_t j;
    int64_t k = 0;
    int64_t buf_size;
    int64_t n_hit;
    int64_t n_mis;
    int64_t min_mis;
    int64_t seed_hit;
    kmer_t kmer;

	min_mis = (int64_t)((1.0 - min_idt) * read->len + 0.5);
	n_hit = 0;

    for (i = read->len - mp->seed_len; i > -(mp->seed_len); i -= mp->seed_len) {
        if (i < 0)
            i = 0;
        kmer.fwd = kmer.rev = 0;
        for (j = 0; j < mp->key_len; ++j) {
            if (read->base[i + j] == 4)
                break;
            kmer.fwd = (kmer.fwd << 2) | (uint64_t)read->base[i + j];
            kmer.rev = (kmer.rev << 2) | (uint64_t)(3^read->base[i + mp->key_len - j - 1]);
        }
        if (j != mp->key_len)
            continue;

        buf_size = mapper_map_seed(&(mp->cmp), &kmer, buf);
		seed_hit = -1;
        for (j = 0; j < buf_size; ++j) {
			n_mis = 0;
            if (buf[j].id > 0) {
                if (buf[j].ofst > mp->cmp.seq[buf[j].id - 1].len - read->len || buf[j].ofst < i) {
					buf[j].id = 0;
                    continue;
				}
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->cmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k])
                        break;
                }
            }
            else {
                if (buf[j].ofst < read->len - mp->key_len) {
					buf[j].id = 0;
                    continue;
				}
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->cmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k]))
                        break;
                }
            }
			if (k == mp->seed_len) {
				if (seed_hit >= 0)
					break;
				seed_hit = j;
			}
			else
				buf[j].id = 0;
		}
		if (j != buf_size || seed_hit < 0)
			continue;

		j = seed_hit;
		n_mis = 0;
		if (buf[j].id > 0) {
			for (k = mp->seed_len; k < read->len; ++k) {
				if (mp->cmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k]) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k != read->len)
				continue;
			for (k = -i; k < 0; ++k) {
				if (mp->cmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k]) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k == i)
				buf[j].ofst -= i;
		}
		else {
			for (k = mp->seed_len; k < read->len; ++k) {
				if (mp->cmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k])) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k != read->len)
				continue;
			for (k = -i; k < 0; ++k) {
				if (mp->cmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k])) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k == i)
				buf[j].ofst += i + mp->key_len - 1;

		}

		if (n_mis < min_mis || n_hit == 0) {
			n_hit = 1;
			min_mis = n_mis;
			buf[0] = buf[j];
		}
		else if (n_hit > 0) {
			if (n_mis == 0)
				return 0;
		}
		else
			++n_hit;


        buf_size = mapper_map_seed(&(mp->bmp), &kmer, buf);
		seed_hit = -1;
        for (j = 0; j < buf_size; ++j) {
			n_mis = 0;
            if (buf[j].id > 0) {
                if (buf[j].ofst > mp->bmp.seq[buf[j].id - 1].len - read->len || buf[j].ofst < i) {
					buf[j].id = 0;
                    continue;
				}
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->bmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k])
                        break;
                }
            }
            else {
                if (buf[j].ofst < read->len - mp->key_len) {
					buf[j].id = 0;
                    continue;
				}
                for (k = mp->key_len; k < mp->seed_len; ++k) {
                    if (mp->bmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k]))
                        break;
                }
            }
			if (k == mp->seed_len) {
				if (seed_hit >= 0)
					break;
				seed_hit = j;
			}
			else
				buf[j].id = 0;
		}
		if (j != buf_size || seed_hit < 0)
			continue;

		j = seed_hit;
		n_mis = 0;
		if (buf[j].id > 0) {
			for (k = mp->seed_len; k < read->len; ++k) {
				if (mp->bmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k]) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k != read->len)
				continue;
			for (k = -i; k < 0; ++k) {
				if (mp->bmp.seq[buf[j].id - 1].seq[buf[j].ofst + k] != read->base[i + k]) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k != i)
				continue;

			buf[j].ofst -= i;
			if (mp->bubble_map[buf[j].id - 1].id > 0) {
				buf[j].ofst = mp->bubble_map[buf[j].id - 1].ofst + buf[j].ofst - mp->bmp.seq[buf[j].id - 1].len / 2;
				buf[j].id = mp->bubble_map[buf[j].id - 1].id;
			}
			else {
				buf[j].ofst = mp->bubble_map[buf[j].id - 1].ofst - buf[j].ofst + mp->bmp.seq[buf[j].id - 1].len / 2;
				buf[j].id = mp->bubble_map[buf[j].id - 1].id;
			}
		}
		else {
			for (k = mp->seed_len; k < read->len; ++k) {
				if (mp->bmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k])) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k != read->len)
				continue;
			for (k = -i; k < 0; ++k) {
				if (mp->bmp.seq[-(buf[j].id) - 1].seq[buf[j].ofst + mp->key_len - k - 1] != (3^read->base[i + k])) {
					++n_mis;
					if (n_mis > min_mis)
						break;
				}
			}
			if (k != i)
				continue;

			buf[j].ofst += i + mp->key_len - 1;
			if (mp->bubble_map[-(buf[j].id) - 1].id > 0) {
				buf[j].ofst = mp->bubble_map[-(buf[j].id) - 1].ofst + buf[j].ofst - mp->bmp.seq[-(buf[j].id) - 1].len / 2;
				buf[j].id = -(mp->bubble_map[-(buf[j].id) - 1].id);
			}
			else {
				buf[j].ofst = mp->bubble_map[-(buf[j].id) - 1].ofst - buf[j].ofst + mp->bmp.seq[-(buf[j].id) - 1].len / 2;
				buf[j].id = -(mp->bubble_map[-(buf[j].id) - 1].id);
			}
		}

		if (n_mis < min_mis || n_hit == 0) {
			n_hit = 1;
			min_mis = n_mis;
			buf[0] = buf[j];
		}
		else if (n_hit > 0) {
			if (n_mis == 0)
				return 0;
		}
		else
			++n_hit;
    }
	if (n_hit == 1)
		return (int8_t)((double)min_mis / read->len);

	return -1;
}
