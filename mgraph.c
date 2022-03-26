#include "mgraph.h"

mstra_t g_mstra_deleted;

void mgraph_init(mgraph_t *mg, const double branch_th, const double bubble_th)
{
    memset(mg, 0, sizeof(mgraph_t));
    mg->bubble_th = bubble_th;
    mg->branch_th = branch_th;
}

void mgraph_destroy_table(mgraph_t *mg)
{
    uint64_t i;

    if (mg->jtable != NULL) {
        my_free(mg->jtable);
        mg->jtable = NULL;
        mg->jtable_size = 0;
    }
    mg->jtable = NULL;

    if (mg->lstable != NULL) {
        for (i = 0; i < mg->stable_size; ++i) {
            if (!mstra_is_emp(mg->lstable[i]) && !mstra_is_del(mg->lstable[i]))
                my_free(mg->lstable[i]);
        }
        my_free(mg->lstable);
        mg->lstable = NULL;
        my_free(mg->rstable);
        mg->rstable = NULL;
        mg->stable_size = 0;
    }
    mg->lstable = NULL;
    mg->rstable = NULL;
}

void mgraph_destroy(mgraph_t *mg)
{
    mgraph_destroy_table(mg);
    if (mg->junc_file != NULL)
        fclose(mg->junc_file);
    if (mg->stra_file != NULL)
        fclose(mg->stra_file);
    memset(mg, 0, sizeof(mgraph_t));
}

void mgraph_save(mgraph_t *mg, const mcounter_t *mc, FILE *sorted_key_fp)
{
    uint64_t j;
    uint64_t k_len;
    uint64_t sum;
    uint64_t mask;
    uint16_t *u16_p;
    uint16_t *tmp_u16_p;
    uint16_t *next_u16_p = NULL;
    uint8_t rflags;
    uint8_t lflags;
    uint8_t tmp_flags;
    kmer_t kmer;
    kmer_t st_kmer;
    kmer_t rkmer;
    kmer_t lkmer;
    varr8_t rseq;
    varr8_t lseq;
    varr8_t seq;
    varr64_t stack;
    mjunc_t mjunc;
    mstra_t mstra;

    mg->mem_lim = mc->mem_lim;
    k_len = mg->k_len = mc->k_len;
    mask = mc->k_mask;
    mg->n_junc = mg->n_stra = 0;
    varr8_init(&rseq);
    varr8_init(&lseq);
    varr8_init(&seq);
    varr64_init(&stack);

    fputs("connecting kmers...\n", stderr);

    if (mg->junc_file != NULL)
        fclose(mg->junc_file);
    mg->junc_file = tmpfile_open();
    check_tmpfile(mg->junc_file);

    if (mg->stra_file != NULL)
        fclose(mg->stra_file);
    mg->stra_file = tmpfile_open();
    check_tmpfile(mg->stra_file);

    varr64_resize(&stack, 1);
    varr64_resize(&stack, 0);
    rewind(sorted_key_fp);
    while (fread(stack.ele, sizeof(uint64_t), 1, sorted_key_fp)) {
        ++(stack.size);
        while (stack.size > 0) {
            st_kmer.fwd = varr64_peek(&stack);    
            st_kmer.rev = 0ull;
            for (j = 0; j < k_len; ++j)
                st_kmer.rev = (st_kmer.rev << 2) | (((st_kmer.fwd >> (2*j))&3)^3);
            u16_p = mcounter_ptr(mc, min_u64(st_kmer.fwd, st_kmer.rev));
            if (*u16_p == UINT16_MAX) {
                varr64_pop(&stack);
                continue;
            }

            sum = *u16_p;
            *u16_p = UINT16_MAX;

            lkmer.fwd = st_kmer.fwd >> 2;
            lkmer.rev = (st_kmer.rev << 2) & mask;
            varr8_clear(&lseq);
            lflags = 0;
            for (j = 0; j < 4; ++j) {
                lkmer.fwd = (lkmer.fwd & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1)));
                lkmer.rev = (lkmer.rev & ~3ull) | (uint64_t)(3^j);
                if (mcounter_occ(mc, min_u64(lkmer.fwd, lkmer.rev)))
                    lflags |= 1 << j;
            }
            if (lflags != 0)
                varr8_push(&lseq, flag_base(lflags));

            rkmer.fwd = (st_kmer.fwd << 2) & mask;
            rkmer.rev = st_kmer.rev >> 2;
            varr8_clear(&rseq);
            rflags = 0;
            for (j = 0; j < 4; ++j) {
                rkmer.fwd = (rkmer.fwd & ~3ull) | (uint64_t)j;
                rkmer.rev = (rkmer.rev & (mask >> 2)) | ((uint64_t)(3^j) << (2*(k_len-1)));
                if (mcounter_occ(mc, min_u64(rkmer.fwd, rkmer.rev)))
                    rflags |= 1 << j;
            }
            if (rflags != 0)
                varr8_push(&rseq, flag_base(rflags));

            if ((lseq.size > 0 && lseq.ele[0] == 4) || (rseq.size > 0 && rseq.ele[0] == 4)) {
                varr64_pop(&stack);
                for (j = 0; j < 4; ++j) {
                    if (rflags & (1 << j))
                        varr64_push(&stack, (rkmer.fwd & ~3ull) | (uint64_t)j);
                    if (lflags & (1 << j)) 
                        varr64_push(&stack, (lkmer.fwd & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1))));
                }
                mjunc.key = st_kmer.fwd;
                mjunc.cov = sum;
                mjunc.out = (lflags << 4) | rflags;
                fwrite(&mjunc, sizeof(mjunc_t), 1, mg->junc_file);
                ++mg->n_junc;
                continue;
            }

            if (rseq.size > 0) {
                kmer.fwd = ((st_kmer.fwd << 2) & mask) | (uint64_t)(rseq.ele[0]);
                kmer.rev = (st_kmer.rev >> 2) | ((uint64_t)(3^rseq.ele[0]) << (2*(k_len-1)));
                u16_p = mcounter_ptr(mc, min_u64(kmer.fwd, kmer.rev));
                while (1) {
                    rkmer = kmer;
                    kmer.fwd = kmer.fwd >> 2;
                    kmer.rev = (kmer.rev << 2) & mask;
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer.fwd = (kmer.fwd & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1)));
                        kmer.rev = (kmer.rev & ~3ull) | (uint64_t)(3^j);
                        if (mcounter_occ(mc, min_u64(kmer.fwd, kmer.rev)))
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        rflags = 0 | (1 << rseq.ele[rseq.size-1]);
                        varr8_pop(&rseq);
                        break;
                    }
                    kmer.fwd = ((kmer.fwd << 4) & mask) | (uint64_t)(rseq.ele[rseq.size-1] << 2);
                    kmer.rev = (kmer.rev >> 4) | ((uint64_t)(3^rseq.ele[rseq.size-1]) << (2*(k_len-2)));
                    rflags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer.fwd = (kmer.fwd & ~3ull) | (uint64_t)j;
                        kmer.rev = (kmer.rev & (mask >> 2)) | ((uint64_t)(3^j) << (2*(k_len-1)));
                        tmp_u16_p = mcounter_ptr(mc, min_u64(kmer.fwd, kmer.rev));
                        if (tmp_u16_p != NULL) {
                            rflags |= 1 << j;
                            next_u16_p = tmp_u16_p;
                        }
                    }
                    if (rflags == 0 ) {
                        sum += *u16_p;
                        *u16_p = UINT16_MAX;
                        break;
                    }
                    else if (flag_degree(rflags) > 1) {
                        rflags = 0 | (1 << rseq.ele[rseq.size-1]);
                        varr8_pop(&rseq);
                        break;
                    }
                    varr8_push(&rseq, flag_base(rflags));
                    sum += *u16_p;
                    *u16_p = UINT16_MAX;
                    u16_p = next_u16_p;
                }
            }

            if (lseq.size > 0) {
                kmer.fwd = (st_kmer.fwd >> 2) | ((uint64_t)(lseq.ele[0]) << (2*(k_len-1)));
                kmer.rev = ((st_kmer.rev << 2) & mask) | (uint64_t)(3^lseq.ele[0]);
                u16_p = mcounter_ptr(mc, min_u64(kmer.fwd, kmer.rev));
                while (1) {
                    lkmer = kmer;
                    kmer.fwd = (kmer.fwd << 2) & mask;
                    kmer.rev = kmer.rev >> 2;
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer.fwd = (kmer.fwd & ~3ull) | (uint64_t)j;
                        kmer.rev = (kmer.rev & (mask >> 2)) | ((uint64_t)(3^j) << (2*(k_len-1)));
                        if (mcounter_occ(mc, min_u64(kmer.fwd, kmer.rev)))
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        lflags = 0 | (1 << lseq.ele[lseq.size-1]);
                        varr8_pop(&lseq);
                        break;
                    }
                    kmer.fwd = (kmer.fwd >> 4) | ((uint64_t)(lseq.ele[lseq.size-1]) << (2*(k_len-2)));
                    kmer.rev = ((kmer.rev << 4) & mask) | ((uint64_t)(3^lseq.ele[lseq.size-1]) << 2);
                    lflags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer.fwd = (kmer.fwd & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1)));
                        kmer.rev = (kmer.rev & ~3ull) | (uint64_t)(3^j);
                        tmp_u16_p = mcounter_ptr(mc, min_u64(kmer.fwd, kmer.rev));
                        if (tmp_u16_p != NULL) {
                            lflags |= 1 << j;
                            next_u16_p = tmp_u16_p;
                        }
                    }
                    if (lflags == 0) {
                        sum += *u16_p;
                        *u16_p = UINT16_MAX;
                        break;
                    }
                    else if (flag_degree(lflags) > 1) {
                        lflags = 0 | (1 << lseq.ele[lseq.size-1]);
                        varr8_pop(&lseq);
                        break;
                    }
                    varr8_push(&lseq, flag_base(lflags));
                    sum += *u16_p;
                    *u16_p = UINT16_MAX;
                    u16_p = next_u16_p;
                }
            }

            varr64_pop(&stack);

            if (flag_degree(lflags))
                varr64_push(&stack, (lkmer.fwd & (mask >> 2)) | ((uint64_t)flag_base(lflags) << (2*(k_len-1))));
            if (flag_degree(rflags))
                varr64_push(&stack, (rkmer.fwd & ~3ull) | (uint64_t)flag_base(rflags));

            mstra.len = lseq.size + rseq.size + 1;
            varr8_resize(&seq, (mstra.len+k_len+2)/4);
            mstra.lkey = mstra.rkey = 0;
            for (j = 0; j < lseq.size; ++j) {
                if (j < k_len)
                    mstra.lkey = (mstra.lkey << 2) | lseq.ele[lseq.size-j-1];
                bseq_set(seq.ele, j, lseq.ele[lseq.size-j-1]);
            }
            for (j = 0; j < k_len; ++j) {
                if (lseq.size+j < k_len)
                    mstra.lkey = (mstra.lkey << 2) | ((st_kmer.fwd >> (2*(k_len-j-1)))&3);
                if (j >= rseq.size)
                    mstra.rkey = (mstra.rkey << 2) | ((st_kmer.fwd >> (2*(k_len-j-1)))&3);
                bseq_set(seq.ele, lseq.size+j, (st_kmer.fwd >> (2*(k_len-j-1)))&3);
            }
            for (j = 0; j < rseq.size; ++j) {
                if (j + k_len >= rseq.size)
                    mstra.rkey = (mstra.rkey << 2) | rseq.ele[j];
                bseq_set(seq.ele, lseq.size+k_len+j, rseq.ele[j]);
            }
            mstra.cov = (uint16_t)((double)sum / (double)(mstra.len) + 0.5);
            mstra.out = (lflags << 4) | rflags;
            fwrite(&mstra, sizeof(mstra_t), 1, mg->stra_file);
            fwrite(seq.ele, sizeof(int8_t), (mstra.len+k_len+2)/4, mg->stra_file);
            ++mg->n_stra;
        }
    }

    varr8_destroy(&rseq);
    varr8_destroy(&lseq);
    varr8_destroy(&seq);
    varr64_destroy(&stack);
}

void mgraph_load(mgraph_t *mg)
{
    uint64_t mem_usage;
    uint64_t k_len;
    uint64_t mask;
    mjunc_t mjunc;
    mstra_t mstra;
    mstra_t *sp;

    fputs("loading kmer graph...\n", stderr);

    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;

    if (mg->jtable != NULL)
        mgraph_destroy_table(mg);
    
    fseek(mg->stra_file, 0, SEEK_END);
    mem_usage = ftell(mg->stra_file);

    mg->jtable_size = 1;
    mg->jindex_len = 0;
    while (mg->jtable_size * MAX_LOAD_FACTOR < mg->n_junc) {
        mg->jtable_size *= 2;
        ++mg->jindex_len;
    }

    mg->stable_size = 1;
    mg->sindex_len = 0;
    while (mg->stable_size * MAX_LOAD_FACTOR < mg->n_stra) {
        mg->stable_size *= 2;
        ++mg->sindex_len;
    }

    mem_usage += mg->jtable_size*sizeof(mjunc_t) + 2*mg->stable_size*sizeof(mstra_t *);
    check_mem_usage(mem_usage, mg->mem_lim);

    while (mem_usage + mg->jtable_size*sizeof(mjunc_t) + 2*mg->stable_size*sizeof(mstra_t *) < mg->mem_lim / 4) {
        mem_usage += mg->jtable_size*sizeof(mjunc_t) + 2*mg->stable_size*sizeof(mstra_t *);
        mg->jtable_size *= 2;
        ++mg->jindex_len;
        mg->stable_size *= 2;
        ++mg->sindex_len;
    }

    fprintf(stderr, "JUNCTION_LOAD_FACTOR=%f, ", (double)(mg->n_junc) / (double)mg->jtable_size);
    fprintf(stderr, "STRAIGHT_LOAD_FACTOR=%f, ", (double)(mg->n_stra) / (double)mg->stable_size);
    fprintf(stderr, "MEM_USAGE=%luB\n", mem_usage); 

    mg->jtable = (mjunc_t *)my_calloc(mg->jtable_size, sizeof(mjunc_t));
    mg->lstable = (mstra_t **)my_calloc(mg->stable_size, sizeof(mstra_t *));
    mg->rstable = (mstra_t **)my_calloc(mg->stable_size, sizeof(mstra_t *));

    rewind(mg->junc_file);
    while (fread(&mjunc, sizeof(mjunc_t), 1, mg->junc_file))
        mgraph_mjunc_insert(mg, &mjunc);

    rewind(mg->stra_file);
    while (fread(&mstra, sizeof(mstra_t), 1, mg->stra_file)) {
        sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(mstra.len+k_len+2)/4);
        *(sp) = mstra;
        fread(sp->seq, sizeof(uint8_t), (mstra.len+k_len+2)/4, mg->stra_file);
        mgraph_mstra_insert(mg, sp);
    }
}

mstra_t *mgraph_join_sj(mgraph_t *mg, mstra_t *sp, mjunc_t *jp)
{
    uint64_t i;
    mstra_t *new_sp;

    new_sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(sp->len+mg->k_len+3)/4);
    new_sp->lkey = sp->lkey;
    new_sp->rkey = jp->key;
    new_sp->len = sp->len + 1;
    new_sp->cov = (double)(sp->cov*(sp->len) + jp->cov) / (double)(new_sp->len) + 0.5;
    new_sp->out = (sp->out & 0xf0) | (jp->out & 0x0f);
    for (i = 0; i < sp->len + mg->k_len - 1; ++i)
        bseq_set(new_sp->seq, i, bseq_get(sp->seq, i));
    bseq_set(new_sp->seq, sp->len + mg->k_len - 1, (jp->key)&3);
    mgraph_mstra_delete(mg, sp);
    mjunc_delete(jp);
    mgraph_mstra_insert(mg, new_sp);
    return new_sp;
}

mstra_t *mgraph_join_js(mgraph_t *mg, mjunc_t *jp, mstra_t *sp)
{
    uint64_t i;
    mstra_t *new_sp;

    new_sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(sp->len+mg->k_len+3)/4);
    new_sp->lkey = jp->key;
    new_sp->rkey = sp->rkey;
    new_sp->len = sp->len + 1;
    new_sp->cov  = (double)(jp->cov + sp->cov*(sp->len)) / (double)(new_sp->len) + 0.5;
    new_sp->out = (jp->out & 0xf0) | (sp->out & 0x0f);
    bseq_set(new_sp->seq, 0, (jp->key >> (2*(mg->k_len - 1)))&3);
    for (i = 0; i < sp->len + mg->k_len - 1; ++i)
        bseq_set(new_sp->seq, i+1, bseq_get(sp->seq, i));
    mgraph_mstra_delete(mg, sp);
    mjunc_delete(jp);
    mgraph_mstra_insert(mg, new_sp);
    return new_sp;
}

mstra_t *mgraph_join_jj(mgraph_t *mg, mjunc_t *ljp, mjunc_t *rjp)
{
    uint64_t i;
    mstra_t *new_sp;

    if (ljp == rjp)
        return mgraph_change_j2s(mg, ljp);
    new_sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(mg->k_len+4)/4);
    new_sp->lkey = ljp->key;
    new_sp->rkey = rjp->key;
    new_sp->len = 2;
    new_sp->cov = (double)(ljp->cov + rjp->cov) / 2.0 + 0.5;
    new_sp->out = (ljp->out & 0xf0) | (rjp->out & 0x0f);
    bseq_set(new_sp->seq, 0, (ljp->key >> (2*(mg->k_len - 1)))&3);
    for (i = 0; i < mg->k_len; ++i)
        bseq_set(new_sp->seq, i+1, (rjp->key >> (2*(mg->k_len+i-1)))&3);
    mjunc_delete(ljp);
    mjunc_delete(rjp);
    mgraph_mstra_insert(mg, new_sp);
    return new_sp;
}

mstra_t *mgraph_change_j2s(mgraph_t *mg, mjunc_t *jp)
{
    uint64_t i;
    mstra_t *new_sp;

    new_sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(mg->k_len+3)/4);
    new_sp->lkey = new_sp->rkey = jp->key;
    new_sp->len = 1;
    new_sp->cov = jp->cov;
    new_sp->out = jp->out;
    for (i = 0; i < mg->k_len; ++i)
        bseq_set(new_sp->seq, i, (jp->key >> (2*(mg->k_len - i - 1))) & 3);
    mjunc_delete(jp);
    mgraph_mstra_insert(mg, new_sp);
    return new_sp;
}

mstra_t *mgraph_join_ss(mgraph_t *mg, mstra_t *lsp, mstra_t *rsp)
{
    uint64_t i;
    mstra_t *new_sp;

    if (lsp == rsp)
        return lsp;
    new_sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(lsp->len+rsp->len+mg->k_len+2)/4);
    new_sp->lkey = lsp->lkey;
    new_sp->rkey = rsp->rkey;
    new_sp->len = lsp->len + rsp->len;
    new_sp->cov = (double)(lsp->cov*(lsp->len) + rsp->cov*(rsp->len)) / (double)(new_sp->len) + 0.5;
    new_sp->out = (lsp->out & 0xf0) | (rsp->out & 0x0f);
    for (i = 0; i < lsp->len + mg->k_len - 1; ++i)
        bseq_set(new_sp->seq, i, bseq_get(lsp->seq, i));
    for (i = 0; i < rsp->len; ++i)
        bseq_set(new_sp->seq, i + lsp->len + mg->k_len - 1, bseq_get(rsp->seq, i+mg->k_len-1));
    mgraph_mstra_delete(mg, lsp);
    mgraph_mstra_delete(mg, rsp);
    mgraph_mstra_insert(mg, new_sp);
    return new_sp;
}

void mgraph_show_contig(const mgraph_t *mg, const double cov_ratio, const char *out_name)
{
    uint64_t i;
    uint64_t j;
    uint64_t n_seq = 0;
    uint64_t *buf;
    uint64_t buf_size;
    mstra_t *sp;
    FILE *out;

    out = fopen(out_name, "w");
    buf_size = mgraph_sorted_prefix_list(mg, &buf);

    n_seq = 0;
    for (i = 0; i < buf_size ; ++i) {
        sp = mgraph_mstra_find_pre(mg, buf[i]);
        if (sp == NULL)
            continue;
        ++n_seq;
        fprintf(out, ">seq%lu_len%lu_cov%u\n", n_seq, sp->len + mg->k_len - 1, (uint16_t)(sp->cov * cov_ratio + 0.5));
        for (j = 0; j < sp->len + mg->k_len - 1; ++j) {
            putc(num2ascii(bseq_get(sp->seq, j)), out);
            if ((j+1)%LINE_LENGTH == 0)
                putc('\n', out);
        }
        if (j % LINE_LENGTH != 0)
            putc('\n', out);
    }

    fclose(out);
    my_free(buf);
}

void mgraph_extract_read(const mgraph_t *mg, FILE **read_fpp)
{
    uint64_t i;
    uint64_t k_len;
    uint64_t k_mask;
    kmer_t kmer;
    seq_t seq;
    seq_t tmp_seq;
    FILE *new_read_fp;

    k_len = mg->k_len;
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

            if (mgraph_mjunc_find(mg, kmer.fwd) != NULL) {
                seq_write(&seq, new_read_fp);
                break;
            }
            else if (mgraph_mjunc_find(mg, kmer.rev) != NULL) {
                seq_write(&seq, new_read_fp);
                break;
            }
        }
    }

    fclose(*read_fpp);
    *read_fpp = new_read_fp;
}

void mgraph_delete_clear_junc(const mgraph_t *mg, const uint64_t len)
{
    uint64_t i;
    uint64_t j;
    uint64_t k_len;
    uint64_t k_mask;
    uint64_t lkey;
    uint64_t rkey;
    mjunc_t *jp;
    mstra_t *sp;

    k_len = mg->k_len;
    k_mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    
    for (i = 0; i < mg->jtable_size; ++i) {
        jp = &(mg->jtable[i]);
        if (mjunc_is_emp(*jp) || mjunc_is_del(*jp))
            continue;
        lkey = (jp->key >> 2);
        rkey = (jp->key << 2) & k_mask;
        for (j = 0; j < 4; ++j) {
            if (jp->out & (1 << (j+4))) {
                lkey = (lkey & (k_mask >> 2)) | (j << (2*(k_len-1)));
                if ((sp = mgraph_mstra_find_suf(mg, lkey)) == NULL || sp->len < len)
                    break;
            }
            if (jp->out & (1 << j)) {
                rkey = (rkey & ~3ull) | j;
                if ((sp = mgraph_mstra_find_pre(mg, rkey)) == NULL || sp->len < len)
                    break;
            }
        }

        if (j == 4)
            mjunc_delete(jp);
    }
}

uint64_t mgraph_cut_branch(mgraph_t *mg)
{
    uint64_t i;
    uint64_t j;
    uint64_t k_len;
    uint64_t mask;
    uint64_t key;
    uint64_t max;
    uint64_t n_delete;
    uint64_t *buf;
    uint64_t buf_size;
    mjunc_t *tmp_jp;
    mjunc_t *jp;
    mstra_t *sp;
    mstra_t *tmp_sp;

    n_delete = 0;
    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    buf_size = mgraph_sorted_prefix_list(mg, &buf);
    for (i = 0; i < buf_size; ++i) {
        sp = mgraph_mstra_find_pre(mg, buf[i]);
        if (sp == NULL || sp->len > k_len)
            continue;

        if (sp->out >> 4) {
            if (sp->out & 0xf)
                continue;
            key = (sp->lkey >> 2) | ((uint64_t)flag_base(sp->out >> 4) << (2*(k_len-1)));
            jp = mgraph_mjunc_find(mg, key);
            if (jp == NULL)
                continue;

            key = (jp->key << 2) & mask;
            max = 0;
            for (j = 0; j < 4; ++j) {
                if (!(jp->out & (1 << j)))
                    continue;
                key = (key & ~3ull) | j;
                tmp_sp = mgraph_mstra_find_pre(mg, key);
                if (tmp_sp != sp && tmp_sp != NULL && tmp_sp->cov > max) {
                    max = tmp_sp->cov;
                    continue;
                }

                tmp_jp = mgraph_mjunc_find(mg, key);
                if (tmp_jp != NULL && tmp_jp->cov > max)
                    max = tmp_jp->cov;
            }
            if (sp->cov > max * mg->branch_th)
                continue;

            jp->out ^= 1 << bseq_get(sp->seq, k_len - 1);
        }
        else if (sp->out & 0xf) {
            key = ((sp->rkey << 2) & mask) | (uint64_t)flag_base(sp->out & 0xf);
            jp = mgraph_mjunc_find(mg, key);
            if (jp == NULL)
                continue;

            key = jp->key >> 2;
            max = 0;
            for (j = 0; j < 4; ++j) {
                if (!(jp->out & (1 << (j+4))))
                    continue;
                key = (key & (mask >> 2)) | (j << (2*(k_len-1)));
                tmp_sp = mgraph_mstra_find_suf(mg, key);
                if (tmp_sp != sp && tmp_sp != NULL && tmp_sp->cov > max) {
                    max = tmp_sp->cov;
                    continue;
                }

                tmp_jp = mgraph_mjunc_find(mg, key);
                if (tmp_jp != NULL && tmp_jp->cov > max)
                    max = tmp_jp->cov;
            }
            if (sp->cov > max * mg->branch_th)
                continue;

            jp->out ^= 1 << (4 + bseq_get(sp->seq, sp->len - 1));
        }
        else 
            continue;

        mgraph_mstra_delete(mg, sp);
        ++n_delete;
    }

    my_free(buf);

    return n_delete;
}

void mgraph_join_nodes(mgraph_t *mg)
{
    uint64_t i;
    uint64_t k_len;
    uint64_t mask;
    uint64_t key;
    uint64_t *buf;
    uint64_t buf_size;
    mjunc_t *jp;
    mjunc_t *ljp;
    mjunc_t *rjp;
    mstra_t *sp;
    mstra_t *lsp;
    mstra_t *rsp;

    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;

    buf_size = mgraph_sorted_junction_list(mg, &buf);
    for (i = 0; i < buf_size; ++i) {
        jp = mgraph_mjunc_find(mg, buf[i]);
        if (jp == NULL || flag_degree(jp->out >> 4) > 1 || flag_degree(jp->out & 0xf) > 1)
            continue;

        if (flag_degree(jp->out >> 4)) {
            key = (jp->key >> 2) | ((uint64_t)flag_base(jp->out >> 4) << (2*(k_len-1)));
            lsp = mgraph_mstra_find_suf(mg, key);
            if (lsp == NULL) {
                ljp = mgraph_mjunc_find(mg, key);
                if (ljp == NULL || flag_degree(ljp->out >> 4) > 1 || flag_degree(ljp->out & 0xf) > 1)
                    ljp = NULL;
            }
        }
        else {
            lsp = NULL;
            ljp = NULL;
        }

        if (lsp != NULL)
            lsp = mgraph_join_sj(mg, lsp, jp);
        else if (ljp != NULL)
            lsp = mgraph_join_jj(mg, ljp, jp);
        else
            lsp = mgraph_change_j2s(mg, jp);

        if (flag_degree(lsp->out & 0xf)) {
            key = ((lsp->rkey << 2) & mask) | (uint64_t)flag_base(lsp->out & 0xf);
            rsp = mgraph_mstra_find_pre(mg, key);
            if (rsp == NULL) {
                rjp = mgraph_mjunc_find(mg, key);
                if (rjp == NULL || flag_degree(rjp->out >> 4) > 1 || flag_degree(rjp->out & 0xf) > 1)
                    rjp = NULL;
            }
        }
        else {
            rsp = NULL;
            rjp = NULL;
        }

        if (rsp != NULL)
            mgraph_join_ss(mg, lsp, rsp);
        else if (rjp != NULL)
            mgraph_join_sj(mg, lsp, rjp);
    }
    my_free(buf);

    buf_size = mgraph_sorted_prefix_list(mg, &buf);
    for (i = 0; i < buf_size ; ++i) {
        sp = mgraph_mstra_find_pre(mg, buf[i]);
        if (sp == NULL)
            continue;

        if (flag_degree(sp->out >> 4)) {
            key = (sp->lkey >> 2) | ((uint64_t)flag_base(sp->out >> 4) << (2*(k_len-1)));
            lsp = mgraph_mstra_find_suf(mg, key);
            if (lsp != NULL)
                sp = mgraph_join_ss(mg, lsp, sp);
        }

        if (flag_degree(sp->out & 0xf)) {
            key = ((sp->rkey << 2) & mask) | (uint64_t)flag_base(sp->out & 0xf);
            rsp = mgraph_mstra_find_pre(mg, key);
            if (rsp != NULL)
                mgraph_join_ss(mg, sp, rsp);
        }
    }
    my_free(buf);
}

uint64_t mgraph_crush_bubble(mgraph_t *mg, const double ave_cov, FILE **bubble_fpp)
{
    uint64_t n_crush;
    uint64_t i;
    uint64_t j;
    uint64_t k;
    uint64_t m;
    uint64_t n;
    uint64_t k_len;
    uint64_t mask;
    uint64_t key;
    uint64_t tmp64;
    uint64_t al_th;
    uint64_t *work;
    uint64_t wsize;
    uint64_t al_len1;
    uint64_t al_len2;
    uint64_t clust[4];
    uint64_t *buf;
    uint64_t buf_size;
    uint16_t cov_th;
    mjunc_t *rt_jp;
    mjunc_t *jp[4];
    mstra_t *sp[4];
    mstra_t *max_sp[3];

    n_crush = 0;

    if ((uint64_t)(ave_cov * 2.0 + 0.5) < UINT16_MAX)
        cov_th = (uint16_t)(ave_cov * 2.0 + 0.5);
    else
        cov_th = UINT16_MAX;

    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    wsize = k_len;
    work = (uint64_t *)my_malloc(wsize * 2 * sizeof(uint64_t));;
    buf_size = mgraph_sorted_junction_list(mg, &buf);

    for (i = 0; i < buf_size ; ++i) {
        rt_jp = mgraph_mjunc_find(mg, buf[i]);
        if (rt_jp == NULL || flag_degree(rt_jp->out & 0x0f) < 2)
            continue;

        key = (rt_jp->key << 2) & mask;
        for (j = 0; j < 4; ++j) {
            clust[j] = j;
            sp[j] = NULL;
            jp[j] = NULL;
            if (!(rt_jp->out & (1 << j)))
                continue;
            key = (key & ~3ull) | j;
            sp[j] = mgraph_mstra_find_pre(mg, key);
            if (sp[j] != NULL && (sp[j]->out & 0x0f))
                jp[j] = mgraph_mjunc_find(mg, ((sp[j]->rkey<<2) & mask) | flag_base(sp[j]->out&0x0f));
        }

        for (j = 0; j < 3; ++j) {
            if (sp[j] == NULL)
                continue;
            for (k = j+1; k < 4; ++k) {
                if (sp[k] == NULL || jp[j] == NULL || jp[j] != jp[k] || sp[j]->cov + sp[k]->cov > cov_th)
                    continue;
                al_th = sp[j]->len > sp[k]->len ? (sp[j]->len+k_len-1)*mg->bubble_th+0.5 : (sp[k]->len+k_len-1)*mg->bubble_th+0.5;
                if (sp[j]->len +1 <= k_len || sp[k]->len + 1 <= k_len) {
                    if (abs64((int64_t)sp[j]->len - (int64_t)sp[k]->len) > al_th)
                        continue;

                    if (clust[k] == k)
                        clust[k] = clust[j];
                    else
                        clust[j] = clust[k];
                    continue;
                }
                al_len1 = sp[j]->len - k_len + 1;
                al_len2 = sp[k]->len - k_len + 1;
                if (al_len2 + 1 > wsize) {
                    wsize = al_len2 + 1;
                    work = (uint64_t *)my_realloc(work, wsize*2*sizeof(uint64_t));
                }
                for (n = 0; n < al_len2+1; ++n)
                    work[n] = n;
                for (m = 0; m < al_len1; ++m) {
                    work[(m&1)*wsize] = m;
                    work[(~m&1)*wsize] = m+1;
                    for (n = 0; n < al_len2; ++n) {
                        if (bseq_get(sp[j]->seq, m+k_len-1) == bseq_get(sp[k]->seq, n+k_len-1)) {
                            work[(~m&1)*wsize + n+1] = work[(m&1)*wsize + n];
                            continue;
                        }
                        work[(~m&1)*wsize + n+1] = work[(m&1)*wsize + n] + 1;
                        if (work[(~m&1)*wsize + n+1] > work[(m&1)*wsize + n+1] + 1)
                            work[(~m&1)*wsize + n+1] = work[(m&1)*wsize + n+1] + 1;
                        if (work[(~m&1)*wsize + n+1] > work[(~m&1)*wsize + n] + 1)
                            work[(~m&1)*wsize + n+1] = work[(~m&1)*wsize + n] + 1;
                    }
                }
                if (work[(m&1)*wsize + n] > al_th)
                    continue;

                if (clust[k] == k)
                    clust[k] = clust[j];
                else
                    clust[j] = clust[k];
            }
        }

        for (j = 0; j < 3; ++j) {
            max_sp[j] = NULL;
            for (k = j; k < 4; ++k) {
                if (clust[k] != j || sp[k] == NULL)
                    continue;
                if (max_sp[j] == NULL || sp[k]->cov > max_sp[j]->cov)
                    max_sp[j] = sp[k];
            }
        }

        for (j = 0; j < 3; ++j) {
            if (max_sp[j] == NULL)
                continue;
            for (k = j; k < 4; ++k) {
                if (clust[k] != j || sp[k] == NULL || sp[k] == max_sp[j])
                    continue;
                tmp64 = max_sp[j]->cov*(max_sp[j]->len) + sp[k]->cov*(sp[k]->len);
                tmp64 = (double)tmp64 / (double)(max_sp[j]->len) + 0.5;
                max_sp[j]->cov = tmp64 < UINT16_MAX ? tmp64 : UINT16_MAX-1;
                rt_jp->out ^= 1 << k;
                jp[k]->out ^= 1 << (((sp[k]->rkey >> 2*(k_len-1)) & 3)+4);

                if (bubble_fpp != NULL) {
                    tmp64 = sp[k]->len + k_len - 1;
                    fwrite(&tmp64, sizeof(uint64_t), 1, *bubble_fpp);
                    fwrite(&(sp[k]->cov), sizeof(uint16_t), 1, *bubble_fpp);
                    fwrite(sp[k]->seq, sizeof(uint8_t), (tmp64 + 3)/4, *bubble_fpp);
                }

                mgraph_mstra_delete(mg, sp[k]);
                ++n_crush;
            }
        }
    }

    my_free(buf);
    my_free(work);

    return n_crush;
}

FILE *mgraph_merge(mgraph_t *mg, mgraph_t *ext_mg)
{
    uint64_t i;
    uint64_t j;
    uint64_t k_len;
    uint64_t mask;
    uint64_t len;
    uint8_t out;
    FILE *contig_fp;
    varr64_t stack;
    kmer_t kmer;
    uint8_t tmp_seq[MAX_READ_LEN] = {0};
    mstra_t *sp;
    mjunc_t *jp;

    fprintf(stderr, "merging kmer graphs...\n");

    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    varr64_init(&stack);

    contig_fp = mgraph_save_contig(mg, mg->k_len);
    fseek(contig_fp, 0, SEEK_END);
    for (i = 0; i < mg->jtable_size; ++i) {
        jp = &(mg->jtable[i]);
        if (mjunc_is_emp(*jp) || mjunc_is_del(*jp))
            continue;
        for (j = 0; j < k_len; ++j)
            bseq_set(tmp_seq, j, (jp->key >> (2 * (k_len - j - 1))) & 3ull);
        fwrite(&k_len, sizeof(uint64_t), 1, contig_fp);
        fwrite(&(jp->cov), sizeof(uint16_t), 1, contig_fp);
        fwrite(tmp_seq, sizeof(uint8_t), (k_len + 3)/4, contig_fp);
    }

    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp))
            continue;
        if (!(sp->out >> 4)) {
            kmer.fwd = sp->lkey >> 2;
            for (j = 0; j < 4; ++j)
                varr64_push(&stack, (kmer.fwd & (mask >> 2)) | (j << (2*(k_len-1))));
        }
        if (!(sp->out & 0xf)) {
            kmer.fwd = (sp->rkey << 2) & mask;
            for (j = 0; j < 4; ++j)
                varr64_push(&stack, (kmer.fwd & ~3ull) | j);
        }

        while (stack.size > 0) {
            kmer.fwd = varr64_peek(&stack);
            varr64_pop(&stack);
            kmer.rev = 0ull;
            for (j = 0; j < k_len; ++j)
                kmer.rev = (kmer.rev << 2) | (((kmer.fwd >> (2*j))&3)^3);

            if ((sp = mgraph_mstra_find_pre(ext_mg, kmer.fwd)) != NULL ||
                (sp = mgraph_mstra_find_pre(ext_mg, kmer.rev)) != NULL ||
                (sp = mgraph_mstra_find_suf(ext_mg, kmer.fwd)) != NULL ||
                (sp = mgraph_mstra_find_suf(ext_mg, kmer.rev)) != NULL) {

                if (sp->cov == UINT16_MAX)
                    continue;
                len = sp->len + k_len - 1;
                fwrite(&len, sizeof(uint64_t), 1, contig_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, contig_fp);
                fwrite(sp->seq, sizeof(uint8_t), (len + 3)/4, contig_fp);
                sp->cov = UINT16_MAX;
                out = sp->out;
                kmer.fwd = sp->lkey >> 2;
                kmer.rev = (sp->rkey << 2) & mask;
            }
            else if ((jp = mgraph_mjunc_find(ext_mg, kmer.fwd)) != NULL ||
                (jp = mgraph_mjunc_find(ext_mg, kmer.rev)) != NULL ||
                (jp = mgraph_mjunc_find(ext_mg, kmer.fwd)) != NULL ||
                (jp = mgraph_mjunc_find(ext_mg, kmer.rev)) != NULL) {

                if (jp->cov == UINT16_MAX)
                    continue;
                for (j = 0; j < k_len; ++j)
                    bseq_set(tmp_seq, j, (kmer.fwd >> (2 * (k_len - j - 1))) & 3ull);
                fwrite(&k_len, sizeof(uint64_t), 1, contig_fp);
                fwrite(&(jp->cov), sizeof(uint16_t), 1, contig_fp);
                fwrite(tmp_seq, sizeof(uint8_t), (k_len + 3)/4, contig_fp);
                jp->cov = UINT16_MAX;
                out = jp->out;
                kmer.fwd = jp->key >> 2;
                kmer.rev = (jp->key << 2) & mask;
            }
            else
                continue;

            for (j = 0; j < 4; ++j) {
                if (out & (1 << (j+4))) {
                    kmer.fwd = (kmer.fwd & (mask >> 2)) | (j << (2*(k_len-1)));
                    varr64_push(&stack, kmer.fwd);
                }

                if (out & (1 << j)) {
                    kmer.rev = (kmer.rev & ~3ull) | j;
                    varr64_push(&stack, kmer.rev);
                }
            }
        }
    }

    varr64_destroy(&stack);
    mgraph_destroy_table(ext_mg);
    mgraph_destroy_table(mg);

    return contig_fp;
}

FILE *mgraph_save_contig(const mgraph_t *mg, const uint64_t next_k)
{
    uint64_t i;
    uint64_t j;
    uint64_t len;
    uint64_t ex_len;
    uint64_t dif;
    uint64_t key;
    uint64_t k_mask;
    FILE *contig_fp;
    mstra_t *sp;
    mstra_t *tmp_sp;
    mjunc_t *jp;
    varr8_t seq;

    k_mask = mg->k_len < 32 ? ~(~0ull << (2*mg->k_len)) : ~0ull;
    dif = next_k - mg->k_len;
    varr8_init(&seq);
    contig_fp = tmpfile_open();
    check_tmpfile(contig_fp);

    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp))
            continue;

        len = 0;
        if (sp->out >> 4) {
            key = (sp->lkey >> 2) | (flag_base(sp->out >> 4) << (2*(mg->k_len - 1)));
            jp = mgraph_mjunc_find(mg, key);
            if (jp != NULL) {
                if (flag_degree(jp->out >> 4) == 1) {
                    key = (key >> 2) | (flag_base(jp->out >> 4) << (2*(mg->k_len - 1)));
                    tmp_sp = mgraph_mstra_find_suf(mg, key);
                    if (tmp_sp != NULL) {
                        ex_len = min_u64(tmp_sp->len, dif);
                        varr8_resize(&seq, (len + ex_len + 3) / 4);
                        for (j = 0; j < ex_len; ++j) {
                            bseq_set(seq.ele, j, bseq_get(tmp_sp->seq, tmp_sp->len - ex_len + j));
                        }
                        len += ex_len;
                    }
                }
                varr8_resize(&seq, len / 4 + 1);
                bseq_set(seq.ele, len, jp->key >> (2*(mg->k_len - 1)));
                ++len;
            }
        }

        varr8_resize(&seq, (len + sp->len + mg->k_len + 2) / 4);
        for (j = 0; j < sp->len + mg->k_len - 1; ++j)
            bseq_set(seq.ele, len + j, bseq_get(sp->seq, j));
        len += sp->len + mg->k_len - 1;

        if (sp->out & 0xf) {
            key = ((sp->rkey << 2) & k_mask) | flag_base(sp->out & 0xf);
            jp = mgraph_mjunc_find(mg, key);
            if (jp != NULL) {
                varr8_resize(&seq, len / 4 + 1);
                bseq_set(seq.ele, len, jp->key & 3ull);
                ++len;
                if (flag_degree(jp->out & 0xf) == 1) {
                    key = ((key << 2) & k_mask) | flag_base(jp->out & 0xf);
                    tmp_sp = mgraph_mstra_find_pre(mg, key);
                    if (tmp_sp != NULL) {
                        ex_len = min_u64(tmp_sp->len, dif);
                        varr8_resize(&seq, (len + ex_len + 3) / 4);
                        for (j = 0; j < ex_len; ++j) {
                            bseq_set(seq.ele, len + j, bseq_get(tmp_sp->seq, mg->k_len - 1 + j));
                        }
                        len += ex_len;
                    }
                }
            }
        }

        fwrite(&len, sizeof(uint64_t), 1, contig_fp);
        fwrite(&(sp->cov), sizeof(uint16_t), 1, contig_fp);
        fwrite(seq.ele, sizeof(uint8_t), (len + 3)/4, contig_fp);
    }

    varr8_destroy(&seq);

    return contig_fp;
}

void mgraph_cut_branch_iterative(mgraph_t *mg)
{
    uint64_t n_delete;
    uint64_t total_delete;

    fputs("removing branches...\n", stderr);
    fprintf(stderr, "BRANCH_DELETE_THRESHOLD=%f\n", mg->branch_th);
    total_delete = 0;
    do {
        n_delete = mgraph_cut_branch(mg);
        fprintf(stderr, "NUM_CUT=%lu\n", n_delete);
        mgraph_join_nodes(mg);
        total_delete += n_delete;
    } while (n_delete > 0);
    fprintf(stderr, "TOTAL_NUM_CUT=%lu\n", total_delete);
}

FILE *mgraph_crush_bubble_iterative(mgraph_t *mg, const double ave_cov)
{
    uint64_t n_crush;
    uint64_t total_crush;
    FILE *bubble_fp;

    fputs("removing bubbles...\n", stderr);
    fprintf(stderr, "BUBBLE_IDENTITY_THRESHOLD=%f\n", mg->bubble_th);

    bubble_fp = tmpfile_open();
    check_tmpfile(bubble_fp);

    total_crush = 0;
    do {
        n_crush = mgraph_crush_bubble(mg, ave_cov, &bubble_fp);
        fprintf(stderr, "NUM_CRUSH=%lu\n", n_crush);
        mgraph_join_nodes(mg);
        total_crush += n_crush;
    } while (n_crush > 0);
    fprintf(stderr, "TOAL_NUM_CRUSH=%lu\n", total_crush);

    return bubble_fp;
}

uint64_t mgraph_sorted_prefix_list(const mgraph_t *mg, uint64_t **buf)
{
    uint64_t i;
    uint64_t buf_size;
    mstra_t *sp;

    buf_size = 0;
    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp))
            continue;
        ++buf_size;
    }
    *buf = (uint64_t *)my_malloc(buf_size * sizeof(uint64_t));
    buf_size = 0;
    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp))
            continue;
        (*buf)[buf_size] = sp->lkey;
        ++buf_size;
    }
    qsort(*buf, buf_size, sizeof(uint64_t), cmp_kmer_uniq);
    
    return buf_size;
}

uint64_t mgraph_sorted_junction_list(const mgraph_t *mg, uint64_t **buf)
{
    uint64_t i;
    uint64_t buf_size;
    mjunc_t *jp;

    buf_size = 0;
    for (i = 0; i < mg->jtable_size; ++i) {
        jp = &(mg->jtable[i]);
        if (mjunc_is_emp(*jp) || mjunc_is_del(*jp))
            continue;
        ++buf_size;
    }
    *buf = (uint64_t *)my_malloc(buf_size * sizeof(uint64_t));

    buf_size = 0;
    for (i = 0; i < mg->jtable_size; ++i) {
        jp = &(mg->jtable[i]);
        if (mjunc_is_emp(*jp) || mjunc_is_del(*jp))
            continue;
        (*buf)[buf_size] = jp->key;
        ++buf_size;
    }
    qsort(*buf, buf_size, sizeof(uint64_t), cmp_kmer_uniq);

    return buf_size;
}

void mgraph_save_edge_kmer(const mgraph_t *mg, mcounter_t *mc, const uint64_t next_k)
{
    uint64_t i;
    uint64_t j;
    uint64_t k_len;
    uint64_t mask;
    uint64_t dif;
    kmer_t lkmer;
    kmer_t rkmer;
    mstra_t *sp;

    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2 * k_len)) : ~0ull;
    dif = next_k - mg->k_len;
    if (mc->kmer_fp != NULL)
        fclose(mc->kmer_fp);
    mc->kmer_fp = tmpfile_open();
    check_tmpfile(mc->kmer_fp);

    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp))
            continue;

        if (sp->len < dif * 2) {
            lkmer.fwd = lkmer.rev = 0;
            for (j = 0; j < k_len - 1; ++j) {
                lkmer.fwd = (lkmer.fwd << 2) | (uint64_t)bseq_get(sp->seq, j);
                lkmer.rev = (lkmer.rev | (uint64_t)(3^bseq_get(sp->seq, k_len - j - 2))) << 2;
            }
            for (j = 0; j < sp->len; ++j) {
                lkmer.fwd = ((lkmer.fwd << 2) & mask) | (uint64_t)bseq_get(sp->seq, j + k_len - 1);
                lkmer.rev = (lkmer.rev >> 2) | ((uint64_t)(3^bseq_get(sp->seq, j + k_len - 1)) << (2*(k_len - 1)));
                if (lkmer.fwd < lkmer.rev)
                    fwrite(&(lkmer.fwd), sizeof(uint64_t), 1, mc->kmer_fp);
                else
                    fwrite(&(lkmer.rev), sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, mc->kmer_fp);
            }
        }
        else {
            lkmer.fwd = lkmer.rev = rkmer.fwd = rkmer.rev = 0;
            for (j = 0; j < k_len - 1; ++j) {
                lkmer.fwd = (lkmer.fwd << 2) | (uint64_t)bseq_get(sp->seq, j);
                lkmer.rev = (lkmer.rev | (uint64_t)(3^bseq_get(sp->seq, k_len - j - 2))) << 2;
                rkmer.fwd = (rkmer.fwd << 2) | (uint64_t)bseq_get(sp->seq, sp->len - dif + j);
                rkmer.rev = (rkmer.rev | (uint64_t)(3^bseq_get(sp->seq, sp->len - dif + k_len - j - 2))) << 2;
            }
            for (j = 0; j < dif; ++j) {
                lkmer.fwd = ((lkmer.fwd << 2) & mask) | (uint64_t)bseq_get(sp->seq, j + k_len - 1);
                lkmer.rev = (lkmer.rev >> 2) | ((uint64_t)(3^bseq_get(sp->seq, j + k_len - 1)) << (2*(k_len - 1)));
                if (lkmer.fwd < lkmer.rev)
                    fwrite(&(lkmer.fwd), sizeof(uint64_t), 1, mc->kmer_fp);
                else
                    fwrite(&(lkmer.rev), sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, mc->kmer_fp);
                rkmer.fwd = ((rkmer.fwd << 2) & mask) | (uint64_t)bseq_get(sp->seq, j + sp->len - dif + k_len - 1);
                rkmer.rev = (rkmer.rev >> 2) | ((uint64_t)(3^bseq_get(sp->seq, j + sp->len - dif + k_len - 1)) << (2*(k_len - 1)));
                if (rkmer.fwd < rkmer.rev)
                    fwrite(&(rkmer.fwd), sizeof(uint64_t), 1, mc->kmer_fp);
                else
                    fwrite(&(rkmer.rev), sizeof(uint64_t), 1, mc->kmer_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, mc->kmer_fp);
            }
        }
    }
}

void mgraph_make_graph_quick(mgraph_t *mg, const mcounter_t *mc, const uint64_t min_cov)
{
    uint64_t i;
    uint64_t j;
    uint64_t k_len;
    uint64_t sum;
    uint64_t mask;
    uint16_t *u16_p;
    uint16_t *tmp_u16_p;
    uint16_t *next_u16_p = NULL;
    uint8_t rflags;
    uint8_t lflags;
    uint8_t tmp_flags;
    uint64_t kmer;
    uint64_t st_kmer;
    uint64_t rkmer;
    uint64_t lkmer;
    varr8_t rseq;
    varr8_t lseq;
    varr8_t seq;
    varr64_t stack;
    mjunc_t mjunc;
    mstra_t mstra;
    mstra_t *sp;

    k_len = mg->k_len = mc->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    varr8_init(&rseq);
    varr8_init(&lseq);
    varr8_init(&seq);
    varr64_init(&stack);

    for (i = j =  0; i < mc->table_size; ++i) {
        if (mc->val_table[i] >= min_cov)
            ++j;
    }

    if ((double)mg->jtable_size * MAX_LOAD_FACTOR < j) {
        mg->jtable_size = mg->stable_size = 2;
        mg->jindex_len = mg->sindex_len = 1;
        while ((double)mg->jtable_size * MAX_LOAD_FACTOR < j) {
            mg->jtable_size *= 2;
            mg->stable_size *= 2;
            ++mg->jindex_len;
            ++mg->sindex_len;
        }
        mg->jtable = (mjunc_t *)my_realloc(mg->jtable, mg->jtable_size * sizeof(mjunc_t));
        mg->lstable = (mstra_t **)my_realloc(mg->lstable, mg->stable_size * sizeof(mstra_t *));
        mg->rstable = (mstra_t **)my_realloc(mg->rstable, mg->stable_size * sizeof(mstra_t *));
    }

    memset(mg->jtable, 0, mg->jtable_size * sizeof(mjunc_t));
    memset(mg->lstable, 0, mg->stable_size * sizeof(mstra_t *));
    memset(mg->rstable, 0, mg->stable_size * sizeof(mstra_t *));

    for (i = 0; i < mc->table_size; ++i) {
        if (mc->val_table[i] < min_cov)
            continue;
        varr64_push(&stack, mc->key_table[i]);
        while (stack.size > 0) {
            st_kmer = varr64_peek(&stack);    
            u16_p = mcounter_ptr(mc, st_kmer);
            if (*u16_p == UINT16_MAX) {
                varr64_pop(&stack);
                continue;
            }

            sum = *u16_p;
            *u16_p = UINT16_MAX;

            lkmer = st_kmer >> 2;
            varr8_clear(&lseq);
            lflags = 0;
            for (j = 0; j < 4; ++j) {
                lkmer = (lkmer & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1)));
                if (mcounter_occ(mc, lkmer) >= min_cov)
                    lflags |= 1 << j;
            }
            if (lflags != 0)
                varr8_push(&lseq, flag_base(lflags));

            rkmer = (st_kmer << 2) & mask;
            varr8_clear(&rseq);
            rflags = 0;
            for (j = 0; j < 4; ++j) {
                rkmer = (rkmer & ~3ull) | (uint64_t)j;
                if (mcounter_occ(mc, rkmer) >= min_cov)
                    rflags |= 1 << j;
            }
            if (rflags != 0)
                varr8_push(&rseq, flag_base(rflags));

            if ((lseq.size > 0 && lseq.ele[0] == 4) || (rseq.size > 0 && rseq.ele[0] == 4)) {
                varr64_pop(&stack);
                for (j = 0; j < 4; ++j) {
                    if (rflags & (1 << j))
                        varr64_push(&stack, (rkmer & ~3ull) | (uint64_t)j);
                    if (lflags & (1 << j)) 
                        varr64_push(&stack, (lkmer & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1))));
                }
                mjunc.key = st_kmer;
                mjunc.cov = sum;
                mjunc.out = (lflags << 4) | rflags;
                mgraph_mjunc_insert(mg, &mjunc);
                continue;
            }

            if (rseq.size > 0) {
                kmer = ((st_kmer << 2) & mask) | (uint64_t)(rseq.ele[0]);
                u16_p = mcounter_ptr(mc, kmer);
                while (1) {
                    rkmer = kmer;
                    kmer = kmer >> 2;
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer = (kmer & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1)));
                        if (mcounter_occ(mc, kmer) >= min_cov)
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        rflags = 0 | (1 << rseq.ele[rseq.size-1]);
                        varr8_pop(&rseq);
                        break;
                    }
                    kmer = ((kmer << 4) & mask) | (uint64_t)(rseq.ele[rseq.size-1] << 2);
                    rflags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer = (kmer & ~3ull) | (uint64_t)j;
                        tmp_u16_p = mcounter_ptr(mc, kmer);
                        if (tmp_u16_p != NULL && *tmp_u16_p >= min_cov) {
                            rflags |= 1 << j;
                            next_u16_p = tmp_u16_p;
                        }
                    }
                    if (rflags == 0 ) {
                        sum += *u16_p;
                        *u16_p = UINT16_MAX;
                        break;
                    }
                    else if (flag_degree(rflags) > 1) {
                        rflags = 0 | (1 << rseq.ele[rseq.size-1]);
                        varr8_pop(&rseq);
                        break;
                    }
                    varr8_push(&rseq, flag_base(rflags));
                    sum += *u16_p;
                    *u16_p = UINT16_MAX;
                    u16_p = next_u16_p;
                }
            }

            if (lseq.size > 0) {
                kmer = (st_kmer >> 2) | ((uint64_t)(lseq.ele[0]) << (2*(k_len-1)));
                u16_p = mcounter_ptr(mc, kmer);
                while (1) {
                    lkmer = kmer;
                    kmer = (kmer << 2) & mask;
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer = (kmer & ~3ull) | (uint64_t)j;
                        if (mcounter_occ(mc, kmer) >= min_cov)
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        lflags = 0 | (1 << lseq.ele[lseq.size-1]);
                        varr8_pop(&lseq);
                        break;
                    }
                    kmer = (kmer >> 4) | ((uint64_t)(lseq.ele[lseq.size-1]) << (2*(k_len-2)));
                    lflags = 0;
                    for (j = 0; j < 4; ++j) {
                        kmer = (kmer & (mask >> 2)) | ((uint64_t)j << (2*(k_len-1)));
                        tmp_u16_p = mcounter_ptr(mc, kmer);
                        if (tmp_u16_p != NULL && *tmp_u16_p >= min_cov) {
                            lflags |= 1 << j;
                            next_u16_p = tmp_u16_p;
                        }
                    }
                    if (lflags == 0) {
                        sum += *u16_p;
                        *u16_p = UINT16_MAX;
                        break;
                    }
                    else if (flag_degree(lflags) > 1) {
                        lflags = 0 | (1 << lseq.ele[lseq.size-1]);
                        varr8_pop(&lseq);
                        break;
                    }
                    varr8_push(&lseq, flag_base(lflags));
                    sum += *u16_p;
                    *u16_p = UINT16_MAX;
                    u16_p = next_u16_p;
                }
            }

            varr64_pop(&stack);

            if (flag_degree(lflags))
                varr64_push(&stack, (lkmer & (mask >> 2)) | ((uint64_t)flag_base(lflags) << (2*(k_len-1))));
            if (flag_degree(rflags))
                varr64_push(&stack, (rkmer & ~3ull) | (uint64_t)flag_base(rflags));

            mstra.len = lseq.size + rseq.size + 1;
            varr8_resize(&seq, (mstra.len+k_len+2)/4);
            mstra.lkey = mstra.rkey = 0;
            for (j = 0; j < lseq.size; ++j) {
                if (j < k_len)
                    mstra.lkey = (mstra.lkey << 2) | lseq.ele[lseq.size-j-1];
                bseq_set(seq.ele, j, lseq.ele[lseq.size-j-1]);
            }
            for (j = 0; j < k_len; ++j) {
                if (lseq.size+j < k_len)
                    mstra.lkey = (mstra.lkey << 2) | ((st_kmer >> (2*(k_len-j-1)))&3);
                if (j >= rseq.size)
                    mstra.rkey = (mstra.rkey << 2) | ((st_kmer >> (2*(k_len-j-1)))&3);
                bseq_set(seq.ele, lseq.size+j, (st_kmer >> (2*(k_len-j-1)))&3);
            }
            for (j = 0; j < rseq.size; ++j) {
                if (j + k_len >= rseq.size)
                    mstra.rkey = (mstra.rkey << 2) | rseq.ele[j];
                bseq_set(seq.ele, lseq.size+k_len+j, rseq.ele[j]);
            }
            mstra.cov = (uint16_t)((double)sum / (double)(mstra.len) + 0.5);
            mstra.out = (lflags << 4) | rflags;

            sp = (mstra_t *)my_malloc(sizeof(mstra_t) + sizeof(uint8_t)*(mstra.len+k_len+2)/4);
            *(sp) = mstra;
            memcpy(sp->seq, seq.ele, (mstra.len+k_len+2)/4 * sizeof(uint8_t));
            mgraph_mstra_insert(mg, sp);
        }
    }

    varr8_destroy(&rseq);
    varr8_destroy(&lseq);
    varr8_destroy(&seq);
    varr64_destroy(&stack);
}

void mgraph_straight_stats(const mgraph_t *mg, uint64_t *len_cutoff, uint64_t *cov_cutoff, double *ave_cov)
{
    uint64_t i;
    uint64_t tmp;
    uint64_t sum;
    uint64_t num;
    mstra_t *sp;
    varr64_t len_distr;
    varr64_t cov_distr;

    varr64_init(&len_distr);
    varr64_init(&cov_distr);
    varr64_resize(&cov_distr, UINT16_MAX);
    memset(cov_distr.ele, 0, cov_distr.size * sizeof(uint64_t));

    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp))
            continue;

        if (sp->len + 1 > len_distr.size) {
            tmp = len_distr.size;
            varr64_resize(&len_distr, sp->len + 1);
            memset(&(len_distr.ele[tmp]), 0, (len_distr.size - tmp) * sizeof(uint64_t));
        }
        ++len_distr.ele[sp->len];
        cov_distr.ele[sp->cov] += sp->len;
    }

    *cov_cutoff = get_left_minimal(cov_distr.ele, cov_distr.size - 1);
    *len_cutoff = get_left_minimal(len_distr.ele, len_distr.size - 1);

    sum = num = 0;
    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp) || (sp->cov < *cov_cutoff && sp->len < *len_cutoff))
            continue;
        sum += sp->cov * sp->len;
        num += sp->len;
    }

    *ave_cov = (double)sum / num;

    varr64_destroy(&len_distr);
    varr64_destroy(&cov_distr);
}

void mgraph_delete_erroneous_straight_node(mgraph_t *mg, const uint64_t len_cutoff, const uint64_t cov_cutoff)
{
    uint64_t i;
    uint64_t n_delete;
    uint64_t k_len;
    uint64_t mask;
    uint64_t key;
    mstra_t *sp;
    mjunc_t *jp;

    fputs("removing erroneous nodes...\n", stderr);

    k_len = mg->k_len;
    mask = k_len < 32 ? ~(~0ull << (2*k_len)) : ~0ull;
    n_delete = 0;

    for (i = 0; i < mg->stable_size; ++i) {
        sp = mg->lstable[i];
        if (mstra_is_emp(sp) || mstra_is_del(sp) || sp->cov >= cov_cutoff || sp->len >= len_cutoff)
            continue;

        if (sp->out << 4) {
            key = (sp->lkey >> 2) | ((uint64_t)flag_base(sp->out >> 4) << (2*(k_len-1)));
            jp = mgraph_mjunc_find(mg, key);
            if (jp != NULL)
                jp->out ^= 1 << bseq_get(sp->seq, k_len - 1);
        }

        if (sp->out & 0xf) {
            key = ((sp->rkey << 2) & mask) | (uint64_t)flag_base(sp->out & 0xf);
            jp = mgraph_mjunc_find(mg, key);
            if (jp != NULL)
                jp->out ^= 1 << (4 + bseq_get(sp->seq, sp->len - 1));
        }

        mgraph_mstra_delete(mg, sp);
        ++n_delete;
    }

    mgraph_join_nodes(mg);

    fprintf(stderr, "NUM_DELETE=%lu\n", n_delete);
}
