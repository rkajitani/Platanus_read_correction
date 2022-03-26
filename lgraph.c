#include "lgraph.h"

lstra_t g_lstra_deleted;

void lgraph_init(lgraph_t *lg, const double branch_th, const double bubble_th)
{
    memset(lg, 0, sizeof(lgraph_t));
    lg->bubble_th = bubble_th;
    lg->branch_th = branch_th;
}

void lgraph_destroy_table(lgraph_t *lg)
{
    uint64_t i;

    if (lg->jval_table != NULL) {
        my_free(lg->jval_table);
        my_free(lg->jkey_table);
        lg->jtable_size = 0;
    }
    lg->jval_table = NULL;

    if (lg->lsval_table != NULL) {
        for (i = 0; i < lg->stable_size; ++i) {
            if (!lstra_is_emp(lg->lsval_table[i]) && !lstra_is_del(lg->lsval_table[i]))
                my_free(lg->lsval_table[i]);
        }
        my_free(lg->lsval_table);
        my_free(lg->lskey_table);
        my_free(lg->rsval_table);
        my_free(lg->rskey_table);
        lg->stable_size = 0;
    }
    lg->lsval_table = NULL;
    lg->rsval_table = NULL;
}

void lgraph_destroy(lgraph_t *lg)
{
    lgraph_destroy_table(lg);
    if (lg->junc_file != NULL)
        fclose(lg->junc_file);
    if (lg->stra_file != NULL)
        fclose(lg->stra_file);
    memset(lg, 0, sizeof(lgraph_t));
}

void lgraph_save(lgraph_t *lg, const lcounter_t *lc, FILE *sorted_key_fp)
{
    uint64_t j;
    uint64_t l_len;
    uint64_t bl;
    uint64_t sum;
    lmer_t fwd;
    lmer_t rev;
    lmer_t st_fwd;
    lmer_t st_rev;
    lmer_t lfwd;
    lmer_t lrev;
    lmer_t rfwd;
    lmer_t rrev;
    uint16_t *u16_p;
    uint16_t *tmp_u16_p;
    uint16_t *next_u16_p = NULL;
    uint8_t rflags;
    uint8_t lflags;
    uint8_t tmp_flags;
    varr8_t rseq;
    varr8_t lseq;
    varr8_t seq;
    varr64_t stack;
    ljunc_t ljunc;
    lstra_t lstra;

    lg->mem_lim = lc->mem_lim;
    l_len = lg->l_len = lc->l_len;
    bl = lg->bl = lc->bl;
    lg->n_junc = lg->n_stra = 0;
    varr8_init(&rseq);
    varr8_init(&lseq);
    varr8_init(&seq);
    varr64_init(&stack);

    fputs("connecting kmers...\n", stderr);

    if (lg->junc_file != NULL)
        fclose(lg->junc_file);
    lg->junc_file = tmpfile_open();
    check_tmpfile(lg->junc_file);

    if (lg->stra_file != NULL)
        fclose(lg->stra_file);
    lg->stra_file = tmpfile_open();
    check_tmpfile(lg->stra_file);

    varr64_resize(&stack, bl);
    varr64_resize(&stack, 0);
    rewind(sorted_key_fp);
    while (fread(stack.ele, sizeof(uint64_t), bl, sorted_key_fp)) {
        stack.size += bl;

        while (stack.size > 0) {
            lmer_cpy(st_fwd, &(stack.ele[stack.size - bl]), l_len);    
            for (j = 0; j < l_len; ++j)
                lmer_set(st_rev, l_len-j-1, 3^lmer_get(st_fwd, j));
            u16_p = lcounter_ptr(lc, lmer_cmp(st_fwd, st_rev, l_len) < 0 ? st_fwd : st_rev);
            if (*u16_p == UINT16_MAX) {
                varr64_resize(&stack, stack.size - bl);
                continue;
            }
            sum = *u16_p;
            *u16_p = UINT16_MAX;

            lmer_cpy(lfwd, st_fwd, l_len);
            lmer_cpy(lrev, st_rev, l_len);
            lmer_rshift(lfwd, l_len);
            lmer_lshift(lrev, l_len);
            varr8_clear(&lseq);
            lflags = 0;
            for (j = 0; j < 4; ++j) {
                lmer_set(lfwd, 0, j);
                lmer_set(lrev, l_len-1, 3^j);
                if (lcounter_occ(lc, lmer_cmp(lfwd, lrev, l_len) < 0 ? lfwd : lrev))
                    lflags |= 1 << j;
            }
            if (lflags != 0)
                varr8_push(&lseq, flag_base(lflags));

            lmer_cpy(rfwd, st_fwd, l_len);
            lmer_cpy(rrev, st_rev, l_len);
            lmer_lshift(rfwd, l_len);
            lmer_rshift(rrev, l_len);
            varr8_clear(&rseq);
            rflags = 0;
            for (j = 0; j < 4; ++j) {
                lmer_set(rfwd, l_len-1, j);
                lmer_set(rrev, 0, 3^j);
                if (lcounter_occ(lc, lmer_cmp(rfwd, rrev, l_len) < 0 ? rfwd : rrev))
                    rflags |= 1 << j;
            }
            if (rflags != 0)
                varr8_push(&rseq, flag_base(rflags));

            if ((lseq.size > 0 && lseq.ele[0] == 4) || (rseq.size > 0 && rseq.ele[0] == 4)) {
                varr64_resize(&stack, stack.size - bl);
                for (j = 0; j < 4; ++j) {
                    if (rflags & (1 << j)) {
                        varr64_resize(&stack, stack.size + bl);
                        lmer_set(rfwd, l_len-1, j);
                        lmer_cpy(&(stack.ele[stack.size-bl]), rfwd, l_len);
                    }
                    if (lflags & (1 << j)) {
                        varr64_resize(&stack, stack.size + bl);
                        lmer_set(lfwd, 0, j);
                        lmer_cpy(&(stack.ele[stack.size-bl]), lfwd, l_len);
                    }
                }
                ljunc.cov = sum;
                ljunc.out = (lflags << 4) | rflags;
                fwrite(st_fwd, sizeof(uint64_t), bl, lg->junc_file);
                fwrite(&ljunc, sizeof(ljunc_t), 1, lg->junc_file);
                ++lg->n_junc;
                continue;
            }

            if (rseq.size > 0) {
                lmer_cpy(fwd, st_fwd, l_len);
                lmer_cpy(rev, st_rev, l_len);
                lmer_lshift(fwd, l_len);
                lmer_rshift(rev, l_len);
                lmer_set(fwd, l_len-1, rseq.ele[0]);
                lmer_set(rev, 0, 3^rseq.ele[0]);
                u16_p = lcounter_ptr(lc, lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev);
                while (1) {
                    lmer_cpy(rfwd, fwd, l_len);
                    lmer_cpy(rrev, rev, l_len);
                    lmer_rshift(fwd, l_len);
                    lmer_lshift(rev, l_len);
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, 0, j);
                        lmer_set(rev, l_len-1, 3^j);
                        if (lcounter_occ(lc, lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev))
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        rflags = 0 | (1 << varr8_peek(&rseq));
                        varr8_pop(&rseq);
                        break;
                    }
                    lmer_lshift(fwd, l_len);
                    lmer_lshift(fwd, l_len);
                    lmer_rshift(rev, l_len);
                    lmer_rshift(rev, l_len);
                    lmer_set(fwd, l_len-2, varr8_peek(&rseq));
                    lmer_set(rev, 1, 3^varr8_peek(&rseq));
                    rflags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, l_len-1, j);
                        lmer_set(rev, 0, 3^j);
                        tmp_u16_p = lcounter_ptr(lc, lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev);
                        if (tmp_u16_p != NULL) {
                            rflags |= 1 << j;
                            next_u16_p = tmp_u16_p;
                        }
                    }
                    if (rflags == 0) {
                        sum += *u16_p;
                        *u16_p = UINT16_MAX;
                        break;
                    }
                    else if (flag_degree(rflags) > 1) {
                        rflags = 0 | (1 << varr8_peek(&rseq));
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
                lmer_cpy(fwd, st_fwd, l_len);
                lmer_cpy(rev, st_rev, l_len);
                lmer_rshift(fwd, l_len);
                lmer_lshift(rev, l_len);
                lmer_set(fwd, 0, lseq.ele[0]);
                lmer_set(rev, l_len-1, 3^lseq.ele[0]);
                u16_p = lcounter_ptr(lc, lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev);
                while (1) {
                    lmer_cpy(lfwd, fwd, l_len);
                    lmer_cpy(lrev, rev, l_len);
                    lmer_lshift(fwd, l_len);
                    lmer_rshift(rev, l_len);
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, l_len-1, j);
                        lmer_set(rev, 0, 3^j);
                        if (lcounter_occ(lc, lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev))
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        lflags = 0 | (1 << varr8_peek(&lseq));
                        varr8_pop(&lseq);
                        break;
                    }
                    lmer_rshift(fwd, l_len);
                    lmer_rshift(fwd, l_len);
                    lmer_lshift(rev, l_len);
                    lmer_lshift(rev, l_len);
                    lmer_set(fwd, 1, varr8_peek(&lseq));
                    lmer_set(rev, l_len-2, 3^varr8_peek(&lseq));
                    lflags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, 0, j);
                        lmer_set(rev, l_len-1, 3^j);
                        tmp_u16_p = lcounter_ptr(lc, lmer_cmp(fwd, rev, l_len) < 0 ? fwd : rev);
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
                        lflags = 0 | (1 << varr8_peek(&lseq));
                        varr8_pop(&lseq);
                        break;
                    }
                    varr8_push(&lseq, flag_base(lflags));
                    sum += *u16_p;
                    *u16_p = UINT16_MAX;
                    u16_p = next_u16_p;
                }
            }

            varr64_resize(&stack, stack.size - bl);

            if (flag_degree(lflags)) {
                varr64_resize(&stack, stack.size + bl);
                lmer_set(lfwd, 0, flag_base(lflags));
                lmer_cpy(&(stack.ele[stack.size-bl]), lfwd, l_len);
            }
            if (flag_degree(rflags)) {
                varr64_resize(&stack, stack.size + bl);
                lmer_set(rfwd, l_len-1, flag_base(rflags));
                lmer_cpy(&(stack.ele[stack.size-bl]), rfwd, l_len);
            }

            lstra.len = lseq.size + rseq.size + 1;
            varr8_resize(&seq, (lstra.len+l_len+2)/4);
            for (j = 0; j < lseq.size; ++j)
                bseq_set(seq.ele, j, lseq.ele[lseq.size-j-1]);
            for (j = 0; j < l_len; ++j)
                bseq_set(seq.ele, lseq.size+j, lmer_get(st_fwd, j));
            for (j = 0; j < rseq.size; ++j)
                bseq_set(seq.ele, lseq.size+l_len+j, rseq.ele[j]);
            lstra.cov = (uint16_t)((double)sum / (double)(lstra.len) + 0.5);
            lstra.out = (lflags << 4) | rflags;
            fwrite(&lstra, sizeof(lstra_t), 1, lg->stra_file);
            fwrite(seq.ele, sizeof(uint8_t), (lstra.len+l_len+2)/4, lg->stra_file);
            ++lg->n_stra;
        }
    }

    varr8_destroy(&rseq);
    varr8_destroy(&lseq);
    varr8_destroy(&seq);
    varr64_destroy(&stack);
}

void lgraph_load(lgraph_t *lg)
{
    uint64_t mem_usage;
    uint64_t l_len;
    uint64_t bl;
    lmer_t lkey;
    ljunc_t ljunc;
    lstra_t lstra;
    lstra_t *sp;

    fputs("loading kmer graph...\n", stderr);

    l_len = lg->l_len;
    bl = lg->bl;

    if (lg->jval_table != NULL)
        lgraph_destroy_table(lg);
    
    fseek(lg->stra_file, 0, SEEK_END);
    mem_usage = ftell(lg->stra_file) - 2*lg->n_stra*bl*sizeof(uint64_t);

    lg->jtable_size = 1;
    lg->jindex_len = 0;
    while (lg->jtable_size * MAX_LOAD_FACTOR < lg->n_junc) {
        lg->jtable_size *= 2;
        ++lg->jindex_len;
    }

    lg->stable_size = 1;
    lg->sindex_len = 0;
    while (lg->stable_size * MAX_LOAD_FACTOR < lg->n_stra) {
        lg->stable_size *= 2;
        ++lg->sindex_len;
    }

    mem_usage += lg->jtable_size * (sizeof(ljunc_t) + bl*sizeof(uint64_t));
    mem_usage += lg->stable_size * 2 * (sizeof(lstra_t *) + bl*sizeof(uint64_t));
    check_mem_usage(mem_usage, lg->mem_lim);

    while (1) {
        mem_usage += lg->jtable_size * (sizeof(ljunc_t) + sizeof(uint64_t));
        mem_usage += lg->stable_size * 2 * (sizeof(lstra_t *) + bl*sizeof(uint64_t));
        if (mem_usage >= lg->mem_lim / 4)
            break;
        lg->jtable_size *= 2;
        ++lg->jindex_len;
        lg->stable_size *= 2;
        ++lg->sindex_len;
    }
    mem_usage -= lg->jtable_size * (sizeof(ljunc_t) + sizeof(uint64_t));
    mem_usage -= lg->stable_size * 2 * (sizeof(lstra_t *) + bl*sizeof(uint64_t));

    fprintf(stderr, "JUNCTION_LOAD_FACTOR=%f, ", (double)(lg->n_junc) / (double)lg->jtable_size);
    fprintf(stderr, "STRAIGHT_LOAD_FACTOR=%f, ", (double)(lg->n_stra) / (double)lg->stable_size);
    fprintf(stderr, "MEM_USAGE=%luB\n", mem_usage); 

    lg->jval_table = (ljunc_t *)my_calloc(lg->jtable_size, sizeof(ljunc_t));
    lg->jkey_table = (uint64_t *)my_calloc(lg->jtable_size * bl, sizeof(uint64_t));
    lg->lsval_table = (lstra_t **)my_calloc(lg->stable_size, sizeof(lstra_t *));
    lg->lskey_table = (uint64_t *)my_calloc(lg->stable_size * bl, sizeof(uint64_t));
    lg->rsval_table = (lstra_t **)my_calloc(lg->stable_size, sizeof(lstra_t *));
    lg->rskey_table = (uint64_t *)my_calloc(lg->stable_size * bl, sizeof(uint64_t));

    rewind(lg->junc_file);
    while (fread(lkey, sizeof(uint64_t), bl, lg->junc_file)) {
        fread(&ljunc, sizeof(ljunc_t), 1, lg->junc_file);
        lgraph_ljunc_insert(lg, &ljunc, lkey);
    }

    rewind(lg->stra_file);
    while (fread(&lstra, sizeof(lstra_t), 1, lg->stra_file)) {
        sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(lstra.len+l_len+2)/4);
        *(sp) = lstra;
        fread(sp->seq, sizeof(uint8_t), (lstra.len+l_len+2)/4, lg->stra_file);
        lgraph_lstra_insert(lg, sp);
    }
}

uint64_t lgraph_crush_bubble(lgraph_t *lg, const double ave_cov, FILE **bubble_fpp)
{
    uint64_t n_crush;
    uint64_t i;
    uint64_t j;
    uint64_t k;
    uint64_t m;
    uint64_t n;
    uint64_t l_len;
    uint64_t bl;
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
    lmer_t key;
    lmer_t rkey;
    ljunc_t *rt_jp;
    ljunc_t *jp[4];
    lstra_t *sp[4];
    lstra_t *max_sp[3];

    n_crush = 0;

    if ((uint64_t)(ave_cov * 2.0 + 0.5) < UINT16_MAX)
        cov_th = (uint16_t)(ave_cov * 2.0 + 0.5);
    else
        cov_th = UINT16_MAX;

    l_len = lg->l_len;
    bl = lg->bl;
    wsize = l_len;
    work = (uint64_t *)my_malloc(wsize * 2 * sizeof(uint64_t));;
    buf_size = lgraph_sorted_junction_list(lg, &buf);

    for (i = 0; i < buf_size ; ++i) {
        rt_jp = lgraph_ljunc_find(lg, &(buf[i * bl]));
        if (rt_jp == NULL || flag_degree(rt_jp->out & 0x0f) < 2)
            continue;

        lmer_cpy(key, &(buf[i * bl]), l_len);
        lmer_lshift(key, l_len);
        for (j = 0; j < 4; ++j) {
            clust[j] = j;
            sp[j] = NULL;
            jp[j] = NULL;
            if (!(rt_jp->out & (1 << j)))
                continue;
            lmer_set(key, l_len-1, j);
            sp[j] = lgraph_lstra_find_pre(lg, key);
            if (sp[j] != NULL && (sp[j]->out & 0x0f)) {
                for (k = 0; k < l_len-1; ++k)
                    lmer_set(rkey, k, bseq_get(sp[j]->seq, sp[j]->len+k));
                lmer_set(rkey, l_len-1, flag_base(sp[j]->out & 0x0f));
                jp[j] = lgraph_ljunc_find(lg, rkey);
            }
        }

        for (j = 0; j < 3; ++j) {
            if (sp[j] == NULL)
                continue;
            for (k = j+1; k < 4; ++k) {
                if (sp[k] == NULL || jp[j] == NULL || jp[j] != jp[k] || sp[j]->cov + sp[k]->cov > cov_th)
                    continue;
                al_th = sp[j]->len > sp[k]->len ? (sp[j]->len+l_len-1)*lg->bubble_th+0.5 : (sp[k]->len+l_len-1)*lg->bubble_th+0.5;
                if (sp[j]->len +1 <= l_len || sp[k]->len + 1 <= l_len) {
                    if (abs64((int64_t)sp[j]->len - (int64_t)sp[k]->len) > al_th)
                        continue;
                    if (clust[k] == k)
                        clust[k] = clust[j];
                    else
                        clust[j] = clust[k];
                    continue;
                }
                al_len1 = sp[j]->len - l_len + 1;
                al_len2 = sp[k]->len - l_len + 1;
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
                        if (bseq_get(sp[j]->seq, m+l_len-1) == bseq_get(sp[k]->seq, n+l_len-1)) {
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
                jp[k]->out ^= 1 << (bseq_get(sp[k]->seq, sp[k]->len-1) + 4);

                if (bubble_fpp != NULL) {
                    tmp64 = sp[k]->len + l_len - 1;
                    fwrite(&tmp64, sizeof(uint64_t), 1, *bubble_fpp);
                    fwrite(&(sp[k]->cov), sizeof(uint16_t), 1, *bubble_fpp);
                    fwrite(sp[k]->seq, sizeof(uint8_t), (tmp64 + 3)/4, *bubble_fpp);
                }

                lgraph_lstra_delete(lg, sp[k]);
                ++n_crush;
            }
        }
    }

    my_free(buf);
    my_free(work);
    
    return n_crush;
}

lstra_t *lgraph_join_ss(lgraph_t *lg, lstra_t *lsp, lstra_t *rsp)
{
    uint64_t i;
    lstra_t *new_sp;

    if (lsp == rsp)
        return lsp;
    new_sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(lsp->len+rsp->len+lg->l_len+2)/4);
    new_sp->len = lsp->len + rsp->len;
    new_sp->cov = (double)(lsp->cov*(lsp->len) + rsp->cov*(rsp->len)) / (double)(new_sp->len) + 0.5;
    new_sp->out = (lsp->out & 0xf0) | (rsp->out & 0x0f);
    for (i = 0; i < lsp->len + lg->l_len - 1; ++i)
        bseq_set(new_sp->seq, i, bseq_get(lsp->seq, i));
    for (i = 0; i < rsp->len; ++i)
        bseq_set(new_sp->seq, i + lsp->len + lg->l_len - 1, bseq_get(rsp->seq, i+lg->l_len-1));
    lgraph_lstra_delete(lg, lsp);
    lgraph_lstra_delete(lg, rsp);
    lgraph_lstra_insert(lg, new_sp);
    return new_sp;
}

lstra_t *lgraph_join_sj(lgraph_t *lg, lstra_t *sp, ljunc_t *jp)
{
    uint64_t i;
    lstra_t *new_sp;

    new_sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(sp->len+lg->l_len+3)/4);
    new_sp->len = sp->len + 1;
    new_sp->cov = (double)(sp->cov*(sp->len) + jp->cov) / (double)(new_sp->len) + 0.5;
    new_sp->out = (sp->out & 0xf0) | (jp->out & 0x0f);
    for (i = 0; i < sp->len + lg->l_len - 1; ++i)
        bseq_set(new_sp->seq, i, bseq_get(sp->seq, i));
    i = (uint64_t)(jp - lg->jval_table); 
    bseq_set(new_sp->seq, sp->len + lg->l_len - 1, lmer_get(&(lg->jkey_table[i * lg->bl]), lg->l_len-1));
    lgraph_lstra_delete(lg, sp);
    ljunc_delete(jp);
    lgraph_lstra_insert(lg, new_sp);
    return new_sp;
}

lstra_t *lgraph_join_js(lgraph_t *lg, ljunc_t *jp, lstra_t *sp)
{
    uint64_t i;
    lstra_t *new_sp;

    new_sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(sp->len+lg->l_len+3)/4);
    new_sp->len = sp->len + 1;
    new_sp->cov  = (double)(jp->cov + sp->cov*(sp->len)) / (double)(new_sp->len) + 0.5;
    new_sp->out = (jp->out & 0xf0) | (sp->out & 0x0f);
    i = (uint64_t)(jp - lg->jval_table); 
    bseq_set(new_sp->seq, 0, lmer_get(&(lg->jkey_table[i * lg->bl]), 0));
    for (i = 0; i < sp->len + lg->l_len - 1; ++i)
        bseq_set(new_sp->seq, i+1, bseq_get(sp->seq, i));
    lgraph_lstra_delete(lg, sp);
    ljunc_delete(jp);
    lgraph_lstra_insert(lg, new_sp);
    return new_sp;
}

lstra_t *lgraph_join_jj(lgraph_t *lg, ljunc_t *ljp, ljunc_t *rjp)
{
    uint64_t i;
    lstra_t *new_sp;

    if (ljp == rjp)
        return lgraph_change_j2s(lg, ljp);
    new_sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(lg->l_len+4)/4);
    new_sp->len = 2;
    new_sp->cov = (double)(ljp->cov + rjp->cov) / 2.0 + 0.5;
    new_sp->out = (ljp->out & 0xf0) | (rjp->out & 0x0f);
    bseq_set(new_sp->seq, 0, lmer_get(&(lg->jkey_table[(uint64_t)(ljp-lg->jval_table) * lg->bl]), 0));
    for (i = 0; i < lg->l_len; ++i)
        bseq_set(new_sp->seq, i+1, lmer_get(&(lg->jkey_table[(uint64_t)(rjp-lg->jval_table) * lg->bl]), i));
    ljunc_delete(ljp);
    ljunc_delete(rjp);
    lgraph_lstra_insert(lg, new_sp);
    return new_sp;
}

lstra_t *lgraph_change_j2s(lgraph_t *lg, ljunc_t *jp)
{
    uint64_t i;
    lstra_t *new_sp;

    new_sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(lg->l_len+3)/4);
    new_sp->len = 1;
    new_sp->cov = jp->cov;
    new_sp->out = jp->out;
    for (i = 0; i < lg->l_len; ++i)
        bseq_set(new_sp->seq, i, lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * lg->bl]), i));
    ljunc_delete(jp);
    lgraph_lstra_insert(lg, new_sp);
    return new_sp;
}

void lgraph_lstra_delete(lgraph_t *lg, lstra_t *sp)
{
    uint64_t i;
    uint64_t index;
    uint64_t step;
    lmer_t rkey;
    lmer_t lkey;

    for (i = 0; i < lg->l_len; ++i) {
        lmer_set(lkey, i, bseq_get(sp->seq, i));
        lmer_set(rkey, i, bseq_get(sp->seq, sp->len + i - 1));
    }
        
    index = lmer_hash(lkey, lg->l_len, lg->sindex_len) & (lg->stable_size-1);
    if (!lstra_is_del(lg->lsval_table[index]) && lmer_cmp(lkey, &(lg->lskey_table[index * lg->bl]), lg->l_len) == 0)
        lg->lsval_table[index] = &g_lstra_deleted;
    else {
        step = lmer_rehash(lkey, lg->l_len, lg->sindex_len);
        index = (index + step) & (lg->stable_size-1);
        while (lstra_is_del(lg->lsval_table[index]) || lmer_cmp(lkey, &(lg->lskey_table[index * lg->bl]), lg->l_len) != 0)
            index = (index+step) & (lg->stable_size-1);
        lg->lsval_table[index] = &g_lstra_deleted;
    }

    index = lmer_hash(rkey, lg->l_len, lg->sindex_len) & (lg->stable_size-1);
    if (!lstra_is_del(lg->rsval_table[index]) && lmer_cmp(rkey, &(lg->rskey_table[index * lg->bl]), lg->l_len) == 0)
        lg->rsval_table[index] = &g_lstra_deleted;
    else {
        step = lmer_rehash(rkey, lg->l_len, lg->sindex_len);
        index = (index + step) & (lg->stable_size-1);
        while (lstra_is_del(lg->rsval_table[index]) || lmer_cmp(rkey, &(lg->rskey_table[index * lg->bl]), lg->l_len) != 0)
            index = (index+step) & (lg->stable_size-1);
        lg->rsval_table[index] = &g_lstra_deleted;
    }
    my_free(sp);
}

void lgraph_lstra_insert(lgraph_t *lg, lstra_t *sp)
{
    uint64_t i;
    uint64_t index;
    uint64_t step;
    lmer_t rkey;
    lmer_t lkey;

    for (i = 0; i < lg->l_len; ++i) {
        lmer_set(lkey, i, bseq_get(sp->seq, i));
        lmer_set(rkey, i, bseq_get(sp->seq, sp->len + i - 1));
    }

    index = lmer_hash(lkey, lg->l_len, lg->sindex_len) & (lg->stable_size-1);
    if (lstra_is_emp(lg->lsval_table[index]) || lstra_is_del(lg->lsval_table[index])) {
        lg->lsval_table[index] = sp;
        lmer_cpy(&(lg->lskey_table[index * lg->bl]), lkey, lg->l_len);
    }
    else {
        step = lmer_rehash(lkey, lg->l_len, lg->sindex_len);
        index = (index+step) & (lg->stable_size-1);
        while (!(lstra_is_emp(lg->lsval_table[index]) || lstra_is_del(lg->lsval_table[index])))
            index = (index+step) & (lg->stable_size-1);
        lg->lsval_table[index] = sp;
        lmer_cpy(&(lg->lskey_table[index * lg->bl]), lkey, lg->l_len);
    }

    index = lmer_hash(rkey, lg->l_len, lg->sindex_len) & (lg->stable_size-1);
    if (lstra_is_emp(lg->rsval_table[index]) || lstra_is_del(lg->rsval_table[index])) {
        lg->rsval_table[index] = sp;
        lmer_cpy(&(lg->rskey_table[index * lg->bl]), rkey, lg->l_len);
    }
    else {
        step = lmer_rehash(rkey, lg->l_len, lg->sindex_len);
        index = (index+step) & (lg->stable_size-1);
        while (!(lstra_is_emp(lg->rsval_table[index]) || lstra_is_del(lg->rsval_table[index])))
            index = (index+step) & (lg->stable_size-1);
        lg->rsval_table[index] = sp;
        lmer_cpy(&(lg->rskey_table[index * lg->bl]), rkey, lg->l_len);
    }
}

void lgraph_extract_read(const lgraph_t *lg, FILE **read_fpp)
{
    uint64_t i;
    uint64_t l_len;
    lmer_t fwd;
    lmer_t rev;
    seq_t seq;
    seq_t tmp_seq;
    FILE *new_read_fp;

    l_len = lg->l_len;

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

            if (lgraph_ljunc_find(lg, fwd) != NULL) {
                seq_write(&seq, new_read_fp);
                break;
            }
            else if (lgraph_ljunc_find(lg, rev) != NULL) {
                seq_write(&seq, new_read_fp);
                break;
            }

        }
    }

    fclose(*read_fpp);
    *read_fpp = new_read_fp;
}

void lgraph_delete_clear_junc(const lgraph_t *lg, const uint64_t len)
{
    uint64_t i;
    uint64_t j;
    uint64_t l_len;
    uint64_t bl;
    lmer_t lkey;
    lmer_t rkey;
    ljunc_t *jp;
    lstra_t *sp;

    l_len = lg->l_len;
    bl = lg->bl;
    
    for (i = 0; i < lg->jtable_size; ++i) {
        jp = &(lg->jval_table[i]);
        if (ljunc_is_emp(*jp) || ljunc_is_del(*jp))
            continue;
        lmer_cpy(lkey, &(lg->jkey_table[i * bl]), l_len);
        lmer_cpy(rkey, lkey, l_len);
        lmer_rshift(lkey, l_len);
        lmer_lshift(rkey, l_len);
        for (j = 0; j < 4; ++j) {
            if (jp->out & (1 << (j+4))) {
                lmer_set(lkey, 0, j);
                if ((sp = lgraph_lstra_find_suf(lg, lkey)) == NULL  || sp->len < len)
                    break;
            }
            if (jp->out & (1 << j)) {
                lmer_set(rkey, l_len-1, j);
                if ((sp = lgraph_lstra_find_pre(lg, rkey)) == NULL ||  sp->len < len)
                    break;
            }
        }

        if (j == 4)
            ljunc_delete(jp);
    }
}

uint64_t lgraph_cut_branch(lgraph_t *lg)
{
    uint64_t i;
    uint64_t j;
    uint64_t l_len;
    uint64_t bl;
    uint64_t max;
    uint64_t n_delete;
    uint64_t *buf;
    uint64_t buf_size;
    lmer_t key;
    ljunc_t *tmp_jp;
    ljunc_t *jp;
    lstra_t *sp;
    lstra_t *tmp_sp;

    n_delete = 0;

    l_len = lg->l_len;
    bl = lg->bl;

    buf_size = lgraph_sorted_prefix_list(lg, &buf);
    for (i = 0; i < buf_size ; ++i) {
        sp = lgraph_lstra_find_pre(lg, &(buf[i * bl]));
        if (sp == NULL || sp->len > l_len)
            continue;

        if (sp->out >> 4) {
            if (sp->out & 0xf)
                continue;
            lmer_cpy(key, &(buf[i * bl]), l_len);
            lmer_rshift(key, l_len);
            lmer_set(key, 0, flag_base(sp->out >> 4));
            jp = lgraph_ljunc_find(lg, key);
            if (jp == NULL)
                continue;

            lmer_lshift(key, l_len);
            max = 0;
            for (j = 0; j < 4; ++j) {
                if (!(jp->out & (1 << j)))
                    continue;
                lmer_set(key, l_len - 1, j);
                tmp_sp = lgraph_lstra_find_pre(lg, key);
                if (tmp_sp != sp && tmp_sp != NULL && tmp_sp->cov > max) {
                    max = tmp_sp->cov;
                    continue;
                }

                tmp_jp = lgraph_ljunc_find(lg, key);
                if (tmp_jp != NULL && tmp_jp->cov > max)
                    max = tmp_jp->cov;
            }
            if (sp->cov > max * lg->branch_th)
                continue;

            jp->out ^= 1 << bseq_get(sp->seq, l_len - 1);
        }
        else if (sp->out & 0xf) {
            for (j = 1; j < l_len; ++j)
                lmer_set(key, j - 1, bseq_get(sp->seq, sp->len + j - 1));
            lmer_set(key, l_len - 1, flag_base(sp->out & 0xf));
            jp = lgraph_ljunc_find(lg, key);
            if (jp == NULL)
                continue;

            lmer_rshift(key, l_len);
            max = 0;
            for (j = 0; j < 4; ++j) {
                if (!(jp->out & (1 << (j+4))))
                    continue;
                lmer_set(key, 0, j);
                tmp_sp = lgraph_lstra_find_suf(lg, key);
                if (tmp_sp != sp && tmp_sp != NULL && tmp_sp->cov > max) {
                    max = tmp_sp->cov;
                    continue;
                }

                tmp_jp = lgraph_ljunc_find(lg, key);
                if (tmp_jp != NULL && tmp_jp->cov > max)
                    max = tmp_jp->cov;
            }
            if (sp->cov > max * lg->branch_th)
                continue;

            jp->out ^= 1 << (4 + bseq_get(sp->seq, sp->len - 1));
        }
        else
            continue;

        lgraph_lstra_delete(lg, sp);
        ++n_delete;
    }
    my_free(buf);

    return n_delete;
}

void lgraph_join_nodes(lgraph_t *lg)
{
    uint64_t i;
    uint64_t j;
    uint64_t l_len;
    uint64_t bl;
    uint64_t *buf;
    uint64_t buf_size;
    lmer_t key;
    ljunc_t *jp;
    ljunc_t *ljp;
    ljunc_t *rjp;
    lstra_t *sp;
    lstra_t *lsp;
    lstra_t *rsp;

    l_len = lg->l_len;
    bl = lg->bl;

    buf_size = lgraph_sorted_junction_list(lg, &buf);
    for (i = 0; i < buf_size; ++i) {
        jp = lgraph_ljunc_find(lg, &(buf[i * bl]));
        if (jp == NULL || flag_degree(jp->out >> 4) > 1 || flag_degree(jp->out & 0xf) > 1)
            continue;

        if (flag_degree(jp->out >> 4)) {
            lmer_cpy(key, &(buf[i * bl]), l_len);
            lmer_rshift(key, l_len);
            lmer_set(key, 0, flag_base(jp->out >> 4));
            lsp = lgraph_lstra_find_suf(lg, key);
            if (lsp == NULL) {
                ljp = lgraph_ljunc_find(lg, key);
                if (ljp == NULL || flag_degree(ljp->out >> 4) > 1 || flag_degree(ljp->out & 0xf) > 1)
                    ljp = NULL;
            }
        }
        else {
            lsp = NULL;
            ljp = NULL;
        }

        if (lsp != NULL)
            lsp = lgraph_join_sj(lg, lsp, jp);
        else if (ljp != NULL)
            lsp = lgraph_join_jj(lg, ljp, jp);
        else
            lsp = lgraph_change_j2s(lg, jp);

        if (flag_degree(lsp->out & 0xf)) {
            for (j = 1; j < l_len; ++j)
                lmer_set(key, j - 1, bseq_get(lsp->seq, lsp->len + j - 1));
            lmer_set(key, l_len - 1, flag_base(lsp->out & 0xf));
            rsp = lgraph_lstra_find_pre(lg, key);
            if (rsp == NULL) {
                rjp = lgraph_ljunc_find(lg, key);
                if (rjp == NULL || flag_degree(rjp->out >> 4) > 1 || flag_degree(rjp->out & 0xf) > 1)
                    rjp = NULL;
            }
        }
        else {
            rsp = NULL;
            rjp = NULL;
        }

        if (rsp != NULL)
            lgraph_join_ss(lg, lsp, rsp);
        else if (rjp != NULL)
            lgraph_join_sj(lg, lsp, rjp);
    }
    my_free(buf);

    buf_size = lgraph_sorted_prefix_list(lg, &buf);
    for (i = 0; i < buf_size ; ++i) {
        sp = lgraph_lstra_find_pre(lg, &(buf[i * bl]));
        if (sp == NULL)
            continue;

        if (flag_degree(sp->out >> 4)) {
            lmer_cpy(key, &(buf[i * bl]), l_len);
            lmer_rshift(key, l_len);
            lmer_set(key, 0, flag_base(sp->out >> 4));
            lsp = lgraph_lstra_find_suf(lg, key);
            if (lsp != NULL)
                sp = lgraph_join_ss(lg, lsp, sp);
        }

        if (flag_degree(sp->out & 0xf)) {
            for (j = 1; j < l_len; ++j)
                lmer_set(key, j - 1, bseq_get(sp->seq, sp->len + j - 1));
            lmer_set(key, l_len - 1, flag_base(sp->out & 0xf));
            rsp = lgraph_lstra_find_pre(lg, key);
            if (rsp != NULL)
                lgraph_join_ss(lg, sp, rsp);
        }
    }
    my_free(buf);
}

FILE *lgraph_merge(lgraph_t *lg, lgraph_t *ext_lg)
{
    uint64_t i;
    uint64_t j;
    uint64_t l_len;
    uint64_t bl;
    uint64_t len;
    uint8_t out;
    FILE *contig_fp;
    varr64_t stack;
    lmer_t fwd;
    lmer_t rev;
    uint8_t tmp_seq[MAX_READ_LEN] = {0};
    lstra_t *sp;
    ljunc_t *jp;

    fprintf(stderr, "merging kmer graphs...\n");

    l_len = lg->l_len;
    bl = lg->bl;
    varr64_init(&stack);

    contig_fp = lgraph_save_contig(lg, lg->l_len);
    fseek(contig_fp, 0, SEEK_END);
    for (i = 0; i < lg->jtable_size; ++i) {
        jp = &(lg->jval_table[i]);
        if (ljunc_is_emp(*jp) || ljunc_is_del(*jp))
            continue;
        for (j = 0; j < l_len; ++j)
            bseq_set(tmp_seq, j, lmer_get(&(lg->jkey_table[i * bl]), j));
        fwrite(&l_len, sizeof(uint64_t), 1, contig_fp);
        fwrite(&(jp->cov), sizeof(uint16_t), 1, contig_fp);
        fwrite(tmp_seq, sizeof(uint8_t), (l_len + 3)/4, contig_fp);
    }

    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (!(lstra_is_emp(sp) || lstra_is_del(sp) || (sp->out >> 4))) {
            lmer_cpy(&(lg->lskey_table[i * bl]), fwd, l_len);
            lmer_rshift(fwd, l_len);
            for (j = 0; j < 4; ++j) {
                lmer_set(fwd, 0, j);
                varr64_resize(&stack, stack.size + bl);
                lmer_cpy(&(stack.ele[stack.size - bl]), fwd, l_len);
            }
        }
        sp = lg->rsval_table[i];
        if (!(lstra_is_emp(sp) || lstra_is_del(sp) || (sp->out & 0xf))) {
            lmer_cpy(&(lg->rskey_table[i * bl]), fwd, l_len);
            lmer_lshift(fwd, l_len);
            for (j = 0; j < 4; ++j) {
                lmer_set(fwd, l_len - 1, j);
                varr64_resize(&stack, stack.size + bl);
                lmer_cpy(&(stack.ele[stack.size - bl]), fwd, l_len);
            }
        }

        while (stack.size > 0) {
            lmer_cpy(fwd, &(stack.ele[stack.size - bl]), l_len);
            varr64_resize(&stack, stack.size - bl);
            for (j = 0; j < l_len; ++j)
                lmer_set(rev, j, 3^lmer_get(fwd, l_len - j - 1));

            if ((sp = lgraph_lstra_find_pre(ext_lg, fwd)) != NULL ||
                (sp = lgraph_lstra_find_pre(ext_lg, rev)) != NULL ||
                (sp = lgraph_lstra_find_suf(ext_lg, fwd)) != NULL ||
                (sp = lgraph_lstra_find_suf(ext_lg, rev)) != NULL) {

                if (sp->cov == UINT16_MAX)
                    continue;
                len = sp->len + l_len - 1;
                fwrite(&len, sizeof(uint64_t), 1, contig_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, contig_fp);
                fwrite(sp->seq, sizeof(uint8_t), (len + 3)/4, contig_fp);
                sp->cov = UINT16_MAX;
                out = sp->out;
                for (j = 0; j < l_len - 1; ++j) {
                    lmer_set(fwd, j + 1, bseq_get(sp->seq, j));
                    lmer_set(rev, j, bseq_get(sp->seq, sp->len + j));
                }
            }
            else if ((jp = lgraph_ljunc_find(ext_lg, fwd)) != NULL ||
                (jp = lgraph_ljunc_find(ext_lg, rev)) != NULL ||
                (jp = lgraph_ljunc_find(ext_lg, fwd)) != NULL ||
                (jp = lgraph_ljunc_find(ext_lg, rev)) != NULL) {

                if (jp->cov == UINT16_MAX)
                    continue;
                for (j = 0; j < l_len; ++j)
                    bseq_set(tmp_seq, j, lmer_get(fwd, j));
                fwrite(&l_len, sizeof(uint64_t), 1, contig_fp);
                fwrite(&(jp->cov), sizeof(uint16_t), 1, contig_fp);
                fwrite(tmp_seq, sizeof(uint8_t), (l_len + 3)/4, contig_fp);
                jp->cov = UINT16_MAX;
                out = jp->out;
                lmer_cpy(rev, fwd, l_len);
                lmer_rshift(fwd, l_len);
                lmer_lshift(rev, l_len);
            }
            else
                continue;

            for (j = 0; j < 4; ++j) {
                if (out & (1 << (j+4))) {
                    lmer_set(fwd, 0, j);
                    varr64_resize(&stack, stack.size + bl);
                    lmer_cpy(&(stack.ele[stack.size - bl]), fwd, l_len);
                }

                if (out & (1 << j)) {
                    lmer_set(rev, l_len - 1, j);
                    varr64_resize(&stack, stack.size + bl);
                    lmer_cpy(&(stack.ele[stack.size - bl]), rev, l_len);
                }
            }
        }
    }

    varr64_destroy(&stack);
    lgraph_destroy_table(ext_lg);
    lgraph_destroy_table(lg);

    return contig_fp;
}

void lgraph_show_contig(const lgraph_t *lg, const double cov_ratio, const char *out_name)
{
    uint64_t i;
    uint64_t j;
    uint64_t n_seq;
    uint64_t buf_size;
    uint64_t *buf;
    lstra_t *sp;
    FILE *out;

    out = fopen(out_name, "w");
    check_file_open(out, out_name);
    buf_size = lgraph_sorted_prefix_list(lg, &buf);

    n_seq = 0;
    for (i = 0; i < buf_size ; ++i) {
        sp = lgraph_lstra_find_pre(lg, &(buf[i * lg->bl]));
        if (sp == NULL)
            continue;
        ++n_seq;
        fprintf(out, ">seq%lu_len%lu_cov%u\n", n_seq, sp->len + lg->l_len - 1, (uint16_t)(sp->cov * cov_ratio + 0.5));
        for (j = 0; j < sp->len + lg->l_len - 1; ++j) {
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

FILE *lgraph_save_contig_simple(const lgraph_t *lg)
{
    uint64_t i;
    uint64_t j;
    FILE *contig_fp;
    lstra_t *sp;
    varr8_t seq;

    varr8_init(&seq);
    contig_fp = tmpfile_open();
    check_tmpfile(contig_fp);

    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp))
            continue;

        varr8_resize(&seq, (sp->len + lg->l_len + 2) / 4);
        for (j = 0; j < sp->len + lg->l_len - 1; ++j)
            bseq_set(seq.ele, j, bseq_get(sp->seq, j));
            
        j = sp->len + lg->l_len - 1;
        fwrite(&j, sizeof(uint64_t), 1, contig_fp);
        fwrite(&(sp->cov), sizeof(uint16_t), 1, contig_fp);
        fwrite(seq.ele, sizeof(uint8_t), seq.size, contig_fp);
    }

    varr8_destroy(&seq);

    return contig_fp;
}

FILE *lgraph_save_contig(const lgraph_t *lg, const uint64_t next_l)
{
    uint64_t i;
    uint64_t j;
    uint64_t len;
    uint64_t ex_len;
    uint64_t dif;
    uint64_t l_len;
    uint64_t bl;
    FILE *contig_fp;
    lmer_t key;
    lstra_t *sp;
    lstra_t *tmp_sp;
    ljunc_t *jp;
    varr8_t seq;

    l_len = lg->l_len;
    bl = lg->bl;
    dif = next_l - l_len;
    varr8_init(&seq);
    contig_fp = tmpfile_open();
    check_tmpfile(contig_fp);

    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp))
            continue;

        len = 0;
        if (sp->out >> 4) {
            lmer_cpy(key, &(lg->lskey_table[i * bl]), l_len);
            lmer_rshift(key, l_len);
            lmer_set(key, 0, flag_base(sp->out >> 4));
            jp = lgraph_ljunc_find(lg, key);
            if (jp != NULL) {
                if (flag_degree(jp->out >> 4) == 1) {
                    lmer_rshift(key, l_len);
                    lmer_set(key, 0, flag_base(jp->out >> 4));
                    tmp_sp = lgraph_lstra_find_suf(lg, key);
                    if (tmp_sp != NULL) {
                        ex_len = min_u64(tmp_sp->len, dif);
                        varr8_resize(&seq, (len + ex_len + 3) / 4);
                        for (j = 0; j < ex_len; ++j)
                            bseq_set(seq.ele, j, bseq_get(tmp_sp->seq, tmp_sp->len - ex_len + j));
                        len += ex_len;
                    }
                }
                varr8_resize(&seq, len / 4 + 1);
                bseq_set(seq.ele, len, lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * bl]), 0));
                ++len;
            }
        }

        varr8_resize(&seq, (len + sp->len + l_len + 2) / 4);
        for (j = 0; j < sp->len + l_len - 1; ++j)
            bseq_set(seq.ele, len + j, bseq_get(sp->seq, j));
        len += sp->len + l_len - 1;

        if (sp->out & 0xf) {
            for (j = 0; j < l_len - 1; ++j)
                lmer_set(key, j, bseq_get(sp->seq, sp->len + j));
            lmer_set(key, l_len - 1, flag_base(sp->out & 0xf));
            jp = lgraph_ljunc_find(lg, key);
            if (jp != NULL) {
                varr8_resize(&seq, len / 4 + 1);
                bseq_set(seq.ele, len, lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * bl]), l_len - 1));
                ++len;
                if (flag_degree(jp->out & 0xf) == 1) {
                    lmer_lshift(key, l_len);
                    lmer_set(key, l_len - 1, flag_base(jp->out & 0xf));
                    tmp_sp = lgraph_lstra_find_pre(lg, key);
                    if (tmp_sp != NULL) {
                        ex_len = min_u64(tmp_sp->len, dif);
                        varr8_resize(&seq, (len + ex_len + 3) / 4);
                        for (j = 0; j < ex_len; ++j)
                            bseq_set(seq.ele, len + j, bseq_get(tmp_sp->seq, l_len - 1 + j));
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

void lgraph_cut_branch_iterative(lgraph_t *lg)
{
    uint64_t n_delete;
    uint64_t total_delete;

    fputs("removing branches...\n", stderr);
    fprintf(stderr, "BRANCH_DELETE_THRESHOLD=%f\n", lg->branch_th);
    total_delete = 0;
    do {
        n_delete = lgraph_cut_branch(lg);
        fprintf(stderr, "NUM_CUT=%lu\n", n_delete);
        lgraph_join_nodes(lg);
        total_delete += n_delete;
    } while (n_delete > 0);
    fprintf(stderr, "TOTAL_NUM_CUT=%lu\n", total_delete);
}

FILE *lgraph_crush_bubble_iterative(lgraph_t *lg, const double ave_cov)
{
    uint64_t n_crush;
    uint64_t total_crush;
    FILE *bubble_fp;

    fputs("removing bubbles...\n", stderr);

    bubble_fp = tmpfile_open();
    check_tmpfile(bubble_fp);

    fprintf(stderr, "BUBBLE_IDENTITY_THRESHOLD=%f\n", lg->bubble_th);
    total_crush = 0;
    do {
        n_crush = lgraph_crush_bubble(lg, ave_cov, &bubble_fp);
        fprintf(stderr, "NUM_CRUSH=%lu\n", n_crush);
        lgraph_join_nodes(lg);
        total_crush += n_crush;
    } while (n_crush > 0);
    fprintf(stderr, "TOTAL_NUM_CRUSH=%lu\n", total_crush);

    return bubble_fp;
}

uint64_t lgraph_sorted_prefix_list(const lgraph_t *lg, uint64_t **buf)
{
    uint64_t i;
    uint64_t buf_size;
    lstra_t *sp;

    buf_size = 0;
    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp))
            continue;
        ++buf_size;
    }
    *buf = (uint64_t *)my_malloc(buf_size * lg->bl * sizeof(uint64_t));

    buf_size = 0;
    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp))
            continue;
        lmer_cpy(&((*buf)[buf_size * lg->bl]), &(lg->lskey_table[i * lg->bl]), lg->l_len);
        ++buf_size;
    }
    g_l_len = lg->l_len;
    qsort(*buf, buf_size, sizeof(uint64_t) * lg->bl, cmp_lmer_uniq);

    return buf_size;
}

uint64_t lgraph_sorted_junction_list(const lgraph_t *lg, uint64_t **buf)
{
    uint64_t i;
    uint64_t buf_size;
    ljunc_t *jp;

    buf_size = 0;
    for (i = 0; i < lg->jtable_size; ++i) {
        jp = &(lg->jval_table[i]);
        if (ljunc_is_emp(*jp) || ljunc_is_del(*jp))
            continue;
        ++buf_size;
    }
    *buf = (uint64_t *)my_malloc(buf_size * lg->bl * sizeof(uint64_t));

    buf_size = 0;
    for (i = 0; i < lg->jtable_size; ++i) {
        jp = &(lg->jval_table[i]);
        if (ljunc_is_emp(*jp) || ljunc_is_del(*jp))
            continue;
        lmer_cpy(&((*buf)[buf_size * lg->bl]), &(lg->jkey_table[i * lg->bl]), lg->l_len);
        ++buf_size;
    }
    g_l_len = lg->l_len;
    qsort(*buf, buf_size, sizeof(uint64_t) * lg->bl, cmp_lmer_uniq);

    return buf_size;
}

void lgraph_save_edge_kmer(const lgraph_t *lg, lcounter_t *lc, const uint64_t next_l)
{
    uint64_t i;
    uint64_t j;
    uint64_t l_len;
    uint64_t dif;
    lmer_t lfwd;
    lmer_t lrev;
    lmer_t rfwd;
    lmer_t rrev;
    lstra_t *sp;

    l_len = lg->l_len;
    dif = next_l - lg->l_len;
    if (lc->kmer_fp != NULL)
        fclose(lc->kmer_fp);
    lc->kmer_fp = tmpfile_open();
    check_tmpfile(lc->kmer_fp);

    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp))
            continue;

        if (sp->len < dif * 2) {
            for (j = 0; j < l_len - 1; ++j) {
                lmer_set(lfwd, j + 1, bseq_get(sp->seq, j));
                lmer_set(lrev, j, 3^bseq_get(sp->seq, l_len - j - 2));
            }
            for (j = 0; j < sp->len; ++j) {
                lmer_lshift(lfwd, l_len);
                lmer_set(lfwd, l_len - 1, bseq_get(sp->seq, j + l_len - 1));
                lmer_rshift(lrev, l_len);
                lmer_set(lrev, 0, 3^bseq_get(sp->seq, j + l_len - 1));
                if (lmer_cmp(lfwd, lrev, l_len) < 0)
                    fwrite(lfwd, sizeof(uint64_t), lg->bl, lc->kmer_fp);
                else
                    fwrite(lrev, sizeof(uint64_t), lg->bl, lc->kmer_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, lc->kmer_fp);
            }
        }
        else {
            for (j = 0; j < l_len - 1; ++j) {
                lmer_set(lfwd, j + 1, bseq_get(sp->seq, j));
                lmer_set(lrev, j, 3^bseq_get(sp->seq, l_len - j - 2));
                lmer_set(rfwd, j + 1, bseq_get(sp->seq, sp->len - dif + j));
                lmer_set(rrev, j, 3^bseq_get(sp->seq, sp->len - dif + l_len - j - 2));
            }
            for (j = 0; j < dif; ++j) {
                lmer_lshift(lfwd, l_len);
                lmer_set(lfwd, l_len - 1, bseq_get(sp->seq, j + l_len - 1));
                lmer_rshift(lrev, l_len);
                lmer_set(lrev, 0, 3^bseq_get(sp->seq, j + l_len - 1));
                if (lmer_cmp(lfwd, lrev, l_len) < 0)
                    fwrite(lfwd, sizeof(uint64_t), lg->bl, lc->kmer_fp);
                else
                    fwrite(lrev, sizeof(uint64_t), lg->bl, lc->kmer_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, lc->kmer_fp);

                lmer_lshift(rfwd, l_len);
                lmer_set(rfwd, l_len - 1, bseq_get(sp->seq, j + sp->len - dif + l_len - 1));
                lmer_rshift(rrev, l_len);
                lmer_set(rrev, 0, 3^bseq_get(sp->seq, j + sp->len - dif + l_len - 1));
                if (lmer_cmp(rfwd, rrev, l_len) < 0)
                    fwrite(rfwd, sizeof(uint64_t), lg->bl, lc->kmer_fp);
                else
                    fwrite(rrev, sizeof(uint64_t), lg->bl, lc->kmer_fp);
                fwrite(&(sp->cov), sizeof(uint16_t), 1, lc->kmer_fp);
            }
        }
    }
}

void lgraph_print_dot(const lgraph_t *lg, const double cov_ratio, const char *out_name)
{
    uint64_t i;
    uint64_t j;
    uint64_t b;
    uint64_t l_len;
    uint64_t bl;
    uint64_t buf_size;
    uint64_t *buf;
    lmer_t llmer;
    lmer_t rlmer;
    lstra_t *sp;
    ljunc_t *jp;
    ljunc_t *tmp_jp;
    FILE *out;
    
    l_len = lg->l_len;
    bl = lg->bl;

    if ((out = fopen(out_name, "w")) == NULL) {
        fprintf(stderr, "cannot open %s\n", out_name);
        my_exit(1);
    }

    fputs("digraph de_Bruijn_graph {\n", out);

    buf_size = lgraph_sorted_prefix_list(lg, &buf);
    for (i = 0; i < buf_size ; ++i) {
        sp = lgraph_lstra_find_pre(lg, &(buf[i * bl]));

        for (j = 0; j < l_len - 1; ++j) {
            lmer_set(llmer, j + 1, bseq_get(sp->seq, j));
            lmer_set(rlmer, j, bseq_get(sp->seq, sp->len + j));
        }

        if (sp->out >> 4) {
            lmer_set(llmer, 0, flag_base(sp->out >> 4));
            jp = lgraph_ljunc_find(lg, llmer);
            if (jp != NULL) {
                for (j = 0; j < l_len; ++j)
                    putc(num2ascii(lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * bl]), j)), out);
                fputs(" -> ", out);
                fprintf(out, "seq%lu_len%lu_cov%lu\n", i + 1, sp->len + l_len - 1, (uint64_t)(sp->cov * cov_ratio + 0.5));
            }
        }

        if (sp->out & 0xf) {
            lmer_set(rlmer, l_len - 1, flag_base(sp->out & 0xf));
            jp = lgraph_ljunc_find(lg, llmer);
            if (jp != NULL) {
                fprintf(out, "seq%lu_len%lu_cov%lu", i + 1, sp->len + l_len - 1, (uint64_t)(sp->cov * cov_ratio + 0.5));
                fputs(" -> ", out);
                for (j = 0; j < l_len; ++j)
                    putc(num2ascii(lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * bl]), j)), out);
                putc('\n', out);
            }
        }
    }
    my_free(buf);

    buf_size = lgraph_sorted_junction_list(lg, &buf);
    for (i = 0; i < buf_size; ++i) {
        jp = lgraph_ljunc_find(lg, &(buf[i * bl]));

        for (j = 0; j < l_len - 1; ++j) {
            lmer_set(llmer, j + 1, lmer_get(&(buf[i * bl]), j));
            lmer_set(rlmer, j, lmer_get(&(buf[i * bl]), j + 1));
        }

        for (b = 0; b < 4; ++b) {
            if (!(jp->out & (1 << (4 + b))))
                continue;
            lmer_set(llmer, 0, b);
            tmp_jp = lgraph_ljunc_find(lg, llmer);
            if (tmp_jp != NULL) {
                for (j = 0; j < l_len; ++j)
                    putc(num2ascii(lmer_get(&(lg->jkey_table[(uint64_t)(tmp_jp - lg->jval_table) * bl]), j)), out);
                fputs(" -> ", out);
                for (j = 0; j < l_len; ++j)
                    putc(num2ascii(lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * bl]), j)), out);
                putc('\n', out);
            }
        }

        for (b = 0; b < 4; ++b) {
            if (!(jp->out & (1 << b)))
                continue;
            lmer_set(rlmer, l_len - 1, b);
            tmp_jp = lgraph_ljunc_find(lg, rlmer);
            if (tmp_jp != NULL) {
                for (j = 0; j < l_len; ++j)
                    putc(num2ascii(lmer_get(&(lg->jkey_table[(uint64_t)(jp - lg->jval_table) * bl]), j)), out);
                fputs(" -> ", out);
                for (j = 0; j < l_len; ++j)
                    putc(num2ascii(lmer_get(&(lg->jkey_table[(uint64_t)(tmp_jp - lg->jval_table) * bl]), j)), out);
                putc('\n', out);
            }
        }
    }
    my_free(buf);

    fputs("}\n", out);
    fclose(out);
}

void lgraph_make_graph_quick(lgraph_t *lg, const lcounter_t *lc, const uint64_t min_cov)
{
    uint64_t i;
    uint64_t j;
    uint64_t l_len;
    uint64_t sum;
    uint64_t bl;
    lmer_t fwd;
    lmer_t st_fwd;
    lmer_t lfwd;
    lmer_t rfwd;
    uint16_t *u16_p;
    uint16_t *tmp_u16_p;
    uint16_t *next_u16_p = NULL;
    uint8_t rflags;
    uint8_t lflags;
    uint8_t tmp_flags;
    varr8_t rseq;
    varr8_t lseq;
    varr8_t seq;
    varr64_t stack;
    ljunc_t ljunc;
    lstra_t lstra;
    lstra_t *sp;

    l_len = lg->l_len = lc->l_len;
    bl = lg->bl = lc->bl;
    varr8_init(&rseq);
    varr8_init(&lseq);
    varr8_init(&seq);
    varr64_init(&stack);

    for (i = j = 0; i < lc->table_size; ++i) {
        if (lc->val_table[i] >= min_cov)
            ++j;
    }

    if ((double)lg->jtable_size * MAX_LOAD_FACTOR < j) {
        lg->jtable_size = lg->stable_size = 2;
        lg->jindex_len = lg->sindex_len = 1;
        while ((double)lg->jtable_size * MAX_LOAD_FACTOR < j) {
            lg->jtable_size *= 2;
            lg->stable_size *= 2;
            ++lg->jindex_len;
            ++lg->sindex_len;
        }
        lg->jval_table = (ljunc_t *)my_realloc(lg->jval_table, lg->jtable_size * sizeof(ljunc_t));
        lg->jkey_table = (uint64_t *)my_realloc(lg->jkey_table, lg->jtable_size * bl * sizeof(uint64_t));
        lg->lsval_table = (lstra_t **)my_realloc(lg->lsval_table, lg->stable_size * sizeof(lstra_t *));
        lg->lskey_table = (uint64_t *)my_realloc(lg->lskey_table, lg->stable_size * bl * sizeof(uint64_t));
        lg->rsval_table = (lstra_t **)my_realloc(lg->rsval_table, lg->stable_size * sizeof(lstra_t *));
        lg->rskey_table = (uint64_t *)my_realloc(lg->rskey_table, lg->stable_size * bl * sizeof(uint64_t));
    }

    memset(lg->jval_table, 0, lg->jtable_size * sizeof(ljunc_t));
    memset(lg->jkey_table, 0, lg->jtable_size * bl * sizeof(uint64_t));
    memset(lg->lsval_table, 0, lg->stable_size * sizeof(lstra_t *));
    memset(lg->lskey_table, 0, lg->stable_size * bl * sizeof(uint64_t));
    memset(lg->rsval_table, 0, lg->stable_size * sizeof(lstra_t *));
    memset(lg->rskey_table, 0, lg->stable_size * bl * sizeof(uint64_t));

    for (i = 0; i < lc->table_size; ++i) {
        if (lc->val_table[i] < min_cov)
            continue;
        varr64_resize(&stack, bl);
        lmer_cpy(&(stack.ele[0]), &(lc->key_table[i * bl]), l_len);
        while (stack.size > 0) {
            lmer_cpy(st_fwd, &(stack.ele[stack.size - bl]), l_len);
            u16_p = lcounter_ptr(lc, st_fwd);
            if (*u16_p == UINT16_MAX) {
                varr64_resize(&stack, stack.size - bl);
                continue;
            }

            sum = *u16_p;
            *u16_p = UINT16_MAX;

            lmer_cpy(lfwd, st_fwd, l_len);
            lmer_rshift(lfwd, l_len);
            varr8_clear(&lseq);
            lflags = 0;
            for (j = 0; j < 4; ++j) {
                lmer_set(lfwd, 0, j);
                if (lcounter_occ(lc, lfwd) >= min_cov)
                    lflags |= 1 << j;
            }
            if (lflags != 0)
                varr8_push(&lseq, flag_base(lflags));

            lmer_cpy(rfwd, st_fwd, l_len);
            lmer_lshift(rfwd, l_len);
            varr8_clear(&rseq);
            rflags = 0;
            for (j = 0; j < 4; ++j) {
                lmer_set(rfwd, l_len-1, j);
                if (lcounter_occ(lc, rfwd) >= min_cov)
                    rflags |= 1 << j;
            }
            if (rflags != 0)
                varr8_push(&rseq, flag_base(rflags));

            if ((lseq.size > 0 && lseq.ele[0] == 4) || (rseq.size > 0 && rseq.ele[0] == 4)) {
                varr64_resize(&stack, stack.size - bl);
                for (j = 0; j < 4; ++j) {
                    if (rflags & (1 << j)) {
                        varr64_resize(&stack, stack.size + bl);
                        lmer_set(rfwd, l_len-1, j);
                        lmer_cpy(&(stack.ele[stack.size-bl]), rfwd, l_len);
                    }
                    if (lflags & (1 << j)) {
                        varr64_resize(&stack, stack.size + bl);
                        lmer_set(lfwd, 0, j);
                        lmer_cpy(&(stack.ele[stack.size-bl]), lfwd, l_len);
                    }
                }
                ljunc.cov = sum;
                ljunc.out = (lflags << 4) | rflags;
                lgraph_ljunc_insert(lg, &ljunc, st_fwd);
                continue;
            }

            if (rseq.size > 0) {
                lmer_cpy(fwd, st_fwd, l_len);
                lmer_lshift(fwd, l_len);
                lmer_set(fwd, l_len-1, rseq.ele[0]);
                u16_p = lcounter_ptr(lc, fwd);
                while (1) {
                    lmer_cpy(rfwd, fwd, l_len);
                    lmer_rshift(fwd, l_len);
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, 0, j);
                        if (lcounter_occ(lc, fwd) >= min_cov)
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        rflags = 0 | (1 << varr8_peek(&rseq));
                        varr8_pop(&rseq);
                        break;
                    }
                    lmer_lshift(fwd, l_len);
                    lmer_lshift(fwd, l_len);
                    lmer_set(fwd, l_len-2, varr8_peek(&rseq));
                    rflags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, l_len-1, j);
                        tmp_u16_p = lcounter_ptr(lc, fwd);
                        if (tmp_u16_p != NULL && *tmp_u16_p >= min_cov) {
                            rflags |= 1 << j;
                            next_u16_p = tmp_u16_p;
                        }
                    }
                    if (rflags == 0) {
                        sum += *u16_p;
                        *u16_p = UINT16_MAX;
                        break;
                    }
                    else if (flag_degree(rflags) > 1) {
                        rflags = 0 | (1 << varr8_peek(&rseq));
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
                lmer_cpy(fwd, st_fwd, l_len);
                lmer_rshift(fwd, l_len);
                lmer_set(fwd, 0, lseq.ele[0]);
                u16_p = lcounter_ptr(lc, fwd);
                while (1) {
                    lmer_cpy(lfwd, fwd, l_len);
                    lmer_lshift(fwd, l_len);
                    tmp_flags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, l_len-1, j);
                        if (lcounter_occ(lc, fwd) >= min_cov)
                            tmp_flags |= 1 << j;
                    }
                    if (flag_degree(tmp_flags) > 1 || *u16_p == UINT16_MAX) {
                        lflags = 0 | (1 << varr8_peek(&lseq));
                        varr8_pop(&lseq);
                        break;
                    }
                    lmer_rshift(fwd, l_len);
                    lmer_rshift(fwd, l_len);
                    lmer_set(fwd, 1, varr8_peek(&lseq));
                    lflags = 0;
                    for (j = 0; j < 4; ++j) {
                        lmer_set(fwd, 0, j);
                        tmp_u16_p = lcounter_ptr(lc, fwd);
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
                        lflags = 0 | (1 << varr8_peek(&lseq));
                        varr8_pop(&lseq);
                        break;
                    }
                    varr8_push(&lseq, flag_base(lflags));
                    sum += *u16_p;
                    *u16_p = UINT16_MAX;
                    u16_p = next_u16_p;
                }
            }

            varr64_resize(&stack, stack.size - bl);

            if (flag_degree(lflags)) {
                varr64_resize(&stack, stack.size + bl);
                lmer_set(lfwd, 0, flag_base(lflags));
                lmer_cpy(&(stack.ele[stack.size-bl]), lfwd, l_len);
            }
            if (flag_degree(rflags)) {
                varr64_resize(&stack, stack.size + bl);
                lmer_set(rfwd, l_len-1, flag_base(rflags));
                lmer_cpy(&(stack.ele[stack.size-bl]), rfwd, l_len);
            }

            lstra.len = lseq.size + rseq.size + 1;
            varr8_resize(&seq, (lstra.len+l_len+2)/4);
            for (j = 0; j < lseq.size; ++j)
                bseq_set(seq.ele, j, lseq.ele[lseq.size-j-1]);
            for (j = 0; j < l_len; ++j)
                bseq_set(seq.ele, lseq.size+j, lmer_get(st_fwd, j));
            for (j = 0; j < rseq.size; ++j)
                bseq_set(seq.ele, lseq.size+l_len+j, rseq.ele[j]);
            lstra.cov = (uint16_t)((double)sum / (double)(lstra.len) + 0.5);
            lstra.out = (lflags << 4) | rflags;

            sp = (lstra_t *)my_malloc(sizeof(lstra_t) + sizeof(uint8_t)*(lstra.len+l_len+2)/4);
            *(sp) = lstra;
            memcpy(sp->seq, seq.ele, (lstra.len+l_len+2)/4 * sizeof(uint8_t));
            lgraph_lstra_insert(lg, sp);
        }
    }

    varr8_destroy(&rseq);
    varr8_destroy(&lseq);
    varr8_destroy(&seq);
    varr64_destroy(&stack);
}

void lgraph_straight_stats(const lgraph_t *lg, uint64_t *len_cutoff, uint64_t *cov_cutoff, double *ave_cov)
{
    uint64_t i;
    uint64_t tmp;
    uint64_t sum;
    uint64_t num;
    lstra_t *sp;
    varr64_t len_distr;
    varr64_t cov_distr;

    varr64_init(&len_distr);
    varr64_init(&cov_distr);
    varr64_resize(&cov_distr, UINT16_MAX);
    memset(cov_distr.ele, 0, cov_distr.size * sizeof(uint64_t));

    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp))
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
    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp) || (sp->cov < *cov_cutoff && sp->len < *len_cutoff))
            continue;
        sum += sp->cov * sp->len;
        num += sp->len;
    }

    *ave_cov = (double)sum / num;

    varr64_destroy(&len_distr);
    varr64_destroy(&cov_distr);
}

void lgraph_delete_erroneous_straight_node(lgraph_t *lg, const uint64_t len_cutoff, const uint64_t cov_cutoff)
{
    uint64_t i;
    uint64_t j;
    uint64_t n_delete;
    uint64_t l_len;
    lmer_t key;
    lstra_t *sp;
    ljunc_t *jp;

    fputs("removing erroneous nodes...\n", stderr);

    l_len = lg->l_len;
    n_delete = 0;

    for (i = 0; i < lg->stable_size; ++i) {
        sp = lg->lsval_table[i];
        if (lstra_is_emp(sp) || lstra_is_del(sp) || sp->cov >= cov_cutoff || sp->len >= len_cutoff)
            continue;

        if (sp->out << 4) {
            lmer_cpy(key, &(lg->lskey_table[i * lg->bl]), l_len);
            lmer_rshift(key, l_len);
            lmer_set(key, 0, flag_base(sp->out >> 4));
            jp = lgraph_ljunc_find(lg, key);
            if (jp != NULL)
                jp->out ^= 1 << bseq_get(sp->seq, l_len - 1);
        }

        if (sp->out & 0xf) {
            for (j = 1; j < l_len; ++j)
                lmer_set(key, j - 1, bseq_get(sp->seq, sp->len + j - 1));
            lmer_set(key, l_len - 1, flag_base(sp->out & 0xf));
            jp = lgraph_ljunc_find(lg, key);
            if (jp != NULL)
                jp->out ^= 1 << (4 + bseq_get(sp->seq, sp->len - 1));
        }

        lgraph_lstra_delete(lg, sp);
        ++n_delete;
    }

    lgraph_join_nodes(lg);

    fprintf(stderr, "NUM_DELETE=%lu\n", n_delete);
}
