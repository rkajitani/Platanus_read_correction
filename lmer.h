#ifndef LMER_H
#define LMER_H

#include "common.h"

typedef uint64_t lmer_t[MAX_READ_LEN / 32];

static inline void lmer_lshift(uint64_t *lmer, const uint64_t len)
{
    uint64_t i;

    for (i = 0; i < (len-1)/32; ++i)
        lmer[i] = (lmer[i] >> 2) | (lmer[i+1] << 62);
    if (len % 32)
        lmer[i] &= ~(~0ull << (len%32*2));
    lmer[i] = lmer[i] >> 2;
}

static inline void lmer_rshift(uint64_t *lmer, uint64_t len)
{
    for (len = (len-1)/32; len > 0; --len)
        lmer[len] = (lmer[len] << 2) | (lmer[len-1] >> 62);
    lmer[0] = lmer[0] << 2;
}

static inline void lmer_set(uint64_t *lmer, const uint64_t pos, const uint64_t val)
{
//    lmer[pos >> 5] = lmer[pos >> 5] & ~(3ull << ((pos & 31) << 1)) | (val << ((pos & 31) << 1));
    lmer[pos/32] = (lmer[pos/32] & ~(3ull << ((pos%32)*2))) | (val << ((pos%32)*2));
}

static inline uint64_t lmer_get(uint64_t *lmer, const uint64_t pos)
{
//    return (lmer[pos >> 5] >> ((pos & 31) << 1)) & 3;
    return (lmer[pos/32] >> ((pos%32)*2)) & 3;
}

static inline int lmer_cmp(uint64_t *lmer1, uint64_t *lmer2, uint64_t len)
{
    if (len % 32) {
        lmer1[(len-1)/32] &= ~(~0ull << (len%32*2));
        lmer2[(len-1)/32] &= ~(~0ull << (len%32*2));
    }
    for (len = (len-1)/32; len > 0; --len)
        if (lmer1[len] != lmer2[len])
            return ((int64_t)lmer1[len] < (int64_t)lmer2[len] ? -1 : 1);
    if (lmer1[len] != lmer2[len])
        return ((int64_t)lmer1[len] < (int64_t)lmer2[len] ? -1 : 1);
    return 0;
}

static inline void lmer_cpy(uint64_t *dst, const uint64_t *src, uint64_t len)
{
    for (len = (len-1)/32; len > 0; --len)
        dst[len] = src[len];
    dst[0] = src[0];
}

static inline uint64_t lmer_hash(uint64_t *key, uint64_t len, const uint64_t index_len)
{
    uint64_t val = 0;

    if (len % 32)
        key[(len-1)/32] &= ~(~0ull << (len%32*2));
    for (len = (len-1)/32; len > 0; --len)
        val += (key[len] + (key[len] >> index_len) + (key[len] >> 2*index_len));
    return val + (key[0] + (key[0] >> index_len) + (key[0] >> 2*index_len));
}

static inline uint64_t lmer_rehash(uint64_t *key, uint64_t len, const uint64_t index_len)
{
    uint64_t val = 0;

    if (len % 32)
        key[(len-1)/32] &= ~(~0ull << (len%32*2));
    for (len = (len-1)/32; len > 0; --len)
        val += (~key[len] ^ (key[len] >> index_len) ^ (key[len] >> 2*index_len));
    return (val + (~key[0] ^ (key[0] >> index_len) ^ (key[0] >> 2*index_len))) | 1ull;
}

#endif
