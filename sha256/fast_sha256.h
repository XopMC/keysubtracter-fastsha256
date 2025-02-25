#ifndef SHA256_H
#define SHA256_H
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#if defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#elif defined(__amd64__) || defined(_M_AMD64)
#include <immintrin.h>
#endif

#define mod(x,y) ((x)-((x)/(y)*(y)))
#define shr32(x,n) ((x) >> (n))
#define rotl32(n,d) (((n) << (d)) | ((n) >> (32 - (d))))
#define rotl64(n,d) (((n) << (d)) | ((n) >> (64 - (d))))
#define rotr64(n,d) (((n) >> (d)) | ((n) << (64 - (d))))
#define S0(x) (rotl32 ((x), 25u) ^ rotl32 ((x), 14u) ^ shr32 ((x),  3u))
#define S1(x) (rotl32 ((x), 15u) ^ rotl32 ((x), 13u) ^ shr32 ((x), 10u))
#define S2(x) (rotl32 ((x), 30u) ^ rotl32 ((x), 19u) ^ rotl32 ((x), 10u))
#define S3(x) (rotl32 ((x), 26u) ^ rotl32 ((x), 21u) ^ rotl32 ((x),  7u))

static uint32_t SWAP256(uint32_t val) {
    return (rotl32(((val) & (uint32_t)0x00FF00FF), (uint32_t)24U) | rotl32(((val) & (uint32_t)0xFF00FF00), (uint32_t)8U));
}


#if defined(__amd64__) || defined(_M_AMD64)
static const uint32_t K64[64] = {
    0x428A2F98,0x71374491,0xB5C0FBCF,0xE9B5DBA5,
    0x3956C25B,0x59F111F1,0x923F82A4,0xAB1C5ED5,
    0xD807AA98,0x12835B01,0x243185BE,0x550C7DC3,
    0x72BE5D74,0x80DEB1FE,0x9BDC06A7,0xC19BF174,
    0xE49B69C1,0xEFBE4786,0x0FC19DC6,0x240CA1CC,
    0x2DE92C6F,0x4A7484AA,0x5CB0A9DC,0x76F988DA,
    0x983E5152,0xA831C66D,0xB00327C8,0xBF597FC7,
    0xC6E00BF3,0xD5A79147,0x06CA6351,0x14292967,
    0x27B70A85,0x2E1B2138,0x4D2C6DFC,0x53380D13,
    0x650A7354,0x766A0ABB,0x81C2C92E,0x92722C85,
    0xA2BFE8A1,0xA81A664B,0xC24B8B70,0xC76C51A3,
    0xD192E819,0xD6990624,0xF40E3585,0x106AA070,
    0x19A4C116,0x1E376C08,0x2748774C,0x34B0BCB5,
    0x391C0CB3,0x4ED8AA4A,0x5B9CCA4F,0x682E6FF3,
    0x748F82EE,0x78A5636F,0x84C87814,0x8CC70208,
    0x90BEFFFA,0xA4506CEB,0xBEF9A3F7,0xC67178F2
};


static inline void sha256_compress_block(__m128i state0, __m128i state1, const uint8_t* block) {
    __m128i W0 = _mm_loadu_si128((const __m128i*)block);
    __m128i W1 = _mm_loadu_si128((const __m128i*)(block + 16));
    __m128i W2 = _mm_loadu_si128((const __m128i*)(block + 32));
    __m128i W3 = _mm_loadu_si128((const __m128i*)(block + 48));
    const __m128i SHUF_MASK = _mm_set_epi64x(0x0C0D0E0F08090A0BULL, 0x0405060700010203ULL);
    W0 = _mm_shuffle_epi8(W0, SHUF_MASK);
    W1 = _mm_shuffle_epi8(W1, SHUF_MASK);
    W2 = _mm_shuffle_epi8(W2, SHUF_MASK);
    W3 = _mm_shuffle_epi8(W3, SHUF_MASK);

    __m128i ABEF = state0;
    __m128i CDGH = state1;

    __m128i MSGV = W0;
    __m128i MSGTMP0 = MSGV;
    MSGV = _mm_add_epi32(MSGV, _mm_load_si128((const __m128i*)(&K64[0])));
    CDGH = _mm_sha256rnds2_epu32(CDGH, ABEF, MSGV);
    MSGV = _mm_shuffle_epi32(MSGV, 0x0E);
    ABEF = _mm_sha256rnds2_epu32(ABEF, CDGH, MSGV);

    MSGV = W1;
    __m128i MSGTMP1 = MSGV;
    MSGV = _mm_add_epi32(MSGV, _mm_load_si128((const __m128i*)(&K64[4])));
    CDGH = _mm_sha256rnds2_epu32(CDGH, ABEF, MSGV);
    MSGV = _mm_shuffle_epi32(MSGV, 0x0E);
    ABEF = _mm_sha256rnds2_epu32(ABEF, CDGH, MSGV);
    MSGTMP0 = _mm_sha256msg1_epu32(MSGTMP0, MSGTMP1);
    MSGV = W2;

    __m128i MSGTMP2 = MSGV;
    MSGV = _mm_add_epi32(MSGV, _mm_load_si128((const __m128i*)(&K64[8])));
    CDGH = _mm_sha256rnds2_epu32(CDGH, ABEF, MSGV);
    MSGV = _mm_shuffle_epi32(MSGV, 0x0E);
    ABEF = _mm_sha256rnds2_epu32(ABEF, CDGH, MSGV);
    MSGTMP1 = _mm_sha256msg1_epu32(MSGTMP1, MSGTMP2);

    MSGV = W3;
    __m128i MSGTMP3 = MSGV;
    MSGV = _mm_add_epi32(MSGV, _mm_load_si128((const __m128i*)(&K64[12])));
    CDGH = _mm_sha256rnds2_epu32(CDGH, ABEF, MSGV);
    MSGTMP0 = _mm_add_epi32(MSGTMP0, _mm_alignr_epi8(MSGTMP3, MSGTMP2, 4));
    MSGTMP0 = _mm_sha256msg2_epu32(MSGTMP0, MSGTMP3);
    MSGV = _mm_shuffle_epi32(MSGV, 0x0E);
    ABEF = _mm_sha256rnds2_epu32(ABEF, CDGH, MSGV);
    MSGTMP2 = _mm_sha256msg1_epu32(MSGTMP2, MSGTMP3);

 
#define SHA256ROUND(msgv, msgtmp0, msgtmp1, msgtmp2, msgtmp3, s0, s1, k) \
    do { \
        msgv = msgtmp0; \
        msgv = _mm_add_epi32(msgv, _mm_load_si128((const __m128i*)(k))); \
        s1 = _mm_sha256rnds2_epu32(s1, s0, msgv); \
        msgtmp1 = _mm_add_epi32(msgtmp1, _mm_alignr_epi8(msgtmp0, msgtmp3, 4)); \
        msgtmp1 = _mm_sha256msg2_epu32(msgtmp1, msgtmp0); \
        msgv = _mm_shuffle_epi32(msgv, 0x0E); \
        s0 = _mm_sha256rnds2_epu32(s0, s1, msgv); \
        msgtmp3 = _mm_sha256msg1_epu32(msgtmp3, msgtmp0); \
    } while(0)

    SHA256ROUND(MSGV, MSGTMP0, MSGTMP1, MSGTMP2, MSGTMP3, ABEF, CDGH, &K64[16]);
    SHA256ROUND(MSGV, MSGTMP1, MSGTMP2, MSGTMP3, MSGTMP0, ABEF, CDGH, &K64[20]);
    SHA256ROUND(MSGV, MSGTMP2, MSGTMP3, MSGTMP0, MSGTMP1, ABEF, CDGH, &K64[24]);
    SHA256ROUND(MSGV, MSGTMP3, MSGTMP0, MSGTMP1, MSGTMP2, ABEF, CDGH, &K64[28]);

    SHA256ROUND(MSGV, MSGTMP0, MSGTMP1, MSGTMP2, MSGTMP3, ABEF, CDGH, &K64[32]);
    SHA256ROUND(MSGV, MSGTMP1, MSGTMP2, MSGTMP3, MSGTMP0, ABEF, CDGH, &K64[36]);
    SHA256ROUND(MSGV, MSGTMP2, MSGTMP3, MSGTMP0, MSGTMP1, ABEF, CDGH, &K64[40]);
    SHA256ROUND(MSGV, MSGTMP3, MSGTMP0, MSGTMP1, MSGTMP2, ABEF, CDGH, &K64[44]);

    SHA256ROUND(MSGV, MSGTMP0, MSGTMP1, MSGTMP2, MSGTMP3, ABEF, CDGH, &K64[48]);
    SHA256ROUND(MSGV, MSGTMP1, MSGTMP2, MSGTMP3, MSGTMP0, ABEF, CDGH, &K64[52]);
    SHA256ROUND(MSGV, MSGTMP2, MSGTMP3, MSGTMP0, MSGTMP1, ABEF, CDGH, &K64[56]);
    MSGV = MSGTMP3;
    MSGV = _mm_add_epi32(MSGV, _mm_load_si128((const __m128i*)(&K64[60])));
    CDGH = _mm_sha256rnds2_epu32(CDGH, ABEF, MSGV);
    MSGV = _mm_shuffle_epi32(MSGV, 0x0E);
    ABEF = _mm_sha256rnds2_epu32(ABEF, CDGH, MSGV);
#undef SHA256ROUND

    state0 = _mm_add_epi32(ABEF, state0);
    state1 = _mm_add_epi32(CDGH, state1);
}

static inline void sha256_x64(const uint8_t* data, size_t len, uint8_t hash[32]) {
    const __m128i ABEF_INIT = _mm_set_epi64x(0x6A09E667BB67AE85ULL, 0x510E527F9B05688CULL);
    const __m128i CDGH_INIT = _mm_set_epi64x(0x3C6EF372A54FF53AULL, 0x1F83D9AB5BE0CD19ULL);
    __m128i state0 = ABEF_INIT;
    __m128i state1 = CDGH_INIT;

    size_t num_blocks = len / 64;
    for (size_t i = 0; i < num_blocks; i++) {
        sha256_compress_block(state0, state1, data + i * 64);
    }

    uint8_t block[128] = { 0 };
    size_t rem = len % 64;
    memcpy(block, data + num_blocks * 64, rem);
    block[rem] = 0x80;
    uint64_t bit_len = len * 8;
    if (rem < 56) {
        block[56] = (uint8_t)(bit_len >> 56);
        block[57] = (uint8_t)(bit_len >> 48);
        block[58] = (uint8_t)(bit_len >> 40);
        block[59] = (uint8_t)(bit_len >> 32);
        block[60] = (uint8_t)(bit_len >> 24);
        block[61] = (uint8_t)(bit_len >> 16);
        block[62] = (uint8_t)(bit_len >> 8);
        block[63] = (uint8_t)(bit_len);
        sha256_compress_block(state0, state1, block);
    }
    else {
        sha256_compress_block(state0, state1, block);
        memset(block + 64, 0, 64);
        block[64 + 56] = (uint8_t)(bit_len >> 56);
        block[64 + 57] = (uint8_t)(bit_len >> 48);
        block[64 + 58] = (uint8_t)(bit_len >> 40);
        block[64 + 59] = (uint8_t)(bit_len >> 32);
        block[64 + 60] = (uint8_t)(bit_len >> 24);
        block[64 + 61] = (uint8_t)(bit_len >> 16);
        block[64 + 62] = (uint8_t)(bit_len >> 8);
        block[64 + 63] = (uint8_t)(bit_len);
        sha256_compress_block(state0, state1, block + 64);
    }

    state0 = _mm_shuffle_epi32(state0, 0x1B);
    state1 = _mm_shuffle_epi32(state1, 0xB1);
    __m128i hash0 = _mm_blend_epi16(state0, state1, 0xF0);
    __m128i hash1 = _mm_alignr_epi8(state1, state0, 8);

    const __m128i rev32_mask = _mm_set_epi8(12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3);
    hash0 = _mm_shuffle_epi8(hash0, rev32_mask);
    hash1 = _mm_shuffle_epi8(hash1, rev32_mask);
    _mm_storeu_si128((__m128i*)hash, hash0);
    _mm_storeu_si128((__m128i*)(hash + 16), hash1);
}


#if defined(__AVX512F__)
static inline void sha256_x64_avx512(const uint8_t* data, size_t len, uint8_t hash[32]) {
    const __m128i ABEF_INIT = _mm_set_epi64x(0x6A09E667BB67AE85ULL, 0x510E527F9B05688CULL);
    const __m128i CDGH_INIT = _mm_set_epi64x(0x3C6EF372A54FF53AULL, 0x1F83D9AB5BE0CD19ULL);
    __m128i state0 = ABEF_INIT;
    __m128i state1 = CDGH_INIT;
    size_t num_blocks = len / 64;
    alignas(64) uint8_t temp[64];

    for (size_t i = 0; i < num_blocks; i++) {
        __m512i block512 = _mm512_loadu_si512((const void*)(data + i * 64));
        _mm512_storeu_si512((void*)temp, block512);
        sha256_compress_block(state0, state1, temp);
    }

    uint8_t block[128] = { 0 };
    size_t rem = len % 64;
    memcpy(block, data + num_blocks * 64, rem);
    block[rem] = 0x80;
    uint64_t bit_len = len * 8;
    if (rem < 56) {
        block[56] = (uint8_t)(bit_len >> 56);
        block[57] = (uint8_t)(bit_len >> 48);
        block[58] = (uint8_t)(bit_len >> 40);
        block[59] = (uint8_t)(bit_len >> 32);
        block[60] = (uint8_t)(bit_len >> 24);
        block[61] = (uint8_t)(bit_len >> 16);
        block[62] = (uint8_t)(bit_len >> 8);
        block[63] = (uint8_t)(bit_len);
        sha256_compress_block(state0, state1, block);
    }
    else {
        sha256_compress_block(state0, state1, block);
        memset(block + 64, 0, 64);
        block[64 + 56] = (uint8_t)(bit_len >> 56);
        block[64 + 57] = (uint8_t)(bit_len >> 48);
        block[64 + 58] = (uint8_t)(bit_len >> 40);
        block[64 + 59] = (uint8_t)(bit_len >> 32);
        block[64 + 60] = (uint8_t)(bit_len >> 24);
        block[64 + 61] = (uint8_t)(bit_len >> 16);
        block[64 + 62] = (uint8_t)(bit_len >> 8);
        block[64 + 63] = (uint8_t)(bit_len);
        sha256_compress_block(state0, state1, block + 64);
    }
    state0 = _mm_shuffle_epi32(state0, 0x1B);
    state1 = _mm_shuffle_epi32(state1, 0xB1);
    __m128i hash0 = _mm_blend_epi16(state0, state1, 0xF0);
    __m128i hash1 = _mm_alignr_epi8(state1, state0, 8);
    const __m128i rev32_mask = _mm_set_epi8(12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3);
    hash0 = _mm_shuffle_epi8(hash0, rev32_mask);
    hash1 = _mm_shuffle_epi8(hash1, rev32_mask);
    _mm_storeu_si128((__m128i*)hash, hash0);
    _mm_storeu_si128((__m128i*)(hash + 16), hash1);
}

void old_electrum_sha256_avx512(const uint8_t* data, size_t len, uint64_t num_iters, uint8_t* hash) {
    uint8_t prepared_input[512] = { 0 };
    memcpy(prepared_input, data, len);
    memcpy(prepared_input + len, data, len);


    sha256_x64_avx512(prepared_input, len + len, hash);
    memcpy(hash + 32, prepared_input, len);
    for (uint64_t i = 1; i < num_iters; i++) {
        sha256_x64_avx512(hash, 32 + len, hash);

    }
}

#endif
#if defined(__AVX2__)
static inline void sha256_x64_avx2(const uint8_t* data, size_t len, uint8_t hash[32]) {
    const __m128i ABEF_INIT = _mm_set_epi64x(0x6A09E667BB67AE85ULL, 0x510E527F9B05688CULL);
    const __m128i CDGH_INIT = _mm_set_epi64x(0x3C6EF372A54FF53AULL, 0x1F83D9AB5BE0CD19ULL);
    __m128i state0 = ABEF_INIT;
    __m128i state1 = CDGH_INIT;
    size_t num_blocks = len / 64;
    alignas(32) uint8_t temp[64];

    for (size_t i = 0; i < num_blocks; i++) {
        __m256i part0 = _mm256_loadu_si256((const __m256i*)(data + i * 64));
        __m256i part1 = _mm256_loadu_si256((const __m256i*)(data + i * 64 + 32));
        _mm256_store_si256((__m256i*)temp, part0);
        _mm256_store_si256((__m256i*)(temp + 32), part1);
        sha256_compress_block(state0, state1, temp);
    }

    uint8_t block[128] = { 0 };
    size_t rem = len % 64;
    memcpy(block, data + num_blocks * 64, rem);
    block[rem] = 0x80;
    uint64_t bit_len = len * 8;
    if (rem < 56) {
        block[56] = (uint8_t)(bit_len >> 56);
        block[57] = (uint8_t)(bit_len >> 48);
        block[58] = (uint8_t)(bit_len >> 40);
        block[59] = (uint8_t)(bit_len >> 32);
        block[60] = (uint8_t)(bit_len >> 24);
        block[61] = (uint8_t)(bit_len >> 16);
        block[62] = (uint8_t)(bit_len >> 8);
        block[63] = (uint8_t)(bit_len);
        sha256_compress_block(state0, state1, block);
    }
    else {
        sha256_compress_block(state0, state1, block);
        memset(block + 64, 0, 64);
        block[64 + 56] = (uint8_t)(bit_len >> 56);
        block[64 + 57] = (uint8_t)(bit_len >> 48);
        block[64 + 58] = (uint8_t)(bit_len >> 40);
        block[64 + 59] = (uint8_t)(bit_len >> 32);
        block[64 + 60] = (uint8_t)(bit_len >> 24);
        block[64 + 61] = (uint8_t)(bit_len >> 16);
        block[64 + 62] = (uint8_t)(bit_len >> 8);
        block[64 + 63] = (uint8_t)(bit_len);
        sha256_compress_block(state0, state1, block + 64);
    }
    state0 = _mm_shuffle_epi32(state0, 0x1B);
    state1 = _mm_shuffle_epi32(state1, 0xB1);
    __m128i hash0 = _mm_blend_epi16(state0, state1, 0xF0);
    __m128i hash1 = _mm_alignr_epi8(state1, state0, 8);
    const __m128i rev32_mask = _mm_set_epi8(12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3);
    hash0 = _mm_shuffle_epi8(hash0, rev32_mask);
    hash1 = _mm_shuffle_epi8(hash1, rev32_mask);
    _mm_storeu_si128((__m128i*)hash, hash0);
    _mm_storeu_si128((__m128i*)(hash + 16), hash1);
}

}
#endif





#elif defined(__aarch64__) || defined(_M_ARM64)
// ========================= ARM =============================

static const uint32_t K64[64] = {
  0x428A2F98,0x71374491,0xB5C0FBCF,0xE9B5DBA5,
  0x3956C25B,0x59F111F1,0x923F82A4,0xAB1C5ED5,
  0xD807AA98,0x12835B01,0x243185BE,0x550C7DC3,
  0x72BE5D74,0x80DEB1FE,0x9BDC06A7,0xC19BF174,
  0xE49B69C1,0xEFBE4786,0x0FC19DC6,0x240CA1CC,
  0x2DE92C6F,0x4A7484AA,0x5CB0A9DC,0x76F988DA,
  0x983E5152,0xA831C66D,0xB00327C8,0xBF597FC7,
  0xC6E00BF3,0xD5A79147,0x06CA6351,0x14292967,
  0x27B70A85,0x2E1B2138,0x4D2C6DFC,0x53380D13,
  0x650A7354,0x766A0ABB,0x81C2C92E,0x92722C85,
  0xA2BFE8A1,0xA81A664B,0xC24B8B70,0xC76C51A3,
  0xD192E819,0xD6990624,0xF40E3585,0x106AA070,
  0x19A4C116,0x1E376C08,0x2748774C,0x34B0BCB5,
  0x391C0CB3,0x4ED8AA4A,0x5B9CCA4F,0x682E6FF3,
  0x748F82EE,0x78A5636F,0x84C87814,0x8CC70208,
  0x90BEFFFA,0xA4506CEB,0xBEF9A3F7,0xC67178F2
};

static inline void sha256_compress_block_arm(uint32x4_t* STATE0, uint32x4_t* STATE1, const uint8_t* block) {
   
    uint32x4_t W0 = vld1q_u32((const uint32_t*)block);
    uint32x4_t W1 = vld1q_u32((const uint32_t*)(block + 16));
    uint32x4_t W2 = vld1q_u32((const uint32_t*)(block + 32));
    uint32x4_t W3 = vld1q_u32((const uint32_t*)(block + 48));
   
    W0 = vrev32q_u32(W0);
    W1 = vrev32q_u32(W1);
    W2 = vrev32q_u32(W2);
    W3 = vrev32q_u32(W3);


    uint32x4_t A = *STATE0;
    uint32x4_t B = *STATE1;
    uint32x4_t STATEV;
    uint32x4_t MSGV, MSGTMP0, MSGTMP1, MSGTMP2, MSGTMP3;

  
    MSGV = vaddq_u32(W0, vld1q_u32(&K64[0]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);

    MSGTMP0 = vsha256su0q_u32(W0, W1);

    
    MSGV = vaddq_u32(W1, vld1q_u32(&K64[4]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);
    
    static const uint32_t hpad0cache[4] = { 0x80000000,0x00000000,0x00000000,0x00000000 };
    static const uint32_t hpad1cache[4] = { 0x00000000,0x00000000,0x00000000,0x00000100 };
    const uint32x4_t HPAD0_CACHE = vld1q_u32(hpad0cache);
    const uint32x4_t HPAD1_CACHE = vld1q_u32(hpad1cache);
    MSGTMP0 = vsha256su1q_u32(MSGTMP0, HPAD0_CACHE, HPAD1_CACHE);
    MSGTMP1 = vsha256su0q_u32(W1, HPAD0_CACHE);

    
    MSGV = vaddq_u32(HPAD0_CACHE, vld1q_u32(&K64[8]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);
    MSGTMP1 = vsha256su1q_u32(MSGTMP1, HPAD1_CACHE, MSGTMP0);
    MSGTMP2 = HPAD0_CACHE;

    
    MSGV = vaddq_u32(HPAD1_CACHE, vld1q_u32(&K64[12]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);
    MSGTMP2 = vsha256su1q_u32(MSGTMP2, MSGTMP0, MSGTMP1);
    MSGTMP3 = vsha256su0q_u32(HPAD1_CACHE, MSGTMP0);

#define SHA256ROUND_ARM(msgv, msgtmp0, msgtmp1, msgtmp2, msgtmp3, statev, state0, state1, kvalue) \
    do { \
        msgv = vaddq_u32(msgtmp0, vld1q_u32(kvalue)); \
        statev = state0; \
        state0 = vsha256hq_u32(state0, state1, msgv); \
        state1 = vsha256h2q_u32(state1, statev, msgv); \
        msgtmp3 = vsha256su1q_u32(msgtmp3, msgtmp1, msgtmp2); \
        msgtmp0 = vsha256su0q_u32(msgtmp0, msgtmp1); \
    } while(0)

   
    SHA256ROUND_ARM(MSGV, MSGTMP0, MSGTMP1, MSGTMP2, MSGTMP3, A, A, B, &K64[16]);
    SHA256ROUND_ARM(MSGV, MSGTMP1, MSGTMP2, MSGTMP3, MSGTMP0, A, A, B, &K64[20]);
    SHA256ROUND_ARM(MSGV, MSGTMP2, MSGTMP3, MSGTMP0, MSGTMP1, A, A, B, &K64[24]);
    SHA256ROUND_ARM(MSGV, MSGTMP3, MSGTMP0, MSGTMP1, MSGTMP2, A, A, B, &K64[28]);


    SHA256ROUND_ARM(MSGV, MSGTMP0, MSGTMP1, MSGTMP2, MSGTMP3, A, A, B, &K64[32]);
    SHA256ROUND_ARM(MSGV, MSGTMP1, MSGTMP2, MSGTMP3, MSGTMP0, A, A, B, &K64[36]);
    SHA256ROUND_ARM(MSGV, MSGTMP2, MSGTMP3, MSGTMP0, MSGTMP1, A, A, B, &K64[40]);
    SHA256ROUND_ARM(MSGV, MSGTMP3, MSGTMP0, MSGTMP1, MSGTMP2, A, A, B, &K64[44]);

 
    MSGV = vaddq_u32(MSGTMP0, vld1q_u32(&K64[48]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);
    MSGTMP3 = vsha256su1q_u32(MSGTMP3, MSGTMP1, MSGTMP2);

   
    MSGV = vaddq_u32(MSGTMP1, vld1q_u32(&K64[52]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);

    
    MSGV = vaddq_u32(MSGTMP2, vld1q_u32(&K64[56]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);

   
    MSGV = vaddq_u32(MSGTMP3, vld1q_u32(&K64[60]));
    STATEV = A;
    A = vsha256hq_u32(A, B, MSGV);
    B = vsha256h2q_u32(B, STATEV, MSGV);

   
    uint32x4_t ABCD_INIT = vld1q_u32((const uint32_t[]) { 0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A });
    uint32x4_t EFGH_INIT = vld1q_u32((const uint32_t[]) { 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19 });
    *STATE0 = vaddq_u32(A, ABCD_INIT);
    *STATE1 = vaddq_u32(B, EFGH_INIT);
}

static inline void sha256_arm_hash(const uint8_t* data, size_t len, uint8_t hash[32]) {
    uint32x4_t state0 = vld1q_u32((const uint32_t[]) { 0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A });
    uint32x4_t state1 = vld1q_u32((const uint32_t[]) { 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19 });

    size_t num_blocks = len / 64;
    for (size_t i = 0; i < num_blocks; i++) {
        sha256_compress_block_arm(&state0, &state1, data + i * 64);
    }

    uint8_t block[128] = { 0 };
    size_t rem = len % 64;
    memcpy(block, data + num_blocks * 64, rem);
    block[rem] = 0x80;
    uint64_t bit_len = len * 8;
    if (rem < 56) {
        block[56] = (uint8_t)(bit_len >> 56);
        block[57] = (uint8_t)(bit_len >> 48);
        block[58] = (uint8_t)(bit_len >> 40);
        block[59] = (uint8_t)(bit_len >> 32);
        block[60] = (uint8_t)(bit_len >> 24);
        block[61] = (uint8_t)(bit_len >> 16);
        block[62] = (uint8_t)(bit_len >> 8);
        block[63] = (uint8_t)(bit_len);
        sha256_compress_block_arm(&state0, &state1, block);
    }
    else {
        sha256_compress_block_arm(&state0, &state1, block);
        memset(block + 64, 0, 64);
        block[64 + 56] = (uint8_t)(bit_len >> 56);
        block[64 + 57] = (uint8_t)(bit_len >> 48);
        block[64 + 58] = (uint8_t)(bit_len >> 40);
        block[64 + 59] = (uint8_t)(bit_len >> 32);
        block[64 + 60] = (uint8_t)(bit_len >> 24);
        block[64 + 61] = (uint8_t)(bit_len >> 16);
        block[64 + 62] = (uint8_t)(bit_len >> 8);
        block[64 + 63] = (uint8_t)(bit_len);
        sha256_compress_block_arm(&state0, &state1, block + 64);
    }

    
    state0 = vrev32q_u32(state0);
    state1 = vrev32q_u32(state1);
    vst1q_u32((uint32_t*)hash, state0);
    vst1q_u32((uint32_t*)(hash + 16), state1);
}


#else
#error "Unsupported architecture"
#endif


#endif