#ifndef C25519_SHA512_H
#define C25519_SHA512_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SHA512_BLOCK_SIZE 128
#define SHA512_DIGEST_SIZE 64

typedef struct sha512_ctx {
    uint64_t length;
    uint64_t state[8];
    size_t curlen;
    uint8_t buf[SHA512_BLOCK_SIZE];
} sha512_ctx;

static inline uint64_t sha512_ror64(uint64_t v, unsigned int n)
{
    return (v >> n) | (v << (64U - n));
}

static inline uint64_t sha512_load64_be(const uint8_t* in)
{
    return ((uint64_t)in[0] << 56) | ((uint64_t)in[1] << 48) | ((uint64_t)in[2] << 40) | ((uint64_t)in[3] << 32) | ((uint64_t)in[4] << 24) | ((uint64_t)in[5] << 16) | ((uint64_t)in[6] << 8) | ((uint64_t)in[7]);
}

static inline void sha512_store64_be(uint8_t out[8], uint64_t v)
{
    out[0] = (uint8_t)(v >> 56);
    out[1] = (uint8_t)(v >> 48);
    out[2] = (uint8_t)(v >> 40);
    out[3] = (uint8_t)(v >> 32);
    out[4] = (uint8_t)(v >> 24);
    out[5] = (uint8_t)(v >> 16);
    out[6] = (uint8_t)(v >> 8);
    out[7] = (uint8_t)(v);
}

static inline uint64_t sha512_ch(uint64_t x, uint64_t y, uint64_t z)
{
    return z ^ (x & (y ^ z));
}

static inline uint64_t sha512_maj(uint64_t x, uint64_t y, uint64_t z)
{
    return ((x | y) & z) | (x & y);
}

static inline uint64_t sha512_sigma0(uint64_t x)
{
    return sha512_ror64(x, 28) ^ sha512_ror64(x, 34) ^ sha512_ror64(x, 39);
}

static inline uint64_t sha512_sigma1(uint64_t x)
{
    return sha512_ror64(x, 14) ^ sha512_ror64(x, 18) ^ sha512_ror64(x, 41);
}

static inline uint64_t sha512_gamma0(uint64_t x)
{
    return sha512_ror64(x, 1) ^ sha512_ror64(x, 8) ^ (x >> 7);
}

static inline uint64_t sha512_gamma1(uint64_t x)
{
    return sha512_ror64(x, 19) ^ sha512_ror64(x, 61) ^ (x >> 6);
}

static const uint64_t sha512_k[80] = {
    0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL,
    0xe9b5dba58189dbbcULL, 0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL,
    0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL, 0xd807aa98a3030242ULL,
    0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
    0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL,
    0xc19bf174cf692694ULL, 0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL,
    0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL, 0x2de92c6f592b0275ULL,
    0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
    0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL,
    0xbf597fc7beef0ee4ULL, 0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL,
    0x06ca6351e003826fULL, 0x142929670a0e6e70ULL, 0x27b70a8546d22ffcULL,
    0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL,
    0x650a73548baf63deULL, 0x766a0abb3c77b2a8ULL, 0x81c2c92e47edaee6ULL,
    0x92722c851482353bULL, 0xa2bfe8a14cf10364ULL, 0xa81a664bbc423001ULL,
    0xc24b8b70d0f89791ULL, 0xc76c51a30654be30ULL, 0xd192e819d6ef5218ULL,
    0xd69906245565a910ULL, 0xf40e35855771202aULL, 0x106aa07032bbd1b8ULL,
    0x19a4c116b8d2d0c8ULL, 0x1e376c085141ab53ULL, 0x2748774cdf8eeb99ULL,
    0x34b0bcb5e19b48a8ULL, 0x391c0cb3c5c95a63ULL, 0x4ed8aa4ae3418acbULL,
    0x5b9cca4f7763e373ULL, 0x682e6ff3d6b2b8a3ULL, 0x748f82ee5defb2fcULL,
    0x78a5636f43172f60ULL, 0x84c87814a1f0ab72ULL, 0x8cc702081a6439ecULL,
    0x90befffa23631e28ULL, 0xa4506cebde82bde9ULL, 0xbef9a3f7b2c67915ULL,
    0xc67178f2e372532bULL, 0xca273eceea26619cULL, 0xd186b8c721c0c207ULL,
    0xeada7dd6cde0eb1eULL, 0xf57d4f7fee6ed178ULL, 0x06f067aa72176fbaULL,
    0x0a637dc5a2c898a6ULL, 0x113f9804bef90daeULL, 0x1b710b35131c471bULL,
    0x28db77f523047d84ULL, 0x32caab7b40c72493ULL, 0x3c9ebe0a15c9bebcULL,
    0x431d67c49c100d4cULL, 0x4cc5d4becb3e42b6ULL, 0x597f299cfc657e2aULL,
    0x5fcb6fab3ad6faecULL, 0x6c44198c4a475817ULL
};

static inline void sha512_transform(sha512_ctx* ctx,
    const uint8_t block[SHA512_BLOCK_SIZE])
{
    uint64_t w[16];
    for (size_t i = 0; i < 16; i++) {
        w[i] = sha512_load64_be(block + (i * 8));
    }

    uint64_t a = ctx->state[0];
    uint64_t b = ctx->state[1];
    uint64_t c = ctx->state[2];
    uint64_t d = ctx->state[3];
    uint64_t e = ctx->state[4];
    uint64_t f = ctx->state[5];
    uint64_t g = ctx->state[6];
    uint64_t h = ctx->state[7];

#define SHA512_ROUND(i, wi) \
    do { \
        uint64_t t0 = h + sha512_sigma1(e) + sha512_ch(e, f, g) \
            + sha512_k[(i)] + (wi); \
        uint64_t t1 = sha512_sigma0(a) + sha512_maj(a, b, c); \
        h = g; \
        g = f; \
        f = e; \
        e = d + t0; \
        d = c; \
        c = b; \
        b = a; \
        a = t0 + t1; \
    } while (0)
#define SHA512_SCHED(i) \
    (w[(i) & 15U] += sha512_gamma1(w[((i) + 14U) & 15U]) \
            + w[((i) + 9U) & 15U] \
            + sha512_gamma0(w[((i) + 1U) & 15U]))

    SHA512_ROUND(0, w[0]);
    SHA512_ROUND(1, w[1]);
    SHA512_ROUND(2, w[2]);
    SHA512_ROUND(3, w[3]);
    SHA512_ROUND(4, w[4]);
    SHA512_ROUND(5, w[5]);
    SHA512_ROUND(6, w[6]);
    SHA512_ROUND(7, w[7]);
    SHA512_ROUND(8, w[8]);
    SHA512_ROUND(9, w[9]);
    SHA512_ROUND(10, w[10]);
    SHA512_ROUND(11, w[11]);
    SHA512_ROUND(12, w[12]);
    SHA512_ROUND(13, w[13]);
    SHA512_ROUND(14, w[14]);
    SHA512_ROUND(15, w[15]);
    SHA512_ROUND(16, SHA512_SCHED(16));
    SHA512_ROUND(17, SHA512_SCHED(17));
    SHA512_ROUND(18, SHA512_SCHED(18));
    SHA512_ROUND(19, SHA512_SCHED(19));
    SHA512_ROUND(20, SHA512_SCHED(20));
    SHA512_ROUND(21, SHA512_SCHED(21));
    SHA512_ROUND(22, SHA512_SCHED(22));
    SHA512_ROUND(23, SHA512_SCHED(23));
    SHA512_ROUND(24, SHA512_SCHED(24));
    SHA512_ROUND(25, SHA512_SCHED(25));
    SHA512_ROUND(26, SHA512_SCHED(26));
    SHA512_ROUND(27, SHA512_SCHED(27));
    SHA512_ROUND(28, SHA512_SCHED(28));
    SHA512_ROUND(29, SHA512_SCHED(29));
    SHA512_ROUND(30, SHA512_SCHED(30));
    SHA512_ROUND(31, SHA512_SCHED(31));
    SHA512_ROUND(32, SHA512_SCHED(32));
    SHA512_ROUND(33, SHA512_SCHED(33));
    SHA512_ROUND(34, SHA512_SCHED(34));
    SHA512_ROUND(35, SHA512_SCHED(35));
    SHA512_ROUND(36, SHA512_SCHED(36));
    SHA512_ROUND(37, SHA512_SCHED(37));
    SHA512_ROUND(38, SHA512_SCHED(38));
    SHA512_ROUND(39, SHA512_SCHED(39));
    SHA512_ROUND(40, SHA512_SCHED(40));
    SHA512_ROUND(41, SHA512_SCHED(41));
    SHA512_ROUND(42, SHA512_SCHED(42));
    SHA512_ROUND(43, SHA512_SCHED(43));
    SHA512_ROUND(44, SHA512_SCHED(44));
    SHA512_ROUND(45, SHA512_SCHED(45));
    SHA512_ROUND(46, SHA512_SCHED(46));
    SHA512_ROUND(47, SHA512_SCHED(47));
    SHA512_ROUND(48, SHA512_SCHED(48));
    SHA512_ROUND(49, SHA512_SCHED(49));
    SHA512_ROUND(50, SHA512_SCHED(50));
    SHA512_ROUND(51, SHA512_SCHED(51));
    SHA512_ROUND(52, SHA512_SCHED(52));
    SHA512_ROUND(53, SHA512_SCHED(53));
    SHA512_ROUND(54, SHA512_SCHED(54));
    SHA512_ROUND(55, SHA512_SCHED(55));
    SHA512_ROUND(56, SHA512_SCHED(56));
    SHA512_ROUND(57, SHA512_SCHED(57));
    SHA512_ROUND(58, SHA512_SCHED(58));
    SHA512_ROUND(59, SHA512_SCHED(59));
    SHA512_ROUND(60, SHA512_SCHED(60));
    SHA512_ROUND(61, SHA512_SCHED(61));
    SHA512_ROUND(62, SHA512_SCHED(62));
    SHA512_ROUND(63, SHA512_SCHED(63));
    SHA512_ROUND(64, SHA512_SCHED(64));
    SHA512_ROUND(65, SHA512_SCHED(65));
    SHA512_ROUND(66, SHA512_SCHED(66));
    SHA512_ROUND(67, SHA512_SCHED(67));
    SHA512_ROUND(68, SHA512_SCHED(68));
    SHA512_ROUND(69, SHA512_SCHED(69));
    SHA512_ROUND(70, SHA512_SCHED(70));
    SHA512_ROUND(71, SHA512_SCHED(71));
    SHA512_ROUND(72, SHA512_SCHED(72));
    SHA512_ROUND(73, SHA512_SCHED(73));
    SHA512_ROUND(74, SHA512_SCHED(74));
    SHA512_ROUND(75, SHA512_SCHED(75));
    SHA512_ROUND(76, SHA512_SCHED(76));
    SHA512_ROUND(77, SHA512_SCHED(77));
    SHA512_ROUND(78, SHA512_SCHED(78));
    SHA512_ROUND(79, SHA512_SCHED(79));

#undef SHA512_SCHED
#undef SHA512_ROUND

    ctx->state[0] += a;
    ctx->state[1] += b;
    ctx->state[2] += c;
    ctx->state[3] += d;
    ctx->state[4] += e;
    ctx->state[5] += f;
    ctx->state[6] += g;
    ctx->state[7] += h;
}

static inline void sha512_init(sha512_ctx* ctx)
{
    ctx->curlen = 0;
    ctx->length = 0;
    ctx->state[0] = 0x6a09e667f3bcc908ULL;
    ctx->state[1] = 0xbb67ae8584caa73bULL;
    ctx->state[2] = 0x3c6ef372fe94f82bULL;
    ctx->state[3] = 0xa54ff53a5f1d36f1ULL;
    ctx->state[4] = 0x510e527fade682d1ULL;
    ctx->state[5] = 0x9b05688c2b3e6c1fULL;
    ctx->state[6] = 0x1f83d9abfb41bd6bULL;
    ctx->state[7] = 0x5be0cd19137e2179ULL;
}

static inline void sha512_update(sha512_ctx* ctx, const void* data,
    size_t len)
{
    if (ctx->curlen > sizeof(ctx->buf)) {
        return;
    }

    const uint8_t* in = (const uint8_t*)data;
    while (len > 0) {
        if (ctx->curlen == 0 && len >= SHA512_BLOCK_SIZE) {
            sha512_transform(ctx, in);
            ctx->length += (uint64_t)SHA512_BLOCK_SIZE * 8U;
            in += SHA512_BLOCK_SIZE;
            len -= SHA512_BLOCK_SIZE;
            continue;
        }

        size_t n = SHA512_BLOCK_SIZE - ctx->curlen;
        if (n > len) {
            n = len;
        }
        memcpy(ctx->buf + ctx->curlen, in, n);
        ctx->curlen += n;
        in += n;
        len -= n;
        if (ctx->curlen == SHA512_BLOCK_SIZE) {
            sha512_transform(ctx, ctx->buf);
            ctx->length += (uint64_t)SHA512_BLOCK_SIZE * 8U;
            ctx->curlen = 0;
        }
    }
}

static inline void sha512_final(sha512_ctx* ctx,
    uint8_t out[SHA512_DIGEST_SIZE])
{
    if (ctx->curlen >= sizeof(ctx->buf)) {
        return;
    }

    ctx->length += (uint64_t)ctx->curlen * 8U;
    ctx->buf[ctx->curlen++] = 0x80;

    if (ctx->curlen > 112) {
        while (ctx->curlen < SHA512_BLOCK_SIZE) {
            ctx->buf[ctx->curlen++] = 0;
        }
        sha512_transform(ctx, ctx->buf);
        ctx->curlen = 0;
    }

    while (ctx->curlen < 120) {
        ctx->buf[ctx->curlen++] = 0;
    }

    sha512_store64_be(ctx->buf + 120, ctx->length);
    sha512_transform(ctx, ctx->buf);

    for (size_t i = 0; i < 8; i++) {
        sha512_store64_be(out + (i * 8), ctx->state[i]);
    }
}

static inline void sha512_sum(const void* data, size_t len,
    uint8_t out[SHA512_DIGEST_SIZE])
{
    sha512_ctx ctx;
    sha512_init(&ctx);
    sha512_update(&ctx, data, len);
    sha512_final(&ctx, out);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* C25519_SHA512_H */
