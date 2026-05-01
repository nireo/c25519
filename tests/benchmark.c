#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../c25519.h"

#define BENCH_MIN_NS 100000000ULL

#if defined(_MSC_VER)
#define BENCH_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define BENCH_NOINLINE __attribute__((noinline))
#else
#define BENCH_NOINLINE
#endif

static volatile uint64_t bench_sink = 0;
static int bench_failed = 0;

static uint64_t now_ns(void)
{
    struct timespec ts;
#if defined(TIME_UTC)
    timespec_get(&ts, TIME_UTC);
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

typedef void (*bench_fn)(void* ctx);

static BENCH_NOINLINE double run_bench(bench_fn fn, void* ctx)
{
    size_t iters = 1;
    uint64_t elapsed_ns = 0;

    fn(ctx);

    for (;;) {
        uint64_t start = now_ns();
        for (size_t i = 0; i < iters; i++) {
            fn(ctx);
        }
        uint64_t end = now_ns();
        elapsed_ns = end - start;
        if (elapsed_ns >= BENCH_MIN_NS || iters >= (1U << 30)) {
            break;
        }
        iters *= 2;
    }

    if (iters == 0 || elapsed_ns == 0) {
        return 0.0;
    }
    return (double)elapsed_ns / (double)iters / 1000.0;
}

static void print_result(const char* name, double us_per_op)
{
    double per_sec = us_per_op > 0.0 ? 1000000.0 / us_per_op : 0.0;
    if (us_per_op >= 1.0) {
        printf("%s: %.0fus (%.0f per second)\n", name, us_per_op, per_sec);
        return;
    }

    double scaled = us_per_op * 1000.0;
    const char* unit = "ns";
    if (scaled < 1.0) {
        scaled *= 1000.0;
        unit = "ps";
    }

    if (scaled >= 100.0) {
        printf("%s: %.0f%s (%.0f per second)\n", name, scaled, unit, per_sec);
        return;
    }
    if (scaled >= 10.0) {
        printf("%s: %.1f%s (%.0f per second)\n", name, scaled, unit, per_sec);
        return;
    }
    printf("%s: %.2f%s (%.0f per second)\n", name, scaled, unit, per_sec);
}

struct seed_ctx {
    ed25519_seed seed;
};

static BENCH_NOINLINE void bench_seed(void* vctx)
{
    struct seed_ctx* ctx = vctx;
    if (ed25519_create_seed(ctx->seed) != 0) {
        bench_failed = 1;
        return;
    }
    bench_sink ^= ctx->seed[0];
}

struct keygen_ctx {
    ed25519_seed seed;
    ed25519_public_key pk;
    ed25519_private_key sk;
};

static BENCH_NOINLINE void bench_keygen(void* vctx)
{
    struct keygen_ctx* ctx = vctx;
    ctx->seed[0]++;
    if (ed25519_keypair_from_seed(&ctx->pk, &ctx->sk, ctx->seed) != 0) {
        bench_failed = 1;
        return;
    }
    bench_sink ^= ctx->pk[0];
}

struct sign_ctx {
    ed25519_private_key sk;
    uint8_t msg[32];
    ed25519_signature sig;
};

struct sha512_bench_ctx {
    uint8_t msg[1024];
    uint8_t digest[SHA512_DIGEST_SIZE];
};

static BENCH_NOINLINE void bench_sha512(void* vctx)
{
    struct sha512_bench_ctx* ctx = vctx;
    sha512_sum(ctx->msg, sizeof(ctx->msg), ctx->digest);
    bench_sink ^= ctx->digest[0];
}

static BENCH_NOINLINE void bench_sign(void* vctx)
{
    struct sign_ctx* ctx = vctx;
    const ed25519_private_key* skp = (const ed25519_private_key*)(const void*)&ctx->sk;
    if (ed25519_sign(ctx->sig, ctx->msg, sizeof(ctx->msg), skp) != 0) {
        bench_failed = 1;
        return;
    }
    bench_sink ^= ctx->sig[0];
}

struct verify_ctx {
    ed25519_public_key pk;
    uint8_t msg[32];
    ed25519_signature sig;
};

static BENCH_NOINLINE void bench_verify(void* vctx)
{
    struct verify_ctx* ctx = vctx;
    const ed25519_public_key* pkp = (const ed25519_public_key*)(const void*)&ctx->pk;
    int ok = ed25519_verify(ctx->sig, ctx->msg, sizeof(ctx->msg), pkp);
    if (ok != 0) {
        bench_failed = 1;
        return;
    }
    bench_sink ^= (uint64_t)(ok == 0);
}

struct scalar_add_ctx {
    scalar a;
    scalar b;
    scalar out;
};

static BENCH_NOINLINE void bench_scalar_add(void* vctx)
{
    struct scalar_add_ctx* ctx = vctx;
    scalar_add(&ctx->out, &ctx->a, &ctx->b);
    bench_sink ^= ctx->out.s[0];
}

struct keyex_ctx {
    scalar a;
    point pub_b;
    uint8_t shared[32];
};

static int scalar_from_seed_uniform(scalar* out, const ed25519_seed seed)
{
    uint8_t digest[SHA512_DIGEST_SIZE];
    sha512_sum(seed, ED25519_SEED_SIZE, digest);
    return scalar_set_uniform_bytes(out, digest, sizeof(digest));
}

static BENCH_NOINLINE void bench_key_exchange(void* vctx)
{
    struct keyex_ctx* ctx = vctx;
    point shared;
    point_scalar_mult(&shared, &ctx->a, &ctx->pub_b);
    point_bytes(ctx->shared, &shared);
    bench_sink ^= ctx->shared[0];
}

int main(void)
{
    struct seed_ctx seed_ctx = { 0 };
    double seed_us = run_bench(bench_seed, &seed_ctx);
    if (bench_failed) {
        fprintf(stderr, "seed generation failed\n");
        return 1;
    }
    print_result("seed generation", seed_us);

    struct keygen_ctx keygen_ctx = { 0 };
    for (size_t i = 0; i < sizeof(keygen_ctx.seed); i++) {
        keygen_ctx.seed[i] = (uint8_t)i;
    }
    double keygen_us = run_bench(bench_keygen, &keygen_ctx);
    if (bench_failed) {
        fprintf(stderr, "key generation failed\n");
        return 1;
    }
    print_result("key generation", keygen_us);

    struct sha512_bench_ctx sha512_ctx = { 0 };
    for (size_t i = 0; i < sizeof(sha512_ctx.msg); i++) {
        sha512_ctx.msg[i] = (uint8_t)i;
    }
    double sha512_us = run_bench(bench_sha512, &sha512_ctx);
    print_result("sha512 (1024 bytes)", sha512_us);

    ed25519_seed sign_seed;
    if (ed25519_create_seed(sign_seed) != 0) {
        fprintf(stderr, "sign seed failed\n");
        return 1;
    }
    ed25519_public_key sign_pk;
    ed25519_private_key sign_sk;
    if (ed25519_keypair_from_seed(&sign_pk, &sign_sk, sign_seed) != 0) {
        fprintf(stderr, "sign keypair failed\n");
        return 1;
    }

    struct sign_ctx sign_ctx = { 0 };
    memcpy(sign_ctx.sk, sign_sk, sizeof(sign_ctx.sk));
    for (size_t i = 0; i < sizeof(sign_ctx.msg); i++) {
        sign_ctx.msg[i] = (uint8_t)(i + 1U);
    }
    double sign_us = run_bench(bench_sign, &sign_ctx);
    if (bench_failed) {
        fprintf(stderr, "sign failed\n");
        return 1;
    }
    print_result("message signing (short message)", sign_us);

    struct verify_ctx verify_ctx = { 0 };
    memcpy(verify_ctx.pk, sign_pk, sizeof(verify_ctx.pk));
    memcpy(verify_ctx.msg, sign_ctx.msg, sizeof(verify_ctx.msg));
    {
        const ed25519_private_key* skp = (const ed25519_private_key*)(const void*)&sign_sk;
        if (ed25519_sign(verify_ctx.sig, verify_ctx.msg, sizeof(verify_ctx.msg),
                skp)
            != 0) {
            fprintf(stderr, "verify seed sign failed\n");
            return 1;
        }
    }
    double verify_us = run_bench(bench_verify, &verify_ctx);
    if (bench_failed) {
        fprintf(stderr, "verify failed\n");
        return 1;
    }
    print_result("message verifying (short message)", verify_us);

    struct scalar_add_ctx scalar_ctx = { 0 };
    uint8_t scalar_buf[64];
    for (size_t i = 0; i < sizeof(scalar_buf); i++) {
        scalar_buf[i] = (uint8_t)(i + 3U);
    }
    if (scalar_set_uniform_bytes(&scalar_ctx.a, scalar_buf, sizeof(scalar_buf)) != 0) {
        fprintf(stderr, "scalar init failed\n");
        return 1;
    }
    for (size_t i = 0; i < sizeof(scalar_buf); i++) {
        scalar_buf[i] = (uint8_t)(0xA5U ^ (uint8_t)i);
    }
    if (scalar_set_uniform_bytes(&scalar_ctx.b, scalar_buf, sizeof(scalar_buf)) != 0) {
        fprintf(stderr, "scalar init failed\n");
        return 1;
    }
    double scalar_us = run_bench(bench_scalar_add, &scalar_ctx);
    if (bench_failed) {
        fprintf(stderr, "scalar add failed\n");
        return 1;
    }
    print_result("scalar addition", scalar_us);

    struct keyex_ctx keyex_ctx = { 0 };
    ed25519_seed seed_a;
    ed25519_seed seed_b;
    if (ed25519_create_seed(seed_a) != 0 || ed25519_create_seed(seed_b) != 0) {
        fprintf(stderr, "key exchange seed failed\n");
        return 1;
    }
    scalar b;
    if (scalar_from_seed_uniform(&keyex_ctx.a, seed_a) != 0 || scalar_from_seed_uniform(&b, seed_b) != 0) {
        fprintf(stderr, "key exchange scalar failed\n");
        return 1;
    }
    point_scalar_base_mult(&keyex_ctx.pub_b, &b);
    double keyex_us = run_bench(bench_key_exchange, &keyex_ctx);
    if (bench_failed) {
        fprintf(stderr, "key exchange failed\n");
        return 1;
    }
    print_result("key exchange", keyex_us);

    if (bench_sink == 0xdeadbeefULL) {
        printf("%llu\n", (unsigned long long)bench_sink);
    }

    return 0;
}
