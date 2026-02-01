#ifndef __C25519_H__
#define __C25519_H__

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "point.h"
#include "scalar.h"
#include "scalarmult.h"
#include "seed.h"
#include "sha512.h"

#define C25519_PUBLIC_KEY_SIZE 32
#define C25519_PRIVATE_KEY_SIZE 64
#define C25519_SIGNATURE_SIZE 64
#define C25519_SEED_SIZE 32

typedef uint8_t ed_public_key[C25519_PUBLIC_KEY_SIZE];
typedef uint8_t ed_private_key[C25519_PRIVATE_KEY_SIZE];
typedef uint8_t ed_seed[C25519_SEED_SIZE];

#define ED25519_PUBLIC_KEY_SIZE C25519_PUBLIC_KEY_SIZE
#define ED25519_PRIVATE_KEY_SIZE C25519_PRIVATE_KEY_SIZE
#define ED25519_SIGNATURE_SIZE C25519_SIGNATURE_SIZE
#define ED25519_SEED_SIZE C25519_SEED_SIZE
#define ED25519_SHA512_SIZE 64

typedef ed_public_key ed25519_public_key;
typedef ed_private_key ed25519_private_key;
typedef ed_seed ed25519_seed;
typedef uint8_t ed25519_signature[ED25519_SIGNATURE_SIZE];

static inline void ed_private_to_public(const ed_private_key* sk, ed_public_key* pk)
{
    memcpy(*pk, (*sk) + 32, 32);
}

static inline void ed_private_seed(const ed_private_key* sk, ed_seed* seed)
{
    memcpy(*seed, *sk, 32);
}

static inline int ed25519_dom2(sha512_ctx* ctx, uint8_t phflag, const uint8_t* ctx_data, size_t ctx_len)
{
    static const uint8_t dom2_prefix[] = "SigEd25519 no Ed25519 collisions";

    if (ctx_len > 255) {
        return -1;
    }

    sha512_update(ctx, dom2_prefix, sizeof(dom2_prefix) - 1);
    uint8_t tag[2];
    tag[0] = phflag;
    tag[1] = (uint8_t)ctx_len;
    sha512_update(ctx, tag, sizeof(tag));
    if (ctx_len > 0) {
        sha512_update(ctx, ctx_data, ctx_len);
    }
    return 0;
}

static inline int ed25519_public_key_from_seed(ed25519_public_key* pk, const ed25519_seed seed)
{
    uint8_t digest[SHA512_DIGEST_SIZE];
    sha512_sum(seed, ED25519_SEED_SIZE, digest);

    scalar s;
    if (scalar_set_bytes_with_clamping(&s, digest, 32) != 0) {
        return -1;
    }

    point A;
    point_scalar_base_mult(&A, &s);
    point_bytes(*pk, &A);
    return 0;
}

static inline int ed25519_keypair_from_seed(ed25519_public_key* pk, ed25519_private_key* sk, const ed25519_seed seed)
{
    if (ed25519_public_key_from_seed(pk, seed) != 0) {
        return -1;
    }
    memcpy(*sk, seed, ED25519_SEED_SIZE);
    memcpy(*sk + ED25519_SEED_SIZE, *pk, ED25519_PUBLIC_KEY_SIZE);
    return 0;
}

static inline int ed25519_keypair(ed25519_public_key* pk, ed25519_private_key* sk)
{
    ed25519_seed seed;
    if (ed25519_create_seed(seed) != 0) {
        return -1;
    }
    return ed25519_keypair_from_seed(pk, sk, seed);
}

static inline void ed25519_public_key_from_private(ed25519_public_key* pk, const ed25519_private_key* sk)
{
    memcpy(*pk, (*sk) + ED25519_SEED_SIZE, ED25519_PUBLIC_KEY_SIZE);
}

static inline int ed25519_sign(ed25519_signature sig, const uint8_t* msg, size_t msg_len, const ed25519_private_key* sk)
{
    const uint8_t* seed = *sk;
    const uint8_t* pub = (*sk) + ED25519_SEED_SIZE;

    uint8_t digest[SHA512_DIGEST_SIZE];
    sha512_sum(seed, ED25519_SEED_SIZE, digest);

    scalar s;
    if (scalar_set_bytes_with_clamping(&s, digest, 32) != 0) {
        return -1;
    }

    uint8_t prefix[32];
    memcpy(prefix, digest + 32, sizeof(prefix));

    sha512_ctx ctx;
    sha512_init(&ctx);
    sha512_update(&ctx, prefix, sizeof(prefix));
    sha512_update(&ctx, msg, msg_len);
    sha512_final(&ctx, digest);

    scalar r;
    if (scalar_set_uniform_bytes(&r, digest, sizeof(digest)) != 0) {
        return -1;
    }

    point R;
    point_scalar_base_mult(&R, &r);
    uint8_t r_bytes[32];
    point_bytes(r_bytes, &R);

    sha512_init(&ctx);
    sha512_update(&ctx, r_bytes, sizeof(r_bytes));
    sha512_update(&ctx, pub, ED25519_PUBLIC_KEY_SIZE);
    sha512_update(&ctx, msg, msg_len);
    sha512_final(&ctx, digest);

    scalar k;
    if (scalar_set_uniform_bytes(&k, digest, sizeof(digest)) != 0) {
        return -1;
    }

    scalar S;
    scalar_multiply_add(&S, &k, &s, &r);

    memcpy(sig, r_bytes, 32);
    scalar_bytes(sig + 32, &S);
    return 0;
}

static inline int ed25519_sign_ph(ed25519_signature sig, const uint8_t* msg, size_t msg_len, const ed25519_private_key* sk, const uint8_t* ctx_data, size_t ctx_len)
{
    if (msg_len != ED25519_SHA512_SIZE) {
        return -1;
    }

    const uint8_t* seed = *sk;
    const uint8_t* pub = (*sk) + ED25519_SEED_SIZE;

    uint8_t digest[SHA512_DIGEST_SIZE];
    sha512_sum(seed, ED25519_SEED_SIZE, digest);

    scalar s;
    if (scalar_set_bytes_with_clamping(&s, digest, 32) != 0) {
        return -1;
    }

    uint8_t prefix[32];
    memcpy(prefix, digest + 32, sizeof(prefix));

    sha512_ctx ctx;
    sha512_init(&ctx);
    if (ed25519_dom2(&ctx, 1, ctx_data, ctx_len) != 0) {
        return -1;
    }
    sha512_update(&ctx, prefix, sizeof(prefix));
    sha512_update(&ctx, msg, msg_len);
    sha512_final(&ctx, digest);

    scalar r;
    if (scalar_set_uniform_bytes(&r, digest, sizeof(digest)) != 0) {
        return -1;
    }

    point R;
    point_scalar_base_mult(&R, &r);
    uint8_t r_bytes[32];
    point_bytes(r_bytes, &R);

    sha512_init(&ctx);
    if (ed25519_dom2(&ctx, 1, ctx_data, ctx_len) != 0) {
        return -1;
    }
    sha512_update(&ctx, r_bytes, sizeof(r_bytes));
    sha512_update(&ctx, pub, ED25519_PUBLIC_KEY_SIZE);
    sha512_update(&ctx, msg, msg_len);
    sha512_final(&ctx, digest);

    scalar k;
    if (scalar_set_uniform_bytes(&k, digest, sizeof(digest)) != 0) {
        return -1;
    }

    scalar S;
    scalar_multiply_add(&S, &k, &s, &r);

    memcpy(sig, r_bytes, 32);
    scalar_bytes(sig + 32, &S);
    return 0;
}

static inline int ed25519_verify(const ed25519_signature sig, const uint8_t* msg, size_t msg_len, const ed25519_public_key* pk)
{
    if (sig[63] & 0xe0) {
        return -1;
    }

    point A;
    if (point_set_bytes(&A, *pk, ED25519_PUBLIC_KEY_SIZE) != 0) {
        return -1;
    }

    sha512_ctx ctx;
    uint8_t digest[SHA512_DIGEST_SIZE];
    sha512_init(&ctx);
    sha512_update(&ctx, sig, 32);
    sha512_update(&ctx, *pk, ED25519_PUBLIC_KEY_SIZE);
    sha512_update(&ctx, msg, msg_len);
    sha512_final(&ctx, digest);

    scalar k;
    if (scalar_set_uniform_bytes(&k, digest, sizeof(digest)) != 0) {
        return -1;
    }

    scalar S;
    if (scalar_set_canonical_bytes(&S, sig + 32, 32) != 0) {
        return -1;
    }

    point minusA;
    point_negate(&minusA, &A);

    point R;
    point_vartime_double_scalar_base_mult(&R, &k, &minusA, &S);

    uint8_t r_check[32];
    point_bytes(r_check, &R);
    if (memcmp(r_check, sig, 32) != 0) {
        return -1;
    }

    return 0;
}

static inline int ed25519_verify_ph(
    const ed25519_signature sig, const uint8_t* msg, size_t msg_len, const ed25519_public_key* pk, const uint8_t* ctx_data, size_t ctx_len)
{
    if (msg_len != ED25519_SHA512_SIZE) {
        return -1;
    }
    if (sig[63] & 0xe0) {
        return -1;
    }

    point A;
    if (point_set_bytes(&A, *pk, ED25519_PUBLIC_KEY_SIZE) != 0) {
        return -1;
    }

    sha512_ctx ctx;
    uint8_t digest[SHA512_DIGEST_SIZE];
    sha512_init(&ctx);
    if (ed25519_dom2(&ctx, 1, ctx_data, ctx_len) != 0) {
        return -1;
    }
    sha512_update(&ctx, sig, 32);
    sha512_update(&ctx, *pk, ED25519_PUBLIC_KEY_SIZE);
    sha512_update(&ctx, msg, msg_len);
    sha512_final(&ctx, digest);

    scalar k;
    if (scalar_set_uniform_bytes(&k, digest, sizeof(digest)) != 0) {
        return -1;
    }

    scalar S;
    if (scalar_set_canonical_bytes(&S, sig + 32, 32) != 0) {
        return -1;
    }

    point minusA;
    point_negate(&minusA, &A);

    point R;
    point_vartime_double_scalar_base_mult(&R, &k, &minusA, &S);

    uint8_t r_check[32];
    point_bytes(r_check, &R);
    if (memcmp(r_check, sig, 32) != 0) {
        return -1;
    }

    return 0;
}

#endif
