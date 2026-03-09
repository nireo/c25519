#ifndef C25519_SCALAR_H
#define C25519_SCALAR_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The fiat-crypto implementation for the Ed25519 scalar field. */
#include "fiat_scalar.h"

typedef struct scalar {
    fiat_25519_scalar_montgomery_domain_field_element s;
} scalar;

static const scalar SCALAR_ZERO = { { 0, 0, 0, 0 } };

/* scalarTwo168 and scalarTwo336 are 2^168 and 2^336 modulo l in Montgomery
 * form. */
static const scalar scalar_two_168 = {
    { 0x5b8ab432eac74798ULL, 0x38afddd6de59d5d7ULL, 0xa2c131b399411b7cULL,
        0x6329a7ed9ce5a30ULL }
};
static const scalar scalar_two_336 = {
    { 0xbd3d108e2b35ecc5ULL, 0x5c3a3718bdf9c90bULL, 0x63aa97a331b4f2eeULL,
        0x3d217f5be65cb5cULL }
};

/* scalar_minus_one_bytes is l - 1 in little endian. */
static const uint8_t scalar_minus_one_bytes[32] = {
    236, 211, 245, 92, 26, 99, 18, 88, 214, 156, 247, 162, 222, 249, 222, 20,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16
};

static inline uint64_t scalar_load64_le(const uint8_t* in)
{
    return ((uint64_t)in[0]) | ((uint64_t)in[1] << 8) | ((uint64_t)in[2] << 16) | ((uint64_t)in[3] << 24) | ((uint64_t)in[4] << 32) | ((uint64_t)in[5] << 40) | ((uint64_t)in[6] << 48) | ((uint64_t)in[7] << 56);
}

static inline uint64_t scalar_subborrow_u64(uint64_t a, uint64_t b,
    uint64_t borrow, uint64_t* out)
{
    uint64_t res = a - b - borrow;
    uint64_t borrow_out = (a < b) | (borrow & (a == b));
    *out = res;
    return borrow_out;
}

static inline int scalar_is_reduced_bytes(const uint8_t* s, size_t len)
{
    if (len != 32) {
        return 0;
    }

    uint64_t s0 = scalar_load64_le(s);
    uint64_t s1 = scalar_load64_le(s + 8);
    uint64_t s2 = scalar_load64_le(s + 16);
    uint64_t s3 = scalar_load64_le(s + 24);

    uint64_t l0 = scalar_load64_le(scalar_minus_one_bytes);
    uint64_t l1 = scalar_load64_le(scalar_minus_one_bytes + 8);
    uint64_t l2 = scalar_load64_le(scalar_minus_one_bytes + 16);
    uint64_t l3 = scalar_load64_le(scalar_minus_one_bytes + 24);

    uint64_t out;
    uint64_t borrow = 0;
    borrow = scalar_subborrow_u64(l0, s0, borrow, &out);
    borrow = scalar_subborrow_u64(l1, s1, borrow, &out);
    borrow = scalar_subborrow_u64(l2, s2, borrow, &out);
    borrow = scalar_subborrow_u64(l3, s3, borrow, &out);
    return borrow == 0;
}

/* NewScalar returns a new zero scalar. */
static inline scalar scalar_new(void) { return SCALAR_ZERO; }

static inline scalar* scalar_zero(scalar* s)
{
    *s = SCALAR_ZERO;
    return s;
}

static inline scalar* scalar_one(scalar* s)
{
    fiat_25519_scalar_set_one(s->s);
    return s;
}

static inline scalar* scalar_set(scalar* s, const scalar* x)
{
    *s = *x;
    return s;
}

static inline scalar* scalar_add(scalar* s, const scalar* x, const scalar* y)
{
    fiat_25519_scalar_add(s->s, x->s, y->s);
    return s;
}

static inline scalar* scalar_subtract(scalar* s, const scalar* x,
    const scalar* y)
{
    fiat_25519_scalar_sub(s->s, x->s, y->s);
    return s;
}

static inline scalar* scalar_negate(scalar* s, const scalar* x)
{
    fiat_25519_scalar_opp(s->s, x->s);
    return s;
}

static inline scalar* scalar_mul(scalar* s, const scalar* x, const scalar* y)
{
    fiat_25519_scalar_mul(s->s, x->s, y->s);
    return s;
}

static inline scalar* scalar_multiply_add(scalar* s, const scalar* x,
    const scalar* y, const scalar* z)
{
    scalar zcopy = *z;
    scalar_mul(s, x, y);
    return scalar_add(s, s, &zcopy);
}

static inline scalar* scalar_set_short_bytes(scalar* s, const uint8_t* x,
    size_t len)
{
    if (len >= 32) {
        return NULL;
    }
    uint8_t buf[32];
    memset(buf, 0, sizeof(buf));
    if (len > 0) {
        memcpy(buf, x, len);
    }
    fiat_25519_scalar_non_montgomery_domain_field_element t;
    fiat_25519_scalar_from_bytes(t, buf);
    fiat_25519_scalar_to_montgomery(s->s, t);
    return s;
}

/* SetUniformBytes sets s = x mod l for a 64-byte little-endian integer. */
static inline int scalar_set_uniform_bytes(scalar* s, const uint8_t* x,
    size_t len)
{
    if (len != 64) {
        return -1;
    }
    scalar_set_short_bytes(s, x, 21);
    scalar t;
    scalar_set_short_bytes(&t, x + 21, 21);
    scalar_mul(&t, &t, &scalar_two_168);
    scalar_add(s, s, &t);
    scalar_set_short_bytes(&t, x + 42, 22);
    scalar_mul(&t, &t, &scalar_two_336);
    scalar_add(s, s, &t);
    return 0;
}

/* SetCanonicalBytes sets s = x when x is a canonical 32-byte encoding. */
static inline int scalar_set_canonical_bytes(scalar* s, const uint8_t* x,
    size_t len)
{
    if (len != 32) {
        return -1;
    }
    if (!scalar_is_reduced_bytes(x, len)) {
        return -1;
    }
    scalar tmp;
    fiat_25519_scalar_non_montgomery_domain_field_element t;
    fiat_25519_scalar_from_bytes(t, x);
    fiat_25519_scalar_to_montgomery(tmp.s, t);
    *s = tmp;
    return 0;
}

/* SetBytesWithClamping applies RFC 8032 clamping then reduces modulo l. */
static inline int scalar_set_bytes_with_clamping(scalar* s, const uint8_t* x,
    size_t len)
{
    if (len != 32) {
        return -1;
    }
    uint8_t wide[64];
    memset(wide, 0, sizeof(wide));
    memcpy(wide, x, 32);
    wide[0] &= 248;
    wide[31] &= 63;
    wide[31] |= 64;
    return scalar_set_uniform_bytes(s, wide, sizeof(wide));
}

/* Bytes returns the canonical 32-byte little-endian encoding of s. */
static inline uint8_t* scalar_bytes(uint8_t out[32], const scalar* s)
{
    fiat_25519_scalar_non_montgomery_domain_field_element t;
    fiat_25519_scalar_from_montgomery(t, s->s);
    fiat_25519_scalar_to_bytes(out, t);
    return out;
}

/* Equal returns 1 if s and t are equal, and 0 otherwise. */
static inline int scalar_equal(const scalar* s, const scalar* t)
{
    fiat_25519_scalar_montgomery_domain_field_element diff;
    fiat_25519_scalar_sub(diff, s->s, t->s);
    uint64_t nonzero;
    fiat_25519_scalar_nonzero(&nonzero, diff);
    nonzero |= nonzero >> 32;
    nonzero |= nonzero >> 16;
    nonzero |= nonzero >> 8;
    nonzero |= nonzero >> 4;
    nonzero |= nonzero >> 2;
    nonzero |= nonzero >> 1;
    return (int)((~nonzero) & 1);
}

/* nonAdjacentForm computes a width-w non-adjacent form for this scalar. */
static inline int scalar_non_adjacent_form(int8_t out[256], const scalar* s,
    uint32_t w)
{
    uint8_t b[32];
    scalar_bytes(b, s);
    if (b[31] > 127) {
        return -1;
    }
    if (w < 2 || w > 8) {
        return -1;
    }

    memset(out, 0, 256);
    uint64_t digits[5] = { 0, 0, 0, 0, 0 };
    for (uint32_t i = 0; i < 4; i++) {
        digits[i] = scalar_load64_le(b + i * 8);
    }

    uint64_t width = (uint64_t)1 << w;
    uint64_t window_mask = width - 1;

    uint32_t pos = 0;
    uint64_t carry = 0;
    while (pos < 256) {
        uint32_t index_u64 = pos / 64;
        uint32_t index_bit = pos % 64;
        uint64_t bit_buf;
        if (index_bit < 64 - w) {
            bit_buf = digits[index_u64] >> index_bit;
        } else {
            bit_buf = (digits[index_u64] >> index_bit) | (digits[index_u64 + 1] << (64 - index_bit));
        }

        uint64_t window = carry + (bit_buf & window_mask);

        if ((window & 1) == 0) {
            pos += 1;
            continue;
        }

        if (window < width / 2) {
            carry = 0;
            out[pos] = (int8_t)window;
        } else {
            carry = 1;
            out[pos] = (int8_t)((int64_t)window - (int64_t)width);
        }
        pos += w;
    }
    return 0;
}

static inline int scalar_signed_radix16(int8_t out[64], const scalar* s)
{
    uint8_t b[32];
    scalar_bytes(b, s);
    if (b[31] > 127) {
        return -1;
    }

    for (uint32_t i = 0; i < 32; i++) {
        out[2 * i] = (int8_t)(b[i] & 15);
        out[2 * i + 1] = (int8_t)((b[i] >> 4) & 15);
    }

    for (uint32_t i = 0; i < 63; i++) {
        int carry = (out[i] + 8) >> 4;
        out[i] = (int8_t)(out[i] - (carry << 4));
        out[i + 1] = (int8_t)(out[i + 1] + carry);
    }

    return 0;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* C25519_SCALAR_H */
