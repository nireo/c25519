#ifndef __C25519_FE_H__
#define __C25519_FE_H__

#include <stddef.h>
#include <stdint.h>

#define M52 0xFFFFFFFFFFFFFULL // 52 ones
#define M48 0x0FFFFFFFFFFFFULL // 48 ones

#define MASK51 ((uint64_t)((1ULL << 51) - 1))

// this is a 51-bit representation the benefit of the 5x51 is lazy reduction, such that we can do
// several operations without reducing and then reduce the whole thing. numbers are little-endian
typedef struct {
    uint64_t l0;
    uint64_t l1;
    uint64_t l2;
    uint64_t l3;
    uint64_t l4;
} fe;

typedef struct {
    uint64_t lo;
    uint64_t hi;
} fe_uint128;

static const fe FE_ZERO = { 0, 0, 0, 0, 0 };
static const fe FE_ONE = { 1, 0, 0, 0, 0 };
static const fe FE_SQRT_M1 = { 1718705420411056ULL, 234908883556509ULL,
    2233514472574048ULL, 2117202627021982ULL,
    765476049583133ULL };

// fe_one sets the element 'e' to one.
static inline fe* fe_one(fe* e)
{
    *e = FE_ONE;
    return e;
}

// fe_zero sets the element 'e' to zero.
static inline fe* fe_zero(fe* e)
{
    *e = FE_ZERO;
    return e;
}

// fe_load64_le takes in 8 bytes and converts those into a internal field element representation.
// the 'in' is expected to be little-endian.
static inline uint64_t fe_load64_le(const uint8_t* in)
{
    return ((uint64_t)in[0]) | ((uint64_t)in[1] << 8) | ((uint64_t)in[2] << 16) | ((uint64_t)in[3] << 24) | ((uint64_t)in[4] << 32) | ((uint64_t)in[5] << 40) | ((uint64_t)in[6] << 48) | ((uint64_t)in[7] << 56);
}

// fe_store64_le stores a internal fe as 8 bytes into out in little-endian encoding. Note that the length
// of out should be at least 8.
static inline void fe_store64_le(uint8_t* out, uint64_t v)
{
    out[0] = (uint8_t)(v);
    out[1] = (uint8_t)(v >> 8);
    out[2] = (uint8_t)(v >> 16);
    out[3] = (uint8_t)(v >> 24);
    out[4] = (uint8_t)(v >> 32);
    out[5] = (uint8_t)(v >> 40);
    out[6] = (uint8_t)(v >> 48);
    out[7] = (uint8_t)(v >> 56);
}

static inline fe_uint128 fe_mul_u64(uint64_t a, uint64_t b)
{
    __uint128_t m = (__uint128_t)a * b;
    fe_uint128 r;
    r.lo = (uint64_t)m;
    r.hi = (uint64_t)(m >> 64);
    return r;
}

static inline fe_uint128 fe_add_mul(fe_uint128 v, uint64_t a, uint64_t b)
{
    __uint128_t acc = ((__uint128_t)v.hi << 64) | v.lo;
    acc += (__uint128_t)a * b;
    return (fe_uint128) { (uint64_t)acc, (uint64_t)(acc >> 64) };
}

static inline uint64_t fe_mul19(uint64_t v)
{
    return v + ((v + (v << 3)) << 1);
}

static inline fe_uint128 fe_add_mul19(fe_uint128 v, uint64_t a, uint64_t b)
{
    __uint128_t acc = ((__uint128_t)v.hi << 64) | v.lo;
    acc += (__uint128_t)fe_mul19(a) * b;
    return (fe_uint128) { (uint64_t)acc, (uint64_t)(acc >> 64) };
}

static inline fe_uint128 fe_add_mul38(fe_uint128 v, uint64_t a, uint64_t b)
{
    __uint128_t acc = ((__uint128_t)v.hi << 64) | v.lo;
    acc += (__uint128_t)fe_mul19(a) * (b << 1);
    return (fe_uint128) { (uint64_t)acc, (uint64_t)(acc >> 64) };
}

static inline uint64_t fe_shift_right_51(fe_uint128 a)
{
    return (a.hi << (64 - 51)) | (a.lo >> 51);
}

static inline fe* fe_carry_propagate(fe* v)
{
    uint64_t l0 = v->l0;
    v->l0 = (v->l0 & MASK51) + fe_mul19(v->l4 >> 51);
    v->l4 = (v->l4 & MASK51) + (v->l3 >> 51);
    v->l3 = (v->l3 & MASK51) + (v->l2 >> 51);
    v->l2 = (v->l2 & MASK51) + (v->l1 >> 51);
    v->l1 = (v->l1 & MASK51) + (l0 >> 51);
    return v;
}

static inline fe* reduce(fe* v)
{
    uint64_t c;
    fe_carry_propagate(v);

    c = (v->l0 + 19) >> 51;
    c = (v->l1 + c) >> 51;
    c = (v->l2 + c) >> 51;
    c = (v->l3 + c) >> 51;
    c = (v->l4 + c) >> 51;

    v->l0 += 19 * c;

    v->l1 += v->l0 >> 51;
    v->l0 &= MASK51;
    v->l2 += v->l1 >> 51;
    v->l1 &= MASK51;
    v->l3 += v->l2 >> 51;
    v->l2 &= MASK51;
    v->l4 += v->l3 >> 51;
    v->l3 &= MASK51;
    v->l4 &= MASK51;

    return v;
}

static inline fe* fe_set(fe* v, const fe* a)
{
    *v = *a;
    return v;
}

static inline fe* fe_add(fe* v, const fe* a, const fe* b)
{
    v->l0 = a->l0 + b->l0;
    v->l1 = a->l1 + b->l1;
    v->l2 = a->l2 + b->l2;
    v->l3 = a->l3 + b->l3;
    v->l4 = a->l4 + b->l4;
    return fe_carry_propagate(v);
}

static inline fe* fe_sub(fe* v, const fe* a, const fe* b)
{
    v->l0 = (a->l0 + UINT64_C(0xFFFFFFFFFFFDA)) - b->l0;
    v->l1 = (a->l1 + UINT64_C(0xFFFFFFFFFFFFE)) - b->l1;
    v->l2 = (a->l2 + UINT64_C(0xFFFFFFFFFFFFE)) - b->l2;
    v->l3 = (a->l3 + UINT64_C(0xFFFFFFFFFFFFE)) - b->l3;
    v->l4 = (a->l4 + UINT64_C(0xFFFFFFFFFFFFE)) - b->l4;
    return fe_carry_propagate(v);
}

static inline fe* fe_neg(fe* v, const fe* a)
{
    return fe_sub(v, &FE_ZERO, a);
}

static inline void fe_mul_generic(fe* v, const fe* a, const fe* b)
{
    uint64_t a0 = a->l0, a1 = a->l1, a2 = a->l2, a3 = a->l3, a4 = a->l4;
    uint64_t b0 = b->l0, b1 = b->l1, b2 = b->l2, b3 = b->l3, b4 = b->l4;

    fe_uint128 r0 = fe_mul_u64(a0, b0);
    r0 = fe_add_mul19(r0, a1, b4);
    r0 = fe_add_mul19(r0, a2, b3);
    r0 = fe_add_mul19(r0, a3, b2);
    r0 = fe_add_mul19(r0, a4, b1);

    fe_uint128 r1 = fe_mul_u64(a0, b1);
    r1 = fe_add_mul(r1, a1, b0);
    r1 = fe_add_mul19(r1, a2, b4);
    r1 = fe_add_mul19(r1, a3, b3);
    r1 = fe_add_mul19(r1, a4, b2);

    fe_uint128 r2 = fe_mul_u64(a0, b2);
    r2 = fe_add_mul(r2, a1, b1);
    r2 = fe_add_mul(r2, a2, b0);
    r2 = fe_add_mul19(r2, a3, b4);
    r2 = fe_add_mul19(r2, a4, b3);

    fe_uint128 r3 = fe_mul_u64(a0, b3);
    r3 = fe_add_mul(r3, a1, b2);
    r3 = fe_add_mul(r3, a2, b1);
    r3 = fe_add_mul(r3, a3, b0);
    r3 = fe_add_mul19(r3, a4, b4);

    fe_uint128 r4 = fe_mul_u64(a0, b4);
    r4 = fe_add_mul(r4, a1, b3);
    r4 = fe_add_mul(r4, a2, b2);
    r4 = fe_add_mul(r4, a3, b1);
    r4 = fe_add_mul(r4, a4, b0);

    uint64_t c0 = fe_shift_right_51(r0);
    uint64_t c1 = fe_shift_right_51(r1);
    uint64_t c2 = fe_shift_right_51(r2);
    uint64_t c3 = fe_shift_right_51(r3);
    uint64_t c4 = fe_shift_right_51(r4);

    uint64_t rr0 = (r0.lo & MASK51) + fe_mul19(c4);
    uint64_t rr1 = (r1.lo & MASK51) + c0;
    uint64_t rr2 = (r2.lo & MASK51) + c1;
    uint64_t rr3 = (r3.lo & MASK51) + c2;
    uint64_t rr4 = (r4.lo & MASK51) + c3;

    v->l0 = (rr0 & MASK51) + fe_mul19(rr4 >> 51);
    v->l1 = (rr1 & MASK51) + (rr0 >> 51);
    v->l2 = (rr2 & MASK51) + (rr1 >> 51);
    v->l3 = (rr3 & MASK51) + (rr2 >> 51);
    v->l4 = (rr4 & MASK51) + (rr3 >> 51);
}

static inline fe* fe_mul(fe* v, const fe* x, const fe* y)
{
    fe_mul_generic(v, x, y);
    return v;
}

static inline void fe_square_generic(fe* v, const fe* a)
{
    uint64_t l0 = a->l0, l1 = a->l1, l2 = a->l2, l3 = a->l3, l4 = a->l4;

    fe_uint128 r0 = fe_mul_u64(l0, l0);
    r0 = fe_add_mul38(r0, l1, l4);
    r0 = fe_add_mul38(r0, l2, l3);

    fe_uint128 r1 = fe_mul_u64(l0 * 2, l1);
    r1 = fe_add_mul38(r1, l2, l4);
    r1 = fe_add_mul19(r1, l3, l3);

    fe_uint128 r2 = fe_mul_u64(l0 * 2, l2);
    r2 = fe_add_mul(r2, l1, l1);
    r2 = fe_add_mul38(r2, l3, l4);

    fe_uint128 r3 = fe_mul_u64(l0 * 2, l3);
    r3 = fe_add_mul(r3, l1 * 2, l2);
    r3 = fe_add_mul19(r3, l4, l4);

    fe_uint128 r4 = fe_mul_u64(l0 * 2, l4);
    r4 = fe_add_mul(r4, l1 * 2, l3);
    r4 = fe_add_mul(r4, l2, l2);

    uint64_t c0 = fe_shift_right_51(r0);
    uint64_t c1 = fe_shift_right_51(r1);
    uint64_t c2 = fe_shift_right_51(r2);
    uint64_t c3 = fe_shift_right_51(r3);
    uint64_t c4 = fe_shift_right_51(r4);

    uint64_t rr0 = (r0.lo & MASK51) + fe_mul19(c4);
    uint64_t rr1 = (r1.lo & MASK51) + c0;
    uint64_t rr2 = (r2.lo & MASK51) + c1;
    uint64_t rr3 = (r3.lo & MASK51) + c2;
    uint64_t rr4 = (r4.lo & MASK51) + c3;

    v->l0 = (rr0 & MASK51) + fe_mul19(rr4 >> 51);
    v->l1 = (rr1 & MASK51) + (rr0 >> 51);
    v->l2 = (rr2 & MASK51) + (rr1 >> 51);
    v->l3 = (rr3 & MASK51) + (rr2 >> 51);
    v->l4 = (rr4 & MASK51) + (rr3 >> 51);
}

static inline fe* fe_square(fe* v, const fe* x)
{
    fe_square_generic(v, x);
    return v;
}

static inline void fe_mul51(uint64_t a, uint32_t b, uint64_t* lo, uint64_t* hi)
{
    __uint128_t m = (__uint128_t)a * b;
    *lo = (uint64_t)m & MASK51;
    *hi = (uint64_t)(m >> 51);
}

static inline fe* fe_mul32(fe* v, const fe* x, uint32_t y)
{
    uint64_t x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi, x4lo, x4hi;
    fe_mul51(x->l0, y, &x0lo, &x0hi);
    fe_mul51(x->l1, y, &x1lo, &x1hi);
    fe_mul51(x->l2, y, &x2lo, &x2hi);
    fe_mul51(x->l3, y, &x3lo, &x3hi);
    fe_mul51(x->l4, y, &x4lo, &x4hi);

    v->l0 = x0lo + 19 * x4hi;
    v->l1 = x1lo + x0hi;
    v->l2 = x2lo + x1hi;
    v->l3 = x3lo + x2hi;
    v->l4 = x4lo + x3hi;
    return v;
}

static inline int fe_set_bytes(fe* v, const uint8_t* x, size_t len)
{
    if (len != 32) {
        return -1;
    }

    fe t;
    t.l0 = fe_load64_le(x + 0) & MASK51;
    t.l1 = (fe_load64_le(x + 6) >> 3) & MASK51;
    t.l2 = (fe_load64_le(x + 12) >> 6) & MASK51;
    t.l3 = (fe_load64_le(x + 19) >> 1) & MASK51;
    t.l4 = (fe_load64_le(x + 24) >> 12) & MASK51;
    *v = t;
    return 0;
}

static inline uint8_t* fe_bytes(uint8_t out[32], const fe* v)
{
    fe t = *v;
    reduce(&t);

    uint64_t u0 = (t.l1 << 51) | t.l0;
    uint64_t u1 = (t.l2 << 38) | (t.l1 >> 13);
    uint64_t u2 = (t.l3 << 25) | (t.l2 >> 26);
    uint64_t u3 = (t.l4 << 12) | (t.l3 >> 39);

    fe_store64_le(out + 0, u0);
    fe_store64_le(out + 8, u1);
    fe_store64_le(out + 16, u2);
    fe_store64_le(out + 24, u3);

    return out;
}

static inline int fe_ct_is_zero_u32(uint32_t x)
{
    return (int)(1U ^ ((x | (uint32_t)(0U - x)) >> 31));
}

static inline int fe_equal(const fe* a, const fe* b)
{
    uint8_t ba[32];
    uint8_t bb[32];
    fe_bytes(ba, a);
    fe_bytes(bb, b);
    uint32_t diff = 0;
    for (size_t i = 0; i < 32; i++) {
        diff |= (uint32_t)(ba[i] ^ bb[i]);
    }
    return fe_ct_is_zero_u32(diff);
}

static inline uint64_t fe_mask64(int cond)
{
    return (uint64_t)0 - (uint64_t)(cond & 1);
}

static inline fe* fe_select(fe* v, const fe* a, const fe* b, int cond)
{
    uint64_t m = fe_mask64(cond);
    v->l0 = (m & a->l0) | (~m & b->l0);
    v->l1 = (m & a->l1) | (~m & b->l1);
    v->l2 = (m & a->l2) | (~m & b->l2);
    v->l3 = (m & a->l3) | (~m & b->l3);
    v->l4 = (m & a->l4) | (~m & b->l4);
    return v;
}

static inline void fe_swap(fe* v, fe* u, int cond)
{
    uint64_t m = fe_mask64(cond);
    uint64_t t = m & (v->l0 ^ u->l0);
    v->l0 ^= t;
    u->l0 ^= t;
    t = m & (v->l1 ^ u->l1);
    v->l1 ^= t;
    u->l1 ^= t;
    t = m & (v->l2 ^ u->l2);
    v->l2 ^= t;
    u->l2 ^= t;
    t = m & (v->l3 ^ u->l3);
    v->l3 ^= t;
    u->l3 ^= t;
    t = m & (v->l4 ^ u->l4);
    v->l4 ^= t;
    u->l4 ^= t;
}

static inline int fe_is_negative(const fe* v)
{
    uint8_t out[32];
    fe_bytes(out, v);
    return (int)(out[0] & 1U);
}

static inline fe* fe_abs(fe* v, const fe* u)
{
    fe neg;
    fe_neg(&neg, u);
    return fe_select(v, &neg, u, fe_is_negative(u));
}

static inline fe* fe_pow22523(fe* v, const fe* x)
{
    fe t0, t1, t2;

    fe_square(&t0, x);
    fe_square(&t1, &t0);
    fe_square(&t1, &t1);
    fe_mul(&t1, x, &t1);
    fe_mul(&t0, &t0, &t1);
    fe_square(&t0, &t0);
    fe_mul(&t0, &t1, &t0);
    fe_square(&t1, &t0);
    for (int i = 1; i < 5; i++) {
        fe_square(&t1, &t1);
    }
    fe_mul(&t0, &t1, &t0);
    fe_square(&t1, &t0);
    for (int i = 1; i < 10; i++) {
        fe_square(&t1, &t1);
    }
    fe_mul(&t1, &t1, &t0);
    fe_square(&t2, &t1);
    for (int i = 1; i < 20; i++) {
        fe_square(&t2, &t2);
    }
    fe_mul(&t1, &t2, &t1);
    fe_square(&t1, &t1);
    for (int i = 1; i < 10; i++) {
        fe_square(&t1, &t1);
    }
    fe_mul(&t0, &t1, &t0);
    fe_square(&t1, &t0);
    for (int i = 1; i < 50; i++) {
        fe_square(&t1, &t1);
    }
    fe_mul(&t1, &t1, &t0);
    fe_square(&t2, &t1);
    for (int i = 1; i < 100; i++) {
        fe_square(&t2, &t2);
    }
    fe_mul(&t1, &t2, &t1);
    fe_square(&t1, &t1);
    for (int i = 1; i < 50; i++) {
        fe_square(&t1, &t1);
    }
    fe_mul(&t0, &t1, &t0);
    fe_square(&t0, &t0);
    fe_square(&t0, &t0);
    return fe_mul(v, &t0, x);
}

static inline fe* fe_invert(fe* v, const fe* z)
{
    fe z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;

    fe_square(&z2, z);
    fe_square(&t, &z2);
    fe_square(&t, &t);
    fe_mul(&z9, &t, z);
    fe_mul(&z11, &z9, &z2);
    fe_square(&t, &z11);
    fe_mul(&z2_5_0, &t, &z9);

    fe_square(&t, &z2_5_0);
    for (int i = 0; i < 4; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&z2_10_0, &t, &z2_5_0);

    fe_square(&t, &z2_10_0);
    for (int i = 0; i < 9; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&z2_20_0, &t, &z2_10_0);

    fe_square(&t, &z2_20_0);
    for (int i = 0; i < 19; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&t, &t, &z2_20_0);

    fe_square(&t, &t);
    for (int i = 0; i < 9; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&z2_50_0, &t, &z2_10_0);

    fe_square(&t, &z2_50_0);
    for (int i = 0; i < 49; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&z2_100_0, &t, &z2_50_0);

    fe_square(&t, &z2_100_0);
    for (int i = 0; i < 99; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&t, &t, &z2_100_0);

    fe_square(&t, &t);
    for (int i = 0; i < 49; i++) {
        fe_square(&t, &t);
    }
    fe_mul(&t, &t, &z2_50_0);

    fe_square(&t, &t);
    fe_square(&t, &t);
    fe_square(&t, &t);
    fe_square(&t, &t);
    fe_square(&t, &t);

    return fe_mul(v, &t, &z11);
}

static inline int fe_sqrt_ratio(fe* r, const fe* u, const fe* v)
{
    fe t0, v2, uv3, uv7, rr, check, u_neg, r_prime;

    fe_square(&v2, v);
    fe_mul(&t0, &v2, v);
    fe_mul(&uv3, u, &t0);
    fe_square(&t0, &v2);
    fe_mul(&uv7, &uv3, &t0);
    fe_pow22523(&t0, &uv7);
    fe_mul(&rr, &uv3, &t0);

    fe_square(&t0, &rr);
    fe_mul(&check, v, &t0);

    fe_neg(&u_neg, u);
    int correct = fe_equal(&check, u);
    int flipped = fe_equal(&check, &u_neg);
    fe_mul(&t0, &u_neg, &FE_SQRT_M1);
    int flipped_i = fe_equal(&check, &t0);

    fe_mul(&r_prime, &rr, &FE_SQRT_M1);
    fe_select(&rr, &r_prime, &rr, flipped | flipped_i);
    fe_abs(r, &rr);
    return correct | flipped;
}

#endif
