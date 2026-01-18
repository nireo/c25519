#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "../fe.h"

#define TEST_ITERS 1000

static uint64_t rng_state = 0x243f6a8885a308d3ULL;

static uint64_t rng_u64(void)
{
    uint64_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}

static uint64_t mask52(void)
{
    return ((uint64_t)1 << 52) - 1;
}

static fe random_fe(void)
{
    fe x;
    x.l0 = rng_u64() & mask52();
    x.l1 = rng_u64() & mask52();
    x.l2 = rng_u64() & mask52();
    x.l3 = rng_u64() & mask52();
    x.l4 = rng_u64() & mask52();
    return x;
}

static const uint64_t weird_limbs51[] = {
    0, 0, 0, 0,
    1,
    19 - 1,
    19,
    0x2aaaaaaaaaaaaULL,
    0x5555555555555ULL,
    ((uint64_t)1 << 51) - 20,
    ((uint64_t)1 << 51) - 19,
    ((uint64_t)1 << 51) - 1, ((uint64_t)1 << 51) - 1,
    ((uint64_t)1 << 51) - 1, ((uint64_t)1 << 51) - 1,
};

static const uint64_t weird_limbs52[] = {
    0, 0, 0, 0, 0, 0,
    1,
    19 - 1,
    19,
    0x2aaaaaaaaaaaaULL,
    0x5555555555555ULL,
    ((uint64_t)1 << 51) - 20,
    ((uint64_t)1 << 51) - 19,
    ((uint64_t)1 << 51) - 1, ((uint64_t)1 << 51) - 1,
    ((uint64_t)1 << 51) - 1, ((uint64_t)1 << 51) - 1,
    ((uint64_t)1 << 51) - 1, ((uint64_t)1 << 51) - 1,
    (uint64_t)1 << 51,
    ((uint64_t)1 << 51) + 1,
    ((uint64_t)1 << 52) - 19,
    ((uint64_t)1 << 52) - 1,
};

static fe random_weird_fe(void)
{
    fe x;
    size_t i0 = (size_t)(rng_u64() % (sizeof(weird_limbs52) / sizeof(weird_limbs52[0])));
    size_t i1 = (size_t)(rng_u64() % (sizeof(weird_limbs51) / sizeof(weird_limbs51[0])));
    size_t i2 = (size_t)(rng_u64() % (sizeof(weird_limbs51) / sizeof(weird_limbs51[0])));
    size_t i3 = (size_t)(rng_u64() % (sizeof(weird_limbs51) / sizeof(weird_limbs51[0])));
    size_t i4 = (size_t)(rng_u64() % (sizeof(weird_limbs51) / sizeof(weird_limbs51[0])));
    x.l0 = weird_limbs52[i0];
    x.l1 = weird_limbs51[i1];
    x.l2 = weird_limbs51[i2];
    x.l3 = weird_limbs51[i3];
    x.l4 = weird_limbs51[i4];
    return x;
}

static fe random_fe_maybe_weird(void)
{
    if ((rng_u64() & 1ULL) == 0) {
        return random_weird_fe();
    }
    return random_fe();
}

static int bits_len64(uint64_t x)
{
    if (x == 0) {
        return 0;
    }
    return 64 - __builtin_clzll(x);
}

static int is_in_bounds(const fe* x)
{
    return bits_len64(x->l0) <= 52 &&
           bits_len64(x->l1) <= 52 &&
           bits_len64(x->l2) <= 52 &&
           bits_len64(x->l3) <= 52 &&
           bits_len64(x->l4) <= 52;
}

static int decode_hex(const char* s, uint8_t* out, size_t out_len)
{
    size_t n = strlen(s);
    if (n != out_len * 2) {
        return -1;
    }
    for (size_t i = 0; i < out_len; i++) {
        char c0 = s[i * 2];
        char c1 = s[i * 2 + 1];
        uint8_t v0 = (uint8_t)(c0 >= '0' && c0 <= '9' ? c0 - '0' :
                               c0 >= 'a' && c0 <= 'f' ? c0 - 'a' + 10 :
                               c0 - 'A' + 10);
        uint8_t v1 = (uint8_t)(c1 >= '0' && c1 <= '9' ? c1 - '0' :
                               c1 >= 'a' && c1 <= 'f' ? c1 - 'a' + 10 :
                               c1 - 'A' + 10);
        out[i] = (uint8_t)((v0 << 4) | v1);
    }
    return 0;
}

static int test_mul64_to_128(void)
{
    fe_uint128 r = fe_mul_u64(5, 5);
    if (r.lo != 0x19 || r.hi != 0) {
        fprintf(stderr, "mul_u64 lo-range failed: %llu + %llu*(2^64)\n",
                (unsigned long long)r.lo, (unsigned long long)r.hi);
        return 0;
    }

    uint64_t a = 18014398509481983ULL;
    uint64_t b = 18014398509481983ULL;
    r = fe_mul_u64(a, b);
    if (r.lo != 0xff80000000000001ULL || r.hi != 0xfffffffffffULL) {
        fprintf(stderr, "mul_u64 hi-range failed: %llu + %llu*(2^64)\n",
                (unsigned long long)r.lo, (unsigned long long)r.hi);
        return 0;
    }

    a = 1125899906842661ULL;
    b = 2097155ULL;
    r = fe_mul_u64(a, b);
    r = fe_add_mul(r, a, b);
    r = fe_add_mul(r, a, b);
    r = fe_add_mul(r, a, b);
    r = fe_add_mul(r, a, b);
    if (r.lo != 16888498990613035ULL || r.hi != 640ULL) {
        fprintf(stderr, "add_mul wrong answer: %llu + %llu*(2^64)\n",
                (unsigned long long)r.lo, (unsigned long long)r.hi);
        return 0;
    }
    return 1;
}

static int test_set_bytes_roundtrip(void)
{
    for (int i = 0; i < TEST_ITERS; i++) {
        uint8_t in[32];
        for (size_t j = 0; j < sizeof(in); j++) {
            in[j] = (uint8_t)rng_u64();
        }
        fe x;
        if (fe_set_bytes(&x, in, sizeof(in)) != 0) {
            fprintf(stderr, "fe_set_bytes failed\n");
            return 0;
        }
        in[31] &= 0x7f;
        uint8_t out[32];
        fe_bytes(out, &x);
        if (memcmp(in, out, sizeof(in)) != 0 || !is_in_bounds(&x)) {
            fprintf(stderr, "set_bytes roundtrip mismatch\n");
            return 0;
        }
    }
    return 1;
}

static int test_multiply_distributes_over_add(void)
{
    for (int i = 0; i < TEST_ITERS; i++) {
        fe x = random_fe_maybe_weird();
        fe y = random_fe_maybe_weird();
        fe z = random_fe_maybe_weird();

        fe t1, t2, t3;
        fe_add(&t1, &x, &y);
        fe_mul(&t1, &t1, &z);

        fe_mul(&t2, &x, &z);
        fe_mul(&t3, &y, &z);
        fe_add(&t2, &t2, &t3);

        if (fe_equal(&t1, &t2) != 1 || !is_in_bounds(&t1) || !is_in_bounds(&t2)) {
            fprintf(stderr, "distributive test failed\n");
            return 0;
        }
    }
    return 1;
}

static int test_equal(void)
{
    fe x = {1, 1, 1, 1, 1};
    fe y = {5, 4, 3, 2, 1};
    if (fe_equal(&x, &x) != 1) {
        fprintf(stderr, "equality failed\n");
        return 0;
    }
    if (fe_equal(&x, &y) != 0) {
        fprintf(stderr, "inequality failed\n");
        return 0;
    }
    return 1;
}

static int test_invert(void)
{
    fe x = {1, 1, 1, 1, 1};
    fe one = FE_ONE;
    fe xinv, r;

    fe_invert(&xinv, &x);
    fe_mul(&r, &x, &xinv);
    reduce(&r);

    if (fe_equal(&r, &one) != 1) {
        fprintf(stderr, "inversion identity failed\n");
        return 0;
    }

    uint8_t bytes[32];
    for (size_t i = 0; i < sizeof(bytes); i++) {
        bytes[i] = (uint8_t)rng_u64();
    }
    fe_set_bytes(&x, bytes, sizeof(bytes));

    fe_invert(&xinv, &x);
    fe_mul(&r, &x, &xinv);
    reduce(&r);

    if (fe_equal(&r, &one) != 1) {
        fprintf(stderr, "random inversion identity failed\n");
        return 0;
    }

    fe_zero(&x);
    fe_invert(&xinv, &x);
    if (fe_equal(&xinv, &FE_ZERO) != 1) {
        fprintf(stderr, "inverting zero did not return zero\n");
        return 0;
    }

    return 1;
}

static int test_select_swap(void)
{
    fe a = {358744748052810ULL, 1691584618240980ULL, 977650209285361ULL,
            1429865912637724ULL, 560044844278676ULL};
    fe b = {84926274344903ULL, 473620666599931ULL, 365590438845504ULL,
            1028470286882429ULL, 2146499180330972ULL};

    fe c, d;
    fe_select(&c, &a, &b, 1);
    fe_select(&d, &a, &b, 0);
    if (fe_equal(&c, &a) != 1 || fe_equal(&d, &b) != 1) {
        fprintf(stderr, "select failed\n");
        return 0;
    }

    fe_swap(&c, &d, 0);
    if (fe_equal(&c, &a) != 1 || fe_equal(&d, &b) != 1) {
        fprintf(stderr, "swap failed (cond 0)\n");
        return 0;
    }

    fe_swap(&c, &d, 1);
    if (fe_equal(&c, &b) != 1 || fe_equal(&d, &a) != 1) {
        fprintf(stderr, "swap failed (cond 1)\n");
        return 0;
    }
    return 1;
}

static int test_mult32(void)
{
    for (int i = 0; i < TEST_ITERS; i++) {
        fe x = random_fe_maybe_weird();
        uint32_t y = (uint32_t)rng_u64();

        fe t1, t2, ty;
        fe_mul32(&t1, &x, y);

        fe_zero(&ty);
        ty.l0 = (uint64_t)y;
        fe_mul(&t2, &x, &ty);

        if (fe_equal(&t1, &t2) != 1 || !is_in_bounds(&t1) || !is_in_bounds(&t2)) {
            fprintf(stderr, "mult32 failed\n");
            return 0;
        }
    }
    return 1;
}

static int test_sqrt_ratio(void)
{
    struct test_vec {
        const char* u;
        const char* v;
        int was_square;
        const char* r;
    } tests[] = {
        {
            "0000000000000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000000000000000000000000000000000000000000000",
            1,
            "0000000000000000000000000000000000000000000000000000000000000000",
        },
        {
            "0000000000000000000000000000000000000000000000000000000000000000",
            "0100000000000000000000000000000000000000000000000000000000000000",
            1,
            "0000000000000000000000000000000000000000000000000000000000000000",
        },
        {
            "0100000000000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000000000000000000000000000000000000000000000",
            0,
            "0000000000000000000000000000000000000000000000000000000000000000",
        },
        {
            "0200000000000000000000000000000000000000000000000000000000000000",
            "0100000000000000000000000000000000000000000000000000000000000000",
            0,
            "3c5ff1b5d8e4113b871bd052f9e7bcd0582804c266ffb2d4f4203eb07fdb7c54",
        },
        {
            "0400000000000000000000000000000000000000000000000000000000000000",
            "0100000000000000000000000000000000000000000000000000000000000000",
            1,
            "0200000000000000000000000000000000000000000000000000000000000000",
        },
        {
            "0100000000000000000000000000000000000000000000000000000000000000",
            "0400000000000000000000000000000000000000000000000000000000000000",
            1,
            "f6ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff3f",
        },
    };

    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
        uint8_t ub[32];
        uint8_t vb[32];
        uint8_t rb[32];
        if (decode_hex(tests[i].u, ub, sizeof(ub)) != 0 ||
            decode_hex(tests[i].v, vb, sizeof(vb)) != 0 ||
            decode_hex(tests[i].r, rb, sizeof(rb)) != 0) {
            fprintf(stderr, "hex decode failed\n");
            return 0;
        }
        fe u, v, want, got;
        fe_set_bytes(&u, ub, sizeof(ub));
        fe_set_bytes(&v, vb, sizeof(vb));
        fe_set_bytes(&want, rb, sizeof(rb));

        int was_square = fe_sqrt_ratio(&got, &u, &v);
        if (fe_equal(&got, &want) != 1 || was_square != tests[i].was_square) {
            fprintf(stderr, "sqrt_ratio failed at %zu\n", i);
            return 0;
        }
    }
    return 1;
}

static int test_square_mul_consistency(void)
{
    for (int i = 0; i < TEST_ITERS; i++) {
        fe a = random_fe_maybe_weird();
        fe t1, t2;
        fe_square(&t1, &a);
        fe_mul(&t2, &a, &a);
        if (fe_equal(&t1, &t2) != 1 || !is_in_bounds(&t1) || !is_in_bounds(&t2)) {
            fprintf(stderr, "square/mul mismatch\n");
            fprintf(stderr, "a:  %llu %llu %llu %llu %llu\n",
                    (unsigned long long)a.l0, (unsigned long long)a.l1,
                    (unsigned long long)a.l2, (unsigned long long)a.l3,
                    (unsigned long long)a.l4);
            fprintf(stderr, "sq: %llu %llu %llu %llu %llu\n",
                    (unsigned long long)t1.l0, (unsigned long long)t1.l1,
                    (unsigned long long)t1.l2, (unsigned long long)t1.l3,
                    (unsigned long long)t1.l4);
            fprintf(stderr, "ml: %llu %llu %llu %llu %llu\n",
                    (unsigned long long)t2.l0, (unsigned long long)t2.l1,
                    (unsigned long long)t2.l2, (unsigned long long)t2.l3,
                    (unsigned long long)t2.l4);
            return 0;
        }
    }
    return 1;
}

int main(void)
{
    struct {
        const char* name;
        int (*fn)(void);
    } tests[] = {
        {"mul64_to_128", test_mul64_to_128},
        {"set_bytes_roundtrip", test_set_bytes_roundtrip},
        {"multiply_distributes_over_add", test_multiply_distributes_over_add},
        {"equal", test_equal},
        {"invert", test_invert},
        {"select_swap", test_select_swap},
        {"mult32", test_mult32},
        {"sqrt_ratio", test_sqrt_ratio},
        {"square_mul_consistency", test_square_mul_consistency},
    };

    int failures = 0;
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
        if (!tests[i].fn()) {
            fprintf(stderr, "[FAIL] %s\n", tests[i].name);
            failures++;
        } else {
            fprintf(stdout, "[OK]   %s\n", tests[i].name);
        }
    }

    if (failures) {
        fprintf(stderr, "%d test(s) failed\n", failures);
        return 1;
    }

    printf("All tests passed.\n");
    return 0;
}
