#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../scalar.h"

static uint64_t rng_state = 0x4d595df4d0f33173ULL;

static uint64_t rng_next_u64(void)
{
    uint64_t x = rng_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng_state = x;
    return x * 2685821657736338717ULL;
}

static void rng_bytes(uint8_t* out, size_t len)
{
    size_t i = 0;
    while (i < len) {
        uint64_t v = rng_next_u64();
        size_t n = len - i;
        if (n > 8) {
            n = 8;
        }
        for (size_t j = 0; j < n; j++) {
            out[i + j] = (uint8_t)(v >> (8 * j));
        }
        i += n;
    }
}

static int scalar_from_reduced_bytes(scalar* out, const uint8_t in[32])
{
    return scalar_set_canonical_bytes(out, in, 32);
}

static void random_reduced_bytes(uint8_t out[32])
{
    rng_bytes(out, 32);
    out[31] &= 0x0f;
}

static int test_scalar_generate(void)
{
    for (size_t i = 0; i < 1000; i++) {
        uint8_t b[32];
        random_reduced_bytes(b);
        scalar sc;
        if (scalar_set_canonical_bytes(&sc, b, sizeof(b)) != 0) {
            fprintf(stderr, "failed to parse reduced scalar\n");
            return 0;
        }
        uint8_t repr[32];
        scalar_bytes(repr, &sc);
        if (!scalar_is_reduced_bytes(repr, sizeof(repr))) {
            fprintf(stderr, "generated unreduced scalar\n");
            return 0;
        }
    }
    return 1;
}

static int test_scalar_set_canonical_bytes(void)
{
    for (size_t i = 0; i < 1000; i++) {
        uint8_t in[32];
        random_reduced_bytes(in);
        scalar sc;
        if (scalar_set_canonical_bytes(&sc, in, sizeof(in)) != 0) {
            return 0;
        }
        uint8_t repr[32];
        scalar_bytes(repr, &sc);
        if (memcmp(in, repr, sizeof(in)) != 0) {
            fprintf(stderr, "round-trip bytes->scalar->bytes failed\n");
            return 0;
        }
        if (!scalar_is_reduced_bytes(repr, sizeof(repr))) {
            fprintf(stderr, "round-trip result not reduced\n");
            return 0;
        }
    }

    for (size_t i = 0; i < 1000; i++) {
        uint8_t in[32];
        random_reduced_bytes(in);
        scalar sc1;
        scalar sc2;
        if (scalar_set_canonical_bytes(&sc1, in, sizeof(in)) != 0) {
            return 0;
        }
        uint8_t repr[32];
        scalar_bytes(repr, &sc1);
        if (scalar_set_canonical_bytes(&sc2, repr, sizeof(repr)) != 0) {
            return 0;
        }
        if (!scalar_equal(&sc1, &sc2)) {
            fprintf(stderr, "scalar->bytes->scalar failed\n");
            return 0;
        }
    }

    scalar sc_one;
    scalar_one(&sc_one);

    scalar sc_minus_one;
    scalar zero;
    scalar_zero(&zero);
    scalar_subtract(&sc_minus_one, &zero, &sc_one);

    uint8_t b[32];

    scalar_bytes(b, &sc_minus_one);
    b[0] += 1;
    {
        scalar s = sc_one;
        if (scalar_set_canonical_bytes(&s, b, sizeof(b)) == 0) {
            fprintf(stderr, "accepted non-canonical scalar (b[0]+1)\n");
            return 0;
        }
        if (!scalar_equal(&s, &sc_one)) {
            fprintf(stderr, "receiver modified on error\n");
            return 0;
        }
    }

    scalar_bytes(b, &sc_minus_one);
    b[31] += 1;
    {
        scalar s = sc_one;
        if (scalar_set_canonical_bytes(&s, b, sizeof(b)) == 0) {
            fprintf(stderr, "accepted non-canonical scalar (b[31]+1)\n");
            return 0;
        }
        if (!scalar_equal(&s, &sc_one)) {
            fprintf(stderr, "receiver modified on error\n");
            return 0;
        }
    }

    scalar_bytes(b, &sc_minus_one);
    b[31] |= 0x80;
    {
        scalar s = sc_one;
        if (scalar_set_canonical_bytes(&s, b, sizeof(b)) == 0) {
            fprintf(stderr, "accepted non-canonical scalar (high bit)\n");
            return 0;
        }
        if (!scalar_equal(&s, &sc_one)) {
            fprintf(stderr, "receiver modified on error\n");
            return 0;
        }
    }

    return 1;
}

static int test_scalar_set_uniform_bytes(void)
{
    for (size_t i = 0; i < 1000; i++) {
        uint8_t lo[32];
        uint8_t wide[64];
        random_reduced_bytes(lo);
        memset(wide, 0, sizeof(wide));
        memcpy(wide, lo, sizeof(lo));

        scalar sc_wide;
        scalar sc_canon;
        if (scalar_set_uniform_bytes(&sc_wide, wide, sizeof(wide)) != 0) {
            fprintf(stderr, "SetUniformBytes failed\n");
            return 0;
        }
        if (scalar_set_canonical_bytes(&sc_canon, lo, sizeof(lo)) != 0) {
            fprintf(stderr, "SetCanonicalBytes failed\n");
            return 0;
        }
        if (!scalar_equal(&sc_wide, &sc_canon)) {
            fprintf(stderr, "SetUniformBytes mismatch for narrow input\n");
            return 0;
        }
        uint8_t repr[32];
        scalar_bytes(repr, &sc_wide);
        if (!scalar_is_reduced_bytes(repr, sizeof(repr))) {
            fprintf(stderr, "SetUniformBytes result not reduced\n");
            return 0;
        }
    }
    return 1;
}

static int test_scalar_set_bytes_with_clamping(void)
{
    for (size_t i = 0; i < 1000; i++) {
        uint8_t in[32];
        rng_bytes(in, sizeof(in));

        scalar sc_clamp;
        if (scalar_set_bytes_with_clamping(&sc_clamp, in, sizeof(in)) != 0) {
            fprintf(stderr, "SetBytesWithClamping failed\n");
            return 0;
        }

        uint8_t wide[64];
        memset(wide, 0, sizeof(wide));
        memcpy(wide, in, sizeof(in));
        wide[0] &= 248;
        wide[31] &= 63;
        wide[31] |= 64;

        scalar sc_ref;
        if (scalar_set_uniform_bytes(&sc_ref, wide, sizeof(wide)) != 0) {
            fprintf(stderr, "SetUniformBytes failed for clamp reference\n");
            return 0;
        }
        if (!scalar_equal(&sc_clamp, &sc_ref)) {
            fprintf(stderr, "SetBytesWithClamping mismatch\n");
            return 0;
        }
    }
    return 1;
}

static int test_scalar_distributive(void)
{
    for (size_t i = 0; i < 1000; i++) {
        uint8_t bx[32];
        uint8_t by[32];
        uint8_t bz[32];
        random_reduced_bytes(bx);
        random_reduced_bytes(by);
        random_reduced_bytes(bz);

        scalar x = SCALAR_ZERO;
        scalar y = SCALAR_ZERO;
        scalar z = SCALAR_ZERO;
        scalar_from_reduced_bytes(&x, bx);
        scalar_from_reduced_bytes(&y, by);
        scalar_from_reduced_bytes(&z, bz);

        scalar t1;
        scalar t2;
        scalar t3;

        scalar_add(&t1, &x, &y);
        scalar_mul(&t1, &t1, &z);

        scalar_mul(&t2, &x, &z);
        scalar_mul(&t3, &y, &z);
        scalar_add(&t2, &t2, &t3);

        if (!scalar_equal(&t1, &t2)) {
            fprintf(stderr, "distributive law failed\n");
            return 0;
        }
        uint8_t repr[32];
        scalar_bytes(repr, &t1);
        if (!scalar_is_reduced_bytes(repr, sizeof(repr))) {
            fprintf(stderr, "distributive result not reduced\n");
            return 0;
        }
    }
    return 1;
}

static int test_scalar_add_like_sub_neg(void)
{
    for (size_t i = 0; i < 1000; i++) {
        uint8_t bx[32];
        uint8_t by[32];
        random_reduced_bytes(bx);
        random_reduced_bytes(by);

        scalar x = SCALAR_ZERO;
        scalar y = SCALAR_ZERO;
        scalar_from_reduced_bytes(&x, bx);
        scalar_from_reduced_bytes(&y, by);

        scalar t1;
        scalar t2;

        scalar_subtract(&t1, &x, &y);
        scalar_negate(&t2, &y);
        scalar_add(&t2, &t2, &x);

        if (!scalar_equal(&t1, &t2)) {
            fprintf(stderr, "x-y != -y+x\n");
            return 0;
        }
        uint8_t repr[32];
        scalar_bytes(repr, &t1);
        if (!scalar_is_reduced_bytes(repr, sizeof(repr))) {
            fprintf(stderr, "sub/neg result not reduced\n");
            return 0;
        }
    }
    return 1;
}

static int test_scalar_non_adjacent_form(void)
{
    static const uint8_t bytes[32] = {
        0x1a,
        0x0e,
        0x97,
        0x8a,
        0x90,
        0xf6,
        0x62,
        0x2d,
        0x37,
        0x47,
        0x02,
        0x3f,
        0x8a,
        0xd8,
        0x26,
        0x4d,
        0xa7,
        0x58,
        0xaa,
        0x1b,
        0x88,
        0xe0,
        0x40,
        0xd1,
        0x58,
        0x9e,
        0x7b,
        0x7f,
        0x23,
        0x76,
        0xef,
        0x09,
    };

    static const int8_t expected[256] = {
        0,
        13,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        7,
        0,
        0,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        -11,
        0,
        0,
        0,
        0,
        3,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        9,
        0,
        0,
        0,
        0,
        -5,
        0,
        0,
        0,
        0,
        0,
        0,
        3,
        0,
        0,
        0,
        0,
        11,
        0,
        0,
        0,
        0,
        11,
        0,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        0,
        -3,
        0,
        0,
        0,
        0,
        9,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        -1,
        0,
        0,
        0,
        0,
        0,
        9,
        0,
        0,
        0,
        0,
        -15,
        0,
        0,
        0,
        0,
        -7,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        0,
        5,
        0,
        0,
        0,
        0,
        13,
        0,
        0,
        0,
        0,
        0,
        -3,
        0,
        0,
        0,
        0,
        -11,
        0,
        0,
        0,
        0,
        -7,
        0,
        0,
        0,
        0,
        -13,
        0,
        0,
        0,
        0,
        11,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        -15,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        7,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        5,
        0,
        0,
        0,
        0,
        0,
        13,
        0,
        0,
        0,
        0,
        0,
        0,
        11,
        0,
        0,
        0,
        0,
        0,
        15,
        0,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        -1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        7,
        0,
        0,
        0,
        0,
        0,
        -15,
        0,
        0,
        0,
        0,
        0,
        15,
        0,
        0,
        0,
        0,
        15,
        0,
        0,
        0,
        0,
        15,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
    };

    scalar s;
    if (scalar_set_canonical_bytes(&s, bytes, sizeof(bytes)) != 0) {
        fprintf(stderr, "failed to set scalar for NAF\n");
        return 0;
    }

    int8_t naf[256];
    if (scalar_non_adjacent_form(naf, &s, 5) != 0) {
        fprintf(stderr, "scalar_non_adjacent_form failed\n");
        return 0;
    }

    for (size_t i = 0; i < 256; i++) {
        if (naf[i] != expected[i]) {
            fprintf(stderr, "NAF mismatch at %zu: got %d, expected %d\n", i, naf[i],
                expected[i]);
            return 0;
        }
    }
    return 1;
}

static int test_scalar_equal(void)
{
    scalar sc_one;
    scalar_one(&sc_one);
    scalar sc_minus_one;
    scalar zero;
    scalar_zero(&zero);
    scalar_subtract(&sc_minus_one, &zero, &sc_one);

    if (scalar_equal(&sc_one, &sc_minus_one)) {
        fprintf(stderr, "scalar_equal incorrectly true\n");
        return 0;
    }
    if (!scalar_equal(&sc_minus_one, &sc_minus_one)) {
        fprintf(stderr, "scalar_equal incorrectly false\n");
        return 0;
    }
    return 1;
}

int main(void)
{
    struct {
        const char* name;
        int (*fn)(void);
    } tests[] = {
        { "scalar_generate", test_scalar_generate },
        { "scalar_set_canonical_bytes", test_scalar_set_canonical_bytes },
        { "scalar_set_uniform_bytes", test_scalar_set_uniform_bytes },
        { "scalar_set_bytes_with_clamping", test_scalar_set_bytes_with_clamping },
        { "scalar_distributive", test_scalar_distributive },
        { "scalar_add_like_sub_neg", test_scalar_add_like_sub_neg },
        { "scalar_non_adjacent_form", test_scalar_non_adjacent_form },
        { "scalar_equal", test_scalar_equal },
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
