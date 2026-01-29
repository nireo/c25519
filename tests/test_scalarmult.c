#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../scalarmult.h"

static int check_on_curve(const point* points, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        const point* p = &points[i];
        fe XX, YY, ZZ, ZZZZ;
        fe_square(&XX, &p->x);
        fe_square(&YY, &p->y);
        fe_square(&ZZ, &p->z);
        fe_square(&ZZZZ, &ZZ);

        fe lhs;
        fe rhs;
        fe_sub(&lhs, &YY, &XX);
        fe_mul(&lhs, &lhs, &ZZ);

        fe_mul(&rhs, &XX, &YY);
        fe_mul(&rhs, &rhs, &d);
        fe_add(&rhs, &rhs, &ZZZZ);

        if (fe_equal(&lhs, &rhs) != 1) {
            fprintf(stderr, "point %zu is not on curve\n", i);
            return 0;
        }

        fe_mul(&lhs, &p->x, &p->y);
        fe_mul(&rhs, &p->z, &p->t);
        if (fe_equal(&lhs, &rhs) != 1) {
            fprintf(stderr, "point %zu has invalid T\n", i);
            return 0;
        }
    }
    return 1;
}

static void random_scalar(scalar* out)
{
    uint8_t wide[64];
    for (size_t i = 0; i < sizeof(wide); i++) {
        wide[i] = (uint8_t)(rand() & 0xff);
    }
    if (scalar_set_uniform_bytes(out, wide, sizeof(wide)) != 0) {
        fprintf(stderr, "failed to set random scalar\n");
        exit(1);
    }
}

static int test_scalar_mult_small_scalars(void)
{
    scalar z;
    scalar_zero(&z);
    point p;
    point B = point_new_generator();
    point_scalar_mult(&p, &z, &B);
    point I = point_new_identity();
    if (point_equal(&I, &p) != 1) {
        fprintf(stderr, "0*B != 0\n");
        return 0;
    }
    if (!check_on_curve(&p, 1)) {
        return 0;
    }

    scalar one;
    scalar_one(&one);
    B = point_new_generator();
    point_scalar_mult(&p, &one, &B);
    if (point_equal(&B, &p) != 1) {
        fprintf(stderr, "1*B != B\n");
        return 0;
    }
    return check_on_curve(&p, 1);
}

static int test_scalar_mult_vs_dalek(void)
{
    static const uint8_t dalek_scalar_bytes[32] = {
        219, 106, 114, 9, 174, 249, 155, 89, 69, 203, 201, 93, 92, 116, 234, 187,
        78, 115, 103, 172, 182, 98, 62, 103, 187, 136, 13, 100, 248, 110, 12, 4
    };
    static const uint8_t dalek_scalar_basepoint_bytes[32] = {
        0xf4, 0xef, 0x7c, 0x0a, 0x34, 0x55, 0x7b, 0x9f, 0x72, 0x3b, 0xb6, 0x1e, 0xf9, 0x46, 0x09, 0x91,
        0x1c, 0xb9, 0xc0, 0x6c, 0x17, 0x28, 0x2d, 0x8b, 0x43, 0x2b, 0x05, 0x18, 0x6a, 0x54, 0x3e, 0x48
    };

    scalar s;
    if (scalar_set_canonical_bytes(&s, dalek_scalar_bytes, sizeof(dalek_scalar_bytes)) != 0) {
        fprintf(stderr, "failed to set dalek scalar\n");
        return 0;
    }
    point expected;
    if (point_set_bytes(&expected, dalek_scalar_basepoint_bytes, sizeof(dalek_scalar_basepoint_bytes)) != 0) {
        fprintf(stderr, "failed to set dalek basepoint\n");
        return 0;
    }

    point p;
    point B = point_new_generator();
    point_scalar_mult(&p, &s, &B);
    if (point_equal(&expected, &p) != 1) {
        fprintf(stderr, "scalar mult does not match dalek\n");
        return 0;
    }
    return check_on_curve(&p, 1);
}

static int test_base_mult_vs_dalek(void)
{
    static const uint8_t dalek_scalar_bytes[32] = {
        219, 106, 114, 9, 174, 249, 155, 89, 69, 203, 201, 93, 92, 116, 234, 187,
        78, 115, 103, 172, 182, 98, 62, 103, 187, 136, 13, 100, 248, 110, 12, 4
    };
    static const uint8_t dalek_scalar_basepoint_bytes[32] = {
        0xf4, 0xef, 0x7c, 0x0a, 0x34, 0x55, 0x7b, 0x9f, 0x72, 0x3b, 0xb6, 0x1e, 0xf9, 0x46, 0x09, 0x91,
        0x1c, 0xb9, 0xc0, 0x6c, 0x17, 0x28, 0x2d, 0x8b, 0x43, 0x2b, 0x05, 0x18, 0x6a, 0x54, 0x3e, 0x48
    };

    scalar s;
    if (scalar_set_canonical_bytes(&s, dalek_scalar_bytes, sizeof(dalek_scalar_bytes)) != 0) {
        fprintf(stderr, "failed to set dalek scalar\n");
        return 0;
    }
    point expected;
    if (point_set_bytes(&expected, dalek_scalar_basepoint_bytes, sizeof(dalek_scalar_basepoint_bytes)) != 0) {
        fprintf(stderr, "failed to set dalek basepoint\n");
        return 0;
    }

    point p;
    point_scalar_base_mult(&p, &s);
    if (point_equal(&expected, &p) != 1) {
        fprintf(stderr, "scalar base mult does not match dalek\n");
        return 0;
    }
    return check_on_curve(&p, 1);
}

static int test_vartime_double_base_mult_vs_dalek(void)
{
    static const uint8_t dalek_scalar_bytes[32] = {
        219, 106, 114, 9, 174, 249, 155, 89, 69, 203, 201, 93, 92, 116, 234, 187,
        78, 115, 103, 172, 182, 98, 62, 103, 187, 136, 13, 100, 248, 110, 12, 4
    };
    static const uint8_t dalek_scalar_basepoint_bytes[32] = {
        0xf4, 0xef, 0x7c, 0x0a, 0x34, 0x55, 0x7b, 0x9f, 0x72, 0x3b, 0xb6, 0x1e, 0xf9, 0x46, 0x09, 0x91,
        0x1c, 0xb9, 0xc0, 0x6c, 0x17, 0x28, 0x2d, 0x8b, 0x43, 0x2b, 0x05, 0x18, 0x6a, 0x54, 0x3e, 0x48
    };

    scalar s;
    scalar zero;
    if (scalar_set_canonical_bytes(&s, dalek_scalar_bytes, sizeof(dalek_scalar_bytes)) != 0) {
        fprintf(stderr, "failed to set dalek scalar\n");
        return 0;
    }
    scalar_zero(&zero);
    point expected;
    if (point_set_bytes(&expected, dalek_scalar_basepoint_bytes, sizeof(dalek_scalar_basepoint_bytes)) != 0) {
        fprintf(stderr, "failed to set dalek basepoint\n");
        return 0;
    }

    point p;
    point B = point_new_generator();
    point_vartime_double_scalar_base_mult(&p, &s, &B, &zero);
    if (point_equal(&expected, &p) != 1) {
        fprintf(stderr, "VarTimeDoubleScalarBaseMult fails with b=0\n");
        return 0;
    }
    if (!check_on_curve(&p, 1)) {
        return 0;
    }

    point_vartime_double_scalar_base_mult(&p, &zero, &B, &s);
    if (point_equal(&expected, &p) != 1) {
        fprintf(stderr, "VarTimeDoubleScalarBaseMult fails with a=0\n");
        return 0;
    }
    return check_on_curve(&p, 1);
}

static int test_scalar_mult_distributes_over_add(void)
{
    point B = point_new_generator();
    for (int i = 0; i < 64; i++) {
        scalar x;
        scalar y;
        random_scalar(&x);
        random_scalar(&y);

        scalar z;
        scalar_add(&z, &x, &y);

        point p;
        point q;
        point r;
        point check;
        point_scalar_mult(&p, &x, &B);
        point_scalar_mult(&q, &y, &B);
        point_scalar_mult(&r, &z, &B);
        point_add(&check, &p, &q);

        if (!check_on_curve((point[]) { p, q, r, check }, 4)) {
            return 0;
        }
        if (point_equal(&check, &r) != 1) {
            fprintf(stderr, "scalar mult does not distribute over add\n");
            return 0;
        }
    }
    return 1;
}

static int test_scalar_mult_non_identity_point(void)
{
    point B = point_new_generator();
    for (int i = 0; i < 64; i++) {
        scalar x;
        random_scalar(&x);
        point p;
        point q;

        point_scalar_mult(&p, &x, &B);
        point_scalar_base_mult(&q, &x);
        if (!check_on_curve((point[]) { p, q }, 2)) {
            return 0;
        }
        if (point_equal(&p, &q) != 1) {
            fprintf(stderr, "scalar mult does not match base mult\n");
            return 0;
        }
    }
    return 1;
}

static int test_basepoint_table_generation(void)
{
    const affine_lookup_table* base_table = basepoint_table();

    proj_p1xp1 tmp1;
    proj_p2 tmp2;
    point tmp3 = point_new_generator();

    affine_lookup_table table[32];
    for (int i = 0; i < 32; i++) {
        affine_lookup_table_from_p3(&table[i], &tmp3);
        if (memcmp(&table[i], &base_table[i], sizeof(affine_lookup_table)) != 0) {
            fprintf(stderr, "basepoint table %d does not match\n", i);
            return 0;
        }

        proj_p2_from_p3(&tmp2, &tmp3);
        for (int j = 0; j < 7; j++) {
            proj_p1xp1_double(&tmp1, &tmp2);
            proj_p2_from_p1x1(&tmp2, &tmp1);
        }
        proj_p1xp1_double(&tmp1, &tmp2);
        point_from_p1x1(&tmp3, &tmp1);
        if (!check_on_curve(&tmp3, 1)) {
            return 0;
        }
    }
    return 1;
}

static int test_scalar_mult_matches_base_mult(void)
{
    point B = point_new_generator();
    for (int i = 0; i < 64; i++) {
        scalar x;
        random_scalar(&x);
        point p;
        point q;

        point_scalar_mult(&p, &x, &B);
        point_scalar_base_mult(&q, &x);
        if (!check_on_curve((point[]) { p, q }, 2)) {
            return 0;
        }
        if (point_equal(&p, &q) != 1) {
            fprintf(stderr, "scalar mult does not match base mult\n");
            return 0;
        }
    }
    return 1;
}

static int test_basepoint_naf_table_generation(void)
{
    naf_lookup_table8 table;
    point B = point_new_generator();
    naf_lookup_table8_from_p3(&table, &B);

    if (memcmp(&table, basepoint_naf_table(), sizeof(naf_lookup_table8)) != 0) {
        fprintf(stderr, "basepoint naf table does not match\n");
        return 0;
    }
    return 1;
}

static int test_vartime_double_base_mult_matches_base_mult(void)
{
    point B = point_new_generator();
    for (int i = 0; i < 64; i++) {
        scalar x;
        scalar y;
        random_scalar(&x);
        random_scalar(&y);

        point p;
        point q1;
        point q2;
        point check;
        point_vartime_double_scalar_base_mult(&p, &x, &B, &y);
        point_scalar_base_mult(&q1, &x);
        point_scalar_base_mult(&q2, &y);
        point_add(&check, &q1, &q2);

        if (!check_on_curve((point[]) { p, q1, q2, check }, 4)) {
            return 0;
        }
        if (point_equal(&p, &check) != 1) {
            fprintf(stderr, "VarTimeDoubleScalarBaseMult does not match base mult\n");
            return 0;
        }
    }
    return 1;
}

int main(void)
{
    srand(1);
    struct {
        const char* name;
        int (*fn)(void);
    } tests[] = {
        { "scalar_mult_small_scalars", test_scalar_mult_small_scalars },
        { "scalar_mult_vs_dalek", test_scalar_mult_vs_dalek },
        { "base_mult_vs_dalek", test_base_mult_vs_dalek },
        { "vartime_double_base_mult_vs_dalek", test_vartime_double_base_mult_vs_dalek },
        { "scalar_mult_distributes_over_add", test_scalar_mult_distributes_over_add },
        { "scalar_mult_non_identity_point", test_scalar_mult_non_identity_point },
        { "basepoint_table_generation", test_basepoint_table_generation },
        { "scalar_mult_matches_base_mult", test_scalar_mult_matches_base_mult },
        { "basepoint_naf_table_generation", test_basepoint_naf_table_generation },
        { "vartime_double_base_mult_matches_base_mult", test_vartime_double_base_mult_matches_base_mult },
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
