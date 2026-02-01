#ifndef C25519_SCALARMULT_H
#define C25519_SCALARMULT_H

#include "scalar.h"
#include "tables.h"

#ifdef __cplusplus
extern "C" {
#endif

/* basepoint_table is a set of 32 affine_lookup_table entries for 256^i * B. */
static inline const affine_lookup_table* basepoint_table(void)
{
    static affine_lookup_table table[32];
    static int initialized = 0;

    if (!initialized) {
        point p = point_new_generator();
        for (int i = 0; i < 32; i++) {
            affine_lookup_table_from_p3(&table[i], &p);
            for (int j = 0; j < 8; j++) {
                point_double(&p, &p);
            }
        }
        initialized = 1;
    }
    return table;
}

/* ScalarBaseMult sets v = x * B, where B is the canonical generator. */
static inline point* point_scalar_base_mult(point* v, const scalar* x)
{
    const affine_lookup_table* table = basepoint_table();
    int8_t digits[64];
    if (scalar_signed_radix16(digits, x) != 0) {
        return point_set_identity(v);
    }

    affine_cached multiple;
    proj_p1xp1 tmp1;
    proj_p2 tmp2;

    point_set_identity(v);
    for (int i = 1; i < 64; i += 2) {
        affine_lookup_table_select_into(&multiple, &table[i / 2], digits[i]);
        proj_p1xp1_add_affine(&tmp1, v, &multiple);
        point_from_p1x1(v, &tmp1);
    }

    proj_p2_from_p3(&tmp2, v);
    proj_p1xp1_double(&tmp1, &tmp2);
    proj_p2_from_p1x1(&tmp2, &tmp1);
    proj_p1xp1_double(&tmp1, &tmp2);
    proj_p2_from_p1x1(&tmp2, &tmp1);
    proj_p1xp1_double(&tmp1, &tmp2);
    proj_p2_from_p1x1(&tmp2, &tmp1);
    proj_p1xp1_double(&tmp1, &tmp2);
    point_from_p1x1(v, &tmp1);

    for (int i = 0; i < 64; i += 2) {
        affine_lookup_table_select_into(&multiple, &table[i / 2], digits[i]);
        proj_p1xp1_add_affine(&tmp1, v, &multiple);
        point_from_p1x1(v, &tmp1);
    }

    return v;
}

/* ScalarMult sets v = x * q. */
static inline point* point_scalar_mult(point* v, const scalar* x, const point* q)
{
    proj_lookup_table table;
    proj_lookup_table_from_p3(&table, q);

    int8_t digits[64];
    if (scalar_signed_radix16(digits, x) != 0) {
        return point_set_identity(v);
    }

    proj_cached multiple;
    proj_p1xp1 tmp1;
    proj_p2 tmp2;

    proj_lookup_table_select_into(&multiple, &table, digits[63]);
    point_set_identity(v);
    proj_p1xp1_add(&tmp1, v, &multiple);

    for (int i = 62; i >= 0; i--) {
        proj_p2_from_p1x1(&tmp2, &tmp1);
        proj_p1xp1_double(&tmp1, &tmp2);
        proj_p2_from_p1x1(&tmp2, &tmp1);
        proj_p1xp1_double(&tmp1, &tmp2);
        proj_p2_from_p1x1(&tmp2, &tmp1);
        proj_p1xp1_double(&tmp1, &tmp2);
        proj_p2_from_p1x1(&tmp2, &tmp1);
        proj_p1xp1_double(&tmp1, &tmp2);
        point_from_p1x1(v, &tmp1);
        proj_lookup_table_select_into(&multiple, &table, digits[i]);
        proj_p1xp1_add(&tmp1, v, &multiple);
    }
    point_from_p1x1(v, &tmp1);
    return v;
}

/* basepoint_naf_table is the naf_lookup_table8 for the basepoint. */
static inline const naf_lookup_table8* basepoint_naf_table(void)
{
    static naf_lookup_table8 table;
    static int initialized = 0;

    if (!initialized) {
        point p = point_new_generator();
        naf_lookup_table8_from_p3(&table, &p);
        initialized = 1;
    }
    return &table;
}

/*
 * VarTimeDoubleScalarBaseMult sets v = a * A + b * B, where B is the
 * canonical generator. Execution time depends on the inputs.
 */
static inline point* point_vartime_double_scalar_base_mult(point* v, const scalar* a, const point* A, const scalar* b)
{
    const naf_lookup_table8* base_table = basepoint_naf_table();
    naf_lookup_table5 a_table;
    naf_lookup_table5_from_p3(&a_table, A);

    int8_t a_naf[256];
    int8_t b_naf[256];
    if (scalar_non_adjacent_form(a_naf, a, 5) != 0 || scalar_non_adjacent_form(b_naf, b, 8) != 0) {
        return point_set_identity(v);
    }

    int i = 255;
    for (; i >= 0; i--) {
        if (a_naf[i] != 0 || b_naf[i] != 0) {
            break;
        }
    }
    if (i < 0) {
        return point_set_identity(v);
    }

    proj_cached mult_a;
    affine_cached mult_b;
    proj_p1xp1 tmp1;
    proj_p2 tmp2;
    proj_p2_zero(&tmp2);

    for (; i >= 0; i--) {
        proj_p1xp1_double(&tmp1, &tmp2);

        if (a_naf[i] > 0) {
            point_from_p1x1(v, &tmp1);
            naf_lookup_table5_select_into(&mult_a, &a_table, a_naf[i]);
            proj_p1xp1_add(&tmp1, v, &mult_a);
        } else if (a_naf[i] < 0) {
            point_from_p1x1(v, &tmp1);
            naf_lookup_table5_select_into(&mult_a, &a_table, (int8_t)(-a_naf[i]));
            proj_p1xp1_sub(&tmp1, v, &mult_a);
        }

        if (b_naf[i] > 0) {
            point_from_p1x1(v, &tmp1);
            naf_lookup_table8_select_into(&mult_b, base_table, b_naf[i]);
            proj_p1xp1_add_affine(&tmp1, v, &mult_b);
        } else if (b_naf[i] < 0) {
            point_from_p1x1(v, &tmp1);
            naf_lookup_table8_select_into(&mult_b, base_table, (int8_t)(-b_naf[i]));
            proj_p1xp1_sub_affine(&tmp1, v, &mult_b);
        }

        proj_p2_from_p1x1(&tmp2, &tmp1);
    }

    point_from_p2(v, &tmp2);
    return v;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* C25519_SCALARMULT_H */
