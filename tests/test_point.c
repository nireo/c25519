#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../point.h"

static int decode_hex(const char* s, uint8_t* out, size_t out_len)
{
    size_t n = strlen(s);
    if (n != out_len * 2) {
        return -1;
    }
    for (size_t i = 0; i < out_len; i++) {
        char c0 = s[i * 2];
        char c1 = s[i * 2 + 1];
        uint8_t v0 = (uint8_t)(c0 >= '0' && c0 <= '9' ? c0 - '0' : c0 >= 'a' && c0 <= 'f' ? c0 - 'a' + 10
                                                                                          : c0 - 'A' + 10);
        uint8_t v1 = (uint8_t)(c1 >= '0' && c1 <= '9' ? c1 - '0' : c1 >= 'a' && c1 <= 'f' ? c1 - 'a' + 10
                                                                                          : c1 - 'A' + 10);
        out[i] = (uint8_t)((v0 << 4) | v1);
    }
    return 0;
}

static void encode_hex(char* out, size_t out_len, const uint8_t* in, size_t in_len)
{
    static const char* hex = "0123456789abcdef";
    if (out_len < in_len * 2 + 1) {
        return;
    }
    for (size_t i = 0; i < in_len; i++) {
        out[i * 2] = hex[in[i] >> 4];
        out[i * 2 + 1] = hex[in[i] & 0x0f];
    }
    out[in_len * 2] = '\0';
}

static int check_on_curve(const Point* points, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        const Point* p = &points[i];
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

static int test_generator(void)
{
    Point B = point_new_generator();

    static const char* exp_x = "1ad5258f602d56c9b2a7259560c72c695cdcd6fd31e2a4c0fe536ecdd3366921";
    static const char* exp_y = "5866666666666666666666666666666666666666666666666666666666666666";

    uint8_t bx[32];
    uint8_t by[32];
    fe_bytes(bx, &B.x);
    fe_bytes(by, &B.y);

    char hx[65];
    char hy[65];
    encode_hex(hx, sizeof(hx), bx, sizeof(bx));
    encode_hex(hy, sizeof(hy), by, sizeof(by));

    if (strcmp(hx, exp_x) != 0) {
        fprintf(stderr, "wrong B.x: got %s, expected %s\n", hx, exp_x);
        return 0;
    }
    if (strcmp(hy, exp_y) != 0) {
        fprintf(stderr, "wrong B.y: got %s, expected %s\n", hy, exp_y);
        return 0;
    }
    if (fe_equal(&B.z, &FE_ONE) != 1) {
        fprintf(stderr, "wrong B.z\n");
        return 0;
    }
    return check_on_curve(&B, 1);
}

static int test_add_sub_neg_on_base_point(void)
{
    Point B = point_new_generator();
    Point I = point_new_identity();
    Point check_lhs;
    Point check_rhs;

    point_add(&check_lhs, &B, &B);
    projP2 tmp_p2;
    projP1xP1 tmp_p1x1;
    projP2_from_p3(&tmp_p2, &B);
    projP1xP1_double(&tmp_p1x1, &tmp_p2);
    point_from_p1x1(&check_rhs, &tmp_p1x1);
    if (point_equal(&check_lhs, &check_rhs) != 1) {
        fprintf(stderr, "B + B != [2]B\n");
        return 0;
    }
    if (!check_on_curve((Point[]) { check_lhs, check_rhs }, 2)) {
        return 0;
    }

    point_subtract(&check_lhs, &B, &B);
    Point Bneg;
    point_negate(&Bneg, &B);
    point_add(&check_rhs, &B, &Bneg);
    if (point_equal(&check_lhs, &check_rhs) != 1) {
        fprintf(stderr, "B - B != B + (-B)\n");
        return 0;
    }
    if (point_equal(&I, &check_lhs) != 1) {
        fprintf(stderr, "B - B != 0\n");
        return 0;
    }
    if (point_equal(&I, &check_rhs) != 1) {
        fprintf(stderr, "B + (-B) != 0\n");
        return 0;
    }
    if (!check_on_curve((Point[]) { check_lhs, check_rhs, Bneg }, 3)) {
        return 0;
    }
    return 1;
}

static int test_invalid_encodings(void)
{
    const char* invalid = "efffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f";
    uint8_t buf[32];
    if (decode_hex(invalid, buf, sizeof(buf)) != 0) {
        fprintf(stderr, "invalid hex in test\n");
        return 0;
    }

    Point p = point_new_generator();
    Point orig = p;

    if (point_set_bytes(&p, buf, sizeof(buf)) == 0) {
        fprintf(stderr, "expected error for invalid point\n");
        return 0;
    }
    if (point_equal(&p, &orig) != 1) {
        fprintf(stderr, "Point modified after invalid decode\n");
        return 0;
    }
    return check_on_curve(&p, 1);
}

static int test_non_canonical_points(void)
{
    struct test_case {
        const char* name;
        const char* encoding;
        const char* canonical;
    };

    static const struct test_case tests[] = {
        { "y=1,sign-",
            "0100000000000000000000000000000000000000000000000000000000000080",
            "0100000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+1,sign-",
            "eeffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0100000000000000000000000000000000000000000000000000000000000000" },
        { "y=p-1,sign-",
            "ecffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "ecffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f" },
        { "y=p,sign+",
            "edffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0000000000000000000000000000000000000000000000000000000000000000" },
        { "y=p,sign-",
            "edffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0000000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+1,sign+",
            "eeffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0100000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+3,sign+",
            "f0ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0300000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+3,sign-",
            "f0ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0300000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+4,sign+",
            "f1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0400000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+4,sign-",
            "f1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0400000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+5,sign+",
            "f2ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0500000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+5,sign-",
            "f2ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0500000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+6,sign+",
            "f3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0600000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+6,sign-",
            "f3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0600000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+9,sign+",
            "f6ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0900000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+9,sign-",
            "f6ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0900000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+10,sign+",
            "f7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0a00000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+10,sign-",
            "f7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0a00000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+14,sign+",
            "fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0e00000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+14,sign-",
            "fbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0e00000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+15,sign+",
            "fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "0f00000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+15,sign-",
            "fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "0f00000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+16,sign+",
            "fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "1000000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+16,sign-",
            "fdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "1000000000000000000000000000000000000000000000000000000000000080" },
        { "y=p+18,sign+",
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f",
            "1200000000000000000000000000000000000000000000000000000000000000" },
        { "y=p+18,sign-",
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
            "1200000000000000000000000000000000000000000000000000000000000080" },
    };

    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
        uint8_t enc[32];
        uint8_t canon[32];
        if (decode_hex(tests[i].encoding, enc, sizeof(enc)) != 0 || decode_hex(tests[i].canonical, canon, sizeof(canon)) != 0) {
            fprintf(stderr, "bad test vector %s\n", tests[i].name);
            return 0;
        }

        Point p1;
        Point p2;
        if (point_set_bytes(&p1, enc, sizeof(enc)) != 0) {
            fprintf(stderr, "error decoding non-canonical point %s\n", tests[i].name);
            return 0;
        }
        if (point_set_bytes(&p2, canon, sizeof(canon)) != 0) {
            fprintf(stderr, "error decoding canonical point %s\n", tests[i].name);
            return 0;
        }
        if (point_equal(&p1, &p2) != 1) {
            fprintf(stderr, "equivalent points not equal %s\n", tests[i].name);
            return 0;
        }

        uint8_t out[32];
        char out_hex[65];
        point_bytes(out, &p1);
        encode_hex(out_hex, sizeof(out_hex), out, sizeof(out));
        if (strcmp(out_hex, tests[i].canonical) != 0) {
            fprintf(stderr, "re-encoding mismatch %s\n", tests[i].name);
            return 0;
        }
        if (!check_on_curve((Point[]) { p1, p2 }, 2)) {
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
        { "generator", test_generator },
        { "add_sub_neg_on_base_point", test_add_sub_neg_on_base_point },
        { "invalid_encodings", test_invalid_encodings },
        { "non_canonical_points", test_non_canonical_points },
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
