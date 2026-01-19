#ifndef __POINT_H__
#define __POINT_H__

#include <stddef.h>
#include <stdint.h>

#include "fe.h"

typedef struct projP1xP1 {
    fe X;
    fe Y;
    fe Z;
    fe T;
} projP1xP1;

typedef struct projP2 {
    fe X;
    fe Y;
    fe Z;
} projP2;

typedef struct Point {
    fe x;
    fe y;
    fe z;
    fe t;
} Point;

typedef struct projCached {
    fe YplusX;
    fe YminusX;
    fe Z;
    fe T2d;
} projCached;

typedef struct affineCached {
    fe YplusX;
    fe YminusX;
    fe T2d;
} affineCached;

static const fe d = {
    929955233495203ULL,
    466365720129213ULL,
    1662059464998953ULL,
    2033849074728123ULL,
    1442794654840575ULL,
};

static const fe d2 = {
    1859910466990425ULL,
    932731440258426ULL,
    1072319116312658ULL,
    1815898335770999ULL,
    633789495995903ULL,
};

static inline projP2* projP2_zero(projP2* v)
{
    fe_zero(&v->X);
    fe_one(&v->Y);
    fe_one(&v->Z);
    return v;
}

static inline projCached* projCached_zero(projCached* v)
{
    fe_one(&v->YplusX);
    fe_one(&v->YminusX);
    fe_one(&v->Z);
    fe_zero(&v->T2d);
    return v;
}

static inline affineCached* affineCached_zero(affineCached* v)
{
    fe_one(&v->YplusX);
    fe_one(&v->YminusX);
    fe_zero(&v->T2d);
    return v;
}

static inline Point* point_set(Point* v, const Point* u)
{
    *v = *u;
    return v;
}

static inline Point* point_set_identity(Point* v)
{
    fe_zero(&v->x);
    fe_one(&v->y);
    fe_one(&v->z);
    fe_zero(&v->t);
    return v;
}

static inline int point_set_bytes(Point* v, const uint8_t* x, size_t len);

static inline Point* point_set_generator(Point* v)
{
    static const uint8_t gen[32] = {
        0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
    };
    Point tmp;
    if (point_set_bytes(&tmp, gen, sizeof(gen)) != 0) {
        return NULL;
    }
    *v = tmp;
    return v;
}

static inline Point point_new_identity(void)
{
    Point p;
    point_set_identity(&p);
    return p;
}

static inline Point point_new_generator(void)
{
    Point p;
    point_set_generator(&p);
    return p;
}

static inline uint8_t* point_bytes(uint8_t out[32], const Point* v)
{
    fe z_inv;
    fe x;
    fe y;
    fe_invert(&z_inv, &v->z);
    fe_mul(&x, &v->x, &z_inv);
    fe_mul(&y, &v->y, &z_inv);
    fe_bytes(out, &y);
    out[31] |= (uint8_t)(fe_is_negative(&x) << 7);
    return out;
}

static inline int point_set_bytes(Point* v, const uint8_t* x, size_t len)
{
    if (len != 32) {
        return -1;
    }

    fe y;
    if (fe_set_bytes(&y, x, len) != 0) {
        return -1;
    }

    fe y2;
    fe_square(&y2, &y);

    fe u;
    fe_sub(&u, &y2, &FE_ONE);

    fe vv;
    fe_mul(&vv, &y2, &d);
    fe_add(&vv, &vv, &FE_ONE);

    fe xx;
    if (!fe_sqrt_ratio(&xx, &u, &vv)) {
        return -1;
    }

    fe xx_neg;
    fe_neg(&xx_neg, &xx);
    fe_select(&xx, &xx_neg, &xx, (int)((x[31] >> 7) & 1U));

    Point tmp;
    fe_set(&tmp.x, &xx);
    fe_set(&tmp.y, &y);
    fe_one(&tmp.z);
    fe_mul(&tmp.t, &xx, &y);

    *v = tmp;
    return 0;
}

static inline projP2* projP2_from_p1x1(projP2* v, const projP1xP1* p)
{
    fe_mul(&v->X, &p->X, &p->T);
    fe_mul(&v->Y, &p->Y, &p->Z);
    fe_mul(&v->Z, &p->Z, &p->T);
    return v;
}

static inline projP2* projP2_from_p3(projP2* v, const Point* p)
{
    fe_set(&v->X, &p->x);
    fe_set(&v->Y, &p->y);
    fe_set(&v->Z, &p->z);
    return v;
}

static inline Point* point_from_p1x1(Point* v, const projP1xP1* p)
{
    fe_mul(&v->x, &p->X, &p->T);
    fe_mul(&v->y, &p->Y, &p->Z);
    fe_mul(&v->z, &p->Z, &p->T);
    fe_mul(&v->t, &p->X, &p->Y);
    return v;
}

static inline Point* point_from_p2(Point* v, const projP2* p)
{
    fe_mul(&v->x, &p->X, &p->Z);
    fe_mul(&v->y, &p->Y, &p->Z);
    fe_square(&v->z, &p->Z);
    fe_mul(&v->t, &p->X, &p->Y);
    return v;
}

static inline projCached* projCached_from_p3(projCached* v, const Point* p)
{
    fe_add(&v->YplusX, &p->y, &p->x);
    fe_sub(&v->YminusX, &p->y, &p->x);
    fe_set(&v->Z, &p->z);
    fe_mul(&v->T2d, &p->t, &d2);
    return v;
}

static inline affineCached* affineCached_from_p3(affineCached* v, const Point* p)
{
    fe_add(&v->YplusX, &p->y, &p->x);
    fe_sub(&v->YminusX, &p->y, &p->x);
    fe_mul(&v->T2d, &p->t, &d2);

    fe inv_z;
    fe_invert(&inv_z, &p->z);
    fe_mul(&v->YplusX, &v->YplusX, &inv_z);
    fe_mul(&v->YminusX, &v->YminusX, &inv_z);
    fe_mul(&v->T2d, &v->T2d, &inv_z);
    return v;
}

static inline projP1xP1* projP1xP1_add(projP1xP1* v, const Point* p, const projCached* q)
{
    fe YplusX;
    fe YminusX;
    fe PP;
    fe MM;
    fe TT2d;
    fe ZZ2;

    fe_add(&YplusX, &p->y, &p->x);
    fe_sub(&YminusX, &p->y, &p->x);

    fe_mul(&PP, &YplusX, &q->YplusX);
    fe_mul(&MM, &YminusX, &q->YminusX);
    fe_mul(&TT2d, &p->t, &q->T2d);
    fe_mul(&ZZ2, &p->z, &q->Z);

    fe_add(&ZZ2, &ZZ2, &ZZ2);

    fe_sub(&v->X, &PP, &MM);
    fe_add(&v->Y, &PP, &MM);
    fe_add(&v->Z, &ZZ2, &TT2d);
    fe_sub(&v->T, &ZZ2, &TT2d);
    return v;
}

static inline projP1xP1* projP1xP1_sub(projP1xP1* v, const Point* p, const projCached* q)
{
    fe YplusX;
    fe YminusX;
    fe PP;
    fe MM;
    fe TT2d;
    fe ZZ2;

    fe_add(&YplusX, &p->y, &p->x);
    fe_sub(&YminusX, &p->y, &p->x);

    fe_mul(&PP, &YplusX, &q->YminusX);
    fe_mul(&MM, &YminusX, &q->YplusX);
    fe_mul(&TT2d, &p->t, &q->T2d);
    fe_mul(&ZZ2, &p->z, &q->Z);

    fe_add(&ZZ2, &ZZ2, &ZZ2);

    fe_sub(&v->X, &PP, &MM);
    fe_add(&v->Y, &PP, &MM);
    fe_sub(&v->Z, &ZZ2, &TT2d);
    fe_add(&v->T, &ZZ2, &TT2d);
    return v;
}

static inline projP1xP1* projP1xP1_double(projP1xP1* v, const projP2* p)
{
    fe A;
    fe B;
    fe C;
    fe D;
    fe E;
    fe F;
    fe G;
    fe H;
    fe t0;

    fe_square(&A, &p->X);
    fe_square(&B, &p->Y);
    fe_square(&C, &p->Z);
    fe_add(&C, &C, &C);
    fe_neg(&D, &A);
    fe_add(&t0, &p->X, &p->Y);
    fe_square(&E, &t0);
    fe_sub(&E, &E, &A);
    fe_sub(&E, &E, &B);
    fe_add(&G, &D, &B);
    fe_sub(&F, &G, &C);
    fe_sub(&H, &D, &B);

    fe_set(&v->X, &E);
    fe_set(&v->Y, &H);
    fe_set(&v->Z, &G);
    fe_set(&v->T, &F);
    return v;
}

static inline Point* point_add(Point* v, const Point* p, const Point* q)
{
    projCached q_cached;
    projP1xP1 result;
    Point tmp;

    projCached_from_p3(&q_cached, q);
    projP1xP1_add(&result, p, &q_cached);
    point_from_p1x1(&tmp, &result);
    *v = tmp;
    return v;
}

static inline Point* point_subtract(Point* v, const Point* p, const Point* q)
{
    projCached q_cached;
    projP1xP1 result;
    Point tmp;

    projCached_from_p3(&q_cached, q);
    projP1xP1_sub(&result, p, &q_cached);
    point_from_p1x1(&tmp, &result);
    *v = tmp;
    return v;
}

static inline Point* point_double(Point* v, const Point* p)
{
    projP2 p2;
    projP1xP1 p1;
    Point tmp;

    projP2_from_p3(&p2, p);
    projP1xP1_double(&p1, &p2);
    point_from_p1x1(&tmp, &p1);
    *v = tmp;
    return v;
}

static inline Point* point_negate(Point* v, const Point* p)
{
    fe_neg(&v->x, &p->x);
    fe_set(&v->y, &p->y);
    fe_set(&v->z, &p->z);
    fe_neg(&v->t, &p->t);
    return v;
}

static inline int point_equal(const Point* p, const Point* q)
{
    fe x1;
    fe x2;
    fe y1;
    fe y2;
    fe_mul(&x1, &p->x, &q->z);
    fe_mul(&x2, &q->x, &p->z);
    fe_mul(&y1, &p->y, &q->z);
    fe_mul(&y2, &q->y, &p->z);
    return fe_equal(&x1, &x2) & fe_equal(&y1, &y2);
}

static inline Point NewIdentityPoint(void)
{
    return point_new_identity();
}

static inline Point NewGeneratorPoint(void)
{
    return point_new_generator();
}

static inline uint8_t* Point_Bytes(uint8_t out[32], const Point* v)
{
    return point_bytes(out, v);
}

static inline Point* Point_Set(Point* v, const Point* u)
{
    return point_set(v, u);
}

static inline Point* Point_Add(Point* v, const Point* p, const Point* q)
{
    return point_add(v, p, q);
}

static inline Point* Point_Subtract(Point* v, const Point* p, const Point* q)
{
    return point_subtract(v, p, q);
}

static inline Point* Point_Double(Point* v, const Point* p)
{
    return point_double(v, p);
}

static inline Point* Point_Negate(Point* v, const Point* p)
{
    return point_negate(v, p);
}

static inline int Point_Equal(const Point* p, const Point* q)
{
    return point_equal(p, q);
}

static inline Point* Point_SetBytes(Point* v, const uint8_t* x, size_t len)
{
    return point_set_bytes(v, x, len) == 0 ? v : NULL;
}

static inline projP2* projP2_Zero(projP2* v)
{
    return projP2_zero(v);
}

static inline projCached* projCached_Zero(projCached* v)
{
    return projCached_zero(v);
}

static inline affineCached* affineCached_Zero(affineCached* v)
{
    return affineCached_zero(v);
}

static inline projP1xP1* projP1xP1_Add(projP1xP1* v, const Point* p, const projCached* q)
{
    return projP1xP1_add(v, p, q);
}

static inline projP1xP1* projP1xP1_Sub(projP1xP1* v, const Point* p, const projCached* q)
{
    return projP1xP1_sub(v, p, q);
}

static inline projP1xP1* projP1xP1_Double(projP1xP1* v, const projP2* p)
{
    return projP1xP1_double(v, p);
}

static inline projP2* projP2_FromP1xP1(projP2* v, const projP1xP1* p)
{
    return projP2_from_p1x1(v, p);
}

static inline projP2* projP2_FromP3(projP2* v, const Point* p)
{
    return projP2_from_p3(v, p);
}

static inline projCached* projCached_FromP3(projCached* v, const Point* p)
{
    return projCached_from_p3(v, p);
}

static inline affineCached* affineCached_FromP3(affineCached* v, const Point* p)
{
    return affineCached_from_p3(v, p);
}

#endif // !__POINT_H__
