#ifndef __TABLES_H__
#define __TABLES_H__

#include "point.h"

typedef struct {
    proj_cached points[8];
} proj_lookup_table;

typedef struct {
    affine_cached points[8];
} affine_lookup_table;

typedef struct {
    proj_cached points[8];
} naf_lookup_table5;

typedef struct {
    affine_cached points[64];
} naf_lookup_table8;

static inline int tables_ct_eq_u8(uint8_t a, uint8_t b)
{
    return fe_ct_is_zero_u32((uint32_t)(a ^ b));
}

static inline proj_cached* proj_cached_select(proj_cached* v, const proj_cached* a, const proj_cached* b, int cond)
{
    fe_select(&v->YplusX, &a->YplusX, &b->YplusX, cond);
    fe_select(&v->YminusX, &a->YminusX, &b->YminusX, cond);
    fe_select(&v->Z, &a->Z, &b->Z, cond);
    fe_select(&v->T2d, &a->T2d, &b->T2d, cond);
    return v;
}

static inline affine_cached* affine_cached_select(affine_cached* v, const affine_cached* a, const affine_cached* b, int cond)
{
    fe_select(&v->YplusX, &a->YplusX, &b->YplusX, cond);
    fe_select(&v->YminusX, &a->YminusX, &b->YminusX, cond);
    fe_select(&v->T2d, &a->T2d, &b->T2d, cond);
    return v;
}

static inline proj_cached* proj_cached_condneg(proj_cached* v, int cond)
{
    fe_swap(&v->YplusX, &v->YminusX, cond);
    fe neg_t2d;
    fe_neg(&neg_t2d, &v->T2d);
    fe_select(&v->T2d, &neg_t2d, &v->T2d, cond);
    return v;
}

static inline affine_cached* affine_cached_condneg(affine_cached* v, int cond)
{
    fe_swap(&v->YplusX, &v->YminusX, cond);
    fe neg_t2d;
    fe_neg(&neg_t2d, &v->T2d);
    fe_select(&v->T2d, &neg_t2d, &v->T2d, cond);
    return v;
}

static inline proj_p1xp1* proj_p1xp1_add_affine(proj_p1xp1* v, const point* p, const affine_cached* q)
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
    fe_set(&ZZ2, &p->z);
    fe_add(&ZZ2, &ZZ2, &ZZ2);

    fe_sub(&v->X, &PP, &MM);
    fe_add(&v->Y, &PP, &MM);
    fe_add(&v->Z, &ZZ2, &TT2d);
    fe_sub(&v->T, &ZZ2, &TT2d);
    return v;
}

static inline proj_p1xp1* proj_p1xp1_sub_affine(proj_p1xp1* v, const point* p, const affine_cached* q)
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
    fe_set(&ZZ2, &p->z);
    fe_add(&ZZ2, &ZZ2, &ZZ2);

    fe_sub(&v->X, &PP, &MM);
    fe_add(&v->Y, &PP, &MM);
    fe_sub(&v->Z, &ZZ2, &TT2d);
    fe_add(&v->T, &ZZ2, &TT2d);
    return v;
}

static inline proj_lookup_table* proj_lookup_table_from_p3(proj_lookup_table* v, const point* q)
{
    proj_cached_from_p3(&v->points[0], q);
    point tmp_p3;
    proj_p1xp1 tmp_p1xp1;
    for (int i = 0; i < 7; i++) {
        proj_p1xp1_add(&tmp_p1xp1, q, &v->points[i]);
        point_from_p1x1(&tmp_p3, &tmp_p1xp1);
        proj_cached_from_p3(&v->points[i + 1], &tmp_p3);
    }
    return v;
}

static inline affine_lookup_table* affine_lookup_table_from_p3(affine_lookup_table* v, const point* q)
{
    affine_cached_from_p3(&v->points[0], q);
    point tmp_p3;
    proj_p1xp1 tmp_p1xp1;
    for (int i = 0; i < 7; i++) {
        proj_p1xp1_add_affine(&tmp_p1xp1, q, &v->points[i]);
        point_from_p1x1(&tmp_p3, &tmp_p1xp1);
        affine_cached_from_p3(&v->points[i + 1], &tmp_p3);
    }
    return v;
}

static inline naf_lookup_table5* naf_lookup_table5_from_p3(naf_lookup_table5* v, const point* q)
{
    proj_cached_from_p3(&v->points[0], q);
    point q2;
    point_add(&q2, q, q);
    point tmp_p3;
    proj_p1xp1 tmp_p1xp1;
    for (int i = 0; i < 7; i++) {
        proj_p1xp1_add(&tmp_p1xp1, &q2, &v->points[i]);
        point_from_p1x1(&tmp_p3, &tmp_p1xp1);
        proj_cached_from_p3(&v->points[i + 1], &tmp_p3);
    }
    return v;
}

static inline naf_lookup_table8* naf_lookup_table8_from_p3(naf_lookup_table8* v, const point* q)
{
    affine_cached_from_p3(&v->points[0], q);
    point q2;
    point_add(&q2, q, q);
    point tmp_p3;
    proj_p1xp1 tmp_p1xp1;
    for (int i = 0; i < 63; i++) {
        proj_p1xp1_add_affine(&tmp_p1xp1, &q2, &v->points[i]);
        point_from_p1x1(&tmp_p3, &tmp_p1xp1);
        affine_cached_from_p3(&v->points[i + 1], &tmp_p3);
    }
    return v;
}

static inline proj_cached* proj_lookup_table_select_into(proj_cached* dest, const proj_lookup_table* v, int8_t x)
{
    int8_t xmask = (int8_t)(x >> 7);
    uint8_t xabs = (uint8_t)((x + xmask) ^ xmask);

    proj_cached_zero(dest);
    for (int j = 1; j <= 8; j++) {
        int cond = tables_ct_eq_u8(xabs, (uint8_t)j);
        proj_cached_select(dest, &v->points[j - 1], dest, cond);
    }
    proj_cached_condneg(dest, (int)(xmask & 1));
    return dest;
}

static inline affine_cached* affine_lookup_table_select_into(affine_cached* dest, const affine_lookup_table* v, int8_t x)
{
    int8_t xmask = (int8_t)(x >> 7);
    uint8_t xabs = (uint8_t)((x + xmask) ^ xmask);

    affine_cached_zero(dest);
    for (int j = 1; j <= 8; j++) {
        int cond = tables_ct_eq_u8(xabs, (uint8_t)j);
        affine_cached_select(dest, &v->points[j - 1], dest, cond);
    }
    affine_cached_condneg(dest, (int)(xmask & 1));
    return dest;
}

static inline proj_cached* naf_lookup_table5_select_into(proj_cached* dest, const naf_lookup_table5* v, int8_t x)
{
    *dest = v->points[(uint8_t)x / 2];
    return dest;
}

static inline affine_cached* naf_lookup_table8_select_into(affine_cached* dest, const naf_lookup_table8* v, int8_t x)
{
    *dest = v->points[(uint8_t)x / 2];
    return dest;
}

#endif
