#include <stdint.h>
#include <stdio.h>

#include "../tables.h"

static int test_proj_lookup_table(void) {
  point b = point_new_generator();
  point identity = point_new_identity();

  proj_lookup_table table;
  proj_lookup_table_from_p3(&table, &b);

  proj_cached tmp1;
  proj_cached tmp2;
  proj_cached tmp3;
  proj_lookup_table_select_into(&tmp1, &table, (int8_t)6);
  proj_lookup_table_select_into(&tmp2, &table, (int8_t)-2);
  proj_lookup_table_select_into(&tmp3, &table, (int8_t)-4);

  proj_p1xp1 acc_p1xp1;
  point acc_p3 = identity;

  proj_p1xp1_add(&acc_p1xp1, &acc_p3, &tmp1);
  point_from_p1x1(&acc_p3, &acc_p1xp1);
  proj_p1xp1_add(&acc_p1xp1, &acc_p3, &tmp2);
  point_from_p1x1(&acc_p3, &acc_p1xp1);
  proj_p1xp1_add(&acc_p1xp1, &acc_p3, &tmp3);
  point_from_p1x1(&acc_p3, &acc_p1xp1);

  if (point_equal(&acc_p3, &identity) != 1) {
    fprintf(stderr, "proj_lookup_table_select_into consistency failed\n");
    return 0;
  }
  return 1;
}

static int test_affine_lookup_table(void) {
  point b = point_new_generator();
  point identity = point_new_identity();

  affine_lookup_table table;
  affine_lookup_table_from_p3(&table, &b);

  affine_cached tmp1;
  affine_cached tmp2;
  affine_cached tmp3;
  affine_lookup_table_select_into(&tmp1, &table, (int8_t)3);
  affine_lookup_table_select_into(&tmp2, &table, (int8_t)-7);
  affine_lookup_table_select_into(&tmp3, &table, (int8_t)4);

  proj_p1xp1 acc_p1xp1;
  point acc_p3 = identity;

  proj_p1xp1_add_affine(&acc_p1xp1, &acc_p3, &tmp1);
  point_from_p1x1(&acc_p3, &acc_p1xp1);
  proj_p1xp1_add_affine(&acc_p1xp1, &acc_p3, &tmp2);
  point_from_p1x1(&acc_p3, &acc_p1xp1);
  proj_p1xp1_add_affine(&acc_p1xp1, &acc_p3, &tmp3);
  point_from_p1x1(&acc_p3, &acc_p1xp1);

  if (point_equal(&acc_p3, &identity) != 1) {
    fprintf(stderr, "affine_lookup_table_select_into consistency failed\n");
    return 0;
  }
  return 1;
}

static int test_naf_lookup_table5(void) {
  point b = point_new_generator();

  naf_lookup_table5 table;
  naf_lookup_table5_from_p3(&table, &b);

  proj_cached tmp1;
  proj_cached tmp2;
  proj_cached tmp3;
  proj_cached tmp4;
  naf_lookup_table5_select_into(&tmp1, &table, (int8_t)9);
  naf_lookup_table5_select_into(&tmp2, &table, (int8_t)11);
  naf_lookup_table5_select_into(&tmp3, &table, (int8_t)7);
  naf_lookup_table5_select_into(&tmp4, &table, (int8_t)13);

  proj_p1xp1 acc_p1xp1;
  point lhs = point_new_identity();
  point rhs = point_new_identity();

  proj_p1xp1_add(&acc_p1xp1, &lhs, &tmp1);
  point_from_p1x1(&lhs, &acc_p1xp1);
  proj_p1xp1_add(&acc_p1xp1, &lhs, &tmp2);
  point_from_p1x1(&lhs, &acc_p1xp1);

  proj_p1xp1_add(&acc_p1xp1, &rhs, &tmp3);
  point_from_p1x1(&rhs, &acc_p1xp1);
  proj_p1xp1_add(&acc_p1xp1, &rhs, &tmp4);
  point_from_p1x1(&rhs, &acc_p1xp1);

  if (point_equal(&lhs, &rhs) != 1) {
    fprintf(stderr, "naf_lookup_table5 consistency failed\n");
    return 0;
  }
  return 1;
}

static int test_naf_lookup_table8(void) {
  point b = point_new_generator();

  naf_lookup_table8 table;
  naf_lookup_table8_from_p3(&table, &b);

  affine_cached tmp1;
  affine_cached tmp2;
  affine_cached tmp3;
  affine_cached tmp4;
  naf_lookup_table8_select_into(&tmp1, &table, (int8_t)49);
  naf_lookup_table8_select_into(&tmp2, &table, (int8_t)11);
  naf_lookup_table8_select_into(&tmp3, &table, (int8_t)35);
  naf_lookup_table8_select_into(&tmp4, &table, (int8_t)25);

  proj_p1xp1 acc_p1xp1;
  point lhs = point_new_identity();
  point rhs = point_new_identity();

  proj_p1xp1_add_affine(&acc_p1xp1, &lhs, &tmp1);
  point_from_p1x1(&lhs, &acc_p1xp1);
  proj_p1xp1_add_affine(&acc_p1xp1, &lhs, &tmp2);
  point_from_p1x1(&lhs, &acc_p1xp1);

  proj_p1xp1_add_affine(&acc_p1xp1, &rhs, &tmp3);
  point_from_p1x1(&rhs, &acc_p1xp1);
  proj_p1xp1_add_affine(&acc_p1xp1, &rhs, &tmp4);
  point_from_p1x1(&rhs, &acc_p1xp1);

  if (point_equal(&lhs, &rhs) != 1) {
    fprintf(stderr, "naf_lookup_table8 consistency failed\n");
    return 0;
  }
  return 1;
}

int main(void) {
  struct {
    const char *name;
    int (*fn)(void);
  } tests[] = {
      {"proj_lookup_table", test_proj_lookup_table},
      {"affine_lookup_table", test_affine_lookup_table},
      {"naf_lookup_table5", test_naf_lookup_table5},
      {"naf_lookup_table8", test_naf_lookup_table8},
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
