#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_math.h"

#define TAST_CASE                                                            \
  {                                                                          \
    0.0, -0.0, 1.1, -1.1, 0.1, -0.1, -S21_INFINITY, S21_INFINITY, S21_NAN,   \
        S21_M_E, S21_M_LOG2E, S21_M_LOG10E, S21_M_LN2, S21_M_LN10, S21_M_PI, \
        S21_M_PI_4, S21_M_1_PI, S21_M_2_PI, S21_M_2_SQRTPI, S21_M_SQRT2,     \
        S21_M_SQRT1_2, 1e-10, -1e-10, 1e10, -1e10, S21_M_PI_2, 1e20, -1e20,  \
        1e-20, -1e-20, -3.0, 3.0, -4.0, 4.0, 1.0, -1.0                       \
  }  // pedago говорит это экстремальные значения

#define TAST_CASE_2                                                          \
  {                                                                          \
    0.0, -0.0, 1.1, -1.1, 0.1, -0.1, -S21_INFINITY, S21_INFINITY, S21_NAN,   \
        S21_M_E, S21_M_LOG2E, S21_M_LOG10E, S21_M_LN2, S21_M_LN10, S21_M_PI, \
        S21_M_PI_4, S21_M_1_PI, S21_M_2_PI, S21_M_2_SQRTPI, S21_M_SQRT2,     \
        S21_M_SQRT1_2                                                        \
  }

#define EPS_TEST 1e-6  // точность наших функций?

void test_for_abs(long int (*test)(int), int (*lib)(int));
void test_for_all_one_arg(long double (*test)(double), double (*lib)(double));
void test_for_all_one_arg_easy(long double (*test)(double),
                               double (*lib)(double));
void test_for_all_one_arg_interval(long double (*test)(double),
                                   double (*lib)(double));
void test_for_two_args(long double (*test)(double, double),
                       double (*lib)(double, double));

void print_result(int flag, long double a, long double b,
                  long double test_case);
void print_result_two_arg(int flag, long double a, long double b,
                          long double test_case_1, long double test_case_2);
void print_test(char *name);

START_TEST(abs_check) {
  print_test("ABS");
  test_for_abs(s21_abs, abs);
}
END_TEST

START_TEST(fabs_check) {
  print_test("FABS");
  test_for_all_one_arg(s21_fabs, fabs);
}
END_TEST

START_TEST(ceil_check) {
  print_test("CEIL");
  test_for_all_one_arg(s21_ceil, ceil);
}
END_TEST

START_TEST(floor_check) {
  print_test("FLOOR");
  test_for_all_one_arg(s21_floor, floor);
}
END_TEST

START_TEST(sin_check) {
  print_test("SIN");
  test_for_all_one_arg_easy(s21_sin, sin);
}
END_TEST

START_TEST(sin2_check) {
  print_test("SIN I");
  test_for_all_one_arg_interval(s21_sin, sin);
}
END_TEST

START_TEST(cos_check) {
  print_test("COS");
  test_for_all_one_arg_easy(s21_cos, cos);
}
END_TEST

START_TEST(cos2_check) {
  print_test("COS I");
  test_for_all_one_arg_interval(s21_cos, cos);
}
END_TEST

START_TEST(tan_check) {
  print_test("TAN");
  test_for_all_one_arg_easy(s21_tan, tan);
}
END_TEST

START_TEST(tan2_check) {
  print_test("TAN I");
  test_for_all_one_arg_interval(s21_tan, tan);
}
END_TEST

START_TEST(acos_check) {
  print_test("ACOS");
  test_for_all_one_arg(s21_acos, acos);
}
END_TEST

START_TEST(acos2_check) {
  print_test("ACOS I");
  test_for_all_one_arg_interval(s21_acos, acos);
}
END_TEST

START_TEST(asin_check) {
  print_test("ASIN");
  test_for_all_one_arg(s21_asin, asin);
}
END_TEST

START_TEST(asin2_check) {
  print_test("ASIN I");
  test_for_all_one_arg_interval(s21_asin, asin);
}
END_TEST

START_TEST(atan_check) {
  print_test("ATAN");
  test_for_all_one_arg(s21_atan, atan);
}
END_TEST

START_TEST(atan2_check) {
  print_test("ATAN I");
  test_for_all_one_arg_interval(s21_atan, atan);
}
END_TEST

START_TEST(exp_check) {
  print_test("EXP");
  test_for_all_one_arg(s21_exp, exp);
}
END_TEST

START_TEST(exp2_check) {
  print_test("EXP I");
  test_for_all_one_arg_interval(s21_exp, exp);
}
END_TEST

START_TEST(fmod_check) {
  print_test("FMOD");
  test_for_two_args(s21_fmod, fmod);
}
END_TEST

START_TEST(log_check) {
  print_test("LOG");
  test_for_all_one_arg(s21_log, log);
}
END_TEST

START_TEST(pow_check) {
  print_test("POW");
  test_for_two_args(s21_pow, pow);
}
END_TEST

START_TEST(sqrt_check) {
  print_test("SQRT");
  test_for_all_one_arg(s21_sqrt, sqrt);
}
END_TEST

int test_rdontos() {
  Suite *s = suite_create("All Part");
  TCase *tc = tcase_create("Test1");
  SRunner *sr = srunner_create(s);
  suite_add_tcase(s, tc);
  tcase_add_test(tc, sin_check);
  tcase_add_test(tc, sin2_check);
  tcase_add_test(tc, cos_check);
  tcase_add_test(tc, cos2_check);
  tcase_add_test(tc, tan_check);
  tcase_add_test(tc, tan2_check);
  tcase_add_test(tc, asin_check);
  tcase_add_test(tc, asin2_check);
  tcase_add_test(tc, acos_check);
  tcase_add_test(tc, acos2_check);
  tcase_add_test(tc, atan_check);
  tcase_add_test(tc, atan2_check);
  tcase_add_test(tc, exp_check);
  tcase_add_test(tc, exp2_check);
  tcase_add_test(tc, fmod_check);
  tcase_add_test(tc, abs_check);
  tcase_add_test(tc, fabs_check);
  tcase_add_test(tc, floor_check);
  tcase_add_test(tc, ceil_check);
  tcase_add_test(tc, log_check);
  tcase_add_test(tc, pow_check);
  tcase_add_test(tc, sqrt_check);
  srunner_run_all(sr, CK_ENV);
  int nf = srunner_ntests_failed(sr);
  srunner_free(sr);
  return nf;
}

int main(void) {
  int nf = 0;
  nf += test_rdontos();
  return 0;
}

void test_for_abs(long int (*test)(int), int (*lib)(int)) {
  int arr[] = TAST_CASE_2;
  int n = sizeof(arr) / sizeof(int);
  for (int i = 0; i < n; i++) {
    int a = lib(arr[i]);
    int b = (int)test(arr[i]);
    int flag = 0;
    if (a == b) flag = 1;

    print_result(flag, a, b, arr[i]);
  }
  printf("\n");
}

void test_for_all_one_arg(long double (*test)(double), double (*lib)(double)) {
  long double arr[] = TAST_CASE;
  int n = sizeof(arr) / sizeof(long double);
  for (int i = 0; i < n; i++) {
    long double a = lib(arr[i]);
    long double b = test(arr[i]);
    int flag = 0;
    if (s21_isnan(a) || s21_isnan(b)) {
      if (s21_isnan(a) && s21_isnan(b)) flag = 1;
    } else if (s21_isinf(a) || s21_isinf(b)) {
      if (s21_isinf(a) && s21_isinf(b)) flag = 1;
    } else {
      if (fabsl(a) < 10.0) {
        if (fabsl(a - b) <= EPS_TEST) flag = 1;
      } else if (fabsl(a - b) <= EPS_TEST * fabsl(a)) {
        flag = 1;
      }
    }
    print_result(flag, a, b, arr[i]);
  }
  printf("\n");
}

void test_for_all_one_arg_easy(long double (*test)(double),
                               double (*lib)(double)) {
  long double arr[] = TAST_CASE_2;
  int n = sizeof(arr) / sizeof(long double);
  for (int i = 0; i < n; i++) {
    long double a = lib(arr[i]);
    long double b = test(arr[i]);
    int flag = 0;
    if (s21_isnan(a) || s21_isnan(b)) {
      if (s21_isnan(a) && s21_isnan(b)) flag = 1;
    } else if (s21_isinf(a) || s21_isinf(b)) {
      if (s21_isinf(a) && s21_isinf(b)) flag = 1;
    } else {
      if (fabsl(a) < 10.0) {
        if (fabsl(a - b) <= EPS_TEST) flag = 1;
      } else if (fabsl(a - b) <= EPS_TEST * fabsl(a)) {
        flag = 1;
      }
    }
    print_result(flag, a, b, arr[i]);
  }
  printf("\n");
}

void test_for_all_one_arg_interval(long double (*test)(double),
                                   double (*lib)(double)) {
  long double d = -10 * S21_M_PI;
  long double dx = 1e-1;
  while (d < 10 * S21_M_PI) {
    long double a = lib(d);
    long double b = test(d);
    int flag = 0;
    if (s21_isnan(a) || s21_isnan(b)) {
      if (s21_isnan(a) && s21_isnan(b)) flag = 1;
    } else if (s21_isinf(a) || s21_isinf(b)) {
      if (s21_isinf(a) && s21_isinf(b)) flag = 1;
    } else {
      if (fabsl(a) < 10.0) {
        if (fabsl(a - b) <= EPS_TEST) {
          flag = 1;
        }
      } else if (fabsl(a - b) <= EPS_TEST * fabsl(a)) {
        flag = 1;
      }
    }
    print_result(flag, a, b, d);
    d += dx;
  }
  printf("\n");
}

void test_for_two_args(long double (*test)(double, double),
                       double (*lib)(double, double)) {
  long double arr_1[] = TAST_CASE;
  long double arr_2[] = TAST_CASE;
  int n = sizeof(arr_1) / sizeof(long double);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      long double a = lib(arr_1[i], arr_2[j]);
      long double b = test(arr_1[i], arr_2[j]);
      int flag = 0;
      if (s21_isnan(a) || s21_isnan(b)) {
        if (s21_isnan(a) && s21_isnan(b)) flag = 1;
      } else if (s21_isinf(a) || s21_isinf(b)) {
        if (s21_isinf(a) && s21_isinf(b)) flag = 1;
      } else {
        if (fabsl(a - b) <= EPS_TEST * fabsl(a)) flag = 1;
      }
      print_result_two_arg(flag, a, b, arr_1[i], arr_2[j]);
    }
  }
  printf("\n");
}

void print_result(int flag, long double a, long double b,
                  long double test_case) {
  if (flag) {
    printf("+");
  } else {
    printf("-\narg = %.16Le\norig = %.16Le\nour  = %.16Le\n", test_case, a, b);
    printf("diff = %Le %Le\n", fabsl(a - b), EPS_TEST * fabsl(a));
    printf("-");
  }
}
void print_test(char *name) { printf("%s:\t", name); }

void print_result_two_arg(int flag, long double a, long double b,
                          long double test_case_1, long double test_case_2) {
  if (flag) {
    printf("+");
  } else {
    printf("\narg 1 %.16Le\narg 2 %.16Le\norig = %.16Le\nour  = %.16Le\n",
           test_case_1, test_case_2, a, b);
    printf("-");
  }
}
