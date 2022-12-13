#include "s21_math.h"

int s21_isnan(double x) { return (x != x); }

int s21_isinf(double x) {
  int result = 0;
  if (x == S21_INFINITY) result = 1;
  if (x == -S21_INFINITY) result = -1;
  return result;
}

long double s21_sin(double x) {
  long double result = 0.0;
  if (s21_isnan(x) || s21_isinf(x)) {
    result = S21_NAN;
  } else {
    if (s21_fabs(x) > S21_M_PI + S21_M_PI) x = s21_fmod(x, S21_M_PI + S21_M_PI);
    long double n = x, sum = 0.L;
    int i = 1;
    do {
      sum += n;
      n *= -1.0 * x * x / ((2 * i) * (2 * i + 1));
      ++i;
    } while ((n > 0.0 ? n : -n) > EPS * EPS);
    result = sum;
  }
  return result;
}

long double s21_cos(double x) { return s21_sin(S21_M_PI_2 + x); }

long double s21_tan(double x) {
  long double result = 0.0;
  if (s21_isnan(x) || s21_isinf(x)) {
    result = S21_NAN;
  } else {
    result = s21_sin(x) / s21_cos(x);
  }

  return result;
}

// dpasty

#define __HI(x) *(1 + (int32_t*)&x)
#define __LO(x) *(int32_t*)&x

long double my_fmod(long double x1, long double y1);

long double s21_ceil(double x) {
  long double result;
  if (s21_isnan(x)) {
    result = S21_NAN;
  } else if (s21_isinf(x)) {
    result = S21_INFINITY;
  } else if (x == 0) {
    result = 0;
  } else if (s21_fmod(x, 1.) == 0.) {
    result = x;
  } else if (x > 0) {
    result = x + 1 - s21_fmod(x, 1.);
  } else {
    result = x - s21_fmod(x, 1.);
  }
  return result;
}

long double s21_fabs(double x) { return (x >= 0) ? x : -x; }

long double s21_floor(double x) {
  long double result;
  if (s21_isnan(x)) {
    result = S21_NAN;
  } else if (s21_isinf(x)) {
    result = S21_INFINITY;
  } else if (x == 0) {
    result = 0;
  } else if (s21_fmod(x, 1.) == 0.) {
    result = x;
  } else if (x > 0) {
    result = x - s21_fmod(x, 1.);
  } else {
    result = x - 1. - s21_fmod(x, 1.);
  }
  return result;
}

long double s21_fmod(double x, double y) {
  long double result;
  if (s21_isnan(x) || s21_isnan(y)) {
    result = S21_NAN;
  } else if (s21_isinf(x)) {
    result = S21_NAN;
  } else if (y == 0) {
    result = S21_NAN;
  } else if (s21_isinf(y)) {
    result = x;
  } else if ((x == 0 && y != 0) || (x / y == 0)) {
    result = 0;
  } else if (x == y) {
    result = 0;
  } else if (s21_fabs(y) > s21_fabs(x)) {
    result = x;
  } else {
    int sgn = 1;
    if (x < 0) {
      sgn = -1;
      x *= sgn;
    }
    result = my_fmod(x, y);
    result *= sgn;
  }
  return result;
}

long int s21_abs(int x) {
  long int result = (x >= 0) ? x : -x;
  return result;
}

#define get_high_part(x) *(1 + (int32_t*)&x)
#define get_low_part(x) *(int32_t*)&x

int32_t s21_ilogb(int32_t h) {
  int32_t res;
  res = (h >> 20) - 1023;
  return res;
}

void align(int32_t* hp, uint32_t* lp) {
  int32_t h = *hp;
  uint32_t l = *lp;
  h = 0x00100000 | (0x000fffff & h);
  *hp = h;
  *lp = l;
}

double convert_back(int32_t hx, uint32_t lx, int32_t iy) {
  double res;
  if ((hx | lx) == 0) {
    res = 0;
  } else {
    while (hx < 0x00100000) { /* normalize x */
      hx = hx + hx + (lx >> 31);
      lx = lx + lx;
      iy -= 1;
    }
    hx = ((hx - 0x00100000) | ((iy + 1023) << 20));
    get_high_part(res) = hx;
    get_low_part(res) = lx;
  }
  return res;
}

int fmod_cycle(int n, int32_t* hxp, uint32_t* lxp, int32_t hy, uint32_t ly) {
  int zero_flag = 0;
  int32_t hz;
  uint32_t lz;
  int32_t hx = *hxp;
  uint32_t lx = *lxp;

  while (n--) {
    hz = hx - hy;
    lz = lx - ly;
    if (lx < ly) hz -= 1;
    if (hz < 0) {
      hx = hx + hx + (lx >> 31);
      lx = lx + lx;
    } else {
      if ((hz | lz) == 0) zero_flag = 1;
      hx = hz + hz + (lz >> 31);
      lx = lz + lz;
    }
  }
  hz = hx - hy;
  lz = lx - ly;
  if (lx < ly) hz -= 1;
  if (hz >= 0) {
    hx = hz;
    lx = lz;
  }
  *hxp = hx;
  *lxp = lx;
  return zero_flag;
}

long double my_fmod(long double x1, long double y1) {
  double x = (double)x1;
  double y = (double)y1;
  int32_t hx, hy, ix, iy;
  uint32_t lx, ly;

  hx = get_high_part(x); /* high word of x */
  lx = get_low_part(x);  /* low  word of x */
  hy = get_high_part(y); /* high word of y */
  ly = get_low_part(y);  /* low  word of y */
  // sx = hx & 0x80000000;  /* sign of x */
  hx &= 0x7fffffff; /* |x| */
  hy &= 0x7fffffff; /* |y| */

  ix = s21_ilogb(hx);
  iy = s21_ilogb(hy);

  align(&hx, &lx);
  align(&hy, &ly);

  if (fmod_cycle(ix - iy, &hx, &lx, hy, ly)) {
    x = 0;
  } else {
    x = convert_back(hx, lx, iy);
  }
  return (long double)x;
}

// mcampfir

long double precision_exp(long double c);

long double s21_pow_mcamp(long double a, int b) {
  long double result = 1;
  for (int i = 0; i < b; i++) {
    result *= a;
  }
  return result;
}

long double s21_exp(double c) {
  long double result = 1;
  if (s21_isnan(c)) {
    result = S21_NAN;
  } else if (s21_isinf(c) == 1) {
    result = S21_INFINITY;
  } else if (s21_isinf(c) == -1) {
    result = 0;
  } else {
    long double x = c;
    int sign = 0;
    if (c < 0) {
      x *= (-1);
      sign++;
    }
    if (x > 645) {
      result = precision_exp(x);
    } else {
      long double fact = 1;
      long double p = 1;

      int i = 1;
      do {
        fact *= i;
        p = (s21_pow_mcamp(x, i) / (fact));
        result += p;
        i++;
      } while (p >= EPS);
    }

    if (sign != 0) result = (1.0 / result);
  }
  return result;
}

long double precision_exp(long double c) {
  long double result = c;
  int j = 0;
  int flag = 0;
  if (c > 2 || c < -2) {
    flag = 1;
    while (result > 2) {
      result /= 2;
      j++;
    }
  }
  if (flag == 1) {
    result = s21_exp(result);
    for (int i = 0; i < j; i++) {
      result = s21_pow_mcamp(result, 2);
    }
  } else {
    result = s21_exp(result);
  }
  return result;
}

// jsera

long double s21_asin(double x) {
  long double result = 0.0;
  long double value = x;

  if (x < 0) value *= -1;

  if (value > 0.525 && value <= 1) {
    result = S21_M_PI_2 - 2 * s21_asin(s21_sqrt((1 - value) / 2));
  } else if (value <= 0.525) {
    long double previous = 0.0;
    result = previous = value;
    for (int i = 1; i <= 57; i++) {
      long double apper_fraction = value * value * (2 * (long double)i - 1) *
                                   (2 * (long double)i - 1) *
                                   (2 * (long double)i);
      long double lower_fraction =
          4 * (long double)i * (long double)i * (2 * (long double)i + 1);

      result = result + apper_fraction / lower_fraction * previous;

      previous = apper_fraction / lower_fraction * previous;
    }
  } else {
    result = S21_NAN;
  }

  if (x < 0) result *= -1;

  return result;
}

long double s21_acos(double x) {
  long double result = S21_M_PI_2 - s21_asin(x);
  return result;
}

long double atan_05(long double x) {
  int series_size = 30;

  long double term = 1 * x;
  long double res = term;
  for (int i = 2; i <= series_size; i++) {
    term = term * (-1) * (2 * i - 3) / (2 * i - 1) * x * x;
    res += term;
  }
  return res;
}

long double s21_atan(double x) {
  long double atan_tbl[] = {S21_ATAN_05, S21_ATAN_10, S21_ATAN_15,
                            S21_ATAN_INF};
  long double res = 0;
  int sgn = 1;

  if (x < 0) {
    sgn = -1;
    x = -x;
  }
  if (0 <= x && x <= 7. / 16) {
    res = atan_05(x);
  } else if (7. / 16 < x && x <= 11. / 16) {
    res = atan_tbl[0] + atan_05((x - 0.5L) / (1 + x / 2));
  } else if (11. / 16 < x && x <= 19. / 16) {
    res = atan_tbl[1] + atan_05((x - 1.L) / (1 + x));
  } else if (19. / 16 < x && x <= 39.L / 16) {
    res = atan_tbl[2] + atan_05((x - 3.L / 2) / (1 + 3 * x / 2));
  } else if (39. / 16 < x) {
    res = atan_tbl[3] + atan_05(-1. / x);
  } else {
    res = atan_05(x);
  }
  return res * sgn;
}

long double pow_x0_yint(long double x, long double y) {
  long double res;
  if (y > 0) {
    if (s21_fabs(s21_fmod(y, 2)) == 1)
      res = x;
    else
      res = +0.0L;
  } else {
    if (s21_fabs(s21_fmod(y, 2)) == 1)
      res = x == -0.0L ? -S21_HUGE_VALL : S21_HUGE_VALL;
    else
      res = +S21_HUGE_VALL;
  }
  return res;
}

long double pow_xinf(long double y) {
  long double res = 0;

  if (s21_fabs(s21_fmod(y, 2)) == 1) {
    if (y < 0)
      res = -0.0L;
    else
      res = -S21_INFINITY;
  } else {
    if (y < 0)
      res = +0.0L;
    else
      res = +S21_INFINITY;
  }
  return res;
}

long double s21_pow(double x, double y) {
  long double res;
  if (x == 1.0L) {
    res = 1.0L;
  } else if (y == 0.0L) {
    res = 1.0L;
  } else if (y != y || x != x) {
    res = S21_NAN;
  } else if ((x == 0.0L || x == -0.0L) && s21_ceil(y) == y) {
    res = pow_x0_yint(x, y);
  } else if (x == -1.0L && (y == S21_INFINITY || y == -S21_INFINITY)) {
    res = 1.0L;
  } else if (s21_fabs(x) < 1 && y == -S21_INFINITY) {
    res = S21_INFINITY;
  } else if (s21_fabs(x) > 1 && y == -S21_INFINITY) {
    res = +0.0L;
  } else if (s21_fabs(x) < 1 && y == S21_INFINITY) {
    res = +0.0L;
  } else if (s21_fabs(x) > 1 && y == S21_INFINITY) {
    res = S21_INFINITY;
  } else if (x == +S21_INFINITY) {
    if (y < 0)
      res = +0.0;
    else
      res = +S21_INFINITY;
  } else if (x == -S21_INFINITY) {
    res = pow_xinf(y);
  } else if (s21_ceil(y) == y && x < 0) {
    if (s21_fabs(s21_fmod(y, 2)) == 1)
      res = -s21_exp(y * s21_log(-x));
    else
      res = s21_exp(y * s21_log(-x));
  } else {
    res = s21_exp(y * s21_log(x));
  }
  return res;
}

long double s21_sqrt(double x) {
  long double res;
  if (x == -S21_INFINITY) {
    res = S21_NAN;
  } else {
    res = s21_pow(x, 0.5);
  }
  return res;
}

long double log_taylor(long double x) {
  long double res = 0;
  long double term = -1;
  for (int i = 1; i < 50; i++) {
    term *= -x;
    res += term / i;
  }
  return res;
}

long double log_taylor2(long double x) {
  long double res = S21_M_LN2;
  long double term = -1;
  for (int i = 1; i < 50; i++) {
    term *= (-x/2);
    res += term / i;
  }
  return res;
}

long double s21_log(double x) {
  int k = 0;
  long double res;

  if (x < 0) {
    res = S21_NAN;
  } else if (x == 0) {
    res = -S21_INFINITY;
  } else if (x == S21_INFINITY) {
    res = S21_INFINITY;
  } else {
    long double buf = x;
    if (buf > 2) {
      do {
        buf = buf / 2;
        k++;
      } while (buf > 2);
    }
    if (buf < 1) {
      do {
        buf = buf * 2;
        k--;
      } while (buf < 1);
    }
    if (buf <= 1.5)
      res = (log_taylor(buf - 1) + k * S21_M_LN2);
    else
      res = (log_taylor2(buf - 2) + k * S21_M_LN2);
  }
  return res;
}
