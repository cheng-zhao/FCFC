/*******************************************************************************
* cspline.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "cspline.h"
#include <stdlib.h>

/*******************************************************************************
  Implementation of the "natural" cubic spline interpolation algorithm.
  ref: https://doi.org/10.5281/zenodo.3611922
  see also: https://arxiv.org/abs/2001.09253

  The original source codes are released under a CC0 license by Haysn Hornbeck.
*******************************************************************************/

/******************************************************************************
Function `cspline_ypp`:
  Compute the second derivative of sample points.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `y`:        y coordinates of the sample points;
  * `n`:        number of sample points.
Return:
  The second derivative of `y`.
******************************************************************************/
double *cspline_ypp(const double *x, const double *y, const size_t n) {
  if (n < 2) return NULL;
  double *cp, *ypp;
  if (!(cp = malloc(n * sizeof(double)))) return NULL;
  if (!(ypp = malloc(n * sizeof(double)))) {
    free(cp); return NULL;
  }

  double newx = x[1];
  double newy = y[1];
  double c = x[1] - x[0];
  double newd = (y[1] - y[0]) / c;

  /* natural condition: second derivative = 0 */
  cp[0] = cp[n - 1] = ypp[0] = ypp[n - 1] = 0;

  /* forward substitution */
  size_t j = 1;
  while (j < n - 1) {
    double oldx = newx;
    double oldy = newy;
    double a = c;
    double oldd = newd;

    newx = x[j + 1];
    newy = y[j + 1];
    c = newx - oldx;
    newd = (newy - oldy) / c;

    double b = (c + a) * 2;
    double invd = 1 / (b - a * cp[j - 1]);
    double d = (newd - oldd) * 6;

    ypp[j] = (d - a * ypp[j - 1]) * invd;
    cp[j] = c * invd;

    j += 1;
  }

  /* backward substitution */
  while (j) {
    j -= 1;
    ypp[j] -= cp[j] * ypp[j + 1];
  }

  free(cp);
  return ypp;
}

/******************************************************************************
Function `cspline_eval`:
  Evaluate the cubic spline interpolation.
  It can be further optimised for uniformly spaced sample points.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `y`:        y coordinates of the sample points;
  * `ypp`:      second derivative of y;
  * `xv`:       x coordinate of the value to be evaluated;
  * `i`:        index of the value to be evaluated.
Return:
  The value for the given x coordinate.
******************************************************************************/
double cspline_eval(const double *x, const double *y, const double *ypp,
    const double xv, const size_t i) {
  size_t j = i + 1;
  double ba = x[j] - x[i];
  double xa = xv - x[i];
  double bx = x[j] - xv;
  double ba2 = ba * ba;

  double lower = xa * y[j] + bx * y[i];
  double c = (xa * xa - ba2) * xa * ypp[j];
  double d = (bx * bx - ba2) * bx * ypp[i];

  /* 1/6 = 0x1.5555555555555p-3 */
  return (lower + 0x1.5555555555555p-3 * (c + d)) / ba;
}

