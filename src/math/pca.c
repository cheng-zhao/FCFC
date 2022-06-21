/*******************************************************************************
* pca.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "pca.h"
#include <math.h>

#define ABS(x)          (((x) >= 0) ? (x) : (-(x)))
#define SIGN(x)         (((x) >= 0) ? 1 : -1)

/*============================================================================*\
                     Functions for symmetric eigen problems
\*============================================================================*/

/*******************************************************************************
  Eigen solver for real symmetric 3x3 matrix, implemented following
  Golub & van Loan, Matrix Computations, 4th Edition
  https://jhupbooks.press.jhu.edu/title/matrix-computations
*******************************************************************************/

/******************************************************************************
Function `householder`:
  Compute factors for the Householder reflection.
  Ref: Section 5.1.3 (Algorithm 5.1.1).
Arguments:
  * `x`:        column of the matrix block to be reflected;
  * `v`:        Householder reflection vector (last component);
  * `b`:        normalization factor for Householder reflection.
******************************************************************************/
static inline void householder(const double *x, double *v, double *b) {
  const double sigma = x[1] * x[1];
  *v = x[1];
  if (sigma == 0 && x[0] >= 0) *b = 0;
  else if (sigma == 0 && x[0] < 0) *b = 2;
  else {
    const double mu = sqrt(x[0] * x[0] + sigma);
    const double u = (x[0] <= 0) ? x[0] - mu : -sigma / (x[0] + mu);
    const double u2 = u * u;
    *b = 2 * u2 / (sigma + u2);
    *v /= u;
  }
}

/******************************************************************************
Function `givens`:
  Compute factors for the Givens rotation.
  Ref: Section 5.1.8 (Algorithm 5.1.3).
Arguments:
  * `a`, `b`:   components of the vector to be rotated;
  * `c`, `s`:   elements of the rotation matrix.
******************************************************************************/
static inline void givens_fac(const double a, const double b, double *c,
    double *s) {
  if (b == 0) {
    *c = 1; *s = 0;
  }
  else {
    if (ABS(b) > ABS(a)) {
      const double tau = -a / b;
      *s = 1 / sqrt(1 + tau * tau);
      *c = *s * tau;
    }
    else {
      const double tau = -b / a;
      *c = 1 / sqrt(1 + tau * tau);
      *s = *c * tau;
    }
  }
}

/******************************************************************************
Function `tridiag`:
  Tridiagonalise a symmetric 3x3 matrix.
  Ref: Section 8.3.1 (Algorithm 8.3.1) and Section 5.1.6.
Arguments:
  * `x`:        the symmetric matrix to be tridiagonalised;
  * `q`:        orthogonal matrix for the transformation.
******************************************************************************/
static void tridiag(double x[static 6], double q[static 9]) {
  double v, b;
  householder(x + 1, &v, &b);
  double p[2];
  p[0] = b * (x[3] + x[4] * v);
  p[1] = b * (x[4] + x[5] * v);
  const double fac = b * (p[0] + p[1] * v) * 0.5;
  p[0] -= fac;
  p[1] -= fac * v;

  x[1] = sqrt(x[1] * x[1] + x[2] * x[2]);
  x[2] = 0;
  x[3] -= p[0] * 2;
  x[4] -= v * p[0] + p[1];
  x[5] -= 2 * v * p[1];

  q[0] = 1;
  q[1] = q[2] = q[3] = q[6] = 0;
  q[4] = 1 - b;
  q[5] = q[7] = -b * v;
  q[8] = 1 - b * v * v;
}

/******************************************************************************
Function `symmetric_qr`:
  The symmetric QR algorithm, for solving eigenvalues and eigenvectors.
  Ref: Section 8.3.5 (Algorithm 8.3.3).
Arguments:
  * `x`:        an input symmetric tridiagonal matrix;
  * `v`:        eigenvectors of the input matrix;
  * `eps`:      threshold for setting matrix elements to 0.
******************************************************************************/
static void symmetric_qr(double x[static 6], double v[static 9],
    const double eps) {
  if (ABS(x[1]) <= (ABS(x[0]) + ABS(x[3])) * eps) x[1] = 0;
  if (ABS(x[4]) <= (ABS(x[3]) + ABS(x[5])) * eps) x[4] = 0;
  while (x[1] != 0 || x[4] != 0) {
    if (x[1] != 0 && x[4] != 0) {       /* reduce the full 3x3 matrix */
      const double d = (x[3] - x[5]) * 0.5;
      const double t = x[4] * x[4];
      const double mu = x[5] - t / (d + SIGN(d) * sqrt(d * d + t));
      double a = x[0] - mu;
      double b = x[1];
      double c, s;
      givens_fac(a, b, &c, &s);

      /* Rotate the tridiagonal matrix: 1st iteration. */
      double c2 = c * c;
      double s2 = s * s;
      double cs = c * s;
      double x2 = c * x[2] - s * x[4];
      x[4] = s * x[2] + c * x[4];
      x[2] = x2;
      double x0 = c2 * x[0] - 2 * cs * x[1] + s2 * x[3];
      double x3 = s2 * x[0] + 2 * cs * x[1] + c2 * x[3];
      x[1] = (c2 - s2) * x[1] + cs * (x[0] - x[3]);
      x[0] = x0;
      x[3] = x3;
      /* Rotate the orthogonal matrix. */
      x0 = c * v[0] - s * v[3]; v[3] = s * v[0] + c * v[3]; v[0] = x0;
      x0 = c * v[1] - s * v[4]; v[4] = s * v[1] + c * v[4]; v[1] = x0;
      x0 = c * v[2] - s * v[5]; v[5] = s * v[2] + c * v[5]; v[2] = x0;

      /* Rotate the tridiagonal matrix: 2nd iteration. */
      a = x[1];
      b = x[2];
      givens_fac(a, b, &c, &s);
      c2 = c * c;
      s2 = s * s;
      cs = c * s;
      double x1 = c * x[1] - s * x[2];
      x[2] = s * x[1] + c * x[2];
      x[1] = x1;
      x3 = c2 * x[3] - 2 * cs * x[4] + s2 * x[5];
      double x5 = s2 * x[3] + 2 * cs * x[4] + c2 * x[5];
      x[4] = (c2 - s2) * x[4] + cs * (x[3] - x[5]);
      x[3] = x3;
      x[5] = x5;
      /* Rotate the orthogonal matrix. */
      x0 = c * v[3] - s * v[6]; v[6] = s * v[3] + c * v[6]; v[3] = x0;
      x0 = c * v[4] - s * v[7]; v[7] = s * v[4] + c * v[7]; v[4] = x0;
      x0 = c * v[5] - s * v[8]; v[8] = s * v[5] + c * v[8]; v[5] = x0;
    }
    else if (x[1] != 0) {               /* reduce the first 2x2 block */
      const double d = (x[0] - x[3]) * 0.5;
      const double t = x[1] * x[1];
      const double mu = x[3] - t / (d + SIGN(d) * sqrt(d * d + t));
      const double a = x[0] - mu;
      const double b = x[1];
      double c, s;
      givens_fac(a, b, &c, &s);
      /* Rotate the tridiagonal matrix. */
      const double c2 = c * c;
      const double s2 = s * s;
      const double cs = c * s;
      double x2 = c * x[2] - s * x[4];
      x[4] = s * x[2] + c * x[4];
      x[2] = x2;
      double x0 = c2 * x[0] - 2 * cs * x[1] + s2 * x[3];
      double x3 = s2 * x[0] + 2 * cs * x[1] + c2 * x[3];
      x[1] = (c2 - s2) * x[1] + cs * (x[0] - x[3]);
      x[0] = x0;
      x[3] = x3;
      /* Rotate the orthogonal matrix. */
      x0 = c * v[0] - s * v[3]; v[3] = s * v[0] + c * v[3]; v[0] = x0;
      x0 = c * v[1] - s * v[4]; v[4] = s * v[1] + c * v[4]; v[1] = x0;
      x0 = c * v[2] - s * v[5]; v[5] = s * v[2] + c * v[5]; v[2] = x0;
    }
    else {                              /* reduce the last 2x2 block */
      const double d = (x[3] - x[5]) * 0.5;
      const double t = x[4] * x[4];
      const double mu = x[5] - t / (d + SIGN(d) * sqrt(d * d + t));
      const double a = x[3] - mu;
      const double b = x[4];
      double c, s;
      givens_fac(a, b, &c, &s);
      /* Rotate the tridiagonal matrix. */
      const double c2 = c * c;
      const double s2 = s * s;
      const double cs = c * s;
      double x1 = c * x[1] - s * x[2];
      x[2] = s * x[1] + c * x[2];
      x[1] = x1;
      double x3 = c2 * x[3] - 2 * cs * x[4] + s2 * x[5];
      double x5 = s2 * x[3] + 2 * cs * x[4] + c2 * x[5];
      x[4] = (c2 - s2) * x[4] + cs * (x[3] - x[5]);
      x[3] = x3;
      x[5] = x5;
      /* Rotate the orthogonal matrix. */
      x1 = c * v[3] - s * v[6]; v[6] = s * v[3] + c * v[6]; v[3] = x1;
      x1 = c * v[4] - s * v[7]; v[7] = s * v[4] + c * v[7]; v[4] = x1;
      x1 = c * v[5] - s * v[8]; v[8] = s * v[5] + c * v[8]; v[5] = x1;
    }
    if (ABS(x[1]) <= (ABS(x[0]) + ABS(x[3])) * eps) x[1] = 0;
    if (ABS(x[4]) <= (ABS(x[3]) + ABS(x[5])) * eps) x[4] = 0;
  }
}

/******************************************************************************
Function `sort_eigen`:
  Sort eigenvalues in descending order, along with the eigenvectors.
Arguments:
  * `x`:        a diagonal matrix;
  * `v`:        eigenvectors of the matrix.
******************************************************************************/
static inline void sort_eigen(double x[static 6], double v[static 9]) {
  if (x[0] < x[3]) {
    if (x[3] < x[5]) {
      double tmp = x[0]; x[0] = x[5]; x[5] = tmp;
      tmp = v[0]; v[0] = v[6]; v[6] = tmp;
      tmp = v[1]; v[1] = v[7]; v[7] = tmp;
      tmp = v[2]; v[2] = v[8]; v[8] = tmp;
      return;
    }
    else {
      double tmp = x[0]; x[0] = x[3]; x[3] = tmp;
      tmp = v[0]; v[0] = v[3]; v[3] = tmp;
      tmp = v[1]; v[1] = v[4]; v[4] = tmp;
      tmp = v[2]; v[2] = v[5]; v[5] = tmp;
    }
  }
  else if (x[0] < x[5]) {
    double tmp = x[0]; x[0] = x[5]; x[5] = tmp;
    tmp = v[0]; v[0] = v[6]; v[6] = tmp;
    tmp = v[1]; v[1] = v[7]; v[7] = tmp;
    tmp = v[2]; v[2] = v[8]; v[8] = tmp;
  }
  if (x[3] < x[5]) {
    double tmp = x[3]; x[3] = x[5]; x[5] = tmp;
    tmp = v[3]; v[3] = v[6]; v[6] = tmp;
    tmp = v[4]; v[4] = v[7]; v[7] = tmp;
    tmp = v[5]; v[5] = v[8]; v[8] = tmp;
  }
}


/*============================================================================*\
                Functions for principal component analysis (PCA)
\*============================================================================*/

/******************************************************************************
Function `cov_mat`:
  Compute the covariance matrix of a point set.
Arguments:
  * `x`:        coordinates of the input points;
  * `n`:        number of input data points;
  * `cov`:      the resulting covariance matrix.
******************************************************************************/
static inline void cov_mat(real *const x[static 3], const size_t n,
    double cov[static 6]) {
  double mean[3];
  mean[0] = mean[1] = mean[2] = 0;
  for (size_t i = 0; i < n; i++) {
    mean[0] += x[0][i];
    mean[1] += x[1][i];
    mean[2] += x[2][i];
  }
  mean[0] /= n;
  mean[1] /= n;
  mean[2] /= n;

  cov[0] = cov[1] = cov[2] = cov[3] = cov[4] = cov[5] = 0;
  for (size_t i = 0; i < n; i++) {
    double d0 = x[0][i] - mean[0];
    double d1 = x[1][i] - mean[1];
    double d2 = x[2][i] - mean[2];
    cov[0] += d0 * d0;
    cov[1] += d0 * d1;
    cov[2] += d0 * d2;
    cov[3] += d1 * d1;
    cov[4] += d1 * d2;
    cov[5] += d2 * d2;
  }
  for (int i = 0; i < 6; i++) cov[i] /= (n - 1);
}

/******************************************************************************
Function `pca_vector`:
  Compute the principal components of a point set.
Arguments:
  * `x`:        coordinates of the input points;
  * `n`:        number of input data points;
  * `eps`:      tolerance for the numerial precision;
  * `v`:        the output principal component vectors.
******************************************************************************/
void pca_vector(real *const x[static 3], const size_t n, const double eps,
    double v[static 9]) {
  /* Compute the covariance matrix. */
  double cov[6];
  cov_mat(x, n, cov);

  /* Tridiagonalise the covariance matrix. */
  tridiag(cov, v);

  /* Diagonalise the covariance matrix, and compute the eigenvectors. */
  symmetric_qr(cov, v, eps);

  /* Sort eigenvectors. */
  sort_eigen(cov, v);
}

