/*******************************************************************************
* benchmark/struct/balltree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "create_data.h"
#include "data_struct.h"
#include "pca.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <math.h>

/*============================================================================*\
                       Definition of the data structures
\*============================================================================*/

#ifdef BENCHMARK_TREE_AS_ARRAY
typedef struct {
  size_t n;                             /* number of data points   */
  real cen[3];                          /* centre of the node      */
  real r;                               /* radius of the node      */
  real *x[3];                           /* pointers to coordinates */
} BLT;
#else
typedef struct balltree_struct {
#if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  real *cen;
#else
  real cen[3];                          /* centre of the node      */
#endif
  real r;                               /* radius of the node      */
  size_t n;                             /* number of data points   */
  real *x[3];                           /* pointers to coordinates */
  struct balltree_struct *left;         /* left child              */
  struct balltree_struct *right;        /* right child             */
} BLT;
#endif

/*============================================================================*\
                           Functions for data sorting
\*============================================================================*/

/******************************************************************************
Function `extreme_points`:
  Find the extreme points along a given direction.
Arguments:
  * `x`:        pointer to coordinates of input points;
  * `istart`:   starting index (included) of data points;
  * `iend`:     ending index (not included) of data points;
  * `p`:        vector for the direction to be examined;
  * `iex`:      indices of the extreme points.
******************************************************************************/
static inline void extreme_points(real **x, const size_t istart,
    const size_t iend, const double *p, size_t iex[2]) {
  real min = DBL_MAX, max = -DBL_MAX;
  iex[0] = iex[1] = 0;
  for (size_t i = istart; i < iend; i++) {
    const real dist = x[0][i] * p[0] + x[1][i] * p[1] + x[2][i] * p[2];
    if (min > dist) {
      min = dist;
      iex[0] = i;
    }
    if (max < dist) {
      max = dist;
      iex[1] = i;
    }
  }
}

/******************************************************************************
Function `dist_plane_compare`:
  Compare the distances to a plane for two points.
Arguments:
  * `x`:        pointer to coordinates of input points;
  * `i`:        index of the first data point;
  * `j`:        index of the second data point;
  * `n`:        normal vector of the plane.
Return:
  The order of the two data points.
******************************************************************************/
static int dist_plane_compare(real **x, const size_t i, const size_t j,
    const double *n) {
  real di = x[0][i] * n[0] + x[1][i] * n[1] + x[2][i] * n[2];
  real dj = x[0][j] * n[0] + x[1][j] * n[1] + x[2][j] * n[2];
  if (di < dj) return -1;
  if (di > dj) return 1;
  return 0;
}


/* Import the function for partition data by the median. */

#define SWAP_DATA(x,i,j) {                                              \
  real _tmp = (x)[0][i]; (x)[0][i] = (x)[0][j]; (x)[0][j] = _tmp;     \
  _tmp = (x)[1][i]; (x)[1][i] = (x)[1][j]; (x)[1][j] = _tmp;            \
  _tmp = (x)[2][i]; (x)[2][i] = (x)[2][j]; (x)[2][j] = _tmp;            \
}

#ifdef QSELECT_COMPARE
  #undef QSELECT_COMPARE
#endif
#ifdef QSELECT_DTYPE
  #undef QSELECT_DTYPE
#endif
#ifdef QSELECT_SWAP
  #undef QSELECT_SWAP
#endif

#define QSELECT_COMPARE(x,i,j,n)        dist_plane_compare(x,i,j,n)
#define QSELECT_DTYPE                   real *
#define QSELECT_SWAP                    SWAP_DATA

#include "qselect.c"

#undef QSELECT_DTYPE
#undef QSELECT_COMPARE
#undef QSELECT_SWAP


/*============================================================================*\
                     Functions for minimum bounding sphere
\*============================================================================*/

/******************************************************************************
Function `minimum_circle`:
  Compute the minimum bounding circle for 3 points.
  Ref: https://realtimecollisiondetection.net/blog/?p=20
Arguments:
  * `x`:        coordinates of the input dataset;
  * `a`:        index of the first point;
  * `b`:        index of the second point;
  * `c`:        index of the third point;
  * `cen`:      centre of the circle;
  * `r2`:       squared radius of the circle.
******************************************************************************/
static inline void minimum_circle(real **x, const size_t a, const size_t b,
    const size_t c, double *cen, double *r2) {
  double AB[3], AC[3];
  AB[0] = x[0][b] - x[0][a];
  AB[1] = x[1][b] - x[1][a];
  AB[2] = x[2][b] - x[2][a];
  AC[0] = x[0][c] - x[0][a];
  AC[1] = x[1][c] - x[1][a];
  AC[2] = x[2][c] - x[2][a];
  const double ABAB = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
  const double ABAC = AB[0] * AC[0] + AB[1] * AC[1] + AB[2] * AC[2];
  const double ACAC = AC[0] * AC[0] + AC[1] * AC[1] + AC[2] * AC[2];
  const double d = 2 * (ABAB * ACAC - ABAC * ABAC);
  /* The 3 points are collinear. */
  if (d > -DBL_TOL && d < DBL_TOL) {
    *r2 = 0;
    for (int i = 0; i < 3; i++) {
      double xmin, xmax;
      if (x[i][a] < x[i][b]) {
        xmin = x[i][a]; xmax = x[i][b];
      }
      else {
        xmin = x[i][b]; xmax = x[i][a];
      }
      if (xmin > x[i][c]) xmin = x[i][c];
      else if (xmax < x[i][c]) xmax = x[i][c];
      cen[i] = (xmin + xmax) * 0.5;
      const double dx = cen[i] - xmin;
      *r2 += dx * dx;
    }
  }
  /* The 3 points form a triangle. */
  else {
    const double s = ACAC * (ABAB - ABAC) / d;
    const double t = ABAB * (ACAC - ABAC) / d;
    if (s <= 0) {
      cen[0] = (x[0][a] + x[0][c]) * 0.5;
      cen[1] = (x[1][a] + x[1][c]) * 0.5;
      cen[2] = (x[2][a] + x[2][c]) * 0.5;
    }
    else if (t <= 0) {
      cen[0] = (x[0][a] + x[0][b]) * 0.5;
      cen[1] = (x[1][a] + x[1][b]) * 0.5;
      cen[2] = (x[2][a] + x[2][b]) * 0.5;
    }
    else if (s + t >= 1) {
      cen[0] = (x[0][b] + x[0][c]) * 0.5;
      cen[1] = (x[1][b] + x[1][c]) * 0.5;
      cen[2] = (x[2][b] + x[2][c]) * 0.5;
      const double dx = cen[0] - x[0][b];
      const double dy = cen[1] - x[1][b];
      const double dz = cen[2] - x[2][b];
      *r2 = dx * dx + dy * dy + dz * dz;
      return;
    }
    else {
      cen[0] = x[0][a] + s * AB[0] + t * AC[0];
      cen[1] = x[1][a] + s * AB[1] + t * AC[1];
      cen[2] = x[2][a] + s * AB[2] + t * AC[2];
    }
    const double dx = cen[0] - x[0][a];
    const double dy = cen[1] - x[1][a];
    const double dz = cen[2] - x[2][a];
    *r2 = dx * dx + dy * dy + dz * dz;
  }
}

/******************************************************************************
Function `minimum_sphere`:
  Compute the minimum bounding sphere for 4 points.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `a`:        index of the first point;
  * `b`:        index of the second point;
  * `c`:        index of the third point;
  * `d`:        index of the fourth point;
  * `cen`:      centre of the sphere;
  * `r2`:       squared radius of the sphere.
******************************************************************************/
static void minimum_sphere(real **x, const size_t a, const size_t b,
    const size_t c, const size_t d, double *cen, double *r2) {
  double AB[3], AC[3], AD[3];
  AB[0] = x[0][b] - x[0][a];
  AB[1] = x[1][b] - x[1][a];
  AB[2] = x[2][b] - x[2][a];
  AC[0] = x[0][c] - x[0][a];
  AC[1] = x[1][c] - x[1][a];
  AC[2] = x[2][c] - x[2][a];
  AD[0] = x[0][d] - x[0][a];
  AD[1] = x[1][d] - x[1][a];
  AD[2] = x[2][d] - x[2][a];
  const double ABAB = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
  const double ABAC = AB[0] * AC[0] + AB[1] * AC[1] + AB[2] * AC[2];
  const double ABAD = AB[0] * AD[0] + AB[1] * AD[1] + AB[2] * AD[2];
  const double ACAC = AC[0] * AC[0] + AC[1] * AC[1] + AC[2] * AC[2];
  const double ACAD = AC[0] * AD[0] + AC[1] * AD[1] + AC[2] * AD[2];
  const double ADAD = AD[0] * AD[0] + AD[1] * AD[1] + AD[2] * AD[2];
  const double v = 2 * (ABAB * (ACAC * ADAD - ACAD * ACAD)
      + ABAC * (ABAD * ACAD - ABAC * ADAD)
      + ABAD * (ABAC * ACAD - ABAD * ACAC));
  /* The 4 points are coplanar. */
  if (v > -DBL_TOL && v < DBL_TOL) {
    /* Create the sphere for the circle first, and then add the last point. */
    minimum_circle(x, a, b, c, cen, r2);
    const double dx = cen[0] - x[0][d];
    const double dy = cen[1] - x[1][d];
    const double dz = cen[2] - x[2][d];
    double dist = dx * dx + dy * dy + dz * dz;
    /* Update the sphere if necessary. */
    if (dist > *r2) {
      dist = sqrt(dist);
      double r = sqrt(*r2);
      r = (r + dist) * 0.5;
      const double fac = r / dist;
      cen[0] = x[0][d] + dx * fac;
      cen[1] = x[1][d] + dy * fac;
      cen[2] = x[2][d] + dz * fac;
      *r2 = r * r;
    }
  }
  else {
    double s = (ACAC * ADAD * (ABAB - ABAC - ABAD)
        + ACAD * (ABAD * ACAC + ABAC * ADAD - ABAB * ACAD)) / v;
    double t = (ABAB * ADAD * (ACAC - ABAC - ACAD)
        + ABAD * (ABAB * ACAD + ABAC * ADAD - ACAC * ABAD)) / v;
    double u = (ABAB * ACAC * (ADAD - ABAD - ACAD)
        + ABAC * (ABAB * ACAD + ACAC * ABAD - ABAC * ADAD)) / v;
    if (s > 0 && t > 0 && u > 0 && s + t + u < 1) {     /* circumsphere */
      cen[0] = x[0][a] + s * AB[0] + t * AC[0] + u * AD[0];
      cen[1] = x[1][a] + s * AB[1] + t * AC[1] + u * AD[1];
      cen[2] = x[2][a] + s * AB[2] + t * AC[2] + u * AD[2];
      const double dx = cen[0] - x[0][a];
      const double dy = cen[1] - x[1][a];
      const double dz = cen[2] - x[2][a];
      *r2 = dx * dx + dy * dy + dz * dz;
    }
    else {
      if (s <= 0) {
        minimum_circle(x, a, c, d, cen, r2);
        const double dx = cen[0] - x[0][b];
        const double dy = cen[1] - x[1][b];
        const double dz = cen[2] - x[2][b];
        if (dx * dx + dy * dy + dz * dz < *r2 + DBL_TOL) return;
      }
      if (t <= 0) {
        minimum_circle(x, a, b, d, cen, r2);
        const double dx = cen[0] - x[0][c];
        const double dy = cen[1] - x[1][c];
        const double dz = cen[2] - x[2][c];
        if (dx * dx + dy * dy + dz * dz < *r2 + DBL_TOL) return;
      }
      if (u <= 0) {
        minimum_circle(x, a, b, c, cen, r2);
        const double dx = cen[0] - x[0][d];
        const double dy = cen[1] - x[1][d];
        const double dz = cen[2] - x[2][d];
        if (dx * dx + dy * dy + dz * dz < *r2 + DBL_TOL) return;
      }
      minimum_circle(x, b, c, d, cen, r2);
    }
  }
}

/*============================================================================*\
                        Functions for tree construction
\*============================================================================*/

#if defined(BENCHMARK_TREE_AS_ARRAY) || defined (BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `balltree_nnode`:
  Compute the upper limit of the number of nodes for a balanced ball tree.
Arguments:
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Upper limit of the number of ball tree nodes.
******************************************************************************/
static inline size_t balltree_nnode(const size_t ndata, const size_t nleaf) {
  if (ndata == 0 || nleaf == 0) return 0;
  size_t size = ndata / nleaf;
  if (size >= BENCHMARK_TREE_MAX_SIZE) return 0;
  /* All the nodes are full. */
  if (size * nleaf == ndata && (size & (size - 1)) == 0) return (size << 1) - 1;

  /* Compute the next power of 2, see
   * http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  size |= size >> 1;
  size |= size >> 2;
  size |= size >> 4;
  size |= size >> 8;
  size |= size >> 16;
  size++;
  size += (size == 0);
  return (size << 1) - 1;
}
#endif

#if defined(BENCHMARK_TREE_AS_ARRAY) || defined (BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `balltree_init`:
  Initialise a node of the ball tree.
Arguments:
  * `root`:     root of the tree;
  * `inode`:    index of the current node;
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  True of the node has children; false otherwise.
******************************************************************************/
static bool balltree_init(BLT *root, const size_t inode, real **x,
    const size_t idx, const size_t ndata, const size_t nleaf) {
  root[inode].n = ndata;
  root[inode].x[0] = x[0] + idx;
  root[inode].x[1] = x[1] + idx;
  root[inode].x[2] = x[2] + idx;
#ifdef BENCHMARK_TREE_PREALLOC
  root[inode].left = root[inode].right = NULL;
#endif

  /* Find the directions with largest variances using PCA. */
  double p[9] = {0,0,0,0,0,0,0,0,0};
  if (ndata < 3)  p[0] = p[4] = p[8] = 1;
  else pca_vector(root[inode].x, ndata, BENCHMARK_PCA_EPSILON, p);

  if (ndata > 4) {
    /* Find the extreme points along the principal direction. */
    size_t iex[2];
    extreme_points(root[inode].x, 0, ndata, p, iex);

    if (iex[0] == iex[1]) {     /* there is only a single point. */
      root[inode].cen[0] = root[inode].x[0][0];
      root[inode].cen[1] = root[inode].x[1][0];
      root[inode].cen[2] = root[inode].x[2][0];
      root[inode].r = 0;
    }
    else {
      /* Move extreme points to the beginning of the dataset. */
      SWAP_DATA(root[inode].x, 0, iex[0]);
      SWAP_DATA(root[inode].x, 1, iex[1]);

      /* Find extreme points along the second principal direction. */
      extreme_points(root[inode].x, 2, ndata, p + 3, iex);

      if (iex[0] == iex[1]) {     /* points are along a line */
        root[inode].cen[0] = (root[inode].x[0][0] + root[inode].x[0][1]) * 0.5;
        root[inode].cen[1] = (root[inode].x[1][0] + root[inode].x[1][1]) * 0.5;
        root[inode].cen[2] = (root[inode].x[2][0] + root[inode].x[2][1]) * 0.5;
        const real dx = root[inode].cen[0] - root[inode].x[0][0];
        const real dy = root[inode].cen[1] - root[inode].x[1][0];
        const real dz = root[inode].cen[2] - root[inode].x[2][0];
        root[inode].r = sqrt(dx * dx + dy * dy + dz * dz);
      }
      else {
        /* Move extreme points to the beginning of the dataset. */
        SWAP_DATA(root[inode].x, 2, iex[0]);
        SWAP_DATA(root[inode].x, 3, iex[1]);

        /* Initialise the bounding sphere with the four extreme points. */
        double r2 = 0;
#ifdef SINGLE_PREC
        double cen[3];
        minimum_sphere(root[inode].x, 0, 1, 2, 3, cen, &r2);
        root[inode].cen[0] = cen[0];
        root[inode].cen[1] = cen[1];
        root[inode].cen[2] = cen[2];
#else
        minimum_sphere(root[inode].x, 0, 1, 2, 3, root[inode].cen, &r2);
#endif
        root[inode].r = sqrt(r2);

        /* Check the rest of the points and enlarge the sphere if necessary. */
        for (size_t i = 4; i < ndata; i++) {
          const real dx = root[inode].cen[0] - root[inode].x[0][i];
          const real dy = root[inode].cen[1] - root[inode].x[1][i];
          const real dz = root[inode].cen[2] - root[inode].x[2][i];
          real dist = dx * dx + dy * dy + dz * dz;
          if (dist > r2) {
            dist = sqrt(dist);
            root[inode].r = (root[inode].r + dist) * 0.5;
            const real fac = root[inode].r / dist;
            root[inode].cen[0] = root[inode].x[0][i] + dx * fac;
            root[inode].cen[1] = root[inode].x[1][i] + dy * fac;
            root[inode].cen[2] = root[inode].x[2][i] + dz * fac;
            r2 = root[inode].r * root[inode].r;
          }
        }
      }
    }
  }
  else if (ndata == 4) {
    double r2 = 0;
#ifdef SINGLE_PREC
    double cen[3];
    minimum_sphere(root[inode].x, 0, 1, 2, 3, cen, &r2);
    root[inode].cen[0] = cen[0];
    root[inode].cen[1] = cen[1];
    root[inode].cen[2] = cen[2];
#else
    minimum_sphere(root[inode].x, 0, 1, 2, 3, root[inode].cen, &r2);
#endif
    root[inode].r = sqrt(r2);
  }
  else if (ndata == 3) {
    double r2 = 0;
#ifdef SINGLE_PREC
    double cen[3];
    minimum_sphere(root[inode].x, 0, 1, 2, 3, cen, &r2);
    root[inode].cen[0] = cen[0];
    root[inode].cen[1] = cen[1];
    root[inode].cen[2] = cen[2];
#else
    minimum_circle(root[inode].x, 0, 1, 2, root[inode].cen, &r2);
#endif
    root[inode].r = sqrt(r2);
  }
  else if (ndata == 2) {
    root[inode].cen[0] = (root[inode].x[0][0] + root[inode].x[0][1]) * 0.5;
    root[inode].cen[1] = (root[inode].x[1][0] + root[inode].x[1][1]) * 0.5;
    root[inode].cen[2] = (root[inode].x[2][0] + root[inode].x[2][1]) * 0.5;
    const real dx = root[inode].cen[0] - root[inode].x[0][0];
    const real dy = root[inode].cen[1] - root[inode].x[1][0];
    const real dz = root[inode].cen[2] - root[inode].x[2][0];
    root[inode].r = sqrt(dx * dx + dy * dy + dz * dz);
  }
  else {
    root[inode].cen[0] = root[inode].x[0][0];
    root[inode].cen[1] = root[inode].x[1][0];
    root[inode].cen[2] = root[inode].x[2][0];
    root[inode].r = 0;
  }

  if (ndata <= nleaf) return false;
  size_t n = ndata >> 1;        /* index of the median point */
  /* Split the data by the median in the direction with largest variance. */
  qselect(root[inode].x, 0, n, ndata, p);
  return true;
}
#else
/******************************************************************************
Function `balltree_init`:
  Initialise a node of the ball tree.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Address of the ball tree node.
******************************************************************************/
static BLT *balltree_init(real **x, const size_t idx, const size_t ndata,
    const size_t nleaf) {
  BLT *node = malloc(sizeof(BLT));
  if (!node) return NULL;

  node->n = ndata;
  node->x[0] = x[0] + idx;
  node->x[1] = x[1] + idx;
  node->x[2] = x[2] + idx;
  node->left = node->right = NULL;

  /* Find the directions with largest variances using PCA. */
  double p[9] = {0,0,0,0,0,0,0,0,0};
  if (ndata < 3)  p[0] = p[4] = p[8] = 1;
  else pca_vector(node->x, ndata, BENCHMARK_PCA_EPSILON, p);

  if (ndata > 4) {
    /* Find the extreme points along the principal direction. */
    size_t iex[2];
    extreme_points(node->x, 0, ndata, p, iex);

    if (iex[0] == iex[1]) {     /* there is only a single point. */
      node->cen[0] = node->x[0][0];
      node->cen[1] = node->x[1][0];
      node->cen[2] = node->x[2][0];
      node->r = 0;
    }
    else {
      /* Move extreme points to the beginning of the dataset. */
      SWAP_DATA(node->x, 0, iex[0]);
      SWAP_DATA(node->x, 1, iex[1]);

      /* Find extreme points along the second principal direction. */
      extreme_points(node->x, 2, ndata, p + 3, iex);

      if (iex[0] == iex[1]) {     /* points are along a line */
        node->cen[0] = (node->x[0][0] + node->x[0][1]) * 0.5;
        node->cen[1] = (node->x[1][0] + node->x[1][1]) * 0.5;
        node->cen[2] = (node->x[2][0] + node->x[2][1]) * 0.5;
        const real dx = node->cen[0] - node->x[0][0];
        const real dy = node->cen[1] - node->x[1][0];
        const real dz = node->cen[2] - node->x[2][0];
        node->r = sqrt(dx * dx + dy * dy + dz * dz);
      }
      else {
        /* Move extreme points to the beginning of the dataset. */
        SWAP_DATA(node->x, 2, iex[0]);
        SWAP_DATA(node->x, 3, iex[1]);

        /* Initialise the bounding sphere with the four extreme points. */
        double r2 = 0;
#ifdef SINGLE_PREC
        double cen[3];
        minimum_sphere(node->x, 0, 1, 2, 3, cen, &r2);
        node->cen[0] = cen[0];
        node->cen[1] = cen[1];
        node->cen[2] = cen[2];
#else
        minimum_sphere(node->x, 0, 1, 2, 3, node->cen, &r2);
#endif
        node->r = sqrt(r2);

        /* Check the rest of the points and enlarge the sphere if necessary. */
        for (size_t i = 4; i < ndata; i++) {
          const real dx = node->cen[0] - node->x[0][i];
          const real dy = node->cen[1] - node->x[1][i];
          const real dz = node->cen[2] - node->x[2][i];
          real dist = dx * dx + dy * dy + dz * dz;
          if (dist > r2) {
            dist = sqrt(dist);
            node->r = (node->r + dist) * 0.5;
            const real fac = node->r / dist;
            node->cen[0] = node->x[0][i] + dx * fac;
            node->cen[1] = node->x[1][i] + dy * fac;
            node->cen[2] = node->x[2][i] + dz * fac;
            r2 = node->r * node->r;
          }
        }
      }
    }
  }
  else if (ndata == 4) {
    double r2 = 0;
#ifdef SINGLE_PREC
    double cen[3];
    minimum_sphere(node->x, 0, 1, 2, 3, cen, &r2);
    node->cen[0] = cen[0];
    node->cen[1] = cen[1];
    node->cen[2] = cen[2];
#else
    minimum_sphere(node->x, 0, 1, 2, 3, node->cen, &r2);
#endif
    node->r = sqrt(r2);
  }
  else if (ndata == 3) {
    double r2 = 0;
#ifdef SINGLE_PREC
    double cen[3];
    minimum_sphere(node->x, 0, 1, 2, 3, cen, &r2);
    node->cen[0] = cen[0];
    node->cen[1] = cen[1];
    node->cen[2] = cen[2];
#else
    minimum_circle(node->x, 0, 1, 2, node->cen, &r2);
#endif
    node->r = sqrt(r2);
  }
  else if (ndata == 2) {
    node->cen[0] = (node->x[0][0] + node->x[0][1]) * 0.5;
    node->cen[1] = (node->x[1][0] + node->x[1][1]) * 0.5;
    node->cen[2] = (node->x[2][0] + node->x[2][1]) * 0.5;
    const real dx = node->cen[0] - node->x[0][0];
    const real dy = node->cen[1] - node->x[1][0];
    const real dz = node->cen[2] - node->x[2][0];
    node->r = sqrt(dx * dx + dy * dy + dz * dz);
  }
  else {
    node->cen[0] = node->x[0][0];
    node->cen[1] = node->x[1][0];
    node->cen[2] = node->x[2][0];
    node->r = 0;
  }

  if (ndata <= nleaf) return node;
  size_t n = ndata >> 1;        /* index of the median point */
  /* Split the data by the median in the direction with largest variance. */
  qselect(node->x, 0, n, ndata, p);
  return node;
}
#endif

#ifdef BENCHMARK_TREE_AS_ARRAY
/******************************************************************************
Function `balltree_build`:
  Construct the ball tree recursively given a data set.
Arguments:
  * `root`:     root of the tree;
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
******************************************************************************/
static BLT* balltree_build(BLT *root, real **x, const size_t idx,
    const size_t ndata, const size_t nleaf) {
  if (!balltree_init(root, 0, x, idx, ndata, nleaf)) return root;

  for (size_t level = 1; ; level++) {
    size_t nnode = 1 << level;          /* number of nodes at this level */
    size_t ntot = (nnode << 1) - 1;     /* total number of nodes so far */

    /* Construct children nodes. */
    bool next = false;
    for (size_t i = nnode - 1; i < ntot; ) {
      size_t pid = (i - 1) >> 1;        /* index of parent node */
      size_t nprev = root[pid].n;
      if (nprev <= nleaf) {             /* the parent is already a leaf */
        root[i++].n = 0;
        root[i++].n = 0;
      }
      else {                            /* construct children nodes */
        size_t n = nprev >> 1;
        next = balltree_init(root, i++, root[pid].x, 0, n, nleaf) || next;
        next = balltree_init(root, i++, root[pid].x, n, nprev - n, nleaf)
            || next;
      }
    }
    if (!next) return root;
  }
}
#elif defined(BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `balltree_build`:
  Construct the ball tree recursively given a data set.
Arguments:
  * `root`:     root of the tree;
  * `inode`:    current index of the node;
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Index of the current node.
******************************************************************************/
static size_t balltree_build(BLT *root, size_t *inode, real **x,
    const size_t idx, const size_t ndata, const size_t nleaf) {
  const size_t this_node = *inode;
  if (balltree_init(root, (*inode)++, x, idx, ndata, nleaf)) {
    BLT *node = root + this_node;
    size_t n = ndata >> 1;
    node->left = root + balltree_build(root, inode, node->x, 0, n, nleaf);
    node->right = root +
        balltree_build(root, inode, node->x, n, ndata - n, nleaf);
  }
  return this_node;
}
#else
/******************************************************************************
Function `balltree_build`:
  Construct the ball tree recursively given a data set.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes;
  * `err`:      error indicator.
Return:
  Root of the constructed tree.
******************************************************************************/
static BLT* balltree_build(real **x, const size_t idx, const size_t ndata,
    const size_t nleaf, int *err) {
  if (*err) return NULL;

  BLT *node = balltree_init(x, idx, ndata, nleaf);
  if (!node) {
    *err = EXIT_FAILURE;
    return NULL;
  }

  if (ndata > nleaf) {
    size_t n = ndata >> 1;      /* index of the median point */
    node->left = balltree_build(node->x, 0, n, nleaf, err);
    node->right = balltree_build(node->x, n, ndata - n, nleaf, err);
  }
  return node;
}
#endif

#if !defined(BENCHMARK_TREE_AS_ARRAY) && !defined(BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `balltree_free`:
  Release the memory allocated for the ball tree.
Arguments:
  * `node`:     pointer to a root of the ball tree.
******************************************************************************/
static void balltree_free(BLT *node) {
  if(!node) return;
  balltree_free(node->left);
  balltree_free(node->right);
  free(node);
}
#endif

#ifdef PRINT_TREE
#ifdef BENCHMARK_TREE_AS_ARRAY
  #error `PRINT_TREE` not implemented for `BENCHMARK_TREE_AS_ARRAY`
#endif
static void balltree_print(BLT *node, const int lv) {
  if (!node) return;
  printf("NODE level: %d", lv);
  if (node->left == NULL) printf(" (leaf)");
  printf("\n  ndata: %zu\n  cen: %g %g %g\n  rad: %g\n", node->n,
      node->cen[0], node->cen[1], node->cen[2], node->r);
  if (node->left == NULL) {
    printf("  data: {");
    for (size_t i = 0; i < node->n; i++) {
      printf("{%g,%g,%g},", node->x[0][i], node->x[1][i], node->x[2][i]);
    }
    printf("\b}\n");
  }
  balltree_print(node->left, lv + 1);
  balltree_print(node->right, lv + 1);
}
#endif

/******************************************************************************
Function `create_balltree`:
  Construct the ball tree for a input dataset.
Arguments:
  * `data`:     the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Root of the constructed tree.
******************************************************************************/
static BLT *create_balltree(DATA *data, const size_t nleaf) {
  if (!data || !data->n) {
    P_ERR("the input dataset is not initialised\n");
    return NULL;
  }

  real *x[3];
  x[0] = data->x;
  x[1] = data->y;
  x[2] = data->z;

#if defined(BENCHMARK_TREE_AS_ARRAY) || defined(BENCHMARK_TREE_PREALLOC)
  const size_t size = balltree_nnode(data->n, nleaf);
  if (!size) {
    P_ERR("failed to precompute the size of the tree\n");
    return NULL;
  }
  BLT *root = calloc(size, sizeof *root);
  if (!root) {
    P_ERR("failed to allocate memory for ball tree\n");
    return NULL;
  }
#endif

#ifdef BENCHMARK_TREE_AS_ARRAY
  balltree_build(root, x, 0, data->n, nleaf);
#elif defined(BENCHMARK_TREE_PREALLOC)
  #if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  /* Allocate memory for bounding box properties. */
  real *vol = NULL;
  if (posix_memalign((void **) &vol, BENCHMARK_MEMALIGN_BYTE,
      size * 3 * sizeof(real))) {
    P_ERR("failed to allocate aligned memory for the volume of nodes\n");
    free(root);
    return NULL;
  }
  for (size_t i = 0; i < size; i++) root[i].cen = vol + i * 3;
  #endif
  size_t inode = 0;
  balltree_build(root, &inode, x, 0, data->n, nleaf);
  /* Reduce memory cost if applicable. */
  if (size > inode) {
    BLT *ini = root;
    BLT *tmp = realloc(root, inode * sizeof(BLT));
    if (tmp) {
      for (size_t i = 0; i < inode; i++) {
        if (tmp[i].left) tmp[i].left = tmp[i].left - ini + tmp;
        if (tmp[i].right) tmp[i].right = tmp[i].right - ini + tmp;
      }
      root = tmp;
    }
  }
#else
  int err = 0;
  BLT *root = balltree_build(x, 0, data->n, nleaf, &err);
  if (err) {
    P_ERR("failed to construct the kd-tree\n");
    balltree_free(root);
    return NULL;
  }
#endif

#ifdef PRINT_TREE
  balltree_print(root, 0);
#endif
  return root;
}


/*============================================================================*\
                     Pair counting functions from templates
\*============================================================================*/

/* Clean all the relevant macros first */
#ifdef BENCHMARK_TREE_TYPE
  #undef BENCHMARK_TREE_TYPE
#endif
#ifdef BENCHMARK_PAIRCOUNT_TYPE
  #undef BENCHMARK_PAIRCOUNT_TYPE
#endif
#ifdef BENCHMARK_BIN_SMIN
  #undef BENCHMARK_BIN_SMIN
#endif

#define BENCHMARK_TREE_TYPE     BENCHMARK_STRUCT_BALLTREE
/*
static inline DUAL_NODE *stack_pop(STACK_DUAL_NODE *s) {
  if (!s->size) return NULL;
  s->size -= 1;
  return s->nodes + s->size;
}
*/

/* paircount_box_smin0 */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_BOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_ZERO
#include "metric_common.c"
#include "metric_balltree.c"
#include "dual_tree.c"

/* paircount_box */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_BOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_NONZERO
#include "metric_common.c"
#include "metric_balltree.c"
#include "dual_tree.c"

/* paircount_smin0 */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_NOBOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_ZERO
#include "metric_common.c"
#include "metric_balltree.c"
#include "dual_tree.c"

/* paircount */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_NOBOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_NONZERO
#include "metric_common.c"
#include "metric_balltree.c"
#include "dual_tree.c"


/*============================================================================*\
                          Interface for pair counting
\*============================================================================*/

/******************************************************************************
Function `paircount_balltree`:
  Count pairs using grids.
Arguments:
  * `data`:     structure for the data catalogue;
  * `rmin`:     minimum separation of interest;
  * `rmax`:     maximum separation of interest;
  * `csize`:    number of grids per side, or number of partilces in leaf nodes;
  * `nnode`:    number of distance evaluations between nodes;
  * `ndist`:    number of distance evaluations performed in total;
  * `npair`:    number of pairs in the separation range of interest.
Return:
  Zero on success; EXIT_FAILURE on error.
******************************************************************************/
int paircount_balltree(DATA *data, const double rmin, const double rmax,
    const int csize,
#ifndef BENCHMARK_TIMING
    size_t *nnode, size_t *ndist,
#endif
    size_t *npair) {
  printf("Constructing the ball tree ...");
  fflush(stdout);

  BLT *tree = create_balltree(data, csize);
  if (!tree) {
    P_ERR("failed to construct ball tree\n");
    return EXIT_FAILURE;
  }

  /* Initialise the stack for dual nodes. */
  STACK_DUAL_NODE stack;
  stack.size = stack.capacity = 0;
#ifdef BENCHMARK_TREE_AS_ARRAY
  stack.root = tree;
  stack.nleaf = csize;
  stack_push(&stack, 0, 0);
#else
  stack_push(&stack, tree, tree);
#endif
  const double r2min = rmin * rmin;
  const double r2max = rmax * rmax;
  *npair = 0;

  /* Initialise configurations. */
  void *conf = NULL;
  if (data->isbox) {
    if (rmin == 0) conf = calloc(1, sizeof(conf_balltree_box_smin0));
    else conf = calloc(1, sizeof(conf_balltree_box));
  }
  else {
    if (rmin == 0) conf = calloc(1, sizeof(conf_balltree_smin0));
    else conf = calloc(1, sizeof(conf_balltree));
  }
  if (!conf) {
    P_ERR("failed to allocate memory for pair counting\n");
#if defined(BENCHMARK_TREE_AS_ARRAY) || defined(BENCHMARK_TREE_PREALLOC)
  #if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
    free(tree->cen);
  #endif
    free(tree);
#else
    balltree_free(tree);
#endif
    return EXIT_FAILURE;
  }

  /* Choose the right pair counting function, and set configurations. */
  void (*count_func) (STACK_DUAL_NODE *, void *) = NULL;
  if (data->isbox) {
    if (rmin == 0) {
      count_func = pair_balltree_box_smin0;
      conf_balltree_box_smin0 *config = (conf_balltree_box_smin0 *) conf;
      config->bsize = data->bsize;
      config->hsize = data->bsize * 0.5;
      config->rmax = rmax;
      config->r2max = r2max;
    }
    else {
      count_func = pair_balltree_box;
      conf_balltree_box *config = (conf_balltree_box *) conf;
      config->bsize = data->bsize;
      config->hsize = data->bsize * 0.5;
      config->rmax = rmax;
      config->r2max = r2max;
      config->rmin = rmin;
      config->r2min = r2min;
    }
  }
  else {
    if (rmin == 0) {
      count_func = pair_balltree_smin0;
      conf_balltree_smin0 *config = (conf_balltree_smin0 *) conf;
      config->rmax = rmax;
      config->r2max = r2max;
    }
    else {
      count_func = pair_balltree;
      conf_balltree *config = (conf_balltree *) conf;
      config->rmax = rmax;
      config->r2max = r2max;
      config->rmin = rmin;
      config->r2min = r2min;
    }
  }

  printf(FMT_DONE "Start pair counting ...");
  fflush(stdout);

#ifdef BENCHMARK_TIMING
  clock_t start = clock(), diff;
#endif
  while (stack.size) count_func(&stack, conf);

#ifdef BENCHMARK_TIMING
  diff = clock() - start;
  double sec = (double) diff / CLOCKS_PER_SEC;
#endif

  /* Retrieve pair counts. */
  if (data->isbox) {
    if (rmin == 0) {
      conf_balltree_box_smin0 *config = (conf_balltree_box_smin0 *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
    else {
      conf_balltree_box *config = (conf_balltree_box *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
  }
  else {
    if (rmin == 0) {
      conf_balltree_smin0 *config = (conf_balltree_smin0 *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
    else {
      conf_balltree *config = (conf_balltree *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
  }

  stack_destroy(&stack);
  free(conf);
#if defined(BENCHMARK_TREE_AS_ARRAY) || defined(BENCHMARK_TREE_PREALLOC)
  #if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  free(tree->cen);
#endif
  free(tree);
#else
  balltree_free(tree);
#endif
  printf(FMT_DONE);

#ifdef BENCHMARK_TIMING
  printf("Time for pair counting: " OFMT_DBL " seconds\n", sec);
#endif
  return 0;
}
