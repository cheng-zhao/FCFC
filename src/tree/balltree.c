/*******************************************************************************
* balltree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "balltree.h"
#include "pca.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
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
static inline void extreme_points(real *const x[static 3], const size_t istart,
    const size_t iend, const double *p, size_t iex[2]) {
  real min = REAL_MAX, max = -REAL_MAX;
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
Function `swap_points`:
  Swap two data points on a tree node given their indices.
Arguments:
  * `node`:     pointer to the tree node;
  * `i`:        index of the first point;
  * `j`:        index of the second point.
******************************************************************************/
static inline void swap_points(BLT *node, const size_t i, const size_t j) {
  real tmp = node->x[0][i]; node->x[0][i] = node->x[0][j]; node->x[0][j] = tmp;
  tmp = node->x[1][i]; node->x[1][i] = node->x[1][j]; node->x[1][j] = tmp;
  tmp = node->x[2][i]; node->x[2][i] = node->x[2][j]; node->x[2][j] = tmp;
#if FCFC_XDIM == 4
  tmp = node->x[3][i]; node->x[3][i] = node->x[3][j]; node->x[3][j] = tmp;
#endif
  if (node->w) {
    tmp = node->w[i]; node->w[i] = node->w[j]; node->w[j] = tmp;
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
static int dist_plane_compare(real *const x[static 3],
    const size_t i, const size_t j, const double *n) {
  real di = x[0][i] * n[0] + x[1][i] * n[1] + x[2][i] * n[2];
  real dj = x[0][j] * n[0] + x[1][j] * n[1] + x[2][j] * n[2];
  if (di < dj) return -1;
  if (di > dj) return 1;
  return 0;
}


/* Import the function for partition data by the median. */

#ifdef QSELECT_COMPARE
  #undef QSELECT_COMPARE
#endif
#ifdef QSELECT_DTYPE
  #undef QSELECT_DTYPE
#endif
#ifdef QSELECT_SWAP
  #undef QSELECT_SWAP
#endif
#ifdef QSELECT_NAME
  #undef QSELECT_NAME
#endif

/* `qselect` for coordinates only. */
typedef struct {
  real *x[FCFC_XDIM];
} DATA_x;

#define QSELECT_COMPARE(d,i,j,n)        dist_plane_compare((d)->x,i,j,n)
#define QSELECT_DTYPE                   DATA_x
#define QSELECT_NAME                    x
#if     FCFC_XDIM == 3
  #define QSELECT_SWAP(d,i,j) {                                              \
    real _t = (d)->x[0][i]; (d)->x[0][i] = (d)->x[0][j]; (d)->x[0][j] = _t;  \
    _t = (d)->x[1][i]; (d)->x[1][i] = (d)->x[1][j]; (d)->x[1][j] = _t;       \
    _t = (d)->x[2][i]; (d)->x[2][i] = (d)->x[2][j]; (d)->x[2][j] = _t;       \
  }
#elif   FCFC_XDIM == 4
  #define QSELECT_SWAP(d,i,j) {                                              \
    real _t = (d)->x[0][i]; (d)->x[0][i] = (d)->x[0][j]; (d)->x[0][j] = _t;  \
    _t = (d)->x[1][i]; (d)->x[1][i] = (d)->x[1][j]; (d)->x[1][j] = _t;       \
    _t = (d)->x[2][i]; (d)->x[2][i] = (d)->x[2][j]; (d)->x[2][j] = _t;       \
    _t = (d)->x[3][i]; (d)->x[3][i] = (d)->x[3][j]; (d)->x[3][j] = _t;       \
  }
#else
  #error unexpected definition of `FCFC_XDIM`
#endif

#include "qselect.c"

/* `qselect` for both coordinates and weights. */
typedef struct {
  real *x[FCFC_XDIM];
  real *w;
} DATA_xw;

#define QSELECT_COMPARE(d,i,j,n)        dist_plane_compare((d)->x,i,j,n)
#define QSELECT_DTYPE                   DATA_xw
#define QSELECT_NAME                    xw
#if     FCFC_XDIM == 3
  #define QSELECT_SWAP(d,i,j) {                                              \
    real _t = (d)->x[0][i]; (d)->x[0][i] = (d)->x[0][j]; (d)->x[0][j] = _t;  \
    _t = (d)->x[1][i]; (d)->x[1][i] = (d)->x[1][j]; (d)->x[1][j] = _t;       \
    _t = (d)->x[2][i]; (d)->x[2][i] = (d)->x[2][j]; (d)->x[2][j] = _t;       \
    _t = (d)->w[i]; (d)->w[i] = (d)->w[j]; (d)->w[j] = _t;                   \
  }
#elif   FCFC_XDIM == 4
  #define QSELECT_SWAP(d,i,j) {                                              \
    real _t = (d)->x[0][i]; (d)->x[0][i] = (d)->x[0][j]; (d)->x[0][j] = _t;  \
    _t = (d)->x[1][i]; (d)->x[1][i] = (d)->x[1][j]; (d)->x[1][j] = _t;       \
    _t = (d)->x[2][i]; (d)->x[2][i] = (d)->x[2][j]; (d)->x[2][j] = _t;       \
    _t = (d)->x[3][i]; (d)->x[3][i] = (d)->x[3][j]; (d)->x[3][j] = _t;       \
    _t = (d)->w[i]; (d)->w[i] = (d)->w[j]; (d)->w[j] = _t;                   \
  }
#else
  #error unexpected definition of `FCFC_XDIM`
#endif

#include "qselect.c"


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
static inline void minimum_circle(real *const x[static 3], const size_t a,
    const size_t b, const size_t c, double *cen, double *r2) {
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
  if (d > -DOUBLE_TOL && d < DOUBLE_TOL) {
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
static void minimum_sphere(real *const x[static 3], const size_t a,
    const size_t b, const size_t c, const size_t d, double *cen, double *r2) {
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
  if (v > -DOUBLE_TOL && v < DOUBLE_TOL) {
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
        if (dx * dx + dy * dy + dz * dz < *r2 + DOUBLE_TOL) return;
      }
      if (t <= 0) {
        minimum_circle(x, a, b, d, cen, r2);
        const double dx = cen[0] - x[0][c];
        const double dy = cen[1] - x[1][c];
        const double dz = cen[2] - x[2][c];
        if (dx * dx + dy * dy + dz * dz < *r2 + DOUBLE_TOL) return;
      }
      if (u <= 0) {
        minimum_circle(x, a, b, c, cen, r2);
        const double dx = cen[0] - x[0][d];
        const double dy = cen[1] - x[1][d];
        const double dz = cen[2] - x[2][d];
        if (dx * dx + dy * dy + dz * dz < *r2 + DOUBLE_TOL) return;
      }
      minimum_circle(x, b, c, d, cen, r2);
    }
  }
}

/*============================================================================*\
                        Functions for tree construction
\*============================================================================*/

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
  if (size >= FCFC_TREE_MAX_SIZE) return 0;
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

/******************************************************************************
Function `balltree_realloc`:
  Reallocate memory for a balltree.
Arguments:
  * `root`:     root of the ball tree;
  * `size`:     the original number of tree nodes;
  * `nnode`:    the expected number of tree nodes after memory reallocation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int balltree_realloc(BLT **root, const size_t size, const size_t nnode) {
  if (SIZE_MAX / 3 < nnode) return FCFC_ERR_ARG;

  /* Record addresses of the original tree. */
  BLT *ini = *root;
  real *vini = ini->cen;
  BLT *tmp = realloc(*root, sizeof(BLT) * nnode);
  if (!tmp)  return FCFC_ERR_MEMORY;

  const size_t num = (size < nnode) ? size : nnode;
  /* Compute the new addresses. */
#ifdef OMP
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < num; i++) {
    if (tmp[i].left) tmp[i].left = tmp[i].left - ini + tmp;
    if (tmp[i].right) tmp[i].right = tmp[i].right - ini + tmp;
  }
  *root = tmp;

  real *vtmp = realloc(vini, sizeof(real) * nnode * 3);
  if (!vtmp)  return FCFC_ERR_MEMORY;

#ifdef OMP
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < nnode; i++) tmp[i].cen = vtmp + i * 3;

  return 0;
}

/******************************************************************************
Function `balltree_init`:
  Initialise a node of the ball tree.
Arguments:
  * `root`:     root of the tree;
  * `inode`:    index of the current node;
  * `x`:        coordinates of the input dataset;
  * `w`:        weights of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `bsize`:    side lengths of the periodic box;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  True of the node has children; false otherwise.
******************************************************************************/
static bool balltree_init(BLT *root, const size_t inode,
    real *x[static FCFC_XDIM], real *w, const size_t idx, const size_t ndata,
#ifdef FCFC_METRIC_PERIODIC
    const real bsize[static 3],
#endif
    const size_t nleaf) {
  root[inode].n = ndata;
  for (int k = 0; k < FCFC_XDIM; k++) root[inode].x[k] = x[k] + idx;
  root[inode].w = (w) ? w + idx : NULL;
  root[inode].left = root[inode].right = NULL;

  /* Find the directions with largest variances using PCA. */
  double p[9] = {0,0,0,0,0,0,0,0,0};
  if (ndata < 3)  p[0] = p[4] = p[8] = 1;
  else pca_vector(root[inode].x, ndata, FCFC_PCA_EPSILON, p);

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
      swap_points(root + inode, 0, iex[0]);
      swap_points(root + inode, 1, iex[1]);

      /* Find extreme points along the second principal direction. */
      extreme_points(root[inode].x, 2, ndata, p + 3, iex);

      if (iex[0] == iex[1]) {     /* points are along a line */
        root[inode].cen[0] = (root[inode].x[0][0] + root[inode].x[0][1]) * 0.5;
        root[inode].cen[1] = (root[inode].x[1][0] + root[inode].x[1][1]) * 0.5;
        root[inode].cen[2] = (root[inode].x[2][0] + root[inode].x[2][1]) * 0.5;
        const real dx = root[inode].cen[0] - root[inode].x[0][0];
        const real dy = root[inode].cen[1] - root[inode].x[1][0];
        const real dz = root[inode].cen[2] - root[inode].x[2][0];
        root[inode].r = sqrt(dx * dx + dy * dy + dz * dz) + REAL_EPS;
      }
      else {
        /* Move extreme points to the beginning of the dataset. */
        swap_points(root + inode, 2, iex[0]);
        swap_points(root + inode, 3, iex[1]);

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
        root[inode].r = sqrt(r2) + REAL_EPS;

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
    root[inode].r = sqrt(r2) + REAL_EPS;
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
    root[inode].r = sqrt(r2) + REAL_EPS;
  }
  else if (ndata == 2) {
    root[inode].cen[0] = (root[inode].x[0][0] + root[inode].x[0][1]) * 0.5;
    root[inode].cen[1] = (root[inode].x[1][0] + root[inode].x[1][1]) * 0.5;
    root[inode].cen[2] = (root[inode].x[2][0] + root[inode].x[2][1]) * 0.5;
    const real dx = root[inode].cen[0] - root[inode].x[0][0];
    const real dy = root[inode].cen[1] - root[inode].x[1][0];
    const real dz = root[inode].cen[2] - root[inode].x[2][0];
    root[inode].r = sqrt(dx * dx + dy * dy + dz * dz) + REAL_EPS;
  }
  else {
    root[inode].cen[0] = root[inode].x[0][0];
    root[inode].cen[1] = root[inode].x[1][0];
    root[inode].cen[2] = root[inode].x[2][0];
    root[inode].r = 0;
  }

  if (ndata <= nleaf
#ifdef FCFC_METRIC_PERIODIC
      && root[inode].r < bsize[0] / 4
      && root[inode].r < bsize[1] / 4
      && root[inode].r < bsize[2] / 4
#endif
      ) return false;
  size_t n = ndata >> 1;        /* index of the median point */
  /* Split the data by the median in the direction with largest variance. */
  if (w) {
    DATA_xw data;
    for (int k = 0; k < FCFC_XDIM; k++) data.x[k] = root[inode].x[k];
    data.w = root[inode].w;
    qselect_xw(&data, 0, n, ndata, p);
  }
  else {
    DATA_x data;
    for (int k = 0; k < FCFC_XDIM; k++) data.x[k] = root[inode].x[k];
    qselect_x(&data, 0, n, ndata, p);
  }
  return true;
}

/******************************************************************************
Function `balltree_build`:
  Construct the ball tree recursively given a data set.
Arguments:
  * `root`:     root of the tree;
  * `inode`:    current index of the node;
  * `x`:        coordinates of the input dataset;
  * `w`:        weights of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `bsize`:    side lengths of the periodic box;
  * `size`:     number of nodes allocated for the tree;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Index of the current node.
******************************************************************************/
static size_t balltree_build(
#ifdef FCFC_METRIC_PERIODIC
    BLT **root,
#else
    BLT *root,
#endif
    size_t *inode, real *x[static FCFC_XDIM], real *w,
    const size_t idx, const size_t ndata,
#ifdef FCFC_METRIC_PERIODIC
    const real bsize[static 3], size_t *size,
#endif
    const size_t nleaf) {
  const size_t this_node = (*inode)++;
#ifdef FCFC_METRIC_PERIODIC
  if (this_node >= *size) {
    if (SIZE_MAX / 2 < *size || balltree_realloc(root, *size, *size << 1)) {
      P_EXT("failed to reallocate memory for the tree\n");
      exit(FCFC_ERR_MEMORY);
    }
    *size <<= 1;
  }
  if (balltree_init(*root, this_node, x, w, idx, ndata, bsize, nleaf)) {
    size_t n = ndata >> 1;
    BLT *node = *root + this_node;
    real *this_x[FCFC_XDIM];
    for (int k = 0; k < FCFC_XDIM; k++) this_x[k] = node->x[k];
    size_t i = balltree_build(root, inode, this_x, node->w,
        0, n, bsize, size, nleaf);
    (*root)[this_node].left = *root + i;
    i = balltree_build(root, inode, this_x, (*root)[this_node].w,
        n, ndata - n, bsize, size, nleaf);
    (*root)[this_node].right = *root + i;
  }
#else
  if (balltree_init(root, this_node, x, w, idx, ndata, nleaf)) {
    BLT *node = root + this_node;
    size_t n = ndata >> 1;
    node->left = root +
        balltree_build(root, inode, node->x, node->w, 0, n, nleaf);
    node->right = root +
        balltree_build(root, inode, node->x, node->w, n, ndata - n, nleaf);
  }
#endif
  return this_node;
}


/*============================================================================*\
                      Functions for requesting tree nodes
\*============================================================================*/

#if defined(MPI) || defined(OMP)
/******************************************************************************
Function `balltree_level`:
  Compute the minimum tree level with at least the specified number of nodes.
Arguments:
  * `nnode`:    number of nodes at the same level.
Return:
  The minimum tree level with at least `nnode` nodes.
******************************************************************************/
static inline int balltree_level(uint32_t nnode) {
  /* Round up the number of nodes to the next power of 2:
  http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  nnode--;
  nnode |= nnode >> 1;
  nnode |= nnode >> 2;
  nnode |= nnode >> 4;
  nnode |= nnode >> 8;
  nnode |= nnode >> 16;
  nnode++;

  /* Compute the log base 2:
  http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightMultLookup */
  static const int MultiplyDeBruijnBitPosition[32] = {
      0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
      31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};

  return MultiplyDeBruijnBitPosition[((uint32_t) (nnode * 0x077CB531U)) >> 27];
}

/******************************************************************************
Function `gather_nodes`:
  Gather all tree nodes at the same level.
Arguments:
  * `tree`:     pointer to the root of the tree;
  * `current`:  current level of the tree;
  * `level`:    level for the request;
  * `nodes`:    array for storing the resulting nodes;
  * `inode`:    current index of the resulting nodes.
******************************************************************************/
static void gather_nodes(const BLT *tree, const int current,
    const int level, const BLT **nodes, size_t *inode) {
  if (!tree) return;
  if (current == level) nodes[(*inode)++] = tree;
  else {
    gather_nodes(tree->left, current + 1, level, nodes, inode);
    gather_nodes(tree->right, current + 1, level, nodes, inode);
  }
}
#endif


/*============================================================================*\
                            Interfaces for ball tree
\*============================================================================*/

/******************************************************************************
Function `create_balltree`:
  Construct the ball tree for a input dataset.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `w`:        weights of the input dataset;
  * `ndata`:    number of input data points;
  * `bsize`:    side lengths of the periodic box;
  * `nleaf`:    maximum number of data points on leaf nodes;
  * `nnode`:    number of tree nodes created successfully.
Return:
  Address of the root node on success; NULL on error.
******************************************************************************/
BLT *create_balltree(real *x[static FCFC_XDIM], real *w, const size_t ndata,
#ifdef FCFC_METRIC_PERIODIC
    const real bsize[static 3],
#endif
    const size_t nleaf, size_t *nnode) {
  if (!x || !ndata) {
    P_ERR("the input dataset is not initialised\n");
    return NULL;
  }
  if (!nleaf) {
    P_ERR("the capacity of leaf nodes must be non-zero\n");
    return NULL;
  }

  size_t size = balltree_nnode(ndata, nleaf);
  if (!size) {
    P_ERR("failed to precompute the size of the tree\n");
    return NULL;
  }
  BLT *root = calloc(size, sizeof *root);
  if (!root) {
    P_ERR("failed to allocate memory for ball tree\n");
    return NULL;
  }

  /* Allocate memory for bounding box properties. */
  real *cen = malloc(sizeof(real) * size * 3);
  if (!cen) {
    P_ERR("failed to allocate memory for the bounding volume of nodes\n");
    free(root);
    return NULL;
  }
#ifdef OMP
  #pragma omp parallel for default(none) shared(size,root,cen)
#endif
  for (size_t i = 0; i < size; i++) root[i].cen = cen + i * 3;

  /* Build the tree. */
  size_t inode = 0;
#ifdef FCFC_METRIC_PERIODIC
  balltree_build(&root, &inode, x, w, 0, ndata, bsize, &size, nleaf);
#else
  balltree_build(root, &inode, x, w, 0, ndata, nleaf);
#endif
  if (inode == 0) {
    P_ERR("failed to create the ball tree\n");
    free(root); free(cen);
    return NULL;
  }

  /* Reduce memory cost if applicable. */
  if (size > inode && balltree_realloc(&root, size, inode)) {
    P_ERR("failed to reallocate memory for the tree\n");
    free(root); free(cen); return NULL;
  }

  *nnode = inode;
  return root;
}

#ifdef MPI
/******************************************************************************
Function `balltree_broadcast`:
  Broadcast the ball tree to all MPI tasks.
Arguments:
  * `root`:     pointer to the root of the ball tree;
  * `nnode`:    total number of tree nodes;
  * `wt`:       indicate whether to broadcast weights;
  * `src`:      the source for the broadcast;
  * `rank`:     ID of MPI task.
******************************************************************************/
void balltree_broadcast(BLT **root, size_t *nnode, const bool wt, const int src,
    const int rank) {
  if (!root || !nnode ||
      (src == rank && (!(*root) || !(*nnode) || !((*root)->n)))) {
    P_ERR("the ball tree is not initialized for broadcasting\n");
    FCFC_QUIT(FCFC_ERR_ARG);
  }
  if (src < 0 || rank < 0) {
    P_ERR("invalid MPI ranks for broadcasting the ball tree\n");
    FCFC_QUIT(FCFC_ERR_ARG);
  }

  /* Broadcast the number of tree nodes and data points. */
  size_t ndata = 0;
  if (rank == src) ndata = (*root)->n;
  MPI_Request req[5 + FCFC_XDIM];
  if (MPI_Ibcast(nnode, 1, FCFC_MPI_SIZE_T, src, MPI_COMM_WORLD, req) ||
      MPI_Ibcast(&ndata, 1, FCFC_MPI_SIZE_T, src, MPI_COMM_WORLD, req + 1) ||
      MPI_Waitall(2, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast the ball tree\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  /* Allocate memory for the addresses and numbers of data points. */
  ptrdiff_t *addr;
  size_t *num;
  real *rad;
  if (!(addr = calloc(*nnode * 3, sizeof(ptrdiff_t))) ||
      !(num = calloc(*nnode, sizeof(size_t))) ||
      !(rad = calloc(*nnode, sizeof(real)))) {
    P_ERR("failed to allocate memory for broadcasting the ball tree\n");
    FCFC_QUIT(FCFC_ERR_MEMORY);
  }

  real *cen, *x[FCFC_XDIM], *w;
  cen = w = NULL;
  for (int k = 0; k < FCFC_XDIM; k++) x[k] = NULL;
  if (rank == src) {
    cen = (*root)->cen;
    for (int k = 0; k < FCFC_XDIM; k++) x[k] = (*root)->x[k];
    w = (*root)->w;

    /* Gather the number of data points and compute memory offsets. */
#ifdef OMP
  #pragma omp parallel for default(none) shared(nnode,num,root,rad,addr)
#endif
    for (size_t i = 0; i < *nnode; i++) {
      num[i] = (*root)[i].n;
      rad[i] = (*root)[i].r;
      size_t j = i * 3;
      addr[j++] = (*root)[i].x[0] - (*root)[0].x[0];
      addr[j++] = ((*root)[i].left) ? (*root)[i].left - (*root) : 0;
      addr[j] = ((*root)[i].right) ? (*root)[i].right - (*root) : 0;
    }
  }
  else {        /* rank != src */
    /* Allocate memory. */
    if (!(*root = malloc(*nnode * sizeof(BLT))) ||
        !(cen = malloc(*nnode * 3 * sizeof(real))) ||
        !(x[0] = malloc(ndata * sizeof(real))) ||
        !(x[1] = malloc(ndata * sizeof(real))) ||
        !(x[2] = malloc(ndata * sizeof(real))) ||
#if FCFC_XDIM == 4
        !(x[3] = malloc(ndata * sizeof(real))) ||
#endif
        (wt && !(w = malloc(ndata * sizeof(real))))) {
      P_ERR("failed to allocate memory for the ball tree\n");
      FCFC_QUIT(FCFC_ERR_MEMORY);
    }
  }

  /* Broadcast information of the tree. */
  if (MPI_Ibcast(cen, *nnode * 3, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req) ||
      MPI_Ibcast(rad, *nnode, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req + 1) ||
      MPI_Ibcast(num, *nnode, FCFC_MPI_SIZE_T, src, MPI_COMM_WORLD, req + 2) ||
      MPI_Ibcast(addr, *nnode * 3, MPI_AINT, src, MPI_COMM_WORLD, req + 3) ||
      MPI_Ibcast(x[0], ndata, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req + 4) ||
      MPI_Ibcast(x[1], ndata, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req + 5) ||
      MPI_Ibcast(x[2], ndata, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req + 6) ||
#if FCFC_XDIM == 4
      MPI_Ibcast(x[3], ndata, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req + 7) ||
#endif
      (wt && (MPI_Ibcast(w, ndata, FCFC_MPI_REAL, src, MPI_COMM_WORLD,
          req + 4 + FCFC_XDIM)
          || MPI_Waitall(5 + FCFC_XDIM, req, MPI_STATUSES_IGNORE))) ||
      (!wt && MPI_Waitall(4 + FCFC_XDIM, req, MPI_STATUSES_IGNORE))) {
    P_ERR("failed to broadcast the ball tree\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  /* Restore the tree structure. */
  if (rank != src) {
#ifdef OMP
  #pragma omp parallel for
#endif
    for (size_t i = 0; i < *nnode; i++) {
      const size_t j = i * 3;
      (*root)[i].cen = cen + j;
      (*root)[i].r = rad[i];
      (*root)[i].n = num[i];
      for (int k = 0; k < FCFC_XDIM; k++) (*root)[i].x[k] = x[k] + addr[j];
      (*root)[i].w = (wt) ? w + addr[j] : NULL;
      (*root)[i].left = (addr[j + 1]) ? (*root) + addr[j + 1] : NULL;
      (*root)[i].right = (addr[j + 2]) ? (*root) + addr[j + 2] : NULL;
    }
  }

  free(addr);
  free(num);
  free(rad);
}
#endif

#if defined(MPI) || defined(OMP)
/******************************************************************************
Function `balltree_get_nodes`:
  Get all tree nodes at the level with at least the specific number of nodes.
Arguments:
  * `root`:     pointer to the root of the ball tree;
  * `nmin`:     the minimum number of nodes to be requested;
  * `nnode`:    the number of requested nodes.
Return:
  Address of the node list on success; NULL on error.
******************************************************************************/
BLT **balltree_get_nodes(const BLT *root, const uint32_t nmin, size_t *nnode) {
  /* Compute the tree level with at least `nmin` nodes. */
  const int level = balltree_level(nmin);
  /* The number of tree nodes at the given level assuming a balanced tree. */
  size_t nmax = 1U << level;

  /* Request tree nodes. */
  BLT **nodes = malloc(nmax * sizeof(BLT *));
  if (!nodes) {
    P_ERR("failed to request tree nodes\n");
    return NULL;
  }
  *nnode = 0;
  gather_nodes(root, 0, level, (const BLT **) nodes, nnode);

  /* Reduce memory cost if applicable. */
  if (*nnode && *nnode < nmax) {
    BLT **tmp = realloc(nodes, *nnode * sizeof(BLT *));
    if (tmp) nodes = tmp;
  }
  return nodes;
}
#endif

/******************************************************************************
Function `balltree_free`:
  Release the memory allocated for the ball tree.
Arguments:
  * `root`:     pointer to a root of the ball tree.
******************************************************************************/
void balltree_free(BLT *root) {
  if(!root) return;
  free(root->cen);
  free(root);
}
