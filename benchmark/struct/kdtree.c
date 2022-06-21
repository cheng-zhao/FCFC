/*******************************************************************************
* benchmark/struct/kdtree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "create_data.h"
#include "data_struct.h"
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>

/*============================================================================*\
                       Definition of the data structures
\*============================================================================*/

#ifdef BENCHMARK_TREE_AS_ARRAY
typedef struct {
  size_t n;                     /* number of data points   */
  real cen[3];                  /* center of the box       */
  real wid[3];                  /* half-width of the box   */
  real *x[3];                   /* pointers to coordinates */
} KDT;
#else
typedef struct kdtree_struct {
#if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  real *cen;                    /* Store bounding information separately */
  real *wid;                    /* for vectorization. */
#else
  real cen[3];                  /* center of the box       */
  real wid[3];                  /* half-width of the box   */
#endif
  size_t n;                     /* number of data points   */
  real *x[3];                   /* pointers to coordinates */
  struct kdtree_struct *left;   /* left child              */
  struct kdtree_struct *right;  /* right child             */
} KDT;
#endif

/*============================================================================*\
                         Functions for data processing
\*============================================================================*/

/******************************************************************************
Function `kth_compare`:
  Compare the k-th coordinate of two data points in composite space.
  (cf. section 5.5 of de Berg et al. 2008, Computing Geometry 3rd Edition)
Arguments:
  * `x`:        pointer to coordinates of input points;
  * `i`:        index of the first data point;
  * `j`:        index of the second data point;
  * `k`:        the direction of the coordinate to be compared.
Return:
  The order of the two data points.
******************************************************************************/
static int kth_compare(real **x, const size_t i, const size_t j,
    int k) {
  if (x[k][i] < x[k][j]) return -1;
  if (x[k][i] > x[k][j]) return 1;
  k = (k + 1) % 3;
  if (x[k][i] < x[k][j]) return -1;
  if (x[k][i] > x[k][j]) return 1;
  k = (k + 1) % 3;
  if (x[k][i] < x[k][j]) return -1;
  if (x[k][i] > x[k][j]) return 1;
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

#define QSELECT_COMPARE(x,i,j,k)        kth_compare(x,i,j,*((int *)(k)))
#define QSELECT_DTYPE                   real *
#define QSELECT_SWAP(x,i,j) {                                           \
  real _tmp = (x)[0][i]; (x)[0][i] = (x)[0][j]; (x)[0][j] = _tmp;     \
  _tmp = (x)[1][i]; (x)[1][i] = (x)[1][j]; (x)[1][j] = _tmp;            \
  _tmp = (x)[2][i]; (x)[2][i] = (x)[2][j]; (x)[2][j] = _tmp;            \
}

#include "qselect.c"

#undef QSELECT_DTYPE
#undef QSELECT_COMPARE
#undef QSELECT_SWAP


/*============================================================================*\
                        Functions for tree construction
\*============================================================================*/

#if defined(BENCHMARK_TREE_AS_ARRAY) || defined (BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `kdtree_nnode`:
  Compute the upper limit of the number of nodes for a balanced k-d tree.
Arguments:
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Upper limit of the number of k-d tree nodes.
******************************************************************************/
static inline size_t kdtree_nnode(const size_t ndata, const size_t nleaf) {
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
Function `kdtree_init`:
  Initialise a node of the k-D tree.
Arguments:
  * `root`:     root of the tree;
  * `i`:        index of the current node;
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  True if the node has children; false otherwise.
******************************************************************************/
static bool kdtree_init(KDT *root, const size_t i, real **x,
    const size_t idx, const size_t ndata, const size_t nleaf) {
  root[i].n = ndata;
  root[i].x[0] = x[0] + idx;
  root[i].x[1] = x[1] + idx;
  root[i].x[2] = x[2] + idx;
#ifdef BENCHMARK_TREE_PREALLOC
  root[i].left = NULL;
  root[i].right = NULL;
#endif

  /* Compute only the corner of `data` if this is a leaf node. */
  if (ndata <= nleaf) {
    for (int k = 0; k < 3; k++) {
      real min, max;
      min = max = root[i].x[k][0];
      for (size_t j = 1; j < ndata; j++) {
        if (min > root[i].x[k][j]) min = root[i].x[k][j];
        if (max < root[i].x[k][j]) max = root[i].x[k][j];
      }
      root[i].cen[k] = (min + max) * 0.5;
      root[i].wid[k] = (max - min) * 0.5;
    }
    return false;
  }

  /* Find the direction with the largest variance, and corners of `data`. */
  int dir = 0;          /* direction with the largest variance */
  real var_max = 0;
  for (int k = 0; k < 3; k++) {
    real mean = 0;
    real min, max;
    min = max = root[i].x[k][0];

    /* Mean and min/max. */
    for (size_t j = 0; j < ndata; j++) {
      mean += root[i].x[k][j];
      if (min > root[i].x[k][j]) min = root[i].x[k][j];
      if (max < root[i].x[k][j]) max = root[i].x[k][j];
    }
    mean /= (real) ndata;
    root[i].cen[k] = (min + max) * 0.5;
    root[i].wid[k] = (max - min) * 0.5;

    real var = 0;
    for (size_t j = 0; j < ndata; j++) {
      real d = root[i].x[k][j] - mean;
      var += d * d;
    }

    if (var > var_max) {
      dir = k;
      var_max = var;
    }
  }

  size_t n = ndata >> 1;        /* index of the median point */
  /* Split the dataset by the median in the direction with largest variance. */
  qselect(root[i].x, 0, n, ndata, &dir);
  return true;
}
#else
/******************************************************************************
Function `kdtree_init`:
  Initialise a node of the k-D tree.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Address of the k-D tree node.
******************************************************************************/
static KDT *kdtree_init(real **x, const size_t idx, const size_t ndata,
    const size_t nleaf) {
  KDT *node = malloc(sizeof(KDT));
  if (!node) return NULL;

  node->n = ndata;
  node->x[0] = x[0] + idx;
  node->x[1] = x[1] + idx;
  node->x[2] = x[2] + idx;
  node->left = node->right = NULL;

  /* Compute only the corner of `data` if this is a leaf node. */
  if (ndata <= nleaf) {
    for (int k = 0; k < 3; k++) {
      real min, max;
      min = max = node->x[k][0];
      for (size_t i = 1; i < ndata; i++) {
        if (min > node->x[k][i]) min = node->x[k][i];
        if (max < node->x[k][i]) max = node->x[k][i];
      }
      node->cen[k] = (min + max) * 0.5;
      node->wid[k] = (max - min) * 0.5;
    }
    return node;
  }

  /* Find the direction with the largest variance, and corners of `data`. */
  int dir = 0;          /* direction with the largest variance */
  real var_max = 0;
  for (int k = 0; k < 3; k++) {
    real mean = 0;
    real min, max;
    min = max = node->x[k][0];

    /* Mean and min/max. */
    for (size_t i = 0; i < ndata; i++) {
      mean += node->x[k][i];
      if (min > node->x[k][i]) min = node->x[k][i];
      if (max < node->x[k][i]) max = node->x[k][i];
    }
    mean /= (real) ndata;
    node->cen[k] = (min + max) * 0.5;
    node->wid[k] = (max - min) * 0.5;

    real var = 0;
    for (size_t i = 0; i < ndata; i++) {
      real d = node->x[k][i] - mean;
      var += d * d;
    }

    if (var > var_max) {
      dir = k;
      var_max = var;
    }
  }

  size_t n = ndata >> 1;        /* index of the median point */
  /* Split the dataset by the median in the direction with largest variance. */
  qselect(node->x, 0, n, ndata, &dir);
  return node;
}
#endif

#ifdef BENCHMARK_TREE_AS_ARRAY
/******************************************************************************
Function `kdtree_build`:
  Construct the k-D tree recursively given a data set.
Arguments:
  * `root`:     root of the tree;
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
******************************************************************************/
static void kdtree_build(KDT *root, real **x, const size_t idx,
    const size_t ndata, const size_t nleaf) {
  if (!kdtree_init(root, 0, x, idx, ndata, nleaf)) return;

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
        next = kdtree_init(root, i++, root[pid].x, 0, n, nleaf) || next;
        next = kdtree_init(root, i++, root[pid].x, n, nprev - n, nleaf) || next;
      }
    }
    if (!next) return;
  }
}
#elif defined(BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `kdtree_build`:
  Construct the k-D tree recursively given a data set.
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
static size_t kdtree_build(KDT *root, size_t *inode, real **x,
    const size_t idx, const size_t ndata, const size_t nleaf) {
  const size_t this_node = *inode;
  if (kdtree_init(root, (*inode)++, x, idx, ndata, nleaf)) {
    KDT *node = root + this_node;
    size_t n = ndata >> 1;
    node->left = root + kdtree_build(root, inode, node->x, 0, n, nleaf);
    node->right = root +
        kdtree_build(root, inode, node->x, n, ndata - n, nleaf);
  }
  return this_node;
}
#else
/******************************************************************************
Function `kdtree_build`:
  Construct the k-D tree recursively given a data set.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes;
  * `err`:      error indicator.
Return:
  Root of the constructed tree.
******************************************************************************/
static KDT* kdtree_build(real **x, const size_t idx, const size_t ndata,
    const size_t nleaf, int *err) {
  if (*err) return NULL;

  KDT *node = kdtree_init(x, idx, ndata, nleaf);
  if (!node) {
    *err = EXIT_FAILURE;
    return NULL;
  }

  if (ndata > nleaf) {
    size_t n = ndata >> 1;      /* index of the median point */
    node->left = kdtree_build(node->x, 0, n, nleaf, err);
    node->right = kdtree_build(node->x, n, ndata - n, nleaf, err);
  }

  return node;
}
#endif

#if !defined(BENCHMARK_TREE_AS_ARRAY) && !defined(BENCHMARK_TREE_PREALLOC)
/******************************************************************************
Function `kdtree_free`:
  Release the memory allocated for the k-D tree.
Arguments:
  * `node`:     pointer to a root of the k-D tree.
******************************************************************************/
static void kdtree_free(KDT *node) {
  if(!node) return;
  kdtree_free(node->left);
  kdtree_free(node->right);
  free(node);
}
#endif

#ifdef PRINT_TREE
#ifdef BENCHMARK_TREE_AS_ARRAY
  #error `PRINT_TREE` not implemented for `BENCHMARK_TREE_AS_ARRAY`
#endif
static void kdtree_print(KDT *node, const int lv) {
  if (!node) return;
  printf("NODE level: %d", lv);
  if (node->left == NULL) printf(" (leaf)");
  printf("\n  ndata: %zu\n  cen: %g %g %g\n  wid: %g %g %g\n", node->n,
      node->cen[0], node->cen[1], node->cen[2],
      node->wid[0], node->wid[1], node->wid[2]);
  if (node->left == NULL) {
    printf("  data:");
    for (size_t i = 0; i < node->n; i++) {
      printf("  (%g,%g,%g)", node->x[0][i], node->x[1][i], node->x[2][i]);
    }
    printf("\n");
  }
  kdtree_print(node->left, lv + 1);
  kdtree_print(node->right, lv + 1);
}
#endif

/******************************************************************************
Function `create_kdtree`:
  Construct the k-D tree for a input dataset.
Arguments:
  * `data`:     the input dataset;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  Root of the constructed tree.
******************************************************************************/
static KDT *create_kdtree(DATA *data, const size_t nleaf) {
  if (!data || !data->n) {
    P_ERR("the input dataset is not initialised\n");
    return NULL;
  }
  if (!nleaf) {
    P_ERR("the size of leaf nodes must be non-zero\n");
    return NULL;
  }

  real *x[3];
  x[0] = data->x;
  x[1] = data->y;
  x[2] = data->z;

#if defined(BENCHMARK_TREE_AS_ARRAY) || defined(BENCHMARK_TREE_PREALLOC)
  const size_t size = kdtree_nnode(data->n, nleaf);
  if (!size) {
    P_ERR("failed to precompute the size of the tree\n");
    return NULL;
  }
  KDT *root = calloc(size, sizeof *root);
  if (!root) {
    P_ERR("failed to allocate memory for k-d tree\n");
    return NULL;
  }
#endif

#ifdef BENCHMARK_TREE_AS_ARRAY
  kdtree_build(root, x, 0, data->n, nleaf);
#elif defined(BENCHMARK_TREE_PREALLOC)
  #if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  /* Allocate memory for bounding box properties. */
  real *vol = NULL;
  if (posix_memalign((void **) &vol, BENCHMARK_MEMALIGN_BYTE,
      size * 6 * sizeof(real))) {
    P_ERR("failed to allocate aligned memory for the volume of nodes\n");
    free(root);
    return NULL;
  }
  for (size_t i = 0; i < size; i++) {
    root[i].cen = vol + i * 6;
    root[i].wid = root[i].cen + 3;
  }
  #endif
  size_t inode = 0;
  kdtree_build(root, &inode, x, 0, data->n, nleaf);
  /* Reduce memory cost if applicable. */
  if (size > inode) {
    KDT *ini = root;
    KDT *tmp = realloc(root, inode * sizeof(KDT));
    if (tmp) {
      for (size_t i = 0; i < inode; i++) {
        if (tmp[i].left) tmp[i].left = (tmp[i].left - ini) + tmp;
        if (tmp[i].right) tmp[i].right = (tmp[i].right - ini) + tmp;
      }
      root = tmp;
    }
  }
#else
  int err = 0;
  KDT *root = kdtree_build(x, 0, data->n, nleaf, &err);
  if (err) {
    P_ERR("failed to construct the k-d tree\n");
    kdtree_free(root);
    return NULL;
  }
#endif

#ifdef PRINT_TREE
  kdtree_print(root, 0);
#endif
  return root;
}


/*============================================================================*\
                        Function for distance evaluation
\*============================================================================*/

/******************************************************************************
Function `squared_distance`:
  Compute the squared Euclidean distance between two points in 3-D space.
Arguments:
  * `x1`:       x coordinate of the first point;
  * `y1`:       y coordinate of the first point;
  * `z1`:       z coordinate of the first point;
  * `x2`:       x coordinate of the second point;
  * `y2`:       y coordinate of the second point;
  * `z2`:       z coordinate of the second point.
Return:
  The squared Euclidean distance of the two points.
******************************************************************************/
/*
static inline real squared_distance(const real x1, const real y1, const real z1,
    const real x2, const real y2, const real z2) {
  register real dx = x1 - x2;
  register real dy = y1 - y2;
  register real dz = z1 - z2;
  return dx * dx + dy * dy + dz * dz;
}
*/


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

#define BENCHMARK_TREE_TYPE     BENCHMARK_STRUCT_KDTREE

/* paircount_box_smin0 */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_BOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_ZERO
#include "metric_common.c"
#include "metric_kdtree.c"
#include "dual_tree.c"

/* paircount_box */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_BOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_NONZERO
#include "metric_common.c"
#include "metric_kdtree.c"
#include "dual_tree.c"

/* paircount_smin0 */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_NOBOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_ZERO
#include "metric_common.c"
#include "metric_kdtree.c"
#include "dual_tree.c"

/* paircount */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_NOBOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_NONZERO
#include "metric_common.c"
#include "metric_kdtree.c"
#include "dual_tree.c"


/*============================================================================*\
                          Interface for pair counting
\*============================================================================*/

/******************************************************************************
Function `paircount_kdtree`:
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
int paircount_kdtree(DATA *data, const double rmin, const double rmax,
    const int csize,
#ifndef BENCHMARK_TIMING
    size_t *nnode, size_t *ndist,
#endif
    size_t *npair) {
  printf("Constructing the kd-tree ...");
  fflush(stdout);

  KDT *tree = create_kdtree(data, csize);
  if (!tree) {
    P_ERR("failed to construct kd-tree\n");
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
    if (rmin == 0) conf = calloc(1, sizeof(conf_kdtree_box_smin0));
    else conf = calloc(1, sizeof(conf_kdtree_box));
  }
  else {
    if (rmin == 0) conf = calloc(1, sizeof(conf_kdtree_smin0));
    else conf = calloc(1, sizeof(conf_kdtree));
  }
  if (!conf) {
    P_ERR("failed to allocate memory for pair counting\n");
#if defined(BENCHMARK_TREE_AS_ARRAY) || defined(BENCHMARK_TREE_PREALLOC)
  #if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
    free(tree->cen);
  #endif
    free(tree);
#else
    kdtree_free(tree);
#endif
    return EXIT_FAILURE;
  }

  /* Choose the right pair counting function, and set configurations. */
  void (*count_func) (STACK_DUAL_NODE *, void *) = NULL;
  if (data->isbox) {
    if (rmin == 0) {
      count_func = pair_kdtree_box_smin0;
      conf_kdtree_box_smin0 *config = (conf_kdtree_box_smin0 *) conf;
      config->bsize = data->bsize;
      config->hsize = data->bsize * 0.5;
      config->r2max = r2max;
    }
    else {
      count_func = pair_kdtree_box;
      conf_kdtree_box *config = (conf_kdtree_box *) conf;
      config->bsize = data->bsize;
      config->hsize = data->bsize * 0.5;
      config->r2max = r2max;
      config->r2min = r2min;
    }
  }
  else {
    if (rmin == 0) {
      count_func = pair_kdtree_smin0;
      conf_kdtree_smin0 *config = (conf_kdtree_smin0 *) conf;
      config->r2max = r2max;
    }
    else {
      count_func = pair_kdtree;
      conf_kdtree *config = (conf_kdtree *) conf;
      config->r2max = r2max;
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
      conf_kdtree_box_smin0 *config = (conf_kdtree_box_smin0 *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
    else {
      conf_kdtree_box *config = (conf_kdtree_box *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
  }
  else {
    if (rmin == 0) {
      conf_kdtree_smin0 *config = (conf_kdtree_smin0 *) conf;
      *npair = config->npair;
#ifndef BENCHMARK_TIMING
      *ndist = config->ndist;
      *nnode = config->nnode;
#endif
    }
    else {
      conf_kdtree *config = (conf_kdtree *) conf;
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
  kdtree_free(tree);
#endif
  printf(FMT_DONE);

#ifdef BENCHMARK_TIMING
  printf("Time for pair counting: " OFMT_DBL " seconds\n", sec);
#endif
  return 0;
}
