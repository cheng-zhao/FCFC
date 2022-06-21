/*******************************************************************************
* kdtree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "kdtree.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef MPI
#include <mpi.h>
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
static int kth_compare(real *const x[static 3], const size_t i, const size_t j,
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
#ifdef QSELECT_NAME
  #undef QSELECT_NAME
#endif

/* `qselect` for coordinates only. */
typedef struct {
  real *x[FCFC_XDIM];
} DATA_x;

#define QSELECT_COMPARE(d,i,j,k)        kth_compare((d)->x,i,j,*((int *)(k)))
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

#define QSELECT_COMPARE(d,i,j,k)        kth_compare((d)->x,i,j,*((int *)(k)))
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
                        Functions for tree construction
\*============================================================================*/

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
Function `kdtree_realloc`:
  Reallocate memory for a k-d tree.
Arguments:
  * `root`:     root of the k-d tree;
  * `size`:     the original number of tree nodes;
  * `nnode`:    the expected number of tree nodes after memory reallocation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int kdtree_realloc(KDT **root, const size_t size, const size_t nnode) {
  if (SIZE_MAX / 3 < nnode) return FCFC_ERR_ARG;

  /* Record addresses of the original tree. */
  KDT *ini = *root;
  real *cini = ini->cen;
  real *wini = ini->wid;
  KDT *tmp = realloc(*root, sizeof(KDT) * nnode);
  if (!tmp) return FCFC_ERR_MEMORY;

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

  real *vtmp = realloc(cini, sizeof(real) * nnode * 3);
  if (!vtmp) return FCFC_ERR_MEMORY;
#ifdef OMP
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < nnode; i++) tmp[i].cen = vtmp + i * 3;

  if (!(vtmp = realloc(wini, sizeof(real) * nnode * 3))) return FCFC_ERR_MEMORY;
#ifdef OMP
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < nnode; i++) tmp[i].wid = vtmp + i * 3;

  return 0;
}

/******************************************************************************
Function `kdtree_init`:
  Initialise a node of the k-d tree.
Arguments:
  * `root`:     root of the tree;
  * `i`:        index of the current node;
  * `x`:        coordinates of the input dataset;
  * `w`:        weights of the input dataset;
  * `idx`:      starting index of the dataset;
  * `ndata`:    number of elements of the input dataset;
  * `bsize`:    side lengths of the periodic box;
  * `nleaf`:    maximum number of data points on leaf nodes.
Return:
  True if the node has children; false otherwise.
******************************************************************************/
static bool kdtree_init(KDT *root, const size_t i, real *x[static FCFC_XDIM],
    real *w, const size_t idx, const size_t ndata,
#ifdef FCFC_METRIC_PERIODIC
    const real bsize[static 3],
#endif
    const size_t nleaf) {
  root[i].n = ndata;
  for (int k = 0; k < FCFC_XDIM; k++) root[i].x[k] = x[k] + idx;
  root[i].w = (w) ? w + idx : NULL;
  root[i].left = root[i].right = NULL;

  /* Compute only the AABB of `data` if this is a leaf node candidate. */
  if (ndata <= nleaf) {
#ifdef FCFC_METRIC_PERIODIC
    bool isleaf = true;
#endif
    for (int k = 0; k < 3; k++) {
      register real min, max;
      min = max = root[i].x[k][0];
      for (size_t j = 1; j < ndata; j++) {
        if (min > root[i].x[k][j]) min = root[i].x[k][j];
        if (max < root[i].x[k][j]) max = root[i].x[k][j];
      }
#ifdef FCFC_METRIC_PERIODIC
      /* Subdivide the node if the extents are larger than half box sizes. */
      if (max - min >= bsize[k] / 2) {
        isleaf = false;
        break;
      }
#endif
      root[i].cen[k] = (min + max) * 0.5;
      root[i].wid[k] = (max - min) * 0.5;
    }
#ifdef FCFC_METRIC_PERIODIC
    if (isleaf) return false;
#else
    return false;
#endif
  }

  /* Find the direction with the largest variance, and corners of `data`. */
  int dir = 0;          /* direction with the largest variance */
  real var_max = 0;
  for (int k = 0; k < 3; k++) {
    register real mean, min, max;
    mean = min = max = root[i].x[k][0];

    /* Mean and min/max. */
    for (size_t j = 1; j < ndata; j++) {
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
  if (w) {
    DATA_xw data;
    for (int k = 0; k < FCFC_XDIM; k++) data.x[k] = root[i].x[k];
    data.w = root[i].w;
    qselect_xw(&data, 0, n, ndata, &dir);
  }
  else {
    DATA_x data;
    for (int k = 0; k < FCFC_XDIM; k++) data.x[k] = root[i].x[k];
    qselect_x(&data, 0, n, ndata, &dir);
  }
  return true;
}

/******************************************************************************
Function `kdtree_build`:
  Construct the k-d tree recursively given a data set.
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
static size_t kdtree_build(
#ifdef FCFC_METRIC_PERIODIC
    KDT **root,
#else
    KDT *root,
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
    if (SIZE_MAX / 2 < *size || kdtree_realloc(root, *size, *size << 1)) {
      P_EXT("failed to reallocate memory for the tree\n");
      exit(FCFC_ERR_MEMORY);
    }
    *size <<= 1;
  }
  if (kdtree_init(*root, this_node, x, w, idx, ndata, bsize, nleaf)) {
    size_t n = ndata >> 1;
    KDT *node = *root + this_node;
    real *this_x[FCFC_XDIM];
    for (int k = 0; k < FCFC_XDIM; k++) this_x[k] = node->x[k];
    size_t i = kdtree_build(root, inode, this_x, node->w,
        0, n, bsize, size, nleaf);
    (*root)[this_node].left = *root + i;
    i = kdtree_build(root, inode, this_x, (*root)[this_node].w,
        n, ndata - n, bsize, size, nleaf);
    (*root)[this_node].right = *root + i;
  }
#else
  if (kdtree_init(root, this_node, x, w, idx, ndata, nleaf)) {
    KDT *node = root + this_node;
    size_t n = ndata >> 1;
    node->left = root +
        kdtree_build(root, inode, node->x, node->w, 0, n, nleaf);
    node->right = root +
        kdtree_build(root, inode, node->x, node->w, n, ndata - n, nleaf);
  }
#endif
  return this_node;
}


/*============================================================================*\
                      Functions for requesting tree nodes
\*============================================================================*/

#if defined(MPI) || defined(OMP)
/******************************************************************************
Function `kdtree_level`:
  Compute the minimum tree level with at least the specified number of nodes.
Arguments:
  * `nnode`:    number of nodes at the same level.
Return:
  The minimum tree level with at least `nnode` nodes.
******************************************************************************/
static inline int kdtree_level(uint32_t nnode) {
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
static void gather_nodes(const KDT *tree, const int current,
    const int level, const KDT **nodes, size_t *inode) {
  if (!tree) return;
  if (current == level) nodes[(*inode)++] = tree;
  else {
    gather_nodes(tree->left, current + 1, level, nodes, inode);
    gather_nodes(tree->right, current + 1, level, nodes, inode);
  }
}
#endif


/*============================================================================*\
                            Interfaces for k-d tree
\*============================================================================*/

/******************************************************************************
Function `create_kdtree`:
  Construct the k-d tree for a input dataset.
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
KDT *create_kdtree(real *x[static FCFC_XDIM], real *w, const size_t ndata,
#ifdef FCFC_METRIC_PERIODIC
    const real bsize[static 3],
#endif
    const size_t nleaf, size_t *nnode) {
  if (!x || !ndata) {
    P_ERR("the input dataset is not initialised\n");
    return NULL;
  }
#ifdef FCFC_METRIC_PERIODIC
  if (!bsize) {
    P_ERR("the periodic box sizes are not specified\n");
    return NULL;
  }
#endif
  if (!nleaf) {
    P_ERR("the capacity of leaf nodes must be non-zero\n");
    return NULL;
  }

  size_t size = kdtree_nnode(ndata, nleaf);
  if (!size) {
    P_ERR("failed to precompute the size of the tree\n");
    return NULL;
  }
  KDT *root = calloc(size, sizeof *root);
  if (!root) {
    P_ERR("failed to allocate memory for k-d tree\n");
    return NULL;
  }

  /* Allocate memory for bounding box properties. */
  real *cen, *wid;
  cen = wid = NULL;
  if (!(cen = malloc(sizeof(real) * size * 3)) ||
      !(wid = malloc(sizeof(real) * size * 3))) {
    P_ERR("failed to allocate memory for the bounding volume of nodes\n");
    if (cen) free(cen);
    free(root);
    return NULL;
  }
#ifdef OMP
  #pragma omp parallel for default(none) shared(size,root,cen,wid)
#endif
  for (size_t i = 0; i < size; i++) {
    root[i].cen = cen + i * 3;
    root[i].wid = wid + i * 3;
  }

  /* Build the tree. */
  size_t inode = 0;
#ifdef FCFC_METRIC_PERIODIC
  kdtree_build(&root, &inode, x, w, 0, ndata, bsize, &size, nleaf);
#else
  kdtree_build(root, &inode, x, w, 0, ndata, nleaf);
#endif
  if (inode == 0) {
    P_ERR("failed to create the k-d tree\n");
    free(root); free(cen); free(wid);
    return NULL;
  }

  /* Reduce memory cost if applicable. */
  if (size > inode && kdtree_realloc(&root, size, inode)) {
    P_ERR("failed to reallocate memory for the tree\n");
    free(root); free(cen); free(wid);
    return NULL;
  }

  *nnode = inode;
  return root;
}

#ifdef MPI
/******************************************************************************
Function `kdtree_broadcast`:
  Broadcast the k-d tree to all MPI tasks.
Arguments:
  * `root`:     pointer to the root of the k-d tree;
  * `nnode`:    total number of tree nodes;
  * `wt`:       indicate whether to broadcast weights;
  * `src`:      the source for the broadcast;
  * `rank`:     ID of MPI task.
******************************************************************************/
void kdtree_broadcast(KDT **root, size_t *nnode, const bool wt, const int src,
    const int rank) {
  if (!root || !nnode ||
      (src == rank && (!(*root) || !(*nnode) || !((*root)->n)))) {
    P_ERR("the k-d tree is not initialized for broadcasting\n");
    FCFC_QUIT(FCFC_ERR_ARG);
  }
  if (src < 0 || rank < 0) {
    P_ERR("invalid MPI ranks for broadcasting the k-d tree\n");
    FCFC_QUIT(FCFC_ERR_ARG);
  }

  /* Broadcast the number of tree nodes and data points. */
  size_t ndata = 0;
  if (rank == src) ndata = (*root)->n;
  MPI_Request req[5 + FCFC_XDIM];
  if (MPI_Ibcast(nnode, 1, FCFC_MPI_SIZE_T, src, MPI_COMM_WORLD, req) ||
      MPI_Ibcast(&ndata, 1, FCFC_MPI_SIZE_T, src, MPI_COMM_WORLD, req + 1) ||
      MPI_Waitall(2, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast the k-d tree\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  /* Allocate memory for the addresses and numbers of data points. */
  ptrdiff_t *addr;
  size_t *num;
  if (!(addr = calloc(*nnode * 3, sizeof(ptrdiff_t))) ||
      !(num = calloc(*nnode, sizeof(size_t)))) {
    P_ERR("failed to allocate memory for broadcasting the k-d tree\n");
    FCFC_QUIT(FCFC_ERR_MEMORY);
  }

  real *cen, *wid, *x[FCFC_XDIM], *w;
  cen = wid = w = NULL;
  for (int k = 0; k < FCFC_XDIM; k++) x[k] = NULL;
  if (rank == src) {
    cen = (*root)->cen;
    wid = (*root)->wid;
    for (int k = 0; k < FCFC_XDIM; k++) x[k] = (*root)->x[k];
    w = (*root)->w;

    /* Gather the number of data points and compute memory offsets. */
#ifdef OMP
  #pragma omp parallel for default(none) shared(nnode,root,num,addr)
#endif
    for (size_t i = 0; i < *nnode; i++) {
      num[i] = (*root)[i].n;
      size_t j = i * 3;
      addr[j++] = (*root)[i].x[0] - (*root)[0].x[0];
      addr[j++] = ((*root)[i].left) ? (*root)[i].left - (*root) : 0;
      addr[j] = ((*root)[i].right) ? (*root)[i].right - (*root) : 0;
    }
  }
  else {        /* rank != src */
    /* Allocate memory. */
    if (!(*root = malloc(*nnode * sizeof(KDT))) ||
        !(cen = malloc(*nnode * 3 * sizeof(real))) ||
        !(wid = malloc(*nnode * 3 * sizeof(real))) ||
        !(x[0] = malloc(ndata * sizeof(real))) ||
        !(x[1] = malloc(ndata * sizeof(real))) ||
        !(x[2] = malloc(ndata * sizeof(real))) ||
#if FCFC_XDIM == 4
        !(x[3] = malloc(ndata * sizeof(real))) ||
#endif
        (wt && !(w = malloc(ndata * sizeof(real))))) {
      P_ERR("failed to allocate memory for the k-d tree\n");
      FCFC_QUIT(FCFC_ERR_MEMORY);
    }
  }

  /* Broadcast information of the tree. */
  if (MPI_Ibcast(cen, *nnode * 3, FCFC_MPI_REAL, src, MPI_COMM_WORLD, req) ||
      MPI_Ibcast(wid, *nnode * 3, FCFC_MPI_REAL, src, MPI_COMM_WORLD,
          req + 1) ||
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
    P_ERR("failed to broadcast the k-d tree\n");
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
      (*root)[i].wid = wid + j;
      (*root)[i].n = num[i];
      for (int k = 0; k < FCFC_XDIM; k++) (*root)[i].x[k] = x[k] + addr[j];
      (*root)[i].w = (wt) ? w + addr[j] : NULL;
      (*root)[i].left = (addr[j + 1]) ? (*root) + addr[j + 1] : NULL;
      (*root)[i].right = (addr[j + 2]) ? (*root) + addr[j + 2] : NULL;
    }
  }

  free(addr);
  free(num);
}
#endif

#if defined(MPI) || defined(OMP)
/******************************************************************************
Function `kdtree_get_nodes`:
  Get all tree nodes at the level with at least the specific number of nodes.
Arguments:
  * `root`:     pointer to the root of the k-d tree;
  * `nmin`:     the minimum number of nodes to be requested;
  * `nnode`:    the number of requested nodes.
Return:
  Address of the node list on success; NULL on error.
******************************************************************************/
KDT **kdtree_get_nodes(const KDT *root, const uint32_t nmin, size_t *nnode) {
  /* Compute the tree level with at least `nmin` nodes. */
  const int level = kdtree_level(nmin);
  /* The number of tree nodes at the given level assuming a balanced tree. */
  size_t nmax = 1U << level;

  /* Request tree nodes. */
  KDT **nodes = malloc(nmax * sizeof(KDT *));
  if (!nodes) {
    P_ERR("failed to request tree nodes\n");
    return NULL;
  }
  *nnode = 0;
  gather_nodes(root, 0, level, (const KDT **) nodes, nnode);

  /* Reduce memory cost if applicable. */
  if (*nnode < nmax) {
    KDT **tmp = realloc(nodes, *nnode * sizeof(KDT *));
    if (tmp) nodes = tmp;
  }
  return nodes;
}
#endif

/******************************************************************************
Function `kdtree_free`:
  Release the memory allocated for the k-d tree.
Arguments:
  * `root`:     pointer to the root of the k-d tree.
******************************************************************************/
void kdtree_free(KDT *root) {
  if(!root) return;
  free(root->cen);
  free(root->wid);
  free(root);
}
