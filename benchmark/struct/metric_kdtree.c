/*******************************************************************************
* benchmark/struct/metric_kdtree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#if defined(BENCHMARK_TREE_TYPE) && defined(BENCHMARK_PAIRCOUNT_TYPE) && \
  defined(BENCHMARK_BIN_SMIN)

#if !defined(BENCHMARK_PAIRCOUNT_NAME) || !defined(BENCHMARK_BIN_SMIN_NAME)
  #error please include `metric_common.c` before `metric_kdtree.c`
#endif

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef BENCHMARK_TREE_NAME
  #undef BENCHMARK_TREE_NAME
#endif
#ifdef BENCHMARK_TREE_DTYPE
  #undef BENCHMARK_TREE_DTYPE
#endif

#if     BENCHMARK_TREE_TYPE == BENCHMARK_STRUCT_KDTREE
  #define BENCHMARK_TREE_DTYPE          KDT
  #define BENCHMARK_TREE_NAME           kdtree
#else
  #error unexpected definition of `BENCHMARK_TREE_TYPE` for `metric_kdtree.c`
#endif

/* Macros for generating datatype and function names. */
#ifndef CONCAT_STRING
  #define CONCAT_STRING(a,b,c,d)        a##_##b##c##d
#endif

#ifndef PRIVATE_NAME
  #define PRIVATE_NAME(a,b,c,d)         CONCAT_STRING(a,b,c,d)
#endif


/*============================================================================*\
                             Macros for node access
\*============================================================================*/

#ifdef BENCHMARK_GET_NODES
  #undef BENCHMARK_GET_NODES
#endif
#ifdef BENCHMARK_LCHILD
  #undef BENCHMARK_LCHILD
#endif
#ifdef BENCHMARK_RCHILD
  #undef BENCHMARK_RCHILD
#endif
#ifdef BENCHMARK_SELF_NODE
  #undef BENCHMARK_SELF_NODE
#endif
#ifdef BENCHMARK_NOT_LEAF
  #undef BENCHMARK_NOT_LEAF
#endif

#ifdef BENCHMARK_TREE_AS_ARRAY
  #define       BENCHMARK_GET_NODES(node1,node2)                        \
    const size_t num_##node1 = nodes->a;                                \
    const size_t num_##node2 = nodes->b;                                \
    BENCHMARK_TREE_DTYPE *(node1) =                                     \
        (BENCHMARK_TREE_DTYPE *) stack->root + num_##node1;             \
    BENCHMARK_TREE_DTYPE *(node2) =                                     \
        (BENCHMARK_TREE_DTYPE *) stack->root + num_##node2;
  #define       BENCHMARK_LCHILD(node)          (((num_##node) << 1) + 1)
  #define       BENCHMARK_RCHILD(node)          (((num_##node) + 1) << 1)
  #define       BENCHMARK_SELF_NODE(node)       (num_##node)
  #define       BENCHMARK_NOT_LEAF(node)        ((node)->n > stack->nleaf)
#else
  #define       BENCHMARK_GET_NODES(node1,node2)                        \
    BENCHMARK_TREE_DTYPE *(node1) = (BENCHMARK_TREE_DTYPE *) nodes->a;  \
    BENCHMARK_TREE_DTYPE *(node2) = (BENCHMARK_TREE_DTYPE *) nodes->b;
  #define       BENCHMARK_LCHILD(node)          (node)->left
  #define       BENCHMARK_RCHILD(node)          (node)->right
  #define       BENCHMARK_SELF_NODE(node)       (node)
  #define       BENCHMARK_NOT_LEAF(node)        (node)->left
#endif

/*============================================================================*\
             Data structure and functions for distance evaluations
\*============================================================================*/

/* Data structure of configuration parameters for pair counting. */
typedef struct {
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
  real bsize;           /* side length of periodic box */
  real hsize;           /* half the side length        */
#endif
  real r2max;           /* maximum squared distance of interest */
#if     BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  real r2min;           /* minimum squared distance of interest */
#endif
  size_t npair;         /* number of pairs */
#ifndef BENCHMARK_TIMING
  size_t nnode;         /* number of visited nodes */
  size_t ndist;         /* number of distance evaluations */
#endif
} PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME);


/******************************************************************************
Function `compare_node_dist_<BENCHMARK_TREE_NAME><BENCHMARK_PAIRCOUNT_NAME>
  <BENCHMARK_BIN_SMIN_NAME>`:
  Compute the minimum/maximum squared distance between two nodes,
  and check if they are outside the distance range of interest.
Arguments:
  * `node1`:    pointer to the first node;
  * `node2`:    pointer to the second node;
  * `conf`:     the structure for pair count configurations;
  * `shift`:    coordinate offsets for periodic boundary conditions.
Return:
  True of the distance between the nodes is not of interest; false otherwise.
******************************************************************************/
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
static inline int PRIVATE_NAME(compare_node_dist, BENCHMARK_TREE_NAME,
    BENCHMARK_PAIRCOUNT_NAME, BENCHMARK_BIN_SMIN_NAME)
    (const BENCHMARK_TREE_DTYPE *node1, const BENCHMARK_TREE_DTYPE *node2,
    const PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME) *conf, real shift[3]) {
  #if           BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  real dx = node1->cen[0] - node2->cen[0];
  if (dx > conf->hsize) {
    shift[0] = -conf->bsize;
    dx -= conf->bsize;
  }
  else if (dx < -conf->hsize) {
    shift[0] = conf->bsize;
    dx += conf->bsize;
  }
  if (dx < 0) dx = -dx;
  real dy = node1->cen[1] - node2->cen[1];
  if (dy > conf->hsize) {
    shift[1] = -conf->bsize;
    dy -= conf->bsize;
  }
  else if (dy < -conf->hsize) {
    shift[1] = conf->bsize;
    dy += conf->bsize;
  }
  if (dy < 0) dy = -dy;
  real dz = node1->cen[2] - node2->cen[2];
  if (dz > conf->hsize) {
    shift[2] = -conf->bsize;
    dz -= conf->bsize;
  }
  else if (dz < -conf->hsize) {
    shift[2] = conf->bsize;
    dz += conf->bsize;
  }
  if (dz < 0) dz = -dz;

  const real wx = node1->wid[0] + node2->wid[0];
  const real wy = node1->wid[1] + node2->wid[1];
  const real wz = node1->wid[2] + node2->wid[2];

  register real sum = 0;
  register real d = dx - wx;
  if (d > 0) sum += d * d;
  d = dy - wy;
  if (d > 0) sum += d * d;
  d = dz - wz;
  if (d > 0) sum += d * d;

  if (sum >= conf->r2max) return 1;
    #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  const real sum_min = sum;
    #endif
  sum = 0;
  d = dx + wx;
  sum += d * d;
  d = dy + wy;
  sum += d * d;
  d = dz + wz;
  sum += d * d;
  #else      /* BENCHMARK_SIMD  !=  BENCHMARK_SIMD_NONE */
  /* Load only 3 elements to the vector. */
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
      #ifdef SINGLE_PREC
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0,
      0, 0, 0, 0);
      #else       /* double precision */
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFFFFFFFFFFULL,
      0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0);
      #endif
  vec_r x1 = SIMD_MSKLOAD_REAL(node1->cen, mask3);
  vec_r x2 = SIMD_MSKLOAD_REAL(node2->cen, mask3);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
  const vec_msk mask3 = 7;
  vec_r x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->cen);
  vec_r x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->cen);
    #endif
  vec_r vd = SIMD_SUB_REAL(x1, x2);
  /* Periodic wrapping. */
  const vec_r upper = SIMD_SET1_REAL(conf->hsize);
  const vec_r lower = SIMD_SET1_REAL(-conf->hsize);
  vec_msk mask = SIMD_CMP_REAL(vd, upper, _CMP_GT_OQ);
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
  x1 = SIMD_AND_REAL(mask, SIMD_SET1_REAL(-conf->bsize));
  mask = SIMD_CMP_REAL(vd, lower, _CMP_LT_OQ);
  x2 = SIMD_AND_REAL(mask, SIMD_SET1_REAL(conf->bsize));
  x1 = SIMD_OR_REAL(x1, x2);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
  x1 = SIMD_MSK0_MOV_REAL(mask, SIMD_SET1_REAL(-conf->bsize));
  mask = SIMD_CMP_REAL(vd, lower, _CMP_LT_OQ);
  x1 = SIMD_MSK_MOV_REAL(x1, mask, SIMD_SET1_REAL(conf->bsize));
    #endif
  SIMD_MSKSTORE_REAL(shift, mask3, x1);
  vd = SIMD_ADD_REAL(vd, x1);
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
  vd = SIMD_ANDNOT_REAL(SIMD_SET1_REAL(-0.0), vd);      /* absolute value */
  x1 = SIMD_MSKLOAD_REAL(node1->wid, mask3);
  x2 = SIMD_MSKLOAD_REAL(node2->wid, mask3);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
  vd = SIMD_ABS_REAL(vd);
  x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->wid);
  x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->wid);
    #endif
  vec_r vw = SIMD_ADD_REAL(x1, x2);

  /* Minimum distance. */
  vec_r delta = SIMD_SUB_REAL(vd, vw);
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
  mask = SIMD_CMP_REAL(delta, SIMD_ZEROS_REAL(), _CMP_GT_OQ);
  delta = SIMD_AND_REAL(delta, mask);
  delta = SIMD_MUL_REAL(delta, delta);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
  mask = SIMD_CMP_REAL(delta, SIMD_ZEROS_REAL(), _CMP_GT_OQ);
  delta = SIMD_MSK0_MUL_REAL(mask, delta, delta);
    #endif
  real *res = (real *) &delta;
  register real sum = res[0] + res[1] + res[2];
  if (sum >= conf->r2max) return 1;
    #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  const real sum_min = sum;
    #endif

  /* Maximum distance. */
  delta = SIMD_ADD_REAL(vd, vw);
  delta = SIMD_MUL_REAL(delta, delta);
  res = (real *) &delta;
  sum = res[0] + res[1] + res[2];
  #endif     /* BENCHMARK_SIMD */

  #ifndef BENCHMARK_TIMING
  conf->nnode++;
  #endif
  #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_ZERO
  if (sum < conf->r2max) return -1;
  #else
  if (sum < conf->r2min) return 1;
  if (sum_min >= conf->r2min && sum < conf->r2max) return -1;
  #endif
  return 0;
}
#else   /* BENCHMARK_PAIRCOUNT_TYPE  !=  BENCHMARK_PAIRCOUNT_BOX */
static inline int PRIVATE_NAME(compare_node_dist, BENCHMARK_TREE_NAME,
    BENCHMARK_PAIRCOUNT_NAME, BENCHMARK_BIN_SMIN_NAME)
    (const BENCHMARK_TREE_DTYPE *node1, const BENCHMARK_TREE_DTYPE *node2,
    PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME) *conf) {
  #if           BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  const real dx = (node1->cen[0] > node2->cen[0]) ?
      node1->cen[0] - node2->cen[0] : node2->cen[0] - node1->cen[0];
  const real dy = (node1->cen[1] > node2->cen[1]) ?
      node1->cen[1] - node2->cen[1] : node2->cen[1] - node1->cen[1];
  const real dz = (node1->cen[2] > node2->cen[2]) ?
      node1->cen[2] - node2->cen[2] : node2->cen[2] - node1->cen[2];
  const real wx = node1->wid[0] + node2->wid[0];
  const real wy = node1->wid[1] + node2->wid[1];
  const real wz = node1->wid[2] + node2->wid[2];

  register real sum = 0;
  register real d = dx - wx;
  if (d > 0) sum += d * d;
  d = dy - wy;
  if (d > 0) sum += d * d;
  d = dz - wz;
  if (d > 0) sum += d * d;

  if (sum >= conf->r2max) return 1;
    #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  const real sum_min = sum;
    #endif
  sum = 0;
  d = dx + wx;
  sum += d * d;
  d = dy + wy;
  sum += d * d;
  d = dz + wz;
  sum += d * d;
  #else      /* BENCHMARK_SIMD  !=  BENCHMARK_SIMD_NONE */
  /* Load only 3 elements to the vector. */
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
      #ifdef SINGLE_PREC
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0,
      0, 0, 0, 0);
      #else       /* double precision */
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFFFFFFFFFFULL,
      0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0);
      #endif
  vec_r x1 = SIMD_MSKLOAD_REAL(node1->cen, mask3);
  vec_r x2 = SIMD_MSKLOAD_REAL(node2->cen, mask3);
  vec_r vd = SIMD_SUB_REAL(x1, x2);
  vd = SIMD_ANDNOT_REAL(SIMD_SET1_REAL(-0.0), vd);      /* absolute value */
  x1 = SIMD_MSKLOAD_REAL(node1->wid, mask3);
  x2 = SIMD_MSKLOAD_REAL(node2->wid, mask3);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
  const vec_msk mask3 = 7;
  vec_r x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->cen);
  vec_r x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->cen);
  vec_r vd = SIMD_SUB_REAL(x1, x2);
  vd = SIMD_ABS_REAL(vd);
  x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->wid);
  x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->wid);
    #endif
  vec_r vw = SIMD_ADD_REAL(x1, x2);

  /* Minimum distance. */
  vec_r delta = SIMD_SUB_REAL(vd, vw);
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
  vec_msk mask = SIMD_CMP_REAL(delta, SIMD_ZEROS_REAL(), _CMP_GT_OQ);
  delta = SIMD_AND_REAL(delta, mask);
  delta = SIMD_MUL_REAL(delta, delta);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
  vec_msk mask = SIMD_CMP_REAL(delta, SIMD_ZEROS_REAL(), _CMP_GT_OQ);
  delta = SIMD_MSK0_MUL_REAL(mask, delta, delta);
    #endif
  real *res = (real *) &delta;
  register real sum = res[0] + res[1] + res[2];
  if (sum >= conf->r2max) return 1;
    #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  const real sum_min = sum;
    #endif

  /* Maximum distance. */
  delta = SIMD_ADD_REAL(vd, vw);
  delta = SIMD_MUL_REAL(delta, delta);
  res = (real *) &delta;
  sum = res[0] + res[1] + res[2];
  #endif     /* BENCHMARK_SIMD */

  #ifndef BENCHMARK_TIMING
  conf->nnode++;
  #endif
  #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_ZERO
  if (sum < conf->r2max) return -1;
  #else
  if (sum < conf->r2min) return 1;
  if (sum_min >= conf->r2min && sum < conf->r2max) return -1;
  #endif
  return 0;
}
#endif

#endif
