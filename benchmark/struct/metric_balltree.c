/*******************************************************************************
* benchmark/struct/metric_balltree.c: this file is part of the FCFC program.

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

#if     BENCHMARK_TREE_TYPE == BENCHMARK_STRUCT_BALLTREE
  #define BENCHMARK_TREE_DTYPE          BLT
  #define BENCHMARK_TREE_NAME           balltree
#else
  #error unexpected definition of `BENCHMARK_TREE_TYPE` for `metric_balltree.c`
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
  real rmax;            /* maximum distance of interest */
  real r2max;           /* maximum squared distance of interest */
#if     BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  real rmin;            /* minimum distance of interest */
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
    PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME) *conf, real shift[3]) {
  register real sum = 0;
  #if           BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  register real d = node1->cen[0] - node2->cen[0];
  if (d > conf->hsize) shift[0] = -conf->bsize;
  else if (d < -conf->hsize) shift[0] = conf->bsize;
  d += shift[0];
  sum += d * d;
  d = node1->cen[1] - node2->cen[1];
  if (d > conf->hsize) shift[1] = -conf->bsize;
  else if (d < -conf->hsize) shift[1] = conf->bsize;
  d += shift[1];
  sum += d * d;
  d = node1->cen[2] - node2->cen[2];
  if (d > conf->hsize) shift[2] = -conf->bsize;
  else if (d < -conf->hsize) shift[2] = conf->bsize;
  d += shift[2];
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

  /* Distance between spherical centres. */
  vd = SIMD_MUL_REAL(vd, vd);
  real *res = (real *) &vd;
  sum = res[0] + res[1] + res[2];
  #endif

  #ifndef BENCHMARK_TIMING
  conf->nnode++;
  #endif
  register real sr = node1->r + node2->r;
  register real dtmp = conf->rmax + sr;
  if (sum >= dtmp * dtmp) return 1;

  #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_ZERO
  dtmp = conf->rmax - sr;
  if (dtmp > 0 && sum < dtmp * dtmp) return -1;
  #else
  dtmp = conf->rmin - sr;
  if (dtmp > 0 && sum < dtmp * dtmp) return 1;
  dtmp = conf->rmax - sr;
  if (dtmp > 0 && sum < dtmp * dtmp) {
    dtmp = conf->rmin + sr;
    if (sum >= dtmp * dtmp) return -1;
  }
  #endif
  return 0;
}
#else
static inline int PRIVATE_NAME(compare_node_dist, BENCHMARK_TREE_NAME,
    BENCHMARK_PAIRCOUNT_NAME, BENCHMARK_BIN_SMIN_NAME)
    (const BENCHMARK_TREE_DTYPE *node1, const BENCHMARK_TREE_DTYPE *node2,
    PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME) *conf) {
  #if           BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  register real sum = 0;
  register real d = node1->cen[0] - node2->cen[0];
  sum += d * d;
  d = node1->cen[1] - node2->cen[1];
  sum += d * d;
  d = node1->cen[2] - node2->cen[2];
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
  vd = SIMD_MUL_REAL(vd, vd);
  real *res = (real *) &vd;
  register real sum = res[0] + res[1] + res[2];
  #endif     /* BENCHMARK_SIMD */

  #ifndef BENCHMARK_TIMING
  conf->nnode++;
  #endif
  register real sr = node1->r + node2->r;
  register real dcmp = conf->rmax + sr;
  if (sum >= dcmp * dcmp) return 1;

  #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_ZERO
  dcmp = conf->rmax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return -1;
  #else
  dcmp = conf->rmin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return 1;
  dcmp = conf->rmax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) {
    dcmp = conf->rmin + sr;
    if (sum >= dcmp * dcmp) return -1;
  }
  #endif
  return 0;
}
#endif

#endif
