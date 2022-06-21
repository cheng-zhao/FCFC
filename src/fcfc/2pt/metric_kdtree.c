/*******************************************************************************
* 2pt/metric_kdtree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#if defined(FCFC_TREE_TYPE) && defined(FCFC_CNT_TYPE) && \
  defined(FCFC_BIN_TYPE) && defined(FCFC_TAB_TYPE) && \
  defined(FCFC_STAB_WIDTH) && defined(FCFC_PTAB_WIDTH) && \
  defined(FCFC_BIN_SMIN) && defined(FCFC_BIN_PMIN) && defined(FCFC_CNT_WT)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#if !defined(FCFC_TREE_NAME) || !defined(FCFC_TREE_DTYPE) ||            \
  !defined(FCFC_CNT_NAME) || !defined(FCFC_BIN_NAME) ||                 \
  !defined(FCFC_TAB_NAME) || !defined(FCFC_SWIDTH_NAME) ||              \
  !defined(FCFC_STAB_DTYPE) || !defined(FCFC_PWIDTH_NAME) ||            \
  !defined(FCFC_PTAB_DTYPE) || !defined(FCFC_SMIN_NAME) ||              \
  !defined(FCFC_PMIN_NAME) || !defined(FCFC_STAB_MASK) ||               \
  !defined(FCFC_PTAB_MASK) || !defined(FCFC_WT_NAME)
  #error please include `metric_kdtree.c` in `dual_tree.c`
#endif

/*============================================================================*\
                     Functions for checking node distances
\*============================================================================*/

/******************************************************************************
Function `compare_dual_node_<IDENTIFIERS>`:
  Compute the minimum/maximum squared distance between two nodes,
  and check if they are outside the distance range of interest.
Arguments:
  * `node1`:    pointer to the first node;
  * `node2`:    pointer to the second node;
  * `conf`:     the structure for pair count configurations.
Return:
  Indicate whether the separation range is inside, is outside, or intersects
  with the distance range of interest.
******************************************************************************/
static inline fcfc_node_dist_e PRIVATE_NAME(compare_dual_node) (
    const FCFC_TREE_DTYPE *node1, const FCFC_TREE_DTYPE *node2,
    const PRIVATE_NAME(conf) *conf) {
#if             FCFC_SIMD  ==  FCFC_SIMD_NONE
  register real dx = REAL_ABS(node1->cen[0] - node2->cen[0]);
  register real dy = REAL_ABS(node1->cen[1] - node2->cen[1]);
  register real dz = REAL_ABS(node1->cen[2] - node2->cen[2]);
  register real wx = node1->wid[0] + node2->wid[0];
  register real wy = node1->wid[1] + node2->wid[1];
  register real wz = node1->wid[2] + node2->wid[2];

  /* Check the minimum separation. */
  register real sum = 0;
  register real d = dx - wx;
  if (d > 0) sum += d * d;
  d = dy - wy;
  if (d > 0) sum += d * d;
  d = dz - wz;
  if (d > 0) sum += d * d;

  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  if (sum >= conf->ts2max) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  if (sum >= conf->s2max) return FCFC_NODE_DIST_OUT;
  #endif
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const real sum_min = sum;
  #endif

  /* Check the maximum separation. */
  sum = 0;
  d = dx + wx;
  sum += d * d;
  d = dy + wy;
  sum += d * d;
  d = dz + wz;
  sum += d * d;
#else        /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
  /* Load only 3 elements to the vector. */
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
    #ifdef      SINGLE_PREC
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0,
      0, 0, 0, 0);
    #else    /* double precision */
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFFFFFFFFFFULL,
      0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0);
    #endif
  vec_r x1 = SIMD_MSKLOAD_REAL(node1->cen, mask3);
  vec_r x2 = SIMD_MSKLOAD_REAL(node2->cen, mask3);
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  const vec_msk mask3 = 7;              /* 0b0111 */
  vec_r x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->cen);
  vec_r x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->cen);
  #endif
  vec_r vd = SIMD_ABS_REAL(SIMD_SUB_REAL(x1, x2));
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  x1 = SIMD_MSKLOAD_REAL(node1->wid, mask3);
  x2 = SIMD_MSKLOAD_REAL(node2->wid, mask3);
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->wid);
  x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->wid);
  #endif
  vec_r vw = SIMD_ADD_REAL(x1, x2);

  /* Check the minimum separation. */
  vec_r delta = SIMD_SUB_REAL(vd, vw);
  vec_msk mask = SIMD_CMP_REAL(delta, SIMD_ZEROS_REAL(), _CMP_GT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  delta = SIMD_AND_REAL(delta, mask);
  delta = SIMD_MUL_REAL(delta, delta);
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  delta = SIMD_MSK0_MUL_REAL(mask, delta, delta);
  #endif
  real *res = (real *) &delta;
  register real sum = res[0] + res[1] + res[2];

  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  if (sum >= conf->ts2max) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  if (sum >= conf->s2max) return FCFC_NODE_DIST_OUT;
  #endif
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const real sum_min = sum;
  #endif

  /* Check the maximum separation. */
  delta = SIMD_ADD_REAL(vd, vw);
  delta = SIMD_MUL_REAL(delta, delta);
  res = (real *) &delta;
  sum = res[0] + res[1] + res[2];
#endif       /* FCFC_SIMD */

#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO  &&             \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (sum < conf->ts2min) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (sum < conf->s2min) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (sum < conf->p2min) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  if (sum < conf->ms2max) return FCFC_NODE_DIST_IN;
  #endif
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
  if (sum < conf->s2max) return FCFC_NODE_DIST_IN;
  #else
  if (sum < conf->s2min) return FCFC_NODE_DIST_OUT;
  if (sum_min >= conf->s2min && sum < conf->s2max) return FCFC_NODE_DIST_IN;
  #endif
#endif       /* FCFC_BIN_TYPE */

  return FCFC_NODE_DIST_X;
}

/******************************************************************************
Function `compare_single_node_<IDENTIFIERS>`:
  Compute the maximum (squared) separation between points in a single node,
  and check if it is not smaller than the minimum distance of interest.
Arguments:
  * `node`:     pointer to the node;
  * `conf`:     the structure for pair count configurations.
Return:
  Indicate whether the separation range is inside, is outside, or intersects
  with the distance range of interest.
******************************************************************************/
static inline fcfc_node_dist_e PRIVATE_NAME(compare_single_node) (
    const FCFC_TREE_DTYPE *node, const PRIVATE_NAME(conf) *conf) {
  register real wx = node->wid[0] * 2;
  register real wy = node->wid[1] * 2;
  register real wz = node->wid[2] * 2;
  register real dist = wx * wx + wy * wy + wz * wz;     /* maximum separation */

#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO  &&             \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist < conf->ts2min) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist < conf->s2min) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist < conf->p2min) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist < conf->ms2max) return FCFC_NODE_DIST_IN;
  #endif
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist < conf->s2min) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist < conf->s2max) return FCFC_NODE_DIST_IN;
  #endif
#endif       /* FCFC_BIN_TYPE */
  return FCFC_NODE_DIST_X;
}

#endif
