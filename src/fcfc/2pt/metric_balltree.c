/*******************************************************************************
* 2pt/metric_balltree.c: this file is part of the FCFC program.

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
  #error please include `metric_balltree.c` in `dual_tree.c`
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
  register real sum = 0;
  register real d = node1->cen[0] - node2->cen[0];
  sum += d * d;
  d = node1->cen[1] - node2->cen[1];
  sum += d * d;
  d = node1->cen[2] - node2->cen[2];
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
  vec_r vd = SIMD_SUB_REAL(x1, x2);
  vd = SIMD_MUL_REAL(vd, vd);
  real *res = (real *) &vd;
  register real sum = res[0] + res[1] + res[2];
#endif       /* FCFC_SIMD */

  /* Check the minimum separation. */
  register real sr = node1->r + node2->r;
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  register real dcmp = conf->tsmax + sr;
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  register real dcmp = conf->smax + sr;
#endif
  if (sum >= dcmp * dcmp) return FCFC_NODE_DIST_OUT;

  /* Check the maximum separation. */
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO  &&             \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  dcmp = conf->tsmin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  dcmp = conf->smin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  dcmp = conf->pmin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  dcmp = conf->msmax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_IN;
  #endif
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
  dcmp = conf->smax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_IN;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO */
  dcmp = conf->smin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  dcmp = conf->smax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) {
    dcmp = conf->smin + sr;
    if (sum >= dcmp * dcmp) return FCFC_NODE_DIST_IN;
  }
  #endif     /* FCFC_BIN_SMIN */
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
  register real diameter = node->r * 2;
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO  &&             \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (diameter < conf->tsmin) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (diameter < conf->smin) return FCFC_NODE_DIST_OUT;
  #elif         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (diameter < conf->pmin) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  if (diameter < conf->msmax) return FCFC_NODE_DIST_IN;
  #endif
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (diameter < conf->smin) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (diameter < conf->smax) return FCFC_NODE_DIST_IN;
  #endif
#endif       /* FCFC_BIN_TYPE */
  return FCFC_NODE_DIST_X;
}

#endif
