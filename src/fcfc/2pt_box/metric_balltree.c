/*******************************************************************************
* 2pt_box/metric_balltree.c: this file is part of the FCFC program.

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
  * `conf`:     the structure for pair count configurations;
  * `shift`:    coordinate offsets for periodic boundary conditions.
Return:
  Indicate whether the separation range is inside, is outside, or intersects
  with the distance range of interest.
******************************************************************************/
static inline fcfc_node_dist_e PRIVATE_NAME(compare_dual_node) (
    const FCFC_TREE_DTYPE *node1, const FCFC_TREE_DTYPE *node2,
    const PRIVATE_NAME(conf) *conf, real shift[3]) {
#if     FCFC_SIMD  ==  FCFC_SIMD_NONE  ||  FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  register real dx = node1->cen[0] - node2->cen[0];
  if (dx > conf->hsize[0]) shift[0] = -conf->bsize[0];
  else if (dx < -conf->hsize[0]) shift[0] = conf->bsize[0];
  dx += shift[0];

  register real dy = node1->cen[1] - node2->cen[1];
  if (dy > conf->hsize[1]) shift[1] = -conf->bsize[1];
  else if (dy < -conf->hsize[1]) shift[1] = conf->bsize[1];
  dy += shift[1];

  register real dz = node1->cen[2] - node2->cen[2];
  if (dz > conf->hsize[2]) shift[2] = -conf->bsize[2];
  else if (dz < -conf->hsize[2]) shift[2] = conf->bsize[2];
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  dz = REAL_ABS(dz + shift[2]);
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  dz += shift[2];
  #endif
#endif
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  /* Check the range of pi. */
  register real sr = node1->r + node2->r;
  if (dz - sr >= conf->pmax) return FCFC_NODE_DIST_OUT;
  register real zmax = dz + sr;
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (zmax < conf->pmin) return FCFC_NODE_DIST_OUT;
  #endif

  /* Check the range of s_perp. */
  register real sum = dx * dx + dy * dy;
  register real dcmp = conf->smax + sr;
  if (sum >= dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
  dcmp = conf->smax - sr;
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dcmp > 0 && sum < dcmp * dcmp &&
      dz - sr >= conf->pmin && dz + sr < conf->pmax) return FCFC_NODE_DIST_IN;
    #else    /* FCFC_BIN_PMIN  == FCFC_BIN_MIN_ZERO */
  if (dcmp > 0 && sum < dcmp * dcmp && dz + sr < conf->pmax) return FCFC_NODE_DIST_IN;
    #endif
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO */
  dcmp = conf->smin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  dcmp = conf->smax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) {
    dcmp = conf->smin + sr;
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
    if (sum >= dcmp * dcmp && dz - sr >= conf->pmin && dz + sr < conf->pmax)
      return FCFC_NODE_DIST_IN;
    #else    /* FCFC_BIN_PMIN  == FCFC_BIN_MIN_ZERO */
    if (sum >= dcmp * dcmp && dz + sr < conf->pmax) return FCFC_NODE_DIST_IN;
    #endif
  }
  #endif

  /* Compare maximum separations with box sizes,
   * to see if individual periodic wrap checks are necessary. */
  if (REAL_ABS(dx) + sr >= conf->hsize[0] ||
      REAL_ABS(dy) + sr >= conf->hsize[1] ||
      zmax >= conf->hsize[2]) return FCFC_NODE_DIST_X_WRAP;
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  #if           FCFC_SIMD  ==  FCFC_SIMD_NONE
  const real sum = dx * dx + dy * dy + dz * dz;
  #else      /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
  /* Load only 3 elements to the vector. */
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
      #ifdef SINGLE_PREC
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0,
      0, 0, 0, 0);
      #else       /* double precision */
  const vec_i mask3 = SIMD_SETR_INT(0xFFFFFFFFFFFFFFFFULL,
      0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0);
      #endif
  vec_r x1 = SIMD_MSKLOAD_REAL(node1->cen, mask3);
  vec_r x2 = SIMD_MSKLOAD_REAL(node2->cen, mask3);
  const vec_r vbsize = SIMD_MSKLOAD_REAL(conf->bsize, mask3);
  const vec_r vhsize = SIMD_MSKLOAD_REAL(conf->hsize, mask3);
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  const vec_msk mask3 = 7;              /* 0b0111 */
  vec_r x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->cen);
  vec_r x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->cen);
  const vec_r vbsize = SIMD_MSK0_LOADU_REAL(mask3, conf->bsize);
  const vec_r vhsize = SIMD_MSK0_LOADU_REAL(mask3, conf->hsize);
    #endif
  vec_r vd = SIMD_SUB_REAL(x1, x2);
  /* Periodic wrapping. */
  vec_msk mask = SIMD_CMP_REAL(vd, vhsize, _CMP_GT_OQ);
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  x1 = SIMD_AND_REAL(mask, SIMD_NEG_REAL(vbsize));
  mask = SIMD_CMP_REAL(vd, SIMD_NEG_REAL(vhsize), _CMP_LT_OQ);
  x2 = SIMD_AND_REAL(mask, vbsize);
  x1 = SIMD_OR_REAL(x1, x2);
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  x1 = SIMD_MSK0_MOV_REAL(mask, SIMD_NEG_REAL(vbsize));
  mask = SIMD_CMP_REAL(vd, SIMD_NEG_REAL(vhsize), _CMP_LT_OQ);
  x1 = SIMD_MSK_MOV_REAL(x1, mask, vbsize);
    #endif
  SIMD_MSKSTORE_REAL(shift, mask3, x1);
  x1 = SIMD_ADD_REAL(vd, x1);           /* for (dx,dy,dz) */

  /* Distance between spherical centres. */
  vd = SIMD_MUL_REAL(x1, x1);
  real *res = (real *) &vd;
  const real sum = res[0] + res[1] + res[2];
  #endif     /* FCFC_SIMD */

  register real sr = node1->r + node2->r;
  register real dcmp = conf->smax + sr;
  if (sum >= dcmp * dcmp) return FCFC_NODE_DIST_OUT;

  #if   FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
  dcmp = conf->smax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_IN;
  #else
  dcmp = conf->smin - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) return FCFC_NODE_DIST_OUT;
  dcmp = conf->smax - sr;
  if (dcmp > 0 && sum < dcmp * dcmp) {
    dcmp = conf->smin + sr;
    if (sum >= dcmp * dcmp) return FCFC_NODE_DIST_IN;
  }
  #endif

  /* Compare maximum separations with box sizes,
   * to see if individual periodic wrap checks are necessary. */
  #if   FCFC_SIMD  ==  FCFC_SIMD_NONE
  if (REAL_ABS(dx) + sr >= conf->hsize[0] ||
      REAL_ABS(dy) + sr >= conf->hsize[1] ||
      REAL_ABS(dz) + sr >= conf->hsize[2]) return FCFC_NODE_DIST_X_WRAP;
  #else
  x1 = SIMD_ADD_REAL(SIMD_ABS_REAL(x1), SIMD_SET1_REAL(sr));
  mask = SIMD_CMP_REAL(x1, vhsize, _CMP_GE_OQ);
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  if (SIMD_MVMSK_REAL(mask) & 7) return FCFC_NODE_DIST_X_WRAP;
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  if (mask & mask3) return FCFC_NODE_DIST_X_WRAP;
    #endif
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
  register real diameter = node->r * 2;
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (diameter < conf->pmin) return FCFC_NODE_DIST_OUT;
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (diameter < conf->smin) return FCFC_NODE_DIST_OUT;
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&                \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO
  if (diameter < conf->smax && diameter < conf->pmax) return FCFC_NODE_DIST_IN;
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
