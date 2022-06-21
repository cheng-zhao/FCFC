/*******************************************************************************
* 2pt_box/metric_kdtree.c: this file is part of the FCFC program.

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
  dx = REAL_ABS(dx + shift[0]);

  register real dy = node1->cen[1] - node2->cen[1];
  if (dy > conf->hsize[1]) shift[1] = -conf->bsize[1];
  else if (dy < -conf->hsize[1]) shift[1] = conf->bsize[1];
  dy = REAL_ABS(dy + shift[1]);

  register real dz = node1->cen[2] - node2->cen[2];
  if (dz > conf->hsize[2]) shift[2] = -conf->bsize[2];
  else if (dz < -conf->hsize[2]) shift[2] = conf->bsize[2];
  dz = REAL_ABS(dz + shift[2]);
#endif
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  /* Check the minimum pi. */
  const real wz = node1->wid[2] + node2->wid[2];
  if (dz - wz >= conf->pmax) return FCFC_NODE_DIST_OUT;
  /* Check the maximum pi. */
  register real zmax = dz + wz;
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (zmax < conf->pmin) return FCFC_NODE_DIST_OUT;
  #endif

  /* Check the minimum s_perp. */
  const real wx = node1->wid[0] + node2->wid[0];
  const real wy = node1->wid[1] + node2->wid[1];
  register real sum = 0;
  register real d = dx - wx;
  if (d > 0) sum += d * d;
  d = dy - wy;
  if (d > 0) sum += d * d;
  if (sum >= conf->s2max) return FCFC_NODE_DIST_OUT;
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const real sum_min = sum;
  #endif

  /* Check the maximum s_perp. */
  sum = 0;
  register real xmax = dx + wx;
  sum += xmax * xmax;
  register real ymax = dy + wy;
  sum += ymax * ymax;
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (sum < conf->s2max && dz - wz >= conf->pmin && zmax < conf->pmax)
    return FCFC_NODE_DIST_IN;
    #else    /* FCFC_BIN_PMIN  == FCFC_BIN_MIN_ZERO */
  if (sum < conf->s2max && zmax < conf->pmax) return FCFC_NODE_DIST_IN;
    #endif
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO */
  if (sum < conf->s2min) return FCFC_NODE_DIST_OUT;
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (sum_min >= conf->s2min && sum < conf->s2max &&
      dz - wz >= conf->pmin && zmax < conf->pmax) return FCFC_NODE_DIST_IN;
    #else    /* FCFC_BIN_PMIN  == FCFC_BIN_MIN_ZERO */
  if (sum_min >= conf->s2min && sum < conf->s2max &&
      zmax < conf->pmax) return FCFC_NODE_DIST_IN;
    #endif
  #endif

  /* Compare maximum separations with box sizes,
   * to see if individual periodic wrap checks are necessary. */
  if (xmax >= conf->hsize[0] || ymax >= conf->hsize[1] ||
      zmax >= conf->hsize[2]) return FCFC_NODE_DIST_X_WRAP;
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  #if   FCFC_SIMD  ==  FCFC_SIMD_NONE
  /* Check the minimum separation. */
  register real wx = node1->wid[0] + node2->wid[0];
  register real wy = node1->wid[1] + node2->wid[1];
  register real wz = node1->wid[2] + node2->wid[2];

  register real sum = 0;
  register real d = dx - wx;
  if (d > 0) sum += d * d;
  d = dy - wy;
  if (d > 0) sum += d * d;
  d = dz - wz;
  if (d > 0) sum += d * d;

  if (sum >= conf->s2max) return FCFC_NODE_DIST_OUT;
    #if   FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const real sum_min = sum;
    #endif

  /* Check the maximum separation. */
  sum = 0;
  wx += dx;             /* maximum separation in x direction */
  sum += wx * wx;
  wy += dy;             /* maximum separation in y direction */
  sum += wy * wy;
  wz += dz;             /* maximum separation in z direction */
  sum += wz * wz;
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
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
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
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  x1 = SIMD_MSK0_MOV_REAL(mask, SIMD_NEG_REAL(vbsize));
  mask = SIMD_CMP_REAL(vd, SIMD_NEG_REAL(vhsize), _CMP_LT_OQ);
  x1 = SIMD_MSK_MOV_REAL(x1, mask, vbsize);
    #endif
  SIMD_MSKSTORE_REAL(shift, mask3, x1);
  vd = SIMD_ABS_REAL(SIMD_ADD_REAL(vd, x1));
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  x1 = SIMD_MSKLOAD_REAL(node1->wid, mask3);
  x2 = SIMD_MSKLOAD_REAL(node2->wid, mask3);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  x1 = SIMD_MSK0_LOADU_REAL(mask3, node1->wid);
  x2 = SIMD_MSK0_LOADU_REAL(mask3, node2->wid);
    #endif
  vec_r vw = SIMD_ADD_REAL(x1, x2);

  /* Check the minimum separation. */
  vec_r delta = SIMD_SUB_REAL(vd, vw);
  mask = SIMD_CMP_REAL(delta, SIMD_ZEROS_REAL(), _CMP_GT_OQ);
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  delta = SIMD_AND_REAL(delta, mask);
  delta = SIMD_MUL_REAL(delta, delta);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  delta = SIMD_MSK0_MUL_REAL(mask, delta, delta);
    #endif
  real *res = (real *) &delta;
  register real sum = res[0] + res[1] + res[2];
  if (sum >= conf->s2max) return FCFC_NODE_DIST_OUT;
    #if   FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const real sum_min = sum;
    #endif

  /* Check the maximum separation. */
  vw = SIMD_ADD_REAL(vd, vw);
  delta = SIMD_MUL_REAL(vw, vw);
  res = (real *) &delta;
  sum = res[0] + res[1] + res[2];
  #endif     /* FCFC_SIMD */

  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
  if (sum < conf->s2max) return FCFC_NODE_DIST_IN;
  #else
  if (sum < conf->s2min) return FCFC_NODE_DIST_OUT;
  if (sum_min >= conf->s2min && sum < conf->s2max) return FCFC_NODE_DIST_IN;
  #endif

  /* Compare maximum separations with box sizes,
   * to see if individual periodic wrap checks are necessary. */
  #if           FCFC_SIMD  ==  FCFC_SIMD_NONE
  if (wx >= conf->hsize[0] || wy >= conf->hsize[1] || wz >= conf->hsize[2])
    return FCFC_NODE_DIST_X_WRAP;
  #else      /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
  mask = SIMD_CMP_REAL(vw, vhsize, _CMP_GE_OQ);
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
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO  ||             \
               (FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&                \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO)
  register real wz = node->wid[2] * 2;
  #endif
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (wz < conf->pmin) return FCFC_NODE_DIST_OUT;
  #endif

  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO ||              \
               (FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&                \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO)
  register real wx = node->wid[0] * 2;
  register real wy = node->wid[1] * 2;
  register real dist = wx * wx + wy * wy;
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist < conf->s2min) return FCFC_NODE_DIST_OUT;
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO  &&                \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO
  if (dist < conf->s2max && wz < conf->pmax) return FCFC_NODE_DIST_IN;
  #endif
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  register real wx = node->wid[0] * 2;
  register real wy = node->wid[1] * 2;
  register real wz = node->wid[2] * 2;
  register real dist = wx * wx + wy * wy + wz * wz;
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist < conf->s2min) return FCFC_NODE_DIST_OUT;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist < conf->s2max) return FCFC_NODE_DIST_IN;
  #endif
#endif       /* FCFC_BIN_TYPE */
  return FCFC_NODE_DIST_X;
}

#endif
