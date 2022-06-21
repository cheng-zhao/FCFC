/*******************************************************************************
* 2pt_box/metric_common.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#if defined(FCFC_TREE_TYPE) && defined(FCFC_CNT_TYPE) &&                \
  defined(FCFC_BIN_TYPE) && defined(FCFC_TAB_TYPE) &&                   \
  defined(FCFC_STAB_WIDTH) && defined(FCFC_PTAB_WIDTH) &&               \
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
  #error please include `metric_common.c` in `dual_tree.c`
#endif

/*============================================================================*\
                Shortcuts for table lookup and histogram update
\*============================================================================*/

#ifdef FCFC_LOOKUP_HYBRID
  #undef FCFC_LOOKUP_HYBRID
#endif
#ifdef FCFC_UPDATE_HIST
  #undef FCFC_UPDATE_HIST
#endif
#ifdef FCFC_UPDATE_PAIRCNT
  #undef FCFC_UPDATE_PAIRCNT
#endif

/******************************************************************************
Macro `FCFC_LOOKUP_HYBRID`:
  Lookup the (squared) distance bin from a hybrid table.
Arguments:
  * `idx`:      the distance bin index read from the table;
  * `ntab`:     number of entries in the lookup table;
  * `dist`:     the (squared) distance to be examined;
  * `dbin`:     edges of distance bins.
******************************************************************************/
#if             FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_INT
  #define       FCFC_LOOKUP_HYBRID(idx,ntab,dist,dbin)
#else        /* FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID */
  #define       FCFC_LOOKUP_HYBRID(idx,ntab,dist,dbin)                  \
    if (idx >= (ntab)) {                                                \
      idx -= (ntab);                                                    \
      while (idx != 0 && (dist) < (dbin)[idx]) idx--;                   \
    }
#endif

/******************************************************************************
Macro `FCFC_UPDATE_HIST`:
  Add a value to a histogram bin.
Arguments:
  * `idx`:      the histogram bin to be updated;
  * `val`:      value to be added to the histogram bin.
******************************************************************************/
#if             FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  #define       FCFC_UPDATE_HIST(idx,val)                               \
    ((COUNT *) cnt)[(idx)].i += 1
#else        /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  #define       FCFC_UPDATE_HIST(idx,val)                               \
    ((COUNT *) cnt)[(idx)].d += (val)
#endif

/******************************************************************************
Macro `FCFC_UPDATE_PAIRCNT`
  Update pair counts with table lookup.
Arguments:
  * `sidx`:     index of squared distance (or s_perp) in the lookup table;
  * `stab`:     the lookup table for squared distance (or s_perp);
  * `nstab`:    number of entries in the squared s (or s_perp) lookup table;
  * `dist`:     the (squared) distance to be examined;
  * `sbin`:     edges of squared distance (or s_perp) bins;
  * `midx`:     index of squared mu in the lookup table;
  * `mtab`:     the lookup table for squared mu;
  * `pidx`:     index of pi (or squared mu) in the lookup table;
  * `ptab`:     the lookup table for pi;
  * `nptab`:    number of entries in the pi lookup table;
  * `pi`:       the pi to be examined;
  * `pbin`:     edges of pi bins;
  * `wt`:       the weight to be added to the pair count.
******************************************************************************/
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
  #define       FCFC_UPDATE_PAIRCNT(sidx,stab,nstab,dist,sbin,          \
                    midx,mtab,pidx,ptab,nptab,pi,pbin,wt)               \
    int _sid = (stab)[(sidx)];                                          \
    FCFC_LOOKUP_HYBRID(_sid, nstab, dist, sbin);                        \
    FCFC_UPDATE_HIST(_sid, wt)
#elif           FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  #define       FCFC_UPDATE_PAIRCNT(sidx,stab,nstab,dist,sbin,          \
                    midx,mtab,pidx,ptab,nptab,pi,pbin,wt)               \
    int _sid = (stab)[(sidx)];                                          \
    FCFC_LOOKUP_HYBRID(_sid, nstab, dist, sbin);                        \
    int _mid = (mtab)[(midx)];                                          \
    FCFC_UPDATE_HIST(_sid + _mid * (nstab), wt)
#else        /* FCFC_BIN_TYPE  ==  FCFC_BIN_SPI */
  #define       FCFC_UPDATE_PAIRCNT(sidx,stab,nstab,dist,sbin,          \
                    midx,mtab,pidx,ptab,nptab,pi,pbin,wt)               \
    int _sid = (stab)[(sidx)];                                          \
    FCFC_LOOKUP_HYBRID(_sid, nstab, dist, sbin);                        \
    int _pid = (ptab)[(pidx)];                                          \
    FCFC_LOOKUP_HYBRID(_pid, nptab, pi, pbin);                          \
    FCFC_UPDATE_HIST(_sid + _pid * (nstab), wt)
#endif


/*============================================================================*\
            Functions for distance evaluation and pair count update
\*============================================================================*/

/******************************************************************************
Function `compute_dist_hist_scalar_<IDENTIFIERS>`:
  Scalar code for computing distances and updating pair counts.
Arguments:
  * `x1,y1,z1`: coordinates of the first point;
  * `x2,y2,z2`: coordinates of the second point;
  * `w1`:       weight of the first point;
  * `w2`:       weight of the second point;
  * `soff`:     offset for squared distance (or s_perp) lookup;
  * `poff`:     offset for pi lookup;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(compute_dist_hist_scalar) (
    const real x1, const real y1, const real z1,
    const real x2, const real y2, const real z2,
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const real w1, const real w2,
#endif
#if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    const int soff,
#endif
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
    const int poff,
#endif
    const PRIVATE_NAME(conf) *conf, void *cnt) {
  /* Compute (squared) distances. */
  register real dx = x1 - x2;
  register real dy = y1 - y2;
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  register real pi = REAL_ABS(z1 - z2);                 /* pi, a.k.a. s_par */
  /* Check the pi range. */
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (pi >= conf->pmax || pi < conf->pmin) return;
  #else      /* FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  if (pi >= conf->pmax) return;
  #endif
  register real dist = dx * dx + dy * dy;               /* s_perp squared */
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  register real dz = z1 - z2;
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  dz *= dz;
  register real dist = dx * dx + dy * dy + dz;
  #else      /* FCFC_BIN_TYPE  ==  FCFC_BIN_ISO */
  register real dist = dx * dx + dy * dy + dz * dz;
  #endif
#endif
  /* Check the separation (or s_perp) range. */
#if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist >= conf->s2max || dist < conf->s2min) return;
#else        /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist >= conf->s2max) return;
#endif

  /* Compute mu and lookup the mu bin. */
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  int pidx = (dist < REAL_EPS) ? 0 : (dz / dist) * conf->nmu2;
  #ifdef        WITH_MU_ONE
  if (pidx >= conf->nmu2) pidx = conf->nmu2 - 1;
  #else
  if (pidx >= conf->nmu2) return;
  #endif
  pidx = conf->mutab[pidx];
#endif

  /* Lookup the squared separation (or s_perp) bin. */
#if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  int sidx = (int) dist - soff;
#else        /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  int sidx = (int) dist;
#endif
  sidx = conf->stab[sidx];
  FCFC_LOOKUP_HYBRID(sidx, conf->ns, dist, conf->s2bin);

  /* Lookup the pi bin. */
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  int pidx = (int) pi - poff;
  #else      /* FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  int pidx = (int) pi;
  #endif
  pidx = conf->ptab[pidx];
  FCFC_LOOKUP_HYBRID(pidx, conf->np, pi, conf->pbin);
#endif

  /* Increment pair count. */
#if             FCFC_SIMD  <  FCFC_SIMD_AVX512
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
  FCFC_UPDATE_HIST(sidx, w1 * w2);
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_ISO */
  FCFC_UPDATE_HIST(sidx + pidx * conf->ns, w1 * w2);
  #endif
#else        /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  ((int64_t *) cnt)[sidx] += 1;
    #else    /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  ((double *) cnt)[sidx] += w1 * w2;
    #endif
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_ISO */
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  ((int64_t *) cnt)[sidx + pidx * conf->ns] += 1;
    #else    /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  ((double *) cnt)[sidx + pidx * conf->ns] += w1 * w2;
    #endif
  #endif
#endif
}

/******************************************************************************
Function `compute_dist_wrap_hist_scalar_<IDENTIFIERS>`:
  Scalar code for computing distances and updating pair counts,
  with periodic wrap for every point.
Arguments:
  * `x1,y1,z1`: coordinates of the first point;
  * `x2,y2,z2`: coordinates of the second point;
  * `w1`:       weight of the first point;
  * `w2`:       weight of the second point;
  * `soff`:     offset for squared distance (or s_perp) lookup;
  * `poff`:     offset for pi lookup;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(compute_dist_wrap_hist_scalar) (
    const real x1, const real y1, const real z1,
    const real x2, const real y2, const real z2,
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const real w1, const real w2,
#endif
#if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    const int soff,
#endif
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
    const int poff,
#endif
    const PRIVATE_NAME(conf) *conf, void *cnt) {
  /* Compute (squared) distances. */
  register real dx = x1 - x2;
  if (dx > conf->hsize[0]) dx -= conf->bsize[0];
  else if (dx < -conf->hsize[0]) dx += conf->bsize[0];
  register real dy = y1 - y2;
  if (dy > conf->hsize[1]) dy -= conf->bsize[1];
  else if (dy < -conf->hsize[1]) dy += conf->bsize[1];
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  register real pi = z1 - z2;
  if (pi > conf->hsize[2]) pi -= conf->bsize[2];
  else if (pi < -conf->hsize[2]) pi += conf->bsize[2];
  pi = REAL_ABS(pi);                                    /* pi, a.k.a. s_par */
  /* Check the pi range. */
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (pi >= conf->pmax || pi < conf->pmin) return;
  #else      /* FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  if (pi >= conf->pmax) return;
  #endif
  register real dist = dx * dx + dy * dy;               /* s_perp squared */
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  register real dz = z1 - z2;
  if (dz > conf->hsize[2]) dz -= conf->bsize[2];
  else if (dz < -conf->hsize[2]) dz += conf->bsize[2];
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  dz *= dz;
  register real dist = dx * dx + dy * dy + dz;
  #else      /* FCFC_BIN_TYPE  ==  FCFC_BIN_ISO */
  register real dist = dx * dx + dy * dy + dz * dz;
  #endif
#endif
  /* Check the separation (or s_perp) range. */
#if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist >= conf->s2max || dist < conf->s2min) return;
#else        /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist >= conf->s2max) return;
#endif

  /* Compute mu and lookup the mu bin. */
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  int pidx = (dist < REAL_EPS) ? 0 : (dz / dist) * conf->nmu2;
  #ifdef        WITH_MU_ONE
  if (pidx >= conf->nmu2) pidx = conf->nmu2 - 1;
  #else
  if (pidx >= conf->nmu2) return;
  #endif
  pidx = conf->mutab[pidx];
#endif

  /* Lookup the squared separation (or s_perp) bin. */
#if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  int sidx = (int) dist - soff;
#else        /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  int sidx = (int) dist;
#endif
  sidx = conf->stab[sidx];
  FCFC_LOOKUP_HYBRID(sidx, conf->ns, dist, conf->s2bin);

  /* Lookup the pi bin. */
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #if           FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  int pidx = (int) pi - poff;
  #else      /* FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  int pidx = (int) pi;
  #endif
  pidx = conf->ptab[pidx];
  FCFC_LOOKUP_HYBRID(pidx, conf->np, pi, conf->pbin);
#endif

  /* Increment pair count. */
#if             FCFC_SIMD  <  FCFC_SIMD_AVX512
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
  FCFC_UPDATE_HIST(sidx, w1 * w2);
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_ISO */
  FCFC_UPDATE_HIST(sidx + pidx * conf->ns, w1 * w2);
  #endif
#else        /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  ((int64_t *) cnt)[sidx] += 1;
    #else    /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  ((double *) cnt)[sidx] += w1 * w2;
    #endif
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_ISO */
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  ((int64_t *) cnt)[sidx + pidx * conf->ns] += 1;
    #else    /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  ((double *) cnt)[sidx + pidx * conf->ns] += w1 * w2;
    #endif
  #endif
#endif
}

#if             FCFC_SIMD  !=  FCFC_SIMD_NONE
/******************************************************************************
Function `compute_dist_vector_<IDENTIFIERS>`:
  Vector code for computing distances and check them with the ranges.
Arguments:
  * `x1,y1,z1`: coordinates of the first points;
  * `x2,y2,z2`: coordinates of the second points;
  * `vs2min`:   minimum squared distance (or s_perp) of interest;
  * `vs2max`:   maximum squared distance (or s_perp) of interest;
  * `vsoff`:    offset for squared distance (or s_perp) lookup;
  * `vpmin`:    minimum pi of interest;
  * `vpmax`:    maximum pi of interest;
  * `vpoff`:    offset for pi lookup;
  * `vnmu2`:    the squared number of mu bins;
  * `vs2`:      the evaluated squared distance (or s_perp);
  * `vp`:       the evaluated pi;
  * `vsid`:     indices of squared distance bins in the lookup table;
  * `vpid`:     indices of pi (or squared mu) bins in the lookup table;
  * `msk`:      mask indicating whether the distances are of interest.
******************************************************************************/
static inline void PRIVATE_NAME(compute_dist_vector) (
    const vec_r x1, const vec_r y1, const vec_r z1,
    vec_r x2, vec_r y2, vec_r z2,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    const vec_r vs2min,
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
    const vec_r vsoff,
    #else
    const vec_r2i vsoff,
    #endif
  #endif
    const vec_r vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
    const vec_r vpmin,
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
    const vec_r vpoff,
      #else
    const vec_r2i vpoff,
      #endif
    #endif
    const vec_r vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
    const vec_r vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
    vec_r *vs2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    vec_r *vp,
    #endif
  #endif
    vec_r2i *vsid,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    vec_r2i *vpid,
  #endif
    vec_msk *msk) {
  /* Compute (squared) distances. */
  x2 = SIMD_SUB_REAL(x1, x2);
  y2 = SIMD_SUB_REAL(y1, y2);
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  z2 = SIMD_ABS_REAL(SIMD_SUB_REAL(z1, z2));            /* pi, a.k.a. s_par */
  x2 = SIMD_MUL_REAL(x2, x2);
    #ifdef      FCFC_SIMD_FMA
  x2 = SIMD_FMADD_REAL(y2, y2, x2);                     /* s_perp squared */
    #else
  y2 = SIMD_MUL_REAL(y2, y2);
  x2 = SIMD_ADD_REAL(x2, y2);                           /* s_perp squared */
    #endif
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  z2 = SIMD_SUB_REAL(z1, z2);
  z2 = SIMD_MUL_REAL(z2, z2);                           /* pi squared */
    #ifdef      FCFC_SIMD_FMA
  x2 = SIMD_FMADD_REAL(x2, x2, z2);
  x2 = SIMD_FMADD_REAL(y2, y2, x2);                     /* s squared */
    #else
  x2 = SIMD_MUL_REAL(x2, x2);
  x2 = SIMD_ADD_REAL(x2, z2);
  y2 = SIMD_MUL_REAL(y2, y2);
  x2 = SIMD_ADD_REAL(x2, y2);                           /* s squared */
    #endif
  #endif

  /* Check distance ranges. */
  vec_msk mask = SIMD_CMP_REAL(x2, vs2max, _CMP_LT_OQ);
  #if          (FCFC_SIMD  <  FCFC_SIMD_AVX512)  &&                     \
               (FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  ||                     \
                FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO)
  vec_msk mask1;
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  mask1 = SIMD_CMP_REAL(x2, vs2min, _CMP_GE_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  mask = SIMD_MSK_CMP_REAL(mask, x2, vs2min, _CMP_GE_OQ);
    #endif
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  mask1 = SIMD_CMP_REAL(z2, vpmax, _CMP_LT_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  mask1 = SIMD_CMP_REAL(z2, vpmin, _CMP_GE_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
      #endif
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  mask = SIMD_MSK_CMP_REAL(mask, z2, vpmax, _CMP_LT_OQ);
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  mask = SIMD_MSK_CMP_REAL(mask, z2, vpmin, _CMP_GE_OQ);
      #endif
    #endif
  #endif

  /* Compute mu and the index of the lookup table. */
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  vec_msk mask2 = SIMD_CMP_REAL(x2, SIMD_SET1_REAL(REAL_EPS), _CMP_GE_OQ);
  z2 = SIMD_MUL_REAL(z2, vnmu2);
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  z2 = SIMD_DIV_REAL(z2, x2);
  z2 = SIMD_AND_REAL(z2, mask2);
  z2 = SIMD_FLOOR_REAL(z2);
  mask2 = SIMD_CMP_REAL(z2, vnmu2, _CMP_LT_OQ);
      #ifdef    WITH_MU_ONE
  z2 = SIMD_SUB_REAL(z2, SIMD_ANDNOT_REAL(mask2, SIMD_SET1_REAL(1.0)));
      #else
  mask = SIMD_AND_REAL(mask, mask2);
      #endif
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  z2 = SIMD_MSK0_DIV_ROUND_REAL(mask2, z2, x2,
      _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
      #ifdef    WITH_MU_ONE
  mask2 = SIMD_MSK_CMP_REAL(mask, z2, vnmu2, _CMP_GE_OQ);
  z2 = SIMD_MSK_SUB_REAL(z2, mask2, z2, SIMD_SET1_REAL(1.0));
      #else
  mask = SIMD_MSK_CMP_REAL(mask, z2, vnmu2, _CMP_LT_OQ);
      #endif
    #endif
  vec_r2i vpidx = SIMD_R2I(z2);
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  /* Find the index of the squared distance (or s_perp) lookup table. */
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  vec_r x2 = SIMD_SUB_REAL(SIMD_FLOOR_REAL(x2), vsoff);
  vec_r2i vsidx = SIMD_R2I(x2);
    #else
  vec_r2i vsidx = SIMD_R2I(x2);
  vsidx = SIMD_SUB_R2I(vsidx, vsoff);
    #endif
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  vec_r2i vsidx = SIMD_R2I(x2);
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  /* Find the index of the pi lookup table. */
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  vec_r z2 = SIMD_SUB_REAL(SIMD_FLOOR_REAL(z2), vpoff);
  vec_r2i vpidx = SIMD_R2I(z2);
      #else
  vec_r2i vpidx = SIMD_R2I(z2);
  vpidx = SIMD_SUB_R2I(vpidx, vpoff);
      #endif
    #else    /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  vec_r2i vpidx = SIMD_R2I(z2);
    #endif
  #endif

  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  *vs2 = x2;
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  *vp = z2;
    #endif
  #endif
  *vsid = vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  *vpid = vpidx;
  #endif
  *msk = mask;
}

/******************************************************************************
Function `compute_dist_wrap_vector_<IDENTIFIERS>`:
  Vector code for computing distances and check them with the ranges,
  with periodic wrap for every point.
Arguments:
  * `x1,y1,z1`: coordinates of the first points;
  * `x2,y2,z2`: coordinates of the second points;
  * `hx,hy,hz`: half box size on each direction;
  * `bx,by,bz`: box size on each direction;
  * `vs2min`:   minimum squared distance (or s_perp) of interest;
  * `vs2max`:   maximum squared distance (or s_perp) of interest;
  * `vsoff`:    offset for squared distance (or s_perp) lookup;
  * `vpmin`:    minimum pi of interest;
  * `vpmax`:    maximum pi of interest;
  * `vpoff`:    offset for pi lookup;
  * `vnmu2`:    the squared number of mu bins;
  * `vs2`:      the evaluated squared distance (or s_perp);
  * `vp`:       the evaluated pi;
  * `vsid`:     indices of squared distance bins in the lookup table;
  * `vpid`:     indices of pi (or squared mu) bins in the lookup table;
  * `msk`:      mask indicating whether the distances are of interest.
******************************************************************************/
static inline void PRIVATE_NAME(compute_dist_wrap_vector) (
    const vec_r x1, const vec_r y1, const vec_r z1,
    vec_r x2, vec_r y2, vec_r z2,
    const vec_r hx, const vec_r hy, const vec_r hz,
    const vec_r bx, const vec_r by, const vec_r bz,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    const vec_r vs2min,
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
    const vec_r vsoff,
    #else
    const vec_r2i vsoff,
    #endif
  #endif
    const vec_r vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
    const vec_r vpmin,
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
    const vec_r vpoff,
      #else
    const vec_r2i vpoff,
      #endif
    #endif
    const vec_r vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
    const vec_r vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
    vec_r *vs2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    vec_r *vp,
    #endif
  #endif
    vec_r2i *vsid,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    vec_r2i *vpid,
  #endif
    vec_msk *msk) {
  /* Compute (squared) distances. */
  x2 = SIMD_SUB_REAL(x1, x2);
  vec_msk mask = SIMD_CMP_REAL(x2, hx, _CMP_GT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  x2 = SIMD_SUB_REAL(x2, SIMD_AND_REAL(mask, bx));
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  x2 = SIMD_MSK_SUB_REAL(x2, mask, x2, bx);
  #endif
  mask = SIMD_CMP_REAL(x2, SIMD_NEG_REAL(hx), _CMP_LT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  x2 = SIMD_ADD_REAL(x2, SIMD_AND_REAL(mask, bx));
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  x2 = SIMD_MSK_ADD_REAL(x2, mask, x2, bx);
  #endif

  y2 = SIMD_SUB_REAL(y1, y2);
  mask = SIMD_CMP_REAL(y2, hy, _CMP_GT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  y2 = SIMD_SUB_REAL(y2, SIMD_AND_REAL(mask, by));
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  y2 = SIMD_MSK_SUB_REAL(y2, mask, y2, by);
  #endif
  mask = SIMD_CMP_REAL(y2, SIMD_NEG_REAL(hy), _CMP_LT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  y2 = SIMD_ADD_REAL(y2, SIMD_AND_REAL(mask, by));
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  y2 = SIMD_MSK_ADD_REAL(y2, mask, y2, by);
  #endif

  z2 = SIMD_SUB_REAL(z1, z2);
  mask = SIMD_CMP_REAL(z2, hz, _CMP_GT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  z2 = SIMD_SUB_REAL(z2, SIMD_AND_REAL(mask, bz));
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  z2 = SIMD_MSK_SUB_REAL(z2, mask, z2, bz);
  #endif
  mask = SIMD_CMP_REAL(z2, SIMD_NEG_REAL(hz), _CMP_LT_OQ);
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  z2 = SIMD_ADD_REAL(z2, SIMD_AND_REAL(mask, bz));
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  z2 = SIMD_MSK_ADD_REAL(z2, mask, z2, bz);
  #endif

  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  z2 = SIMD_ABS_REAL(z2);                               /* pi, a.k.a. s_par */
  x2 = SIMD_MUL_REAL(x2, x2);
    #ifdef      FCFC_SIMD_FMA
  x2 = SIMD_FMADD_REAL(y2, y2, x2);                     /* s_perp squared */
    #else
  y2 = SIMD_MUL_REAL(y2, y2);
  x2 = SIMD_ADD_REAL(x2, y2);                           /* s_perp squared */
    #endif
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
  z2 = SIMD_MUL_REAL(z2, z2);                           /* pi squared */
    #ifdef      FCFC_SIMD_FMA
  x2 = SIMD_FMADD_REAL(x2, x2, z2);
  x2 = SIMD_FMADD_REAL(y2, y2, x2);                     /* s squared */
    #else
  x2 = SIMD_MUL_REAL(x2, x2);
  x2 = SIMD_ADD_REAL(x2, z2);
  y2 = SIMD_MUL_REAL(y2, y2);
  x2 = SIMD_ADD_REAL(x2, y2);                           /* s squared */
    #endif
  #endif

  /* Check distance ranges. */
  mask = SIMD_CMP_REAL(x2, vs2max, _CMP_LT_OQ);
  #if          (FCFC_SIMD  <  FCFC_SIMD_AVX512)  &&                     \
               (FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  ||                     \
                FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO)
  vec_msk mask1;
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  mask1 = SIMD_CMP_REAL(x2, vs2min, _CMP_GE_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  mask = SIMD_MSK_CMP_REAL(mask, x2, vs2min, _CMP_GE_OQ);
    #endif
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  mask1 = SIMD_CMP_REAL(z2, vpmax, _CMP_LT_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  mask1 = SIMD_CMP_REAL(z2, vpmin, _CMP_GE_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
      #endif
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  mask = SIMD_MSK_CMP_REAL(mask, z2, vpmax, _CMP_LT_OQ);
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  mask = SIMD_MSK_CMP_REAL(mask, z2, vpmin, _CMP_GE_OQ);
      #endif
    #endif
  #endif

  /* Compute mu and the index of the lookup table. */
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  vec_msk mask2 = SIMD_CMP_REAL(x2, SIMD_SET1_REAL(REAL_EPS), _CMP_GE_OQ);
  z2 = SIMD_MUL_REAL(z2, vnmu2);
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  z2 = SIMD_DIV_REAL(z2, x2);
  z2 = SIMD_AND_REAL(z2, mask2);
  z2 = SIMD_FLOOR_REAL(z2);
  mask2 = SIMD_CMP_REAL(z2, vnmu2, _CMP_LT_OQ);
      #ifdef    WITH_MU_ONE
  z2 = SIMD_SUB_REAL(z2, SIMD_ANDNOT_REAL(mask2, SIMD_SET1_REAL(1.0)));
      #else
  mask = SIMD_AND_REAL(mask, mask2);
      #endif
    #else    /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
  z2 = SIMD_MSK0_DIV_ROUND_REAL(mask2, z2, x2,
      _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
      #ifdef    WITH_MU_ONE
  mask2 = SIMD_MSK_CMP_REAL(mask, z2, vnmu2, _CMP_GE_OQ);
  z2 = SIMD_MSK_SUB_REAL(z2, mask2, z2, SIMD_SET1_REAL(1.0));
      #else
  mask = SIMD_MSK_CMP_REAL(mask, z2, vnmu2, _CMP_LT_OQ);
      #endif
    #endif
  vec_r2i vpidx = SIMD_R2I(z2);
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  /* Find the index of the squared distance (or s_perp) lookup table. */
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  vec_r x2 = SIMD_SUB_REAL(SIMD_FLOOR_REAL(x2), vsoff);
  vec_r2i vsidx = SIMD_R2I(x2);
    #else
  vec_r2i vsidx = SIMD_R2I(x2);
  vsidx = SIMD_SUB_R2I(vsidx, vsoff);
    #endif
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  vec_r2i vsidx = SIMD_R2I(x2);
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  /* Find the index of the pi lookup table. */
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  vec_r z2 = SIMD_SUB_REAL(SIMD_FLOOR_REAL(z2), vpoff);
  vec_r2i vpidx = SIMD_R2I(z2);
      #else
  vec_r2i vpidx = SIMD_R2I(z2);
  vpidx = SIMD_SUB_R2I(vpidx, vpoff);
      #endif
    #else    /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  vec_r2i vpidx = SIMD_R2I(z2);
    #endif
  #endif

  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  *vs2 = x2;
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  *vp = z2;
    #endif
  #endif
  *vsid = vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  *vpid = vpidx;
  #endif
  *msk = mask;
}

/******************************************************************************
Function `update_hist_vector_<IDENTIFIERS>`:
  Vector code for update pair counts.
Arguments:
  * `vsidx`:    indices of squared distance bins in the lookup table;
  * `vpidx`:    indices of pi (or squared mu) in the lookup table;
  * `vwt`:      weights to be added to the pair counts;
  * `vdist`:    the squared distances (or s_perp);
  * `vpi`:      the pi values;
  * `vns`:      number of separation (or s_perp) bins;
  * `vnp`:      number of pi bins;
  * `ioff`:     index offsets of private histograms for avoid conflicts;
  * `mask`:     mask indicating whether the distances are of interest;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(update_hist_vector) (
    vec_r2i vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    vec_r2i vpidx,
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    vec_r vwt,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
    vec_r vdist,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    vec_r vpi,
    #endif
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
    int imask,
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    const vec_i vns,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    const vec_i vnp,
    #endif
    const vec_i ioff, vec_msk mask,
  #endif
    const PRIVATE_NAME(conf) *conf, void *cnt) {
  /* Scalar code for histogram update. */
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
  int32_t *sidx = (int32_t *) &vsidx;   /* squared s (or s_perp) bin index */
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  int32_t *pidx = (int32_t *) &vpidx;   /* pi (or squared mu) bin index */
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
  real *wt = (real *) &vwt;
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  real *dist = (real *) &vdist;
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  real *pi = (real *) &vpi;
      #endif
    #endif
  if (imask & 0x1) {
    FCFC_UPDATE_PAIRCNT(sidx[0], conf->stab, conf->ns, dist[0],
        conf->s2bin, pidx[0], conf->mutab, pidx[0], conf->ptab,
        conf->np, pi[0], conf->pbin, wt[0]);
  }
  if (imask & 0x2) {
    FCFC_UPDATE_PAIRCNT(sidx[1], conf->stab, conf->ns, dist[1],
        conf->s2bin, pidx[1], conf->mutab, pidx[1], conf->ptab,
        conf->np, pi[1], conf->pbin, wt[1]);
  }
  if (imask & 0x4) {
    FCFC_UPDATE_PAIRCNT(sidx[2], conf->stab, conf->ns, dist[2],
        conf->s2bin, pidx[2], conf->mutab, pidx[2], conf->ptab,
        conf->np, pi[2], conf->pbin, wt[2]);
  }
  if (imask & 0x8) {
    FCFC_UPDATE_PAIRCNT(sidx[3], conf->stab, conf->ns, dist[3],
        conf->s2bin, pidx[3], conf->mutab, pidx[3], conf->ptab,
        conf->np, pi[3], conf->pbin, wt[3]);
  }
    #ifdef    SINGLE_PREC
  if (imask & 0x10) {
    FCFC_UPDATE_PAIRCNT(sidx[4], conf->stab, conf->ns, dist[4],
        conf->s2bin, pidx[4], conf->mutab, pidx[4], conf->ptab,
        conf->np, pi[4], conf->pbin, wt[4]);
  }
  if (imask & 0x20) {
    FCFC_UPDATE_PAIRCNT(sidx[5], conf->stab, conf->ns, dist[5],
        conf->s2bin, pidx[5], conf->mutab, pidx[5], conf->ptab,
        conf->np, pi[5], conf->pbin, wt[5]);
  }
  if (imask & 0x40) {
    FCFC_UPDATE_PAIRCNT(sidx[6], conf->stab, conf->ns, dist[6],
        conf->s2bin, pidx[6], conf->mutab, pidx[6], conf->ptab,
        conf->np, pi[6], conf->pbin, wt[6]);
  }
  if (imask & 0x80) {
    FCFC_UPDATE_PAIRCNT(sidx[7], conf->stab, conf->ns, dist[7],
        conf->s2bin, pidx[7], conf->mutab, pidx[7], conf->ptab,
        conf->np, pi[7], conf->pbin, wt[7]);
  }
    #endif
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  /* Gather squared distance (or s_perp) bins from the lookup table. */
  vec_i vidx = SIMD_MSK_R2I_GATHER_INT(ioff, mask, vsidx, conf->stab,
      sizeof(FCFC_STAB_DTYPE));         /* the first argument is not used */
  vidx = SIMD_AND_INT(vidx, SIMD_SET1_INT(FCFC_STAB_MASK));
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  vec_msk imask = SIMD_MSK_CMPGE_INT(mask, vidx, vns);
  if (imask != 0) {
    vidx = SIMD_MSK_SUB_INT(vidx, imask, vidx, vns);
    for (;;) {
      imask = SIMD_MSK_CMPNEQ_INT(imask, vidx, SIMD_ZEROS_INT());
      vec_r ref = SIMD_MSK_GATHER_REAL(vdist, imask, vidx, conf->s2bin,
          sizeof(real));
      imask = SIMD_MSK_CMP_REAL(imask, vdist, ref, _CMP_LT_OQ);
      if (imask == 0) break;
      vidx = SIMD_MSK_SUB_INT(vidx, imask, vidx, SIMD_SET1_INT(1));
    }
  }
    #endif
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  /* Gather squared mu bins from the lookup table. */
  vec_i vpid = SIMD_MSK_R2I_GATHER_INT(ioff, mask, vpidx, conf->mutab, 1);
  vpid = SIMD_AND_INT(vpid, SIMD_SET1_INT(0xFF));
    #elif       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  /* Gather pi bins from the lookup table. */
  vec_i vpid = SIMD_MSK_R2I_GATHER_INT(ioff, mask, vpidx, conf->ptab,
      sizeof(FCFC_PTAB_DTYPE));         /* the first argument is not used */
  vpid = SIMD_AND_INT(vpid, SIMD_SET1_INT(FCFC_PTAB_MASK));
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  imask = SIMD_MSK_CMPGE_INT(mask, vpid, vnp);
  if (imask != 0) {
    vpid = SIMD_MSK_SUB_INT(vpid, imask, vpid, vnp);
    for (;;) {
      imask = SIMD_MSK_CMPNEQ_INT(imask, vpid, SIMD_ZEROS_INT());
      vec_r ref = SIMD_MSK_GATHER_REAL(vpi, imask, vpid, conf->pbin,
          sizeof(real));
      imask = SIMD_MSK_CMP_REAL(imask, vpi, ref, _CMP_LT_OQ);
      if (imask == 0) break;
      vpid = SIMD_MSK_SUB_INT(vpid, imask, vpid, SIMD_SET1_INT(1));
    }
  }
      #endif
    #endif   /* FCFC_BIN_TYPE */
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  /* Use 32-bit multiplication when AVX512DQ is unavailable. */
  vpid = SIMD_MUL_INT(vpid, vns);
  vidx = SIMD_ADD_INT(vidx, vpid);
    #endif
  /* Shift indices for private histograms. */
  vidx = SIMD_ADD_INT(vidx, ioff);

  /* Update the pair counts. */
    #ifdef      SINGLE_PREC
  /* Process the lower half lanes. */
  __m256i hidx = _mm512_castsi512_si256(vidx);
  __mmask8 hmask = mask & 0xFF;
      #if       FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  vec_i vcnt = SIMD_MSK_GATHER_I64(vidx, hmask, hidx, cnt, 8);
  vcnt = SIMD_ADD_I64(vcnt, SIMD_SET1_I64(1));
  SIMD_MSK_SCATTER_I64(cnt, hmask, hidx, vcnt, 8);
      #else  /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  vec_r64 hwt = _mm512_cvt_roundps_pd(_mm512_castps512_ps256(vwt),
      _MM_FROUND_NO_EXC);
  vec_r64 vcnt = SIMD_MSK_GATHER_R64(hwt, hmask, hidx, cnt, 8);
  vcnt = SIMD_ADD_R64(vcnt, hwt);
  SIMD_MSK_SCATTER_R64(cnt, hmask, hidx, vcnt, 8);
      #endif
  /* Process the upper half lanes. */
  hidx = _mm512_extracti64x4_epi64(vidx, 1);
  hmask = mask >> 8;
      #if       FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  vcnt = SIMD_MSK_GATHER_I64(vidx, hmask, hidx, cnt, 8);
  vcnt = SIMD_ADD_I64(vcnt, SIMD_SET1_I64(1));
  SIMD_MSK_SCATTER_I64(cnt, hmask, hidx, vcnt, 8);
      #else  /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
        #ifdef  FCFC_SIMD_AVX512DQ
  hwt = _mm512_cvt_roundps_pd(_mm512_extractf32x8_ps(vwt, 1),
      _MM_FROUND_NO_EXC);
        #else
  hwt = _mm512_cvt_roundps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(
      _mm512_castps_pd(vwt), 1)), _MM_FROUND_NO_EXC);
        #endif
  vcnt = SIMD_MSK_GATHER_R64(hwt, hmask, hidx, cnt, 8);
  vcnt = SIMD_ADD_R64(vcnt, hwt);
  SIMD_MSK_SCATTER_R64(cnt, hmask, hidx, vcnt, 8);
      #endif
    #else    /* DOUBLE_PREC */
      #if       FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  vec_i vcnt = SIMD_MSK_GATHER_I64(vidx, mask, vidx, cnt, 8);
  vcnt = SIMD_ADD_I64(vcnt, SIMD_SET1_I64(1));
  SIMD_MSK_SCATTER_I64(cnt, mask, vidx, vcnt, 8);
      #else  /* FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT */
  vec_r vcnt = SIMD_MSK_GATHER_R64(vwt, mask, vidx, cnt, 8);
  vcnt = SIMD_ADD_REAL(vcnt, vwt);
  SIMD_MSK_SCATTER_R64(cnt, mask, vidx, vcnt, 8);
      #endif
    #endif
  #endif     /* FCFC_SIMD */
}
#endif


/*============================================================================*\
                   Function for counting pairs from two nodes
\*============================================================================*/

/******************************************************************************
Function `count_dual_node_<IDENTIFIERS>`:
  Count the number of pairs in the separation range of interest,
  for data points on two tree nodes.
Arguments:
  * `a`:        coordinates of points on the first node;
  * `b`:        coordinates of points on the second node;
  * `wa`:       weights of points on the first node;
  * `wb`:       weights of points on the second node;
  * `na`:       number of points on the first node;
  * `nb`:       number of points on the second node;
  * `shift`:    coordinate offsets for periodic boundary conditions;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(count_dual_node) (
    real *const a[static 3], real *const b[static 3],
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const real *wa, const real *wb,
#endif
    const size_t na, const size_t nb, const real shift[3],
    const PRIVATE_NAME(conf) *conf, void *cnt) {
  /* Scalar code. */
#if             FCFC_SIMD  ==  FCFC_SIMD_NONE
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const int soff = conf->s2min;
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&  \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  const int poff = conf->pmin;
  #endif
  for (size_t i = 0; i < na; i++) {
    real aa[3];
    aa[0] = a[0][i] + shift[0];
    aa[1] = a[1][i] + shift[1];
    aa[2] = a[2][i] + shift[2];
    for (size_t j = 0; j < nb; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_hist_scalar) (aa[0], aa[1], aa[2],
          b[0][j], b[1][j], b[2][j],
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          wa[i], wb[j],
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          soff,
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          poff,
  #endif
          conf, cnt);
    }
  }
  /* Vector code. */
#else        /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const vec_r vs2min = SIMD_SET1_REAL(conf->s2min);
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vsoff = SIMD_FLOOR_REAL(vs2min);  /* offsets for table lookup */
    #else
  const vec_r2i vsoff = SIMD_R2I(vs2min);
    #endif
  #endif
  const vec_r vs2max = SIMD_SET1_REAL(conf->s2max);
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  const vec_r vpmin = SIMD_SET1_REAL(conf->pmin);
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vpoff = SIMD_FLOOR_REAL(vpmin);
      #else
  const vec_r2i vpoff = SIMD_R2I(vpmin);
      #endif
    #endif
  const vec_r vpmax = SIMD_SET1_REAL(conf->pmax);
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  const vec_r vnmu2 = SIMD_SET1_REAL(conf->nmu2);
  #endif
  #if           FCFC_SIMD  >=  FCFC_SIMD_AVX512
  /* Index offsets of private histograms for avoiding conflicts. */
    #ifdef      SINGLE_PREC
  const vec_i ioff = _mm512_setr_epi32(0, conf->ntot, conf->ntot * 2,
      conf->ntot * 3, conf->ntot * 4, conf->ntot * 5, conf->ntot * 6,
      conf->ntot * 7, 0, conf->ntot, conf->ntot * 2, conf->ntot * 3,
      conf->ntot * 4, conf->ntot * 5, conf->ntot * 6, conf->ntot * 7);
    #else
  const vec_i ioff = _mm512_setr_epi64(0, conf->ntot, conf->ntot * 2,
      conf->ntot * 3, conf->ntot * 4, conf->ntot * 5, conf->ntot * 6,
      conf->ntot * 7);
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  const vec_i vns = SIMD_SET1_INT(conf->ns);
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  const vec_i vnp = SIMD_SET1_INT(conf->np);
    #endif
  #endif
  const vec_r vsx = SIMD_SET1_REAL(shift[0]);
  const vec_r vsy = SIMD_SET1_REAL(shift[1]);
  const vec_r vsz = SIMD_SET1_REAL(shift[2]);
  /* Vectorize the loop for node a. */
  const size_t n1 = FCFC_NUM_MASK & na;
  for (size_t i = 0; i < n1; i += FCFC_NUM_REAL) {
    vec_r x1 = SIMD_LOADU_REAL(a[0] + i);
    x1 = SIMD_ADD_REAL(x1, vsx);
    vec_r y1 = SIMD_LOADU_REAL(a[1] + i);
    y1 = SIMD_ADD_REAL(y1, vsy);
    vec_r z1 = SIMD_LOADU_REAL(a[2] + i);
    z1 = SIMD_ADD_REAL(z1, vsz);
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const vec_r w1 = SIMD_LOADU_REAL(wa + i);
  #endif
    for (size_t j = 0; j < nb; j++) {
      vec_r x2 = SIMD_SET1_REAL(b[0][j]);
      vec_r y2 = SIMD_SET1_REAL(b[1][j]);
      vec_r z2 = SIMD_SET1_REAL(b[2][j]);

      vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
      vec_r2i vpidx;
  #endif
      vec_msk mask;

      /* Compute indices of distance bins. */
      PRIVATE_NAME(compute_dist_vector) (x1, y1, z1, x2, y2, z2,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          vs2min, vsoff,
  #endif
          vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          vpmin, vpoff,
    #endif
          vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
          vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
          &x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
          &z2,
    #endif
  #endif
          &vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          &vpidx,
  #endif
          &mask);

  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      vec_r vwt = SIMD_SET1_REAL(wb[j]);
      vwt = SIMD_MUL_REAL(w1, vwt);
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
      int imask = SIMD_MVMSK_REAL(mask);
  #endif

      /* Lookup tables and increment pair counts. */
      PRIVATE_NAME(update_hist_vector) (vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vpidx,
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          vwt,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
          x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
          z2,
    #endif
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
          imask,
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vns,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
          vnp,
    #endif
          ioff, mask,
  #endif
          conf, cnt);
    }
  }

  if (n1 != na) {               /* there are remainders for node1 */
    /* Vectorize the loop for node b. */
    const size_t n2 = FCFC_NUM_MASK & nb;
    for (size_t j = 0; j < n2; j += FCFC_NUM_REAL) {
      vec_r x1 = SIMD_LOADU_REAL(b[0] + j);
      x1 = SIMD_SUB_REAL(x1, vsx);
      vec_r y1 = SIMD_LOADU_REAL(b[1] + j);
      y1 = SIMD_SUB_REAL(y1, vsy);
      vec_r z1 = SIMD_LOADU_REAL(b[2] + j);
      z1 = SIMD_SUB_REAL(z1, vsz);
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      const vec_r w1 = SIMD_LOADU_REAL(wb + j);
  #endif
      size_t i = n1;
      do {      /* for (size_t i = n1; i < na; i++) without initial check */
        vec_r x2 = SIMD_SET1_REAL(a[0][i]);
        vec_r y2 = SIMD_SET1_REAL(a[1][i]);
        vec_r z2 = SIMD_SET1_REAL(a[2][i]);

        vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
  #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_vector) (x1, y1, z1, x2, y2, z2,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
  #endif
            vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vpmin, vpoff,
    #endif
            vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
            vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            &x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            &z2,
    #endif
  #endif
            &vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            &vpidx,
  #endif
            &mask);

  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        vec_r vwt = SIMD_SET1_REAL(wa[i]);
        vwt = SIMD_MUL_REAL(w1, vwt);
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
        int imask = SIMD_MVMSK_REAL(mask);
  #endif

        /* Lookup tables and increment pair counts. */
        PRIVATE_NAME(update_hist_vector) (vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vpidx,
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
            vwt,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            z2,
    #endif
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
            imask,
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vns,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            vnp,
    #endif
            ioff, mask,
  #endif
            conf, cnt);
        i++;
      } while (i < na);
    }

    /* Deal with remainders for both nodes. */
    if (n2 != nb) {
  #if           FCFC_SIMD_MIN_REM_SIZE  >  0  &&                        \
                FCFC_SIMD_MIN_REM_SIZE  <  FCFC_NUM_REAL
      /* Check if it worth processing the remainders with SIMD. */
      const size_t r1 = n1 + FCFC_NUM_REAL - na;
      const size_t r2 = n2 + FCFC_NUM_REAL - nb;
      if (r1 <= FCFC_NUM_REAL - FCFC_SIMD_MIN_REM_SIZE && r1 < r2) {
        /* Vectorize node a. */
        vec_r x1 = SIMD_LOADU_REAL(a[0] + n1);
        x1 = SIMD_ADD_REAL(x1, vsx);
        vec_r y1 = SIMD_LOADU_REAL(a[1] + n1);
        y1 = SIMD_ADD_REAL(y1, vsy);
        vec_r z1 = SIMD_LOADU_REAL(a[2] + n1);
        z1 = SIMD_ADD_REAL(z1, vsz);
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        const vec_r w1 = SIMD_LOADU_REAL(wa + n1);
    #endif
        size_t j = n2;
        do {    /* for (size_t j = n2; j < nb; j++) without initial check */
          vec_r x2 = SIMD_SET1_REAL(b[0][j]);
          vec_r y2 = SIMD_SET1_REAL(b[1][j]);
          vec_r z2 = SIMD_SET1_REAL(b[2][j]);

          vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r2i vpidx;
    #endif
          vec_msk mask;

          /* Compute indices of distance bins. */
          PRIVATE_NAME(compute_dist_vector) (x1, y1, z1, x2, y2, z2,
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
              vs2min, vsoff,
    #endif
              vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
              vpmin, vpoff,
      #endif
              vpmax,
    #elif       FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
              vnmu2,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              &x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              &z2,
      #endif
    #endif
              &vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              &vpidx,
    #endif
              &mask);

    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          vec_r vwt = SIMD_SET1_REAL(wb[j]);
          vwt = SIMD_MUL_REAL(w1, vwt);
    #endif
          /* Mask objects that are not on this node. */
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
          int imask = SIMD_MVMSK_REAL(mask) & (FCFC_REM_MASK >> r1);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
          mask &= (FCFC_REM_MASK >> r1);
    #endif

          /* Lookup tables and increment pair counts. */
          PRIVATE_NAME(update_hist_vector) (vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vpidx,
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
              vwt,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              z2,
      #endif
    #endif
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
              imask,
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vns,
      #endif
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              vnp,
      #endif
              ioff, mask,
    #endif
              conf, cnt);

          j++;
        } while (j < nb);
      }
      else if (r2 <= FCFC_NUM_REAL - FCFC_SIMD_MIN_REM_SIZE) {
        /* Vectorize node b. */
        vec_r x1 = SIMD_LOADU_REAL(b[0] + n2);
        x1 = SIMD_SUB_REAL(x1, vsx);
        vec_r y1 = SIMD_LOADU_REAL(b[1] + n2);
        y1 = SIMD_SUB_REAL(y1, vsy);
        vec_r z1 = SIMD_LOADU_REAL(b[2] + n2);
        z1 = SIMD_SUB_REAL(z1, vsz);
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        const vec_r w1 = SIMD_LOADU_REAL(wb + n2);
    #endif
        size_t i = n1;
        do {    /* for (size_t i = n1; i < na; i++) without initial check */
          vec_r x2 = SIMD_SET1_REAL(a[0][i]);
          vec_r y2 = SIMD_SET1_REAL(a[1][i]);
          vec_r z2 = SIMD_SET1_REAL(a[2][i]);

          vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r2i vpidx;
    #endif
          vec_msk mask;

          /* Compute indices of distance bins. */
          PRIVATE_NAME(compute_dist_vector) (x1, y1, z1, x2, y2, z2,
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
              vs2min, vsoff,
    #endif
              vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
              vpmin, vpoff,
      #endif
              vpmax,
    #elif       FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
              vnmu2,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              &x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              &z2,
      #endif
    #endif
              &vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              &vpidx,
    #endif
              &mask);

    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          vec_r vwt = SIMD_SET1_REAL(wa[i]);
          vwt = SIMD_MUL_REAL(w1, vwt);
    #endif
          /* Mask objects that are not on this node. */
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
          int imask = SIMD_MVMSK_REAL(mask) & (FCFC_REM_MASK >> r2);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
          mask &= (FCFC_REM_MASK >> r2);
    #endif

          /* Lookup tables and increment pair counts. */
          PRIVATE_NAME(update_hist_vector) (vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vpidx,
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
              vwt,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              z2,
      #endif
    #endif
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
              imask,
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vns,
      #endif
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              vnp,
      #endif
              ioff, mask,
    #endif
              conf, cnt);

          i++;
        } while (i < na);
      }
      else {
        /* Scalar code. */
  #endif     /* FCFC_SIMD_MIN_REM_SIZE */

  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
        const int soff = conf->s2min;
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&  \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
        const int poff = conf->pmin;
  #endif
        size_t i = n1;
        do {    /* for (size_t i = n1; i < na; i++) without initial check */
          real aa[3];
          aa[0] = a[0][i] + shift[0];
          aa[1] = a[1][i] + shift[1];
          aa[2] = a[2][i] + shift[2];

          size_t j = n2;
          do {  /* for (size_t j = n2; j < nb; j++) without initial check */
            /* Compute distances and update pair counts. */
            PRIVATE_NAME(compute_dist_hist_scalar) (aa[0], aa[1], aa[2],
                b[0][j], b[1][j], b[2][j],
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
                wa[i], wb[j],
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
                soff,
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
                poff,
  #endif
                conf, cnt);
            j++;
          } while (j < nb);
          i++;
        } while (i < na);

  #if           FCFC_SIMD_MIN_REM_SIZE  >  0  &&                        \
                FCFC_SIMD_MIN_REM_SIZE  <  FCFC_NUM_REAL
      }
  #endif
    }
  }
#endif
}


/******************************************************************************
Function `count_dual_node_wrap_<IDENTIFIERS>`:
  Count the number of pairs in the separation range of interest,
  for data points on two tree nodes, with periodic wrap for every point.
Arguments:
  * `a`:        coordinates of points on the first node;
  * `b`:        coordinates of points on the second node;
  * `wa`:       weights of points on the first node;
  * `wb`:       weights of points on the second node;
  * `na`:       number of points on the first node;
  * `nb`:       number of points on the second node;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(count_dual_node_wrap) (
    real *const a[static 3], real *const b[static 3],
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const real *wa, const real *wb,
#endif
    const size_t na, const size_t nb,
    const PRIVATE_NAME(conf) *conf, void *cnt) {
  /* Scalar code. */
#if             FCFC_SIMD  ==  FCFC_SIMD_NONE
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const int soff = conf->s2min;
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&  \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  const int poff = conf->pmin;
  #endif
  for (size_t i = 0; i < na; i++) {
    for (size_t j = 0; j < nb; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_wrap_hist_scalar) (a[0][i], a[1][i], a[2][i],
          b[0][j], b[1][j], b[2][j],
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          wa[i], wb[j],
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          soff,
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          poff,
  #endif
          conf, cnt);
    }
  }
  /* Vector code. */
#else        /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
  const vec_r vhsize0 = SIMD_SET1_REAL(conf->hsize[0]);
  const vec_r vhsize1 = SIMD_SET1_REAL(conf->hsize[1]);
  const vec_r vhsize2 = SIMD_SET1_REAL(conf->hsize[2]);
  const vec_r vbsize0 = SIMD_SET1_REAL(conf->bsize[0]);
  const vec_r vbsize1 = SIMD_SET1_REAL(conf->bsize[1]);
  const vec_r vbsize2 = SIMD_SET1_REAL(conf->bsize[2]);
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const vec_r vs2min = SIMD_SET1_REAL(conf->s2min);
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vsoff = SIMD_FLOOR_REAL(vs2min);  /* offsets for table lookup */
    #else
  const vec_r2i vsoff = SIMD_R2I(vs2min);
    #endif
  #endif
  const vec_r vs2max = SIMD_SET1_REAL(conf->s2max);
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  const vec_r vpmin = SIMD_SET1_REAL(conf->pmin);
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vpoff = SIMD_FLOOR_REAL(vpmin);
      #else
  const vec_r2i vpoff = SIMD_R2I(vpmin);
      #endif
    #endif
  const vec_r vpmax = SIMD_SET1_REAL(conf->pmax);
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  const vec_r vnmu2 = SIMD_SET1_REAL(conf->nmu2);
  #endif
  #if           FCFC_SIMD  >=  FCFC_SIMD_AVX512
  /* Index offsets of private histograms for avoiding conflicts. */
    #ifdef      SINGLE_PREC
  const vec_i ioff = _mm512_setr_epi32(0, conf->ntot, conf->ntot * 2,
      conf->ntot * 3, conf->ntot * 4, conf->ntot * 5, conf->ntot * 6,
      conf->ntot * 7, 0, conf->ntot, conf->ntot * 2, conf->ntot * 3,
      conf->ntot * 4, conf->ntot * 5, conf->ntot * 6, conf->ntot * 7);
    #else
  const vec_i ioff = _mm512_setr_epi64(0, conf->ntot, conf->ntot * 2,
      conf->ntot * 3, conf->ntot * 4, conf->ntot * 5, conf->ntot * 6,
      conf->ntot * 7);
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  const vec_i vns = SIMD_SET1_INT(conf->ns);
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  const vec_i vnp = SIMD_SET1_INT(conf->np);
    #endif
  #endif
  /* Vectorize the loop for node a. */
  const size_t n1 = FCFC_NUM_MASK & na;
  for (size_t i = 0; i < n1; i += FCFC_NUM_REAL) {
    const vec_r x1 = SIMD_LOADU_REAL(a[0] + i);
    const vec_r y1 = SIMD_LOADU_REAL(a[1] + i);
    const vec_r z1 = SIMD_LOADU_REAL(a[2] + i);
  #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const vec_r w1 = SIMD_LOADU_REAL(wa + i);
  #endif
    for (size_t j = 0; j < nb; j++) {
      vec_r x2 = SIMD_SET1_REAL(b[0][j]);
      vec_r y2 = SIMD_SET1_REAL(b[1][j]);
      vec_r z2 = SIMD_SET1_REAL(b[2][j]);

      vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
      vec_r2i vpidx;
  #endif
      vec_msk mask;

      /* Compute indices of distance bins. */
      PRIVATE_NAME(compute_dist_wrap_vector) (x1, y1, z1, x2, y2, z2,
          vhsize0, vhsize1, vhsize2, vbsize0, vbsize1, vbsize2,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          vs2min, vsoff,
  #endif
          vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          vpmin, vpoff,
    #endif
          vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
          vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
          &x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
          &z2,
    #endif
  #endif
          &vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          &vpidx,
  #endif
          &mask);

  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      vec_r vwt = SIMD_SET1_REAL(wb[j]);
      vwt = SIMD_MUL_REAL(w1, vwt);
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
      int imask = SIMD_MVMSK_REAL(mask);
  #endif

      /* Lookup tables and increment pair counts. */
      PRIVATE_NAME(update_hist_vector) (vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vpidx,
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          vwt,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
          x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
          z2,
    #endif
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
          imask,
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vns,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
          vnp,
    #endif
          ioff, mask,
  #endif
          conf, cnt);
    }
  }

  if (n1 != na) {               /* there are remainders for node1 */
    /* Vectorize the loop for node b. */
    const size_t n2 = FCFC_NUM_MASK & nb;
    for (size_t j = 0; j < n2; j += FCFC_NUM_REAL) {
      const vec_r x1 = SIMD_LOADU_REAL(b[0] + j);
      const vec_r y1 = SIMD_LOADU_REAL(b[1] + j);
      const vec_r z1 = SIMD_LOADU_REAL(b[2] + j);
  #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      const vec_r w1 = SIMD_LOADU_REAL(wb + j);
  #endif

      size_t i = n1;
      do {      /* for (size_t i = n1; i < na; i++) without initial check */
        vec_r x2 = SIMD_SET1_REAL(a[0][i]);
        vec_r y2 = SIMD_SET1_REAL(a[1][i]);
        vec_r z2 = SIMD_SET1_REAL(a[2][i]);

        vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
  #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_wrap_vector) (x1, y1, z1, x2, y2, z2,
          vhsize0, vhsize1, vhsize2, vbsize0, vbsize1, vbsize2,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
  #endif
            vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vpmin, vpoff,
    #endif
            vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
            vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            &x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            &z2,
    #endif
  #endif
            &vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            &vpidx,
  #endif
            &mask);

  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        vec_r vwt = SIMD_SET1_REAL(wa[i]);
        vwt = SIMD_MUL_REAL(w1, vwt);
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
        int imask = SIMD_MVMSK_REAL(mask);
  #endif

        /* Lookup tables and increment pair counts. */
        PRIVATE_NAME(update_hist_vector) (vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vpidx,
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
            vwt,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            z2,
    #endif
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
            imask,
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vns,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            vnp,
    #endif
            ioff, mask,
  #endif
            conf, cnt);
        i++;
      } while (i < na);
    }

    /* Deal with remainders for both nodes. */
    if (n2 != nb) {
  #if           FCFC_SIMD_MIN_REM_SIZE  >  0  &&                        \
                FCFC_SIMD_MIN_REM_SIZE  <  FCFC_NUM_REAL
      /* Check if it worth processing the remainders with SIMD. */
      const size_t r1 = n1 + FCFC_NUM_REAL - na;
      const size_t r2 = n2 + FCFC_NUM_REAL - nb;
      if (r1 <= FCFC_NUM_REAL - FCFC_SIMD_MIN_REM_SIZE && r1 < r2) {
        /* Vectorize node a. */
        vec_r x1 = SIMD_LOADU_REAL(a[0] + n1);
        vec_r y1 = SIMD_LOADU_REAL(a[1] + n1);
        vec_r z1 = SIMD_LOADU_REAL(a[2] + n1);
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        const vec_r w1 = SIMD_LOADU_REAL(wa + n1);
    #endif
        size_t j = n2;
        do {    /* for (size_t j = n2; j < nb; j++) without initial check */
          vec_r x2 = SIMD_SET1_REAL(b[0][j]);
          vec_r y2 = SIMD_SET1_REAL(b[1][j]);
          vec_r z2 = SIMD_SET1_REAL(b[2][j]);

          vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r2i vpidx;
    #endif
          vec_msk mask;

          /* Compute indices of distance bins. */
          PRIVATE_NAME(compute_dist_wrap_vector) (x1, y1, z1, x2, y2, z2,
            vhsize0, vhsize1, vhsize2, vbsize0, vbsize1, vbsize2,
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
              vs2min, vsoff,
    #endif
              vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
              vpmin, vpoff,
      #endif
              vpmax,
    #elif       FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
              vnmu2,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              &x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              &z2,
      #endif
    #endif
              &vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              &vpidx,
    #endif
              &mask);

    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          vec_r vwt = SIMD_SET1_REAL(wb[j]);
          vwt = SIMD_MUL_REAL(w1, vwt);
    #endif
          /* Mask objects that are not on this node. */
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
          int imask = SIMD_MVMSK_REAL(mask) & (FCFC_REM_MASK >> r1);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
          mask &= (FCFC_REM_MASK >> r1);
    #endif

          /* Lookup tables and increment pair counts. */
          PRIVATE_NAME(update_hist_vector) (vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vpidx,
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
              vwt,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              z2,
      #endif
    #endif
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
              imask,
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vns,
      #endif
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              vnp,
      #endif
              ioff, mask,
    #endif
              conf, cnt);

          j++;
        } while (j < nb);
      }
      else if (r2 <= FCFC_NUM_REAL - FCFC_SIMD_MIN_REM_SIZE) {
        /* Vectorize node b. */
        vec_r x1 = SIMD_LOADU_REAL(b[0] + n2);
        vec_r y1 = SIMD_LOADU_REAL(b[1] + n2);
        vec_r z1 = SIMD_LOADU_REAL(b[2] + n2);
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        const vec_r w1 = SIMD_LOADU_REAL(wb + n2);
    #endif
        size_t i = n1;
        do {    /* for (size_t i = n1; i < na; i++) without initial check */
          vec_r x2 = SIMD_SET1_REAL(a[0][i]);
          vec_r y2 = SIMD_SET1_REAL(a[1][i]);
          vec_r z2 = SIMD_SET1_REAL(a[2][i]);

          vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r2i vpidx;
    #endif
          vec_msk mask;

          /* Compute indices of distance bins. */
          PRIVATE_NAME(compute_dist_wrap_vector) (x1, y1, z1, x2, y2, z2,
            vhsize0, vhsize1, vhsize2, vbsize0, vbsize1, vbsize2,
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
              vs2min, vsoff,
    #endif
              vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
              vpmin, vpoff,
      #endif
              vpmax,
    #elif       FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
              vnmu2,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              &x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              &z2,
      #endif
    #endif
              &vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              &vpidx,
    #endif
              &mask);

    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          vec_r vwt = SIMD_SET1_REAL(wa[i]);
          vwt = SIMD_MUL_REAL(w1, vwt);
    #endif
          /* Mask objects that are not on this node. */
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
          int imask = SIMD_MVMSK_REAL(mask) & (FCFC_REM_MASK >> r2);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
          mask &= (FCFC_REM_MASK >> r2);
    #endif

          /* Lookup tables and increment pair counts. */
          PRIVATE_NAME(update_hist_vector) (vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vpidx,
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
              vwt,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
              x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              z2,
      #endif
    #endif
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
              imask,
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              vns,
      #endif
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
              vnp,
      #endif
              ioff, mask,
    #endif
              conf, cnt);

          i++;
        } while (i < na);
      }
      else {
        /* Scalar code. */
  #endif     /* FCFC_SIMD_MIN_REM_SIZE */

  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
        const int soff = conf->s2min;
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&  \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
        const int poff = conf->pmin;
  #endif
        size_t i = n1;
        do {    /* for (size_t i = n1; i < na; i++) without initial check */
          size_t j = n2;
          do {  /* for (size_t j = n2; j < nb; j++) without initial check */
            /* Compute distances and update pair counts. */
            PRIVATE_NAME(compute_dist_wrap_hist_scalar) (
                a[0][i], a[1][i], a[2][i],
                b[0][j], b[1][j], b[2][j],
  #if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
                wa[i], wb[j],
  #endif
  #if             FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
                soff,
  #endif
  #if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                  FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
                poff,
  #endif
                conf, cnt);
            j++;
          } while (j < nb);
          i++;
        } while (i < na);

  #if           FCFC_SIMD_MIN_REM_SIZE  >  0  &&                        \
                FCFC_SIMD_MIN_REM_SIZE  <  FCFC_NUM_REAL
      }
  #endif
    }
  }
#endif
}


/******************************************************************************
Function `count_single_node_<IDENTIFIERS>`:
  Count the number of pairs in the separation range of interest,
  for data points belonging to a single tree node.
Arguments:
  * `a`:        coordinates of points on the node;
  * `wa`:       weights of points on the node;
  * `na`:       number of points on the node;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(count_single_node) (real *const a[static 3],
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
     const real *wa,
#endif
     const size_t na, const PRIVATE_NAME(conf) *conf, void *cnt) {
  /* Scalar code. */
#if             FCFC_SIMD  ==  FCFC_SIMD_NONE
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const int soff = conf->s2min;
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&  \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  const int poff = conf->pmin;
  #endif
  for (size_t i = 0; i < na - 1; i++) {
    for (size_t j = i + 1; j < na; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_hist_scalar) (a[0][i], a[1][i], a[2][i],
          a[0][j], a[1][j], a[2][j],
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          wa[i], wa[j],
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          soff,
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          poff,
  #endif
          conf, cnt);
    }
  }
  /* Vector code. */
#else        /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  const vec_r vs2min = SIMD_SET1_REAL(conf->s2min);
    #if         defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vsoff = SIMD_FLOOR_REAL(vs2min);  /* offsets for table lookup */
    #else
  const vec_r2i vsoff = SIMD_R2I(vs2min);
    #endif
  #endif
  const vec_r vs2max = SIMD_SET1_REAL(conf->s2max);
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  const vec_r vpmin = SIMD_SET1_REAL(conf->pmin);
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vpoff = SIMD_FLOOR_REAL(vpmin);
      #else
  const vec_r2i vpoff = SIMD_R2I(vpmin);
      #endif
    #endif
  const vec_r vpmax = SIMD_SET1_REAL(conf->pmax);
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  const vec_r vnmu2 = SIMD_SET1_REAL(conf->nmu2);
  #endif
  #if           FCFC_SIMD  >=  FCFC_SIMD_AVX512
  /* Index offsets of private histograms for avoiding conflicts. */
    #ifdef      SINGLE_PREC
  const vec_i ioff = _mm512_setr_epi32(0, conf->ntot, conf->ntot * 2,
      conf->ntot * 3, conf->ntot * 4, conf->ntot * 5, conf->ntot * 6,
      conf->ntot * 7, 0, conf->ntot, conf->ntot * 2, conf->ntot * 3,
      conf->ntot * 4, conf->ntot * 5, conf->ntot * 6, conf->ntot * 7);
    #else
  const vec_i ioff = _mm512_setr_epi64(0, conf->ntot, conf->ntot * 2,
      conf->ntot * 3, conf->ntot * 4, conf->ntot * 5, conf->ntot * 6,
      conf->ntot * 7);
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
  const vec_i vns = SIMD_SET1_INT(conf->ns);
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  const vec_i vnp = SIMD_SET1_INT(conf->np);
    #endif
  #endif
  /* Vectorize the second loop. */
  for (size_t i = 0; i < na - 1; i++) {
    size_t j = i + 1;
    const size_t nb = na - j;
    const size_t n2 = nb & FCFC_NUM_MASK;

    /* Use SIMD only when there are enough objects. */
  #if           FCFC_SIMD_MIN_REM_SIZE  >  0  &&                        \
                FCFC_SIMD_MIN_REM_SIZE  <  FCFC_NUM_REAL
    if (nb >= FCFC_SIMD_MIN_REM_SIZE)
  #else
    if (n2 != 0)
  #endif
    {
      const vec_r x1 = SIMD_SET1_REAL(a[0][i]);
      const vec_r y1 = SIMD_SET1_REAL(a[1][i]);
      const vec_r z1 = SIMD_SET1_REAL(a[2][i]);
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      const vec_r w1 = SIMD_SET1_REAL(wa[i]);
  #endif

      for (; j < n2 + i + 1; j += FCFC_NUM_REAL) {
        vec_r x2 = SIMD_LOADU_REAL(a[0] + j);
        vec_r y2 = SIMD_LOADU_REAL(a[1] + j);
        vec_r z2 = SIMD_LOADU_REAL(a[2] + j);

        vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
  #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_vector) (x1, y1, z1, x2, y2, z2,
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
  #endif
            vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vpmin, vpoff,
    #endif
            vpmax,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
            vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            &x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            &z2,
    #endif
  #endif
            &vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            &vpidx,
  #endif
            &mask);

  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
        vec_r vwt = SIMD_MSKLOAD_REAL(wa + j, SIMD_CAST_R2I(mask));
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
        vec_r vwt = SIMD_MSK0_LOADU_REAL(mask, wa + j);
    #endif
        vwt = SIMD_MUL_REAL(w1, vwt);
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
        int imask = SIMD_MVMSK_REAL(mask);
  #endif

        /* Lookup tables and increment pair counts. */
        PRIVATE_NAME(update_hist_vector) (vsidx,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vpidx,
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
            vwt,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            x2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            z2,
    #endif
  #endif
  #if           FCFC_SIMD  <  FCFC_SIMD_AVX512
            imask,
  #else      /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vns,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            vnp,
    #endif
            ioff, mask,
  #endif
            conf, cnt);
      }

  #if           FCFC_SIMD_MIN_REM_SIZE  >  0  &&                        \
                FCFC_SIMD_MIN_REM_SIZE  <  FCFC_NUM_REAL
      /* Vectorize if there are enough remainders. */
      const size_t r2 = n2 + FCFC_NUM_REAL - nb;
      if (r2 <= FCFC_NUM_REAL - FCFC_SIMD_MIN_REM_SIZE) {
        vec_r x2 = SIMD_LOADU_REAL(a[0] + j);
        vec_r y2 = SIMD_LOADU_REAL(a[1] + j);
        vec_r z2 = SIMD_LOADU_REAL(a[2] + j);

        vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
    #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_vector) (x1, y1, z1, x2, y2, z2,
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
    #endif
            vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vpmin, vpoff,
      #endif
            vpmax,
    #elif       FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
            vnmu2,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            &x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            &z2,
      #endif
    #endif
            &vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            &vpidx,
    #endif
            &mask);

    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      #if       FCFC_SIMD  <  FCFC_SIMD_AVX512
        vec_r vwt = SIMD_MSKLOAD_REAL(wa + j, SIMD_CAST_R2I(mask));
      #else  /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
        vec_r vwt = SIMD_MSK0_LOADU_REAL(mask, wa + j);
      #endif
        vwt = SIMD_MUL_REAL(w1, vwt);
    #endif
        /* Mask objects that are not on this node. */
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
        int imask = SIMD_MVMSK_REAL(mask) & (FCFC_REM_MASK >> r2);
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
        mask &= (FCFC_REM_MASK >> r2);
    #endif

        /* Lookup tables and increment pair counts. */
        PRIVATE_NAME(update_hist_vector) (vsidx,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vpidx,
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
            vwt,
    #endif
    #if         FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
            x2,
      #if       FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            z2,
      #endif
    #endif
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
            imask,
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  ||          \
                FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            vns,
      #endif
      #if       FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID  &&          \
                FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
            vnp,
      #endif
            ioff, mask,
    #endif
            conf, cnt);
        j = na;
      }
  #endif
    }

  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
    const int soff = conf->s2min;
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&  \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
    const int poff = conf->pmin;
  #endif
    /* Scalar code for remainders. */
    for (; j < na; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_hist_scalar) (a[0][i], a[1][i], a[2][i],
          a[0][j], a[1][j], a[2][j],
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          wa[i], wa[j],
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          soff,
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          poff,
  #endif
          conf, cnt);
    }
  }
#endif
}


#undef FCFC_LOOKUP_HYBRID
#undef FCFC_UPDATE_HIST
#undef FCFC_UPDATE_PAIRCNT

#endif
