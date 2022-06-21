/*******************************************************************************
* 2pt/metric_common.c: this file is part of the FCFC program.

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
  * `s1`:       sum of coordinate squares of the first point;
  * `x2,y2,z2`: coordinates of the second point;
  * `s2`:       sum of coordinate squares of the second point;
  * `w1`:       weight of the first point;
  * `w2`:       weight of the second point;
  * `soff`:     offset for squared distance (or s_perp) lookup;
  * `poff`:     offset for pi lookup;
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(compute_dist_hist_scalar) (
    const real x1, const real y1, const real z1,
#if             FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    const real s1,
#endif
    const real x2, const real y2, const real z2,
#if             FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    const real s2,
#endif
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
#if             FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
  register real dx = x1 - x2;
  register real dy = y1 - y2;
  register real dz = z1 - z2;
  register real dist = dx * dx + dy * dy + dz * dz;     /* s squared */
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_ISO */
  register real t = (x1 * x2 + y1 * y2 + z1 * z2) * 2;
  register real s = s1 + s2;
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  register real dist = s - t;                           /* s squared */
  /* Check the separation (or s_perp) range. */
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist >= conf->s2max || dist < conf->s2min) return;
    #else    /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist >= conf->s2max) return;
    #endif
  #endif
  register real d = s1 - s2;
  register real pi = d * d / (s + t);                   /* pi (s_par) squared */
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  /* Check the pi range. */
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  if (pi >= conf->p2max || pi < conf->p2min) return;
    #else    /* FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO */
  if (pi >= conf->p2max) return;
    #endif
  register real dist = s - t - pi;                      /* s_perp squared */
  /* Caveat: `dist` can be negative due to round error,
   *         but this should not break the code.*/
  #endif
#endif       /* FCFC_BIN_TYPE */

#if             FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
#else        /* FCFC_BIN_TYPE  !=  FCFC_BIN_SPI */
#endif

  /* Check the separation (or s_perp) range. */
#if             FCFC_BIN_TYPE  !=  FCFC_BIN_SMU
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  if (dist >= conf->s2max || dist < conf->s2min) return;
  #else      /* FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO */
  if (dist >= conf->s2max) return;
  #endif
#else        /* FCFC_BIN_TYPE  ==  FCFC_BIN_SMU */
  /* Compute mu and lookup the mu bin. */
  int pidx = (dist < REAL_EPS) ? 0 : (pi / dist) * conf->nmu2;
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
  FCFC_LOOKUP_HYBRID(pidx, conf->np, pi, conf->p2bin);
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
  * `s1`:       sum of coordinate squares of the first points;
  * `x2,y2,z2`: coordinates of the second points;
  * `s2`:       sum of coordinate squares of the second points;
  * `vs2min`:   minimum squared distance (or s_perp) of interest;
  * `vs2max`:   maximum squared distance (or s_perp) of interest;
  * `vsoff`:    offset for squared distance (or s_perp) lookup;
  * `vp2min`:   minimum squared pi of interest;
  * `vp2max`:   maximum squared pi of interest;
  * `vpoff`:    offset for pi lookup;
  * `vnmu2`:    the squared number of mu bins;
  * `vs2`:      the evaluated squared distance (or s_perp);
  * `vp2`:      the evaluated squared pi;
  * `vsid`:     indices of squared distance bins in the lookup table;
  * `vpid`:     indices of pi (or squared mu) bins in the lookup table;
  * `msk`:      mask indicating whether the distances are of interest.
******************************************************************************/
static inline void PRIVATE_NAME(compute_dist_vector) (
    const vec_r x1, const vec_r y1, const vec_r z1,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    const vec_r s1,
  #endif
    vec_r x2, vec_r y2, vec_r z2,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    vec_r s2,
  #endif
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
    const vec_r vp2min,
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
    const vec_r vpoff,
      #else
    const vec_r2i vpoff,
      #endif
    #endif
    const vec_r vp2max,
  #elif         FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
    const vec_r vnmu2,
  #endif
  #if           FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
    vec_r *vs2,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    vec_r *vp2,
    #endif
  #endif
    vec_r2i *vsid,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    vec_r2i *vpid,
  #endif
    vec_msk *msk) {
  /* Compute (squared) distances. */
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
  x2 = SIMD_SUB_REAL(x1, x2);
  y2 = SIMD_SUB_REAL(y1, y2);
  z2 = SIMD_SUB_REAL(z1, z2);
  x2 = SIMD_MUL_REAL(x2, x2);
    #ifdef      FCFC_SIMD_FMA
  x2 = SIMD_FMADD_REAL(y2, y2, x2);
  x2 = SIMD_FMADD_REAL(z2, z2, x2);                     /* s squared */
    #else
  y2 = SIMD_MUL_REAL(y2, y2);
  x2 = SIMD_ADD_REAL(x2, y2);
  z2 = SIMD_MUL_REAL(z2, z2);
  x2 = SIMD_ADD_REAL(x2, z2);                           /* s squared */
    #endif
  #else      /* FCFC_BIN_TYPE  !=  FCFC_BIN_ISO */
  x2 = SIMD_MUL_REAL(x1, x2);
    #ifdef      FCFC_SIMD_FMA
  x2 = SIMD_FMADD_REAL(y1, y2, x2);
  x2 = SIMD_FMADD_REAL(z1, z2, x2);
    #else
  y2 = SIMD_MUL_REAL(y1, y2);
  x2 = SIMD_ADD_REAL(x2, y2);
  z2 = SIMD_MUL_REAL(z1, z2);
  x2 = SIMD_ADD_REAL(x2, z2);
    #endif
  x2 = SIMD_MUL_REAL(x2, SIMD_SET1_REAL(2.0));  /* t in scalar code */
  y2 = SIMD_ADD_REAL(s1, s2);                   /* s in scalar code */
  s2 = SIMD_SUB_REAL(s1, s2);                   /* d in scalar code */
  s2 = SIMD_MUL_REAL(s2, s2);
  z2 = SIMD_ADD_REAL(x2, y2);
  z2 = SIMD_DIV_REAL(s2, z2);                           /* pi squared */
  x2 = SIMD_SUB_REAL(y2, x2);                           /* s squared */
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  x2 = SIMD_SUB_REAL(x2, z2);                           /* s_perp squared */
  /* Caveat: `x2` can be negative due to round error,
   *         but this should not break the code.*/
    #endif
  #endif     /* FCFC_BIN_TYPE */

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
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  mask = SIMD_MSK_CMP_REAL(mask, x2, vs2min, _CMP_GE_OQ);
    #endif
  #endif
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_SIMD  <  FCFC_SIMD_AVX512
  mask1 = SIMD_CMP_REAL(z2, vp2max, _CMP_LT_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  mask1 = SIMD_CMP_REAL(z2, vp2min, _CMP_GE_OQ);
  mask = SIMD_AND_REAL(mask, mask1);
      #endif
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
  mask = SIMD_MSK_CMP_REAL(mask, z2, vp2max, _CMP_LT_OQ);
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  mask = SIMD_MSK_CMP_REAL(mask, z2, vp2min, _CMP_GE_OQ);
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
    #else    /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
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
  *vp2 = z2;
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
  * `vpidx`:    indices of squared pi (or mu) in the lookup table;
  * `vwt`:      weights to be added to the pair counts;
  * `vdist`:    the squared distances (or s_perp);
  * `vpi`:      the squared pi;
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
        conf->np, pi[0], conf->p2bin, wt[0]);
  }
  if (imask & 0x2) {
    FCFC_UPDATE_PAIRCNT(sidx[1], conf->stab, conf->ns, dist[1],
        conf->s2bin, pidx[1], conf->mutab, pidx[1], conf->ptab,
        conf->np, pi[1], conf->p2bin, wt[1]);
  }
  if (imask & 0x4) {
    FCFC_UPDATE_PAIRCNT(sidx[2], conf->stab, conf->ns, dist[2],
        conf->s2bin, pidx[2], conf->mutab, pidx[2], conf->ptab,
        conf->np, pi[2], conf->p2bin, wt[2]);
  }
  if (imask & 0x8) {
    FCFC_UPDATE_PAIRCNT(sidx[3], conf->stab, conf->ns, dist[3],
        conf->s2bin, pidx[3], conf->mutab, pidx[3], conf->ptab,
        conf->np, pi[3], conf->p2bin, wt[3]);
  }
    #ifdef    SINGLE_PREC
  if (imask & 0x10) {
    FCFC_UPDATE_PAIRCNT(sidx[4], conf->stab, conf->ns, dist[4],
        conf->s2bin, pidx[4], conf->mutab, pidx[4], conf->ptab,
        conf->np, pi[4], conf->p2bin, wt[4]);
  }
  if (imask & 0x20) {
    FCFC_UPDATE_PAIRCNT(sidx[5], conf->stab, conf->ns, dist[5],
        conf->s2bin, pidx[5], conf->mutab, pidx[5], conf->ptab,
        conf->np, pi[5], conf->p2bin, wt[5]);
  }
  if (imask & 0x40) {
    FCFC_UPDATE_PAIRCNT(sidx[6], conf->stab, conf->ns, dist[6],
        conf->s2bin, pidx[6], conf->mutab, pidx[6], conf->ptab,
        conf->np, pi[6], conf->p2bin, wt[6]);
  }
  if (imask & 0x80) {
    FCFC_UPDATE_PAIRCNT(sidx[7], conf->stab, conf->ns, dist[7],
        conf->s2bin, pidx[7], conf->mutab, pidx[7], conf->ptab,
        conf->np, pi[7], conf->p2bin, wt[7]);
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
      vec_r ref = SIMD_MSK_GATHER_REAL(vpi, imask, vpid, conf->p2bin,
          sizeof(real));
      imask = SIMD_MSK_CMP_REAL(imask, vpi, ref, _CMP_LT_OQ);
      if (imask == 0) break;
      vpid = SIMD_MSK_SUB_INT(vpid, imask, vpid, SIMD_SET1_INT(1));
    }
  }
      #endif
    #endif   /* FCFC_BIN_TYPE */
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
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
  * `conf`:     the structure for pair count configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static inline void PRIVATE_NAME(count_dual_node) (
    real *const a[static 4], real *const b[static 4],
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
  const int poff = conf->p2min;
  #endif
  for (size_t i = 0; i < na; i++) {
    for (size_t j = 0; j < nb; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_hist_scalar) (a[0][i], a[1][i], a[2][i],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          a[3][i],
  #endif
          b[0][j], b[1][j], b[2][j],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          b[3][j],
  #endif
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
  const vec_r vp2min = SIMD_SET1_REAL(conf->p2min);
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vpoff = SIMD_FLOOR_REAL(vp2min);  /* offsets for table lookup */
      #else
  const vec_r2i vpoff = SIMD_R2I(vp2min);
      #endif
    #endif
  const vec_r vp2max = SIMD_SET1_REAL(conf->p2max);
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
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
    const vec_r s1 = SIMD_LOADU_REAL(a[3] + i);
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
    const vec_r w1 = SIMD_LOADU_REAL(wa + i);
  #endif
    for (size_t j = 0; j < nb; j++) {
      vec_r x2 = SIMD_SET1_REAL(b[0][j]);
      vec_r y2 = SIMD_SET1_REAL(b[1][j]);
      vec_r z2 = SIMD_SET1_REAL(b[2][j]);
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
      vec_r s2 = SIMD_SET1_REAL(b[3][j]);
  #endif

      vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
      vec_r2i vpidx;
  #endif
      vec_msk mask;

      /* Compute indices of distance bins. */
      PRIVATE_NAME(compute_dist_vector) (x1, y1, z1,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          s1,
  #endif
          x2, y2, z2,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          s2,
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
          vs2min, vsoff,
  #endif
          vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
          vp2min, vpoff,
    #endif
          vp2max,
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
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
      const vec_r s1 = SIMD_LOADU_REAL(b[3] + j);
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      const vec_r w1 = SIMD_LOADU_REAL(wb + j);
  #endif

      size_t i = n1;
      do {      /* for (size_t i = n1; i < na; i++) without initial check */
        vec_r x2 = SIMD_SET1_REAL(a[0][i]);
        vec_r y2 = SIMD_SET1_REAL(a[1][i]);
        vec_r z2 = SIMD_SET1_REAL(a[2][i]);
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r s2 = SIMD_SET1_REAL(a[3][i]);
  #endif

        vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
  #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_vector) (x1, y1, z1,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            s1,
  #endif
            x2, y2, z2,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            s2,
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
  #endif
            vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vp2min, vpoff,
    #endif
            vp2max,
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
        const vec_r x1 = SIMD_LOADU_REAL(a[0] + n1);
        const vec_r y1 = SIMD_LOADU_REAL(a[1] + n1);
        const vec_r z1 = SIMD_LOADU_REAL(a[2] + n1);
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        const vec_r s1 = SIMD_LOADU_REAL(a[3] + n1);
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        const vec_r w1 = SIMD_LOADU_REAL(wa + n1);
    #endif

        size_t j = n2;
        do {    /* for (size_t j = n2; j < nb; j++) without initial check */
          vec_r x2 = SIMD_SET1_REAL(b[0][j]);
          vec_r y2 = SIMD_SET1_REAL(b[1][j]);
          vec_r z2 = SIMD_SET1_REAL(b[2][j]);
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r s2 = SIMD_SET1_REAL(b[3][j]);
    #endif

          vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r2i vpidx;
    #endif
          vec_msk mask;

          /* Compute indices of distance bins. */
          PRIVATE_NAME(compute_dist_vector) (x1, y1, z1,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              s1,
    #endif
              x2, y2, z2,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              s2,
    #endif
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
              vs2min, vsoff,
    #endif
              vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
              vp2min, vpoff,
      #endif
              vp2max,
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
        const vec_r x1 = SIMD_LOADU_REAL(b[0] + n2);
        const vec_r y1 = SIMD_LOADU_REAL(b[1] + n2);
        const vec_r z1 = SIMD_LOADU_REAL(b[2] + n2);
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        const vec_r s1 = SIMD_LOADU_REAL(b[3] + n2);
    #endif
    #if         FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
        const vec_r w1 = SIMD_LOADU_REAL(wb + n2);
    #endif
        size_t i = n1;
        do {    /* for (size_t i = n1; i < na; i++) without initial check */
          vec_r x2 = SIMD_SET1_REAL(a[0][i]);
          vec_r y2 = SIMD_SET1_REAL(a[1][i]);
          vec_r z2 = SIMD_SET1_REAL(a[2][i]);
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r s2 = SIMD_SET1_REAL(a[3][i]);
    #endif

          vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          vec_r2i vpidx;
    #endif
          vec_msk mask;

          /* Compute indices of distance bins. */
          PRIVATE_NAME(compute_dist_vector) (x1, y1, z1,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              s1,
    #endif
              x2, y2, z2,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
              s2,
    #endif
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
              vs2min, vsoff,
    #endif
              vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
              vp2min, vpoff,
      #endif
              vp2max,
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
        const int poff = conf->p2min;
  #endif
        size_t i = n1;
        do {    /* for (size_t i = n1; i < na; i++) without initial check */
          size_t j = n2;
          do {  /* for (size_t j = n2; j < nb; j++) without initial check */
            /* Compute distances and update pair counts. */
            PRIVATE_NAME(compute_dist_hist_scalar) (a[0][i], a[1][i], a[2][i],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
                a[3][i],
  #endif
                b[0][j], b[1][j], b[2][j],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
                b[3][j],
  #endif
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
static inline void PRIVATE_NAME(count_single_node) (real *const a[static 4],
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
  const int poff = conf->p2min;
  #endif
  for (size_t i = 0; i < na - 1; i++) {
    for (size_t j = i + 1; j < na; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_hist_scalar) (a[0][i], a[1][i], a[2][i],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          a[3][i],
  #endif
          a[0][j], a[1][j], a[2][j],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          a[3][j],
  #endif
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
  const vec_r vp2min = SIMD_SET1_REAL(conf->p2min);
      #if       defined(SINGLE_PREC)  &&  FCFC_SIMD == FCFC_SIMD_AVX
  const vec_r vpoff = SIMD_FLOOR_REAL(vp2min);  /* offsets for table lookup */ 
      #else
  const vec_r2i vpoff = SIMD_R2I(vp2min);
      #endif
    #endif
  const vec_r vp2max = SIMD_SET1_REAL(conf->p2max);
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
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
      const vec_r s1 = SIMD_SET1_REAL(a[3][i]);
  #endif
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
      const vec_r w1 = SIMD_SET1_REAL(wa[i]);
  #endif

      for (; j < n2 + i + 1; j += FCFC_NUM_REAL) {
        vec_r x2 = SIMD_LOADU_REAL(a[0] + j);
        vec_r y2 = SIMD_LOADU_REAL(a[1] + j);
        vec_r z2 = SIMD_LOADU_REAL(a[2] + j);
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r s2 = SIMD_LOADU_REAL(a[3] + j);
  #endif

        vec_r2i vsidx;
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
  #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_vector) (x1, y1, z1,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            s1,
  #endif
            x2, y2, z2,
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            s2,
  #endif
  #if           FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
  #endif
            vs2max,
  #if           FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
    #if         FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vp2min, vpoff,
    #endif
            vp2max,
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
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r s2 = SIMD_LOADU_REAL(a[3] + j);
    #endif

        vec_r2i vsidx;
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
        vec_r2i vpidx;
    #endif
        vec_msk mask;

        /* Compute indices of distance bins. */
        PRIVATE_NAME(compute_dist_vector) (x1, y1, z1,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            s1,
    #endif
            x2, y2, z2,
    #if         FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
            s2,
    #endif
    #if         FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
            vs2min, vsoff,
    #endif
            vs2max,
    #if         FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
      #if       FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
            vp2min, vpoff,
      #endif
            vp2max,
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
    const int poff = conf->p2min;
  #endif
    /* Scalar code for remainders. */
    for (; j < na; j++) {
      /* Compute distances and update pair counts. */
      PRIVATE_NAME(compute_dist_hist_scalar) (a[0][i], a[1][i], a[2][i],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          a[3][i],
  #endif
          a[0][j], a[1][j], a[2][j],
  #if           FCFC_BIN_TYPE  !=  FCFC_BIN_ISO
          a[3][j],
  #endif
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
