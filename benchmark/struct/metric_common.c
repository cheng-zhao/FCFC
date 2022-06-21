/*******************************************************************************
* benchmark/struct/metric_common.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#if defined(BENCHMARK_PAIRCOUNT_TYPE) && defined(BENCHMARK_BIN_SMIN)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef BENCHMARK_PAIRCOUNT_NAME
  #undef BENCHMARK_PAIRCOUNT_NAME
#endif
#ifdef BENCHMARK_BIN_SMIN_NAME
  #undef BENCHMARK_BIN_SMIN_NAME
#endif

#if     BENCHMARK_PAIRCOUNT_TYPE == BENCHMARK_PAIRCOUNT_BOX
  #define BENCHMARK_PAIRCOUNT_NAME      _box
#elif   BENCHMARK_PAIRCOUNT_TYPE == BENCHMARK_PAIRCOUNT_NOBOX
  #define BENCHMARK_PAIRCOUNT_NAME
#else
  #error unexpected definition of `BENCHMARK_PAIRCOUNT_TYPE`
#endif

#if     BENCHMARK_BIN_SMIN == BENCHMARK_BIN_MIN_ZERO
  #define BENCHMARK_BIN_SMIN_NAME       _smin0
#elif   BENCHMARK_BIN_SMIN == BENCHMARK_BIN_MIN_NONZERO
  #define BENCHMARK_BIN_SMIN_NAME
#else
  #error unexpected definition of `BENCHMARK_BIN_SMIN`
#endif

/* Macros for generating datatype and function names. */
#ifndef CONCAT_STRING
  #define CONCAT_STRING(a,b,c,d)        a##_##b##c##d
#endif

#ifndef PRIVATE_NAME
  #define PRIVATE_NAME(a,b,c,d)         CONCAT_STRING(a,b,c,d)
#endif


/*============================================================================*\
                   Function for counting pairs from two nodes
\*============================================================================*/

/******************************************************************************
Function `count_dual_node_<BENCHMARK_PAIRCOUNT_NAME><BENCHMARK_BIN_SMIN_NAME>`:
  Count the number of pairs in the separation range of interest,
  for data points belonging to two tree nodes.
Arguments:
  * `na`:       number of points on the first node;
  * `a`:        coordinates of points on the first node;
  * `nb`:       number of points on the second node;
  * `b`:        coordinates of points on the second node;
  * `shift`:    coordinate offsets for periodic boundary conditions;
  * `r2min`:    minimum squared distance of interest;
  * `r2max`:    maximum squared distance of interest;
  * `npair`:    total number of pairs.
******************************************************************************/
static inline void PRIVATE_NAME(count_dual_node, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME, )
    (const size_t na, real **a, const size_t nb, real **b,
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
    const real shift[3],
#endif
#if     BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
    const real r2min,
#endif
    const real r2max, size_t *npair) {
#if             BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  for (size_t i = 0; i < na; i++) {
  #if   BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
    real aa[3];
    aa[0] = a[0][i] + shift[0];
    aa[1] = a[1][i] + shift[1];
    aa[2] = a[2][i] + shift[2];
  #endif
    for (size_t j = 0; j < nb; j++) {
  #if   BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
      register real dx = aa[0] - b[0][j];
      register real dy = aa[1] - b[1][j];
      register real dz = aa[2] - b[2][j];
  #else
      register real dx = a[0][i] - b[0][j];
      register real dy = a[1][i] - b[1][j];
      register real dz = a[2][i] - b[2][j];
  #endif
      register real dist = dx * dx + dy * dy + dz * dz;

  #if   BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
      if (dist >= r2min && dist < r2max) (*npair)++;
  #else
      if (dist < r2max) (*npair)++;
  #endif
    }
  }
#else        /* BENCHMARK_SIMD  !=  BENCHMARK_SIMD_NONE */
  #if BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
  const vec_r vr2min = SIMD_SET1_REAL(r2min);
  #endif
  const vec_r vr2max = SIMD_SET1_REAL(r2max);
  const size_t n2 = BENCHMARK_NUM_MASK & nb;
  const size_t r2 = n2 + BENCHMARK_NUM_REAL - nb;       /* for movemask */
  if (r2 != BENCHMARK_NUM_REAL) {       /* there are remainders for node2 */
    for (size_t i = 0; i != na; i++) {
  #if BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
      const vec_r x1 = SIMD_SET1_REAL(a[0][i] + shift[0]);
      const vec_r y1 = SIMD_SET1_REAL(a[1][i] + shift[1]);
      const vec_r z1 = SIMD_SET1_REAL(a[2][i] + shift[2]);
  #else
      const vec_r x1 = SIMD_SET1_REAL(a[0][i]);
      const vec_r y1 = SIMD_SET1_REAL(a[1][i]);
      const vec_r z1 = SIMD_SET1_REAL(a[2][i]);
  #endif
      for (size_t j = 0; j != n2; j += BENCHMARK_NUM_REAL) {
        /* Compute squared distances. */
        vec_r x2 = SIMD_LOAD_REAL(b[0] + j);
        vec_r y2 = SIMD_LOAD_REAL(b[1] + j);
        vec_r z2 = SIMD_LOAD_REAL(b[2] + j);
        x2 = SIMD_SUB_REAL(x1, x2);
        y2 = SIMD_SUB_REAL(y1, y2);
        z2 = SIMD_SUB_REAL(z1, z2);
        x2 = SIMD_MUL_REAL(x2, x2);
  #ifdef BENCHMARK_SIMD_FMA
        y2 = SIMD_FMADD_REAL(y2, y2, x2);
        z2 = SIMD_FMADD_REAL(z2, z2, y2);
  #else
        y2 = SIMD_MUL_REAL(y2, y2);
        y2 = SIMD_ADD_REAL(x2, y2);
        z2 = SIMD_MUL_REAL(z2, z2);
        z2 = SIMD_ADD_REAL(y2, z2);
  #endif
        /* Check distance range. */
        vec_msk mask = SIMD_CMP_REAL(z2, vr2max, _CMP_LT_OQ);
  #if BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
        vec_msk mask1 = SIMD_CMP_REAL(z2, vr2min, _CMP_GE_OQ);
        mask = SIMD_AND_REAL(mask, mask1);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
        mask = SIMD_MSK_CMP_REAL(mask, z2, vr2min, _CMP_GE_OQ);
    #endif
  #endif
  #if           BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
        const int imask = SIMD_MVMSK_REAL(mask);
        *npair += BENCHMARK_POPCNT(imask);
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
        *npair += BENCHMARK_POPCNT(mask);
  #endif
      }
      /* Deal with node2 remainders. */
      vec_r x2 = SIMD_LOADU_REAL(b[0] + n2);
      vec_r y2 = SIMD_LOADU_REAL(b[1] + n2);
      vec_r z2 = SIMD_LOADU_REAL(b[2] + n2);
      x2 = SIMD_SUB_REAL(x1, x2);
      y2 = SIMD_SUB_REAL(y1, y2);
      z2 = SIMD_SUB_REAL(z1, z2);
      x2 = SIMD_MUL_REAL(x2, x2);
  #ifdef BENCHMARK_SIMD_FMA
      y2 = SIMD_FMADD_REAL(y2, y2, x2);
      z2 = SIMD_FMADD_REAL(z2, z2, y2);
  #else
      y2 = SIMD_MUL_REAL(y2, y2);
      y2 = SIMD_ADD_REAL(x2, y2);
      z2 = SIMD_MUL_REAL(z2, z2);
      z2 = SIMD_ADD_REAL(y2, z2);
  #endif

      vec_msk mask = SIMD_CMP_REAL(z2, vr2max, _CMP_LT_OQ);
  #if BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
      vec_msk mask1 = SIMD_CMP_REAL(z2, vr2min, _CMP_GE_OQ);
      mask = SIMD_AND_REAL(mask, mask1);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
      mask = SIMD_MSK_CMP_REAL(mask, z2, vr2min, _CMP_GE_OQ);
    #endif
  #endif
  #if           BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
      int imask = SIMD_MVMSK_REAL(mask);
      *npair += BENCHMARK_POPCNT(imask & (BENCHMARK_REM_MASK >> r2));
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
      *npair += BENCHMARK_POPCNT(mask & (BENCHMARK_REM_MASK >> r2));
  #endif
    }
  }
  else {      /* r2 == 0 : no remainders for node2 */
    for (size_t i = 0; i != na; i++) {
  #if BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
      const vec_r x1 = SIMD_SET1_REAL(a[0][i] + shift[0]);
      const vec_r y1 = SIMD_SET1_REAL(a[1][i] + shift[1]);
      const vec_r z1 = SIMD_SET1_REAL(a[2][i] + shift[2]);
  #else
      const vec_r x1 = SIMD_SET1_REAL(a[0][i]);
      const vec_r y1 = SIMD_SET1_REAL(a[1][i]);
      const vec_r z1 = SIMD_SET1_REAL(a[2][i]);
  #endif
      for (size_t j = 0; j != n2; j += BENCHMARK_NUM_REAL) {
        /* Compute squared distances. */
        vec_r x2 = SIMD_LOAD_REAL(b[0] + j);
        vec_r y2 = SIMD_LOAD_REAL(b[1] + j);
        vec_r z2 = SIMD_LOAD_REAL(b[2] + j);
        x2 = SIMD_SUB_REAL(x1, x2);
        y2 = SIMD_SUB_REAL(y1, y2);
        z2 = SIMD_SUB_REAL(z1, z2);
        x2 = SIMD_MUL_REAL(x2, x2);
  #ifdef BENCHMARK_SIMD_FMA
        y2 = SIMD_FMADD_REAL(y2, y2, x2);
        z2 = SIMD_FMADD_REAL(z2, z2, y2);
  #else
        y2 = SIMD_MUL_REAL(y2, y2);
        y2 = SIMD_ADD_REAL(x2, y2);
        z2 = SIMD_MUL_REAL(z2, z2);
        z2 = SIMD_ADD_REAL(y2, z2);
  #endif
        /* Check distance range. */
        vec_msk mask = SIMD_CMP_REAL(z2, vr2max, _CMP_LT_OQ);
  #if BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
    #if         BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
        vec_msk mask1 = SIMD_CMP_REAL(z2, vr2min, _CMP_GE_OQ);
        mask = SIMD_AND_REAL(mask, mask1);
    #else    /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
        mask = SIMD_MSK_CMP_REAL(mask, z2, vr2min, _CMP_GE_OQ);
    #endif
  #endif
  #if           BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
        const int imask = SIMD_MVMSK_REAL(mask);
        *npair += BENCHMARK_POPCNT(imask);
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
        *npair += BENCHMARK_POPCNT(mask);
  #endif
      }
    }
  }
#endif
}

#endif
