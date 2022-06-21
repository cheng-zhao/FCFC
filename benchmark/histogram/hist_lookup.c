/*******************************************************************************
* benchmark/struct/hist_lookup.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/*******************************************************************************
  Update the distance histogram using a lookup table.
*******************************************************************************/

/* Macros for the template functions. */
#if defined(BENCHMARK_HIST_TABLE_WIDTH)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef BENCHMARK_HIST_TABLE_DTYPE
  #undef BENCHMARK_HIST_TABLE_DTYPE
#endif
#ifdef BENCHMARK_HIST_TABLE_MASK
  #undef BENCHMARK_HIST_TABLE_MASK
#endif

#if     BENCHMARK_HIST_TABLE_WIDTH == BENCHMARK_HIST_TABLE_W8
  #define BENCHMARK_HIST_TABLE_DTYPE    uint8_t
  #define BENCHMARK_HIST_TABLE_MASK     0xFF
#elif   BENCHMARK_HIST_TABLE_WIDTH == BENCHMARK_HIST_TABLE_W16
  #define BENCHMARK_HIST_TABLE_DTYPE    uint16_t
  #define BENCHMARK_HIST_TABLE_MASK     0xFFFF
#elif   BENCHMARK_HIST_TABLE_WIDTH == BENCHMARK_HIST_TABLE_W32
  #define BENCHMARK_HIST_TABLE_DTYPE    uint32_t
  #define BENCHMARK_HIST_TABLE_MASK     0xFFFFFFFF
#elif   BENCHMARK_HIST_TABLE_WIDTH == BENCHMARK_HIST_TABLE_W64
  #define BENCHMARK_HIST_TABLE_DTYPE    uint64_t
  #define BENCHMARK_HIST_TABLE_MASK     0xFFFFFFFFFFFFFFFF
#else
  #error unexpected definition of `BENCHMARK_HIST_TABLE_WIDTH`
#endif

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b)             a##_##b
#endif

#ifndef HIST_TABLE_FUNC
  #define HIST_TABLE_FUNC(a,b)          CONCAT_FNAME(a,b)
#endif


/*============================================================================*\
           Functions for distance histogram update with lookup table
\*============================================================================*/

/******************************************************************************
Function `create_table_int_<BENCHMARK_HIST_TABLE_DTYPE>`:
  Create the lookup table for integer squared distance bin edges.
Arguments:
  * `hist`:     the structure for distance histogram.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int HIST_TABLE_FUNC(create_table_int, BENCHMARK_HIST_TABLE_DTYPE)
    (HIST *hist) {
  /* Compute table size and allocate memory. */
  const long d2min = hist->d2[0];
  const long d2max = hist->d2[hist->n];
  const long ntab = d2max - d2min;
#if     BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  BENCHMARK_HIST_TABLE_DTYPE *tab;
  /* Allocate extra 64 bits for safe gathering. */
  size_t tsize = ntab * sizeof(BENCHMARK_HIST_TABLE_DTYPE) + 64;
  /* Round up the size of the sample array for memory alignment. */
  const size_t psize = sizeof(void *) - 1;
  tsize = (tsize + psize) & (~psize);
  if (posix_memalign((void **) &tab, BENCHMARK_MEMALIGN_BYTE, tsize)) {
    P_ERR("failed to allocate aligned memory for the lookup table\n");
    return EXIT_FAILURE;
  }
#else
  BENCHMARK_HIST_TABLE_DTYPE *tab = malloc(ntab * sizeof *tab);
  if (!tab) {
    P_ERR("failed to allocate memory for the lookup table\n");
    return EXIT_FAILURE;
  }
#endif

  /* Fill values. */
  const long istart = d2min;
  int n = 1;
  for (long i = 0; i < ntab; i++) {
    if (i + istart < (long) hist->d2[n]) tab[i] = n - 1;
    else {
      if (++n > hist->n) {
        P_ERR("failed to create the lookup table\n");
        free(tab);
        return EXIT_FAILURE;
      }
      i--;
    }
  }

  hist->tab = tab;
  return 0;
}

/******************************************************************************
Function `create_table_hybrid_<BENCHMARK_HIST_TABLE_DTYPE>`:
  Create the hybrid lookup table for squared distance histogram update.
Arguments:
  * `hist`:     the structure for distance histogram.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int HIST_TABLE_FUNC(create_table_hybrid, BENCHMARK_HIST_TABLE_DTYPE)
    (HIST *hist) {
  /* Compute table size and allocate memory. */
  const long d2min = hist->d2[0];
  const long d2max = ceil(hist->d2[hist->n]);
  const long ntab = d2max - d2min;
#if     BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  BENCHMARK_HIST_TABLE_DTYPE *tab;
  /* Allocate extra 64 bits for safe gathering. */
  size_t tsize = ntab * sizeof(BENCHMARK_HIST_TABLE_DTYPE) + 64;
  /* Round up the size of the sample array for memory alignment. */
  const size_t psize = sizeof(void *) - 1;
  tsize = (tsize + psize) & (~psize);
  if (posix_memalign((void **) &tab, BENCHMARK_MEMALIGN_BYTE, tsize)) {
    P_ERR("failed to allocate aligned memory for the lookup table\n");
    return EXIT_FAILURE;
  }
#else
  BENCHMARK_HIST_TABLE_DTYPE *tab = malloc(ntab * sizeof *tab);
  if (!tab) {
    P_ERR("failed to allocate memory for the lookup table\n");
    return EXIT_FAILURE;
  }
#endif

  /* Fill values. */
  const long istart = d2min;
  int n = 1;
  real edge = hist->d2[n];
  for (long i = 0; i < ntab; i++) {
    if (istart + i < (long) edge) tab[i] = n - 1;
    else {
      while (++n <= hist->n) {
        edge = hist->d2[n];
        if (istart + i < (long) edge) {
          tab[i] = hist->n + n - 1;
          break;
        }
      }
      if (n > hist->n) {
        for (; i < ntab; i++) tab[i] = hist->n * 2;
        break;
      }
    }
  }

  hist->tab = tab;
  return 0;
}


/******************************************************************************
Function `hist_table_int_<BENCHMARK_HIST_TABLE_DTYPE>`:
  Update distance histogram using the lookup table for integer squared
  distance bin edges.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void HIST_TABLE_FUNC(hist_table_int, BENCHMARK_HIST_TABLE_DTYPE) (
    const real *d2, const size_t nsp, HIST *hist) {
  BENCHMARK_HIST_TABLE_DTYPE *tab = (BENCHMARK_HIST_TABLE_DTYPE *) hist->tab;
  const real d2min = hist->d2[0];
  const real d2max = hist->d2[hist->n];
#if     BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  const size_t offset = d2min;
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      const size_t idx = (size_t) d2[n] - offset;
      hist->cnt[tab[idx]]++;
    }
  }
#else
  const vec_r vmin = SIMD_SET1_REAL(d2min);
  const vec_r vmax = SIMD_SET1_REAL(d2max);
  #if   defined(SINGLE_PREC)  &&  BENCHMARK_SIMD == BENCHMARK_SIMD_AVX
  const vec_r offs = _mm256_floor_ps(vmin);
  #else
  const vec_r2i offs = SIMD_FLOOR_REAL(vmin);
  #endif
  #if   BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512
  /* Index offsets of private histograms for avoiding conflicts. */
    #ifdef SINGLE_PREC
  const vec_i ioff = _mm512_setr_epi32(0, hist->n, hist->n * 2, hist->n * 3,
      hist->n * 4, hist->n * 5, hist->n * 6, hist->n * 7,
      0, hist->n, hist->n * 2, hist->n * 3,
      hist->n * 4, hist->n * 5, hist->n * 6, hist->n * 7);
    #else
  const vec_i ioff = _mm512_setr_epi64(0, hist->n, hist->n * 2, hist->n * 3,
      hist->n * 4, hist->n * 5, hist->n * 6, hist->n * 7);
    #endif
  #endif
  const size_t npar = BENCHMARK_NUM_MASK & nsp;
  for (size_t n = 0; n != npar; n += BENCHMARK_NUM_REAL) {
    vec_r x = SIMD_LOAD_REAL(d2 + n);
    vec_msk mask = SIMD_CMP_REAL(x, vmin, _CMP_GE_OQ);
  #if           BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
    vec_msk mask1 = SIMD_CMP_REAL(x, vmax, _CMP_LT_OQ);
    mask = SIMD_AND_REAL(mask, mask1);
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
    mask = SIMD_MSK_CMP_REAL(mask, x, vmax, _CMP_LT_OQ);
  #endif
  #if           defined(SINGLE_PREC)  &&  BENCHMARK_SIMD == BENCHMARK_SIMD_AVX
    vec_r dx = SIMD_SUB_REAL(_mm256_floor_ps(x), offs);
    vec_r2i vidx = _mm256_cvtps_epi32(dx);
  #else
    vec_r2i vidx = SIMD_FLOOR_REAL(x);
    vidx = SIMD_SUB_R2I(vidx, offs);
  #endif
  #if           BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX
    /* Convert the mask to numbers for histogram update. */
    int imask = SIMD_MVMSK_REAL(mask);
    int32_t *idx = (int32_t *) &vidx;
    if (imask & 0x1) hist->cnt[tab[idx[0]]]++;
    if (imask & 0x2) hist->cnt[tab[idx[1]]]++;
    if (imask & 0x4) hist->cnt[tab[idx[2]]]++;
    if (imask & 0x8) hist->cnt[tab[idx[3]]]++;
    #ifdef SINGLE_PREC
    if (imask & 0x10) hist->cnt[tab[idx[4]]]++;
    if (imask & 0x20) hist->cnt[tab[idx[5]]]++;
    if (imask & 0x40) hist->cnt[tab[idx[6]]]++;
    if (imask & 0x80) hist->cnt[tab[idx[7]]]++;
    #endif
  #elif         BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX2
    /* Gather values from the lookup table. */
    #if defined(__GNUC__) && !defined(SINGLE_PREC)
    /* `__int64` is interpreted as `long long int` in gcc. */
    vec_i idx = SIMD_MSK_R2I_GATHER_INT(SIMD_ZEROS_INT(), (long long *) tab,
        vidx, SIMD_CAST_R2I(mask), sizeof(BENCHMARK_HIST_TABLE_DTYPE));
    #else
    vec_i idx = SIMD_MSK_R2I_GATHER_INT(SIMD_ZEROS_INT(), (int_wr *) tab,
        vidx, SIMD_CAST_R2I(mask), sizeof(BENCHMARK_HIST_TABLE_DTYPE));
    #endif
    idx = SIMD_AND_INT(idx, SIMD_SET1_INT(BENCHMARK_HIST_TABLE_MASK));
    mask = SIMD_AND_REAL(mask, SIMD_CAST_I2R(SIMD_SET1_INT(1)));
    int_wr *tid = (int_wr *) &idx;
    uint_wr *cnt = (uint_wr *) &mask;
    hist->cnt[tid[0]] += cnt[0];
    hist->cnt[tid[1]] += cnt[1];
    hist->cnt[tid[2]] += cnt[2];
    hist->cnt[tid[3]] += cnt[3];
    #ifdef SINGLE_PREC
    hist->cnt[tid[4]] += cnt[4];
    hist->cnt[tid[5]] += cnt[5];
    hist->cnt[tid[6]] += cnt[6];
    hist->cnt[tid[7]] += cnt[7];
    #endif
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
    /* Gather values from the lookup table. */
    vec_i idx = SIMD_MSK_R2I_GATHER_INT(SIMD_ZEROS_INT(), mask, vidx, tab,
        sizeof(BENCHMARK_HIST_TABLE_DTYPE));
    idx = SIMD_AND_INT(idx, SIMD_SET1_INT(BENCHMARK_HIST_TABLE_MASK));
    idx = SIMD_ADD_INT(idx, ioff);
    #ifdef      SINGLE_PREC
    /* Deal with higher and lower half lanes separately. */
    __m256i hidx = _mm512_castsi512_si256(idx);
    __mmask8 hmask = mask & 0xFF;
    vec_i cnt = SIMD_MSK_GATHER_I64(SIMD_ZEROS_INT(), hmask, hidx,
        hist->cnt, 8);
    cnt = SIMD_ADD_I64(cnt, SIMD_SET1_I64(1));
    SIMD_MSK_SCATTER_I64(hist->cnt, hmask, hidx, cnt, 8);

    hidx = _mm512_extracti64x4_epi64(idx, 1);
    hmask = mask >> 8;
    cnt = SIMD_MSK_GATHER_I64(SIMD_ZEROS_INT(), hmask, hidx, hist->cnt, 8);
    cnt = SIMD_ADD_I64(cnt, SIMD_SET1_I64(1));
    SIMD_MSK_SCATTER_I64(hist->cnt, hmask, hidx, cnt, 8);
    #else    /* DOUBLE_PREC */
    vec_i cnt = SIMD_MSK_GATHER_I64(SIMD_ZEROS_INT(), mask, idx, hist->cnt, 8);
    cnt = SIMD_ADD_INT(cnt, SIMD_SET1_INT(1));
    SIMD_MSK_SCATTER_I64(hist->cnt, mask, idx, cnt, 8);
    #endif
  #endif
  }
  /* Deal with remainders: scalar code. */
  for (size_t n = npar; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      const size_t idx = (size_t) d2[n] - (size_t) d2min;
      hist->cnt[tab[idx]]++;
    }
  }
#endif
}

/******************************************************************************
Function `hist_table_hybrid_<BENCHMARK_HIST_TABLE_DTYPE>`:
  Update distance histogram using the lookup table for non-integer squared
  distance bin edges.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void HIST_TABLE_FUNC(hist_table_hybrid, BENCHMARK_HIST_TABLE_DTYPE) (
    const real *d2, const size_t nsp, HIST *hist) {
  BENCHMARK_HIST_TABLE_DTYPE *tab = (BENCHMARK_HIST_TABLE_DTYPE *) hist->tab;
  const real d2min = hist->d2[0];
  const real d2max = hist->d2[hist->n];
#if     BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  const int offset = d2min;
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      int idx = tab[(int) d2[n] - offset];
      if (idx < hist->n) hist->cnt[idx]++;
      else {
        idx -= hist->n;
        while (idx != 0 && d2[n] < hist->d2[idx]) idx--;
        hist->cnt[idx]++;
      }
      /* For some reason the following codes are slower. */
      /*
      if (idx >= hist->n) {
        idx -= hist->n;
        while (idx != 0 && d2[n] < hist->d2[idx]) idx--;
      }
      hist->cnt[idx]++;
      */
    }
  }
#else
  const vec_r vmin = SIMD_SET1_REAL(d2min);
  const vec_r vmax = SIMD_SET1_REAL(d2max);
  #if   defined(SINGLE_PREC)  &&  BENCHMARK_SIMD == BENCHMARK_SIMD_AVX
  const vec_r offs = _mm256_floor_ps(vmin);
  #else
  const vec_r2i offs = SIMD_FLOOR_REAL(vmin);
  #endif
  #if   BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512
  /* Index offsets of private histograms for avoiding conflicts. */
    #ifdef SINGLE_PREC
  const vec_i ioff = _mm512_setr_epi32(0, hist->n, hist->n * 2, hist->n * 3,
      hist->n * 4, hist->n * 5, hist->n * 6, hist->n * 7,
      0, hist->n, hist->n * 2, hist->n * 3,
      hist->n * 4, hist->n * 5, hist->n * 6, hist->n * 7);
    #else
  const vec_i ioff = _mm512_setr_epi64(0, hist->n, hist->n * 2, hist->n * 3,
      hist->n * 4, hist->n * 5, hist->n * 6, hist->n * 7);
    #endif
  #endif
  #if   BENCHMARK_SIMD  >=  BENCHMARK_SIMD_AVX2
  const vec_i nbin = SIMD_SET1_INT(hist->n);
  #endif
  const size_t npar = BENCHMARK_NUM_MASK & nsp;
  for (size_t n = 0; n != npar; n += BENCHMARK_NUM_REAL) {
    vec_r x = SIMD_LOAD_REAL(d2 + n);
    vec_msk mask = SIMD_CMP_REAL(x, vmin, _CMP_GE_OQ);
  #if           BENCHMARK_SIMD  <  BENCHMARK_SIMD_AVX512
    vec_msk mask1 = SIMD_CMP_REAL(x, vmax, _CMP_LT_OQ);
    mask = SIMD_AND_REAL(mask, mask1);
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
    mask = SIMD_MSK_CMP_REAL(mask, x, vmax, _CMP_LT_OQ);
  #endif
  #if   defined(SINGLE_PREC)  &&  BENCHMARK_SIMD == BENCHMARK_SIMD_AVX
    vec_r dx = SIMD_SUB_REAL(_mm256_floor_ps(x), offs);
    vec_r2i vidx = _mm256_cvtps_epi32(dx);
  #else
    vec_r2i vidx = SIMD_FLOOR_REAL(x);
    vidx = SIMD_SUB_R2I(vidx, offs);
  #endif
  #if   BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX
    const int imask = SIMD_MVMSK_REAL(mask);
    int32_t *idx = (int32_t *) &vidx;
    real *val = (real *) &x;
    if (imask & 0x1) {
      int tid = tab[idx[0]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[0] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    if (imask & 0x2) {
      int tid = tab[idx[1]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[1] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    if (imask & 0x4) {
      int tid = tab[idx[2]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[2] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    if (imask & 0x8) {
      int tid = tab[idx[3]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[3] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    #ifdef SINGLE_PREC
    if (imask & 0x10) {
      int tid = tab[idx[4]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[4] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    if (imask & 0x20) {
      int tid = tab[idx[5]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[5] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    if (imask & 0x40) {
      int tid = tab[idx[6]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[6] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    if (imask & 0x80) {
      int tid = tab[idx[7]];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && val[7] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
    #endif
  #elif BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX2
    /* Gather values from the lookup table. */
    #if defined(__GNUC__) && !defined(SINGLE_PREC)
    /* `__int64` is interpreted as `long long int` in gcc. */
    vec_i idx = SIMD_MSK_R2I_GATHER_INT(SIMD_ZEROS_INT(), (long long *) tab,
        vidx, SIMD_CAST_R2I(mask), sizeof(BENCHMARK_HIST_TABLE_DTYPE));
    #else
    vec_i idx = SIMD_MSK_R2I_GATHER_INT(SIMD_ZEROS_INT(), (int_wr *) tab,
        vidx, SIMD_CAST_R2I(mask), sizeof(BENCHMARK_HIST_TABLE_DTYPE));
    #endif
    idx = SIMD_AND_INT(idx, SIMD_SET1_INT(BENCHMARK_HIST_TABLE_MASK));
    vec_i imask = SIMD_CMPGT_INT(nbin, idx);
    imask = SIMD_ANDNOT_INT(imask, SIMD_CAST_R2I(mask));
    /* Comparisons for non-matched entries. */
    if (SIMD_MVMSK_INT(imask) != 0) {
      idx = SIMD_SUB_INT(idx, SIMD_AND_INT(imask, nbin));
      for (;;) {
        vec_r ref = SIMD_MSK_GATHER_REAL(x, hist->d2, idx,
            SIMD_CAST_I2R(imask), sizeof(real));
        imask = SIMD_CAST_R2I(SIMD_CMP_REAL(x, ref, _CMP_LT_OQ));
        if (SIMD_MVMSK_INT(imask) == 0) break;
        idx = SIMD_SUB_INT(idx, SIMD_AND_INT(imask, SIMD_SET1_INT(1)));
      }
    }
    /* Convert the mask to numbers for histogram update. */
    mask = SIMD_AND_REAL(mask, SIMD_CAST_I2R(SIMD_SET1_INT(1)));
    int_wr *tid = (int_wr *) &idx;
    uint_wr *cnt = (uint_wr *) &mask;
    hist->cnt[tid[0]] += cnt[0];
    hist->cnt[tid[1]] += cnt[1];
    hist->cnt[tid[2]] += cnt[2];
    hist->cnt[tid[3]] += cnt[3];
    #ifdef SINGLE_PREC
    hist->cnt[tid[4]] += cnt[4];
    hist->cnt[tid[5]] += cnt[5];
    hist->cnt[tid[6]] += cnt[6];
    hist->cnt[tid[7]] += cnt[7];
    #endif
  #else      /* BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512 */
    /* Gather values from the lookup table. */
    vec_i idx = SIMD_MSK_R2I_GATHER_INT(SIMD_ZEROS_INT(), mask, vidx, tab,
        sizeof(BENCHMARK_HIST_TABLE_DTYPE));
    idx = SIMD_AND_INT(idx, SIMD_SET1_INT(BENCHMARK_HIST_TABLE_MASK));
    vec_msk imask = SIMD_MSK_CMPGE_INT(mask, idx, nbin);
    /* Comparisons for non-matched entries. */
    if (imask != 0) {
      idx = SIMD_MSK_SUB_INT(idx, imask, idx, nbin);
      for (;;) {
        vec_r ref = SIMD_MSK_GATHER_REAL(x, imask, idx, hist->d2, sizeof(real));
        imask = SIMD_CMP_REAL(x, ref, _CMP_LT_OQ);
        if (imask == 0) break;
        idx = SIMD_MSK_SUB_INT(idx, imask, idx, SIMD_SET1_INT(1));
      }
    }
    idx = SIMD_ADD_INT(idx, ioff);
    /* Update histogram. */
    #ifdef      SINGLE_PREC
    /* Deal with higher and lower half lanes separately. */
    __m256i hidx = _mm512_castsi512_si256(idx);
    __mmask8 hmask = mask & 0xFF;
    vec_i cnt = SIMD_MSK_GATHER_I64(SIMD_ZEROS_INT(), hmask, hidx,
        hist->cnt, 8);
    cnt = SIMD_ADD_I64(cnt, SIMD_SET1_I64(1));
    SIMD_MSK_SCATTER_I64(hist->cnt, hmask, hidx, cnt, 8);

    hidx = _mm512_extracti64x4_epi64(idx, 1);
    hmask = mask >> 8;
    cnt = SIMD_MSK_GATHER_I64(SIMD_ZEROS_INT(), hmask, hidx, hist->cnt, 8);
    cnt = SIMD_ADD_I64(cnt, SIMD_SET1_I64(1));
    SIMD_MSK_SCATTER_I64(hist->cnt, hmask, hidx, cnt, 8);
    #else    /* DOUBLE_PREC */
    vec_i cnt = SIMD_MSK_GATHER_I64(SIMD_ZEROS_INT(), mask, idx, hist->cnt, 8);
    cnt = SIMD_ADD_INT(cnt, SIMD_SET1_INT(1));
    SIMD_MSK_SCATTER_I64(hist->cnt, mask, idx, cnt, 8);
    #endif
  #endif
  }
  /* Deal with remainders: scalar code. */
  for (size_t n = npar; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      int tid = tab[(int) d2[n] - (int) d2min];
      if (tid < hist->n) hist->cnt[tid]++;
      else {
        tid -= hist->n;
        while (tid != 0 && d2[n] < hist->d2[tid]) tid--;
        hist->cnt[tid]++;
      }
    }
  }
#endif
}


#undef BENCHMARK_HIST_TABLE_WIDTH
#undef BENCHMARK_HIST_TABLE_DTYPE
#undef BENCHMARK_HIST_TABLE_MASK

#endif

