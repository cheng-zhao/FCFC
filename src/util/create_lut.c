/*******************************************************************************
* create_lut.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

/*******************************************************************************
  Create the lookup table for updating separation histogram.
*******************************************************************************/

#ifdef FCFC_LOOKUP_TABLE_WIDTH

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef FCFC_LOOKUP_TABLE_DTYPE
  #undef FCFC_LOOKUP_TABLE_DTYPE
#endif

#if     FCFC_LOOKUP_TABLE_WIDTH == FCFC_LOOKUP_TABLE_W8
  #define FCFC_LOOKUP_TABLE_DTYPE       uint8_t
#elif   FCFC_LOOKUP_TABLE_WIDTH == FCFC_LOOKUP_TABLE_W16
  #define FCFC_LOOKUP_TABLE_DTYPE       uint16_t
#else
  #error unexpected definition of `FCFC_LOOKUP_TABLE_WIDTH`
#endif

/* Macros for generating function names. */
#ifndef CONCAT_TABFUNC_NAME
  #define CONCAT_TABFUNC_NAME(a,b)      a##_##b
#endif

#ifndef TABFUNC_NAME
  #define TABFUNC_NAME(a,b)             CONCAT_TABFUNC_NAME(a,b)
#endif


/*============================================================================*\
                      Functions for creating lookup tables
\*============================================================================*/

/******************************************************************************
Function `create_lut_int_<BENCHMARK_HIST_TABLE_DTYPE>`:
  Create the lookup table for integer bin edges.
Arguments:
  * `bins`:     the array for sorted bin edges;
  * `num`:      the number of bins.
Return:
  Address of the lookup table on success; NULL on error.
******************************************************************************/
static FCFC_LOOKUP_TABLE_DTYPE * TABFUNC_NAME(create_lut_int,
    FCFC_LOOKUP_TABLE_DTYPE) (const real *bins, const int num) {
  /* Compute table size and allocate memory. */
  const long min = bins[0];
  const long max = bins[num];
  const long ntab = max - min;
#if     FCFC_SIMD  <  FCFC_SIMD_AVX512
  FCFC_LOOKUP_TABLE_DTYPE *tab = malloc(ntab * sizeof *tab);
#else
  /* Allocate 8 extra bytes for vector gathering. */
  FCFC_LOOKUP_TABLE_DTYPE *tab = calloc(ntab + sizeof(int64_t) / sizeof *tab,
      sizeof *tab);
#endif
  if (!tab) {
    P_ERR("failed to allocate memory for the lookup table\n");
    return NULL;
  }

  /* Fill values. */
  int n = 1;
  for (long i = 0; i < ntab; i++) {
    if (i + min < bins[n]) tab[i] = n - 1;
    else {
      if (++n > num) {
        P_ERR("failed to create the lookup table\n");
        free(tab);
        return NULL;
      }
      i--;
    }
  }

  return tab;
}

/******************************************************************************
Function `create_lut_hybrid_<BENCHMARK_HIST_TABLE_DTYPE>`:
  Create the hybrid lookup table for arbitrary bin edges.
Arguments:
  * `bins`:     the array for sorted bin edges;
  * `num`:      the number of bins.
Return:
  Address of the lookup table on success; NULL on error.
******************************************************************************/
static FCFC_LOOKUP_TABLE_DTYPE * TABFUNC_NAME(create_lut_hybrid,
    FCFC_LOOKUP_TABLE_DTYPE) (const real *bins, const int num) {
  const long min = bins[0];
  const long max = ceil(bins[num]);
  const long ntab = max - min;
#if     FCFC_SIMD  <  FCFC_SIMD_AVX512
  FCFC_LOOKUP_TABLE_DTYPE *tab = malloc(ntab * sizeof *tab);
#else
  /* Allocate 8 extra bytes for vector gathering. */
  FCFC_LOOKUP_TABLE_DTYPE *tab = calloc(ntab + sizeof(int64_t) / sizeof *tab,
      sizeof *tab);
#endif
  if (!tab) {
    P_ERR("failed to allocate memory for the lookup table\n");
    return NULL;
  }

  /* Fill values. */
  int n = 1;
  long edge = bins[n];
  for (long i = 0; i < ntab; i++) {
    if (i + min < edge) tab[i] = n - 1;
    else {
      while (++n <= num) {
        edge = bins[n];
        if (i + min < edge) {
          tab[i] = num + n - 1;
          break;
        }
      }
      if (n > num) {
        for (; i < ntab; i++) tab[i] = num * 2;
        break;
      }
    }
  }

  return tab;
}


/* Clean the definitions. */
#undef FCFC_LOOKUP_TABLE_WIDTH
#undef FCFC_LOOKUP_TABLE_DTYPE

#endif
