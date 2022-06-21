/*******************************************************************************
* benchmark/struct/fast_math.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/*******************************************************************************
  FastSqrt and FastLog2:
  approximate sqrt and log with polynomials for IEEE floating-point numbers.
*******************************************************************************/

/* Macros for the template functions. */
#if defined(BENCHMARK_HIST_APPROX_POLY)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#if BENCHMARK_HIST_APPROX_POLY < 1 || BENCHMARK_HIST_APPROX_POLY > 5
  #error `BENCHMARK_HIST_APPROX_POLY` must be between 1 and 5
#endif

/* Macros for generating function names. */
#ifndef CONCAT_APPROX_FNAME
  #define CONCAT_APPROX_FNAME(a,b)      a##_poly##b
#endif

#ifndef HIST_APPROX_FUNC
  #define HIST_APPROX_FUNC(a,b)         CONCAT_FNAME(a,b)
#endif


/*============================================================================*\
                   Functions for approximating math functions
\*============================================================================*/

/******************************************************************************
Function `fast_sqrt_poly<BENCHMARK_HIST_APPROX_POLY>`:
  Approximate the square root with polynomials for IEEE floating-point numbers.
Arguments:
  * `x`:        an IEEE floating-point number.
Return:
  The square root of the input number.
******************************************************************************/
static inline real HIST_APPROX_FUNC(fast_sqrt, BENCHMARK_HIST_APPROX_POLY)
    (const real x) {
#ifdef SINGLE_PREC
  union { float f; uint32_t u; } v;
  v.f = x;
  int32_t fexp = (v.u & 0x7F800000) >> 23;
  uint32_t esgn = fexp & 1;
  fexp = ((fexp >> 1) - 63) << 23;
  v.u = (v.u & 0x007FFFFF) | 0x3F800000;
#else
  union { double f; uint64_t u; } v;
  v.f = x;
  int64_t fexp = (v.u & 0x7FF0000000000000ULL) >> 52;
  uint64_t esgn = fexp & 1;
  fexp = ((fexp >> 1) - 511) << 52;
  v.u = (v.u & 0xFFFFFFFFFFFFFULL) | 0x3FF0000000000000;
#endif
  static const real parity[] = {0x1.6a09e667f3bcdp-1, 1};       /* 1/sqrt(2) */
#if     BENCHMARK_HIST_APPROX_POLY  ==  1
  static const real c[] = {0.590162067091, 0.417307599638};
  v.f = parity[esgn] * (c[0] + v.f * c[1]);
#elif   BENCHMARK_HIST_APPROX_POLY  ==  2
  static const real c[] = {0.443450799463, 0.629462265951, -0.072268374800};
  v.f = parity[esgn] * (c[0] + v.f * (c[1] + v.f * c[2]));
#elif   BENCHMARK_HIST_APPROX_POLY  ==  3
  static const real c[] = {0.369956389593, 0.787923233437, -0.182747635081,
      0.024937426665};
  v.f = parity[esgn] * (c[0] + v.f * (c[1] + v.f * (c[2] + v.f * c[3])));
#elif   BENCHMARK_HIST_APPROX_POLY  ==  4
  static const real c[] = {0.323953832031, 0.919711880671, -0.321493515070,
      0.088572551188, -0.010736387782};
  v.f = parity[esgn] * (c[0] + v.f * (c[1] + v.f * (c[2] + v.f * (c[3] +
      v.f * c[4]))));
#elif   BENCHMARK_HIST_APPROX_POLY  ==  5
  static const real c[] = {0.291714273044, 1.034903528189, -0.483647097788,
      0.200993605396, -0.049134792827, 0.0051715622196};
  v.f = parity[esgn] * (c[0] + v.f * (c[1] + v.f * (c[2] + v.f * (c[3] +
      v.f * (c[4] + v.f * c[5])))));
#endif
  v.u += fexp;
  return v.f;
}

/******************************************************************************
Function `fast_log2_poly<BENCHMARK_HIST_APPROX_POLY>`:
  Approximate the base-2 log with polynomials for IEEE floating-point numbers.
Arguments:
  * `x`:        an IEEE floating-point number.
Return:
  The base-2 logarithm of the input number.
******************************************************************************/
static inline real HIST_APPROX_FUNC(fast_log2, BENCHMARK_HIST_APPROX_POLY)
    (const real x) {
#ifdef SINGLE_PREC
  union { float f; uint32_t u; } v;
  v.f = x;
  int32_t fexp = (int32_t) ((v.u & 0x7F800000) >> 23) - 125;
  v.u = (v.u & 0x007FFFFF) | 0x3E800000;
#else
  union { double f; uint64_t u; } v;
  v.f = x;
  int64_t fexp = (int64_t) ((v.u & 0x7FF0000000000000ULL) >> 52) - 1021;
  v.u = (v.u & 0xFFFFFFFFFFFFFULL) | 0x3FD0000000000000ULL;
#endif
#if     BENCHMARK_HIST_APPROX_POLY  ==  1
  static const real c[] = {-2.912271003941, 3.883028005255};
  v.f = fexp + c[0] + v.f * c[1];
#elif   BENCHMARK_HIST_APPROX_POLY  ==  2
  static const real c[] = {-3.646305618839, 7.937962186443, -5.304063780513};
  v.f = fexp + c[0] + v.f * (c[1] + v.f * c[2]);
#elif   BENCHMARK_HIST_APPROX_POLY  ==  3
  static const real c[] = {-4.132236627159, 12.010446070979, -16.334276373956,
      9.688098080131};
  v.f = fexp + c[0] + v.f * (c[1] + v.f * (c[2] + v.f * c[3]));
#elif   BENCHMARK_HIST_APPROX_POLY  ==  4
  static const real c[] = {-4.495757411937, 16.087159747285, -33.128656720129,
      39.837917670150, -19.927310928345};
  v.f = fexp + c[0] + v.f * (c[1] + v.f * (c[2] + v.f * (c[3] + v.f * c[4])));
#elif   BENCHMARK_HIST_APPROX_POLY  ==  5
  static const real c[] = {-4.786121308599, 20.165581939530, -55.692362823516,
      101.321919405564, -102.490527887136, 43.739121207160};
  v.f = fexp + c[0] + v.f * (c[1] + v.f * (c[2] + v.f * (c[3] + v.f *
      (c[4] + v.f * c[5]))));
#endif
  return v.f;
}


/*============================================================================*\
           Functions for histogram update with approximate functions
\*============================================================================*/

/******************************************************************************
Function `hist_fast_sqrt_poly<BENCHMARK_HIST_APPROX_POLY>`:
  Update the distance histogram by truncating the approximate distance.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void HIST_APPROX_FUNC(hist_fast_sqrt, BENCHMARK_HIST_APPROX_POLY)
    (const real *d2, const size_t nsp, HIST *hist) {
  const real d2min = hist->d2[0];
  const real d2max = hist->d2[hist->n];
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      int i = HIST_APPROX_FUNC(fast_sqrt, BENCHMARK_HIST_APPROX_POLY) (d2[n])
          - hist->dmin;
      while (i >= hist->n || d2[n] < hist->d2[i]) i--;
      while (d2[n] >= hist->d2[i + 1]) i++;
      hist->cnt[i]++;
    }
  }
}

/******************************************************************************
Function `hist_fast_log2_poly<BENCHMARK_HIST_APPROX_POLY>`:
  Update the distance histogram by truncating the approximate base-2 logarithm
  squared distance.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void HIST_APPROX_FUNC(hist_fast_log2, BENCHMARK_HIST_APPROX_POLY)
    (const real *d2, const size_t nsp, HIST *hist) {
  const real fac = 0x1.62e42fefa39efp-2 * hist->n /
      log(hist->dmax / hist->dmin);     /* the constant is 0.5 * ln(2) */
  const real d2min = hist->d2[0];
  const real d2max = hist->d2[hist->n];
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      int i = HIST_APPROX_FUNC(fast_log2, BENCHMARK_HIST_APPROX_POLY) (d2[n])
          * fac;
      while (i >= hist->n || d2[n] < hist->d2[i]) i--;
      while (d2[n] >= hist->d2[i + 1]) i++;
      hist->cnt[i]++;
    }
  }
}


#undef BENCHMARK_HIST_APPROX_POLY
#endif

