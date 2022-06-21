/*******************************************************************************
* benchmark/common/mt19937.c: this file is part of the FCFC program.
 
* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "mt19937.h"
#include <math.h>

/*******************************************************************************
  Implementation of the Mersenne Twister 19937 random number generator.
  ref: https://doi.org/10.1145%2F272991.272995

  The initialisation with a seed is performed following the 2002 version:
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html

  This variant has been written from scratch, but produces the same
  sequences as the version provided by the authors of the algorithm:
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*******************************************************************************/

#ifndef MT19937_N
#define N               624
#else
#define N               MT19937_N
#endif
#define M               397
#define K               19937
#define MA              0x9908b0dfUL
#define UPPER_MASK(x)   (0x80000000UL & x)      /* most significant w-r bits */
#define LOWER_MASK(x)   (0x7fffffffUL & x)      /* least significant r bits */
#define NORM            0x1p-32                 /* 2^{-32} */
#define NORM2           0x1.fffffffep-33        /* 1 / (2^{32} + 1) */
#define TWO_PI          0x1.921fb54442d18p+2    /* 2 * M_PI */
#define DEFAULT_SEED    1

#ifndef MT19937_UNROLL
  #define MAGIC(y)      (((y) & 1) ? MA : 0)
#else
/* The following trick is taken from
 * https://github.com/cslarsen/mersenne-twister
 * which is distributed under the modified BSD license.
 * It is supposed to work faster with `-ftree-vectorize`, especially for
 * machines with SSE/AVX. */
  #define MAGIC(y)      (((((int32_t)y) << 31) >> 31) & MA)
#endif

/******************************************************************************
Function `mt19937_seed`:
  Initialise the state with an integer.
Arguments:
  * `state`:    the state to be intialised;
  * `seed`:     a positive integer for the initialisation.
******************************************************************************/
void mt19937_seed(mt19937_t *state, uint64_t seed) {
  int i;
  /* Initialise the statee with LCG.
   * Bit truncation with 0xffffffffUL is not necessary for 32-bit statees. */
  state->mt[0] = seed;  /* validation of seed is done in `mt19937_init' */
  for (i = 1; i < N; i++)
    state->mt[i] = 1812433253UL * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i;

  state->idx = i;
  state->g = HUGE_VAL;
}

/******************************************************************************
Function `mt19937_get`:
  Generate an integer and update the state.
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random integer.
******************************************************************************/
uint64_t mt19937_get(mt19937_t *state) {
  uint32_t y;

  if (state->idx >= N) {         /* generate N words at one time */
    int k;
    for (k = 0; k < N - M; k++) {
      y = UPPER_MASK(state->mt[k]) | LOWER_MASK(state->mt[k+1]);
      state->mt[k] = state->mt[k+M] ^ (y >> 1) ^ MAGIC(y);
    }
    for (; k < N - 1; k++) {
      y = UPPER_MASK(state->mt[k]) | LOWER_MASK(state->mt[k+1]);
      state->mt[k] = state->mt[k+M-N] ^ (y >> 1) ^ MAGIC(y);
    }
    y = UPPER_MASK(state->mt[N-1]) | LOWER_MASK(state->mt[0]);
    state->mt[N-1] = state->mt[M-1] ^ (y >> 1) ^ MAGIC(y);
    state->idx = 0;
  }

  /* tempering */
  y = state->mt[state->idx];

  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  state->idx += 1;
  return y;
}

/******************************************************************************
Function `mt19937_get_double`:
  Generate a double-precision floating-point number from the integer.
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
double mt19937_get_double(mt19937_t *state) {
  return mt19937_get(state) * NORM;
}

/******************************************************************************
Function `mt19937_get_gauss`:
  Generate a Gaussian random number.
Arguments:
  * `state`:    the state for the generator;
  * `cen`:      centre of the Gaussian distribution;
  * `sig`:      sigma of the Gaussian distribution.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
double mt19937_get_gauss(mt19937_t *state, const double cen, const double sig) {
  if (state->g != HUGE_VAL) {
    double v = state->g;
    state->g = HUGE_VAL;
    return v;
  }

  /* Generate a pair of uniform random numbers in (0,1). */
  double x = (mt19937_get(state) + 1) * NORM2;
  double y = (mt19937_get(state) + 1) * NORM2;

  double u = sqrt(-2 * log(x)) * cos(TWO_PI * y) * sig + cen;
  state->g = sqrt(-2 * log(x)) * sin(TWO_PI * y) * sig + cen;

  return u;
}
