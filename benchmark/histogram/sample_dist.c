/*******************************************************************************
* benchmark/histogram/sample_dist.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "sample_dist.h"
#include "mt19937.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*============================================================================*\
                 Functions for generating particle coordinates
\*============================================================================*/

/******************************************************************************
Function `sample_dist`:
  Generate an array for squared distances, with the PDF of sqrt(x).
Arguments:
  * `dist`:     the array for squared distances;
  * `dmax`:     the maximum squared distance;
  * `nsp`:      total number of particles to be generated;
  * `seed`:     random seed.
******************************************************************************/
void sample_dist(real *dist, const real dmax, const size_t nsp,
    const int seed) {
  if (!nsp) return;

  mt19937_t ran;
  mt19937_seed(&ran, seed);

  for (size_t i = 0; i < nsp; i++)
    dist[i] = pow(mt19937_get_double(&ran), 0x1.5555555555555p-1) * dmax;

}

