/*******************************************************************************
* benchmark/struct/create_data.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "create_data.h"
#include "mt19937.h"
#include <stdio.h>
#include <stdlib.h>

/*============================================================================*\
                 Functions for generating particle coordinates
\*============================================================================*/

/******************************************************************************
Function `fill_data_box`:
  Generate particle coordinates in a periodic box for the data catalogue.
Arguments:
  * `data`:     structure for the data catalogue.
******************************************************************************/
static inline void fill_data_box(DATA *data) {
  mt19937_t ran;
  mt19937_seed(&ran, BENCHMARK_RAND_SEED);

  for (size_t i = 0; i < data->n; i++) {
    data->x[i] = mt19937_get_double(&ran) * data->bsize;
    data->y[i] = mt19937_get_double(&ran) * data->bsize;
    data->z[i] = mt19937_get_double(&ran) * data->bsize;
  }
}

/******************************************************************************
Function `fill_data_shell`:
  Generate particle coordinates in a shell for the data catalogue, with
  the inner and outer radii being bsize/2 and bsize respectively.
Arguments:
  * `data`:     structure for the data catalogue.
******************************************************************************/
static inline void fill_data_shell(DATA *data) {
  const real d2max = 1;
  const real d2min = 0.25;
  mt19937_t ran;
  mt19937_seed(&ran, BENCHMARK_RAND_SEED);

  register size_t n = 0;
  while (n < data->n) {
    const real x = mt19937_get_double(&ran);
    const real y = mt19937_get_double(&ran);
    const real z = mt19937_get_double(&ran);
    const real dist = x * x + y * y + z * z;
    if (dist < d2min || dist >= d2max) continue;
    data->x[n] = x * data->bsize;
    data->y[n] = y * data->bsize;
    data->z[n++] = z * data->bsize;
  }
}


/*============================================================================*\
                  Interfaces for generating the data catalogue
\*============================================================================*/

/******************************************************************************
Function `create_data`:
  Generate the data catalogue for benchmarking.
Arguments:
  * `bsize`:    side length of the periodic box;
  * `npar`:     total number of particles to be generated;
  * `isbox`:    indicate whether to keep all particles in the periodic box.
Return:
  Address of the data catalogue on success; NULL on error.
******************************************************************************/
DATA *create_data(const double bsize, const size_t npar, const bool isbox) {
  printf("Generating the data catalogue ...");
  fflush(stdout);

  /* Allocate memory. */
  DATA *data = calloc(1, sizeof(DATA));
  if (!data) {
    P_ERR("failed to allocate memory for the data catalog\n");
    return NULL;
  }

  data->x = data->y = data->z = NULL;
  data->n = npar;
  data->bsize = bsize;
  data->isbox = isbox;
#if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  /* Append null data points for safely loading vectors from tree nodes. */
  size_t num = npar + BENCHMARK_MEMALIGN_BYTE / sizeof(real) - 1;
  /* Round up the size of the data array for memory alignment. */
  size_t size = num * sizeof(real);
  size = (size + BENCHMARK_MEMALIGN_BYTE - 1) &
      (~(BENCHMARK_MEMALIGN_BYTE - 1));
  if (posix_memalign((void **) &data->x, BENCHMARK_MEMALIGN_BYTE, size) ||
      posix_memalign((void **) &data->y, BENCHMARK_MEMALIGN_BYTE, size) ||
      posix_memalign((void **) &data->z, BENCHMARK_MEMALIGN_BYTE, size)) {
    P_ERR("failed to allocate aligned memory for the data catalog\n");
    destroy_data(data);
    return NULL;
  }
  /* Initialse the extra data points. */
  num = size / sizeof(real);
  for (size_t i = npar; i < num; i++) {
    data->x[i] = data->y[i] = data->z[i] = 0;
  }
#else
  if (!(data->x = malloc(npar * sizeof(real))) ||
      !(data->y = malloc(npar * sizeof(real))) ||
      !(data->z = malloc(npar * sizeof(real)))) {
    P_ERR("failed to allocate memory for the data catalog\n");
    destroy_data(data);
    return NULL;
  }
#endif

  /* Generate the data catalogue. */
  if (isbox) fill_data_box(data);
  else fill_data_shell(data);

  printf(FMT_DONE);
  return data;
}

/******************************************************************************
Function `destroy_data`:
  Deconstuct the data catalogue.
Arguments:
  * `data`:     pointer to the data catalogue.
******************************************************************************/
void destroy_data(DATA *data) {
  if (!data) return;
  if (data->x) free(data->x);
  if (data->y) free(data->y);
  if (data->z) free(data->z);
  free(data);
}
