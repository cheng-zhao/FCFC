/*******************************************************************************
* 2pt/cnvt_coord.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include "timsort.h"
#include "cspline.h"
#include "legauss.h"
#include "cnvt_coord.h"
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                     Functions for coordinate interpolation
\*============================================================================*/

/******************************************************************************
Function `read_sample`:
  Read sample points from a file.
Arguments:
  * `fname`:    name of the file for coordinate conversion;
  * `cnvt`:     the structure for coordinate interpolation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_sample(const char *fname, COORD_CNVT *cnvt) {
  if (!fname || !cnvt) return FCFC_ERR_ARG;
  if (read_ascii_table(fname, &cnvt->z, &cnvt->d, &cnvt->nsp))
    return FCFC_ERR_FILE;
  /* Sort the sample points. */
  tim_sort(cnvt->z, cnvt->d, cnvt->nsp);
  return 0;
}

/******************************************************************************
Function `sample_interp`:
  Interpolate the coordinate conversion samples with binary search.
Arguments:
  * `cnvt`:     data structure for storing sample points;
  * `z`:        the redshift to be converted to a comoving distance.
Return:
  The radial comoving distance on success; HUGE_VAL on error.
******************************************************************************/
static inline double sample_interp(const COORD_CNVT *cnvt, const double z) {
  size_t i, l, u;
  i = l = 0;
  u = cnvt->nsp - 1;
  if (z < cnvt->z[l] || z >= cnvt->z[u]) return HUGE_VAL;

  while (l <= u) {
    i = (l + u) >> 1;
    if (cnvt->z[i + 1] <= z) l = i + 1;
    else if (cnvt->z[i] > z) u = i - 1;
    else break;
  }

  return cspline_eval(cnvt->z, cnvt->d, cnvt->ypp, z, i);
}

/******************************************************************************
Function `cnvt_coord_interp`:
  Apply coordinate conversion for a catalog using cubic spline interpolation.
Arguments:
  * `cnvt`:     data structure for storing sample points;
  * `data`:     data strucutre for the coordinates;
  * `ndata`:    number of objects.
******************************************************************************/
static void cnvt_coord_interp(const COORD_CNVT *cnvt, DATA *data,
    const size_t ndata) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    double ra = data[i].x[0] * DEGREE_2_RAD;
    double dec = data[i].x[1] * DEGREE_2_RAD;
    double dist = sample_interp(cnvt, data[i].x[2]);

    data[i].x[0] = dist * cos(dec) * cos(ra);
    data[i].x[1] = dist * cos(dec) * sin(ra);
    data[i].x[2] = dist * sin(dec);
  }
}


/*============================================================================*\
                      Functions for coordinate integration
\*============================================================================*/

/******************************************************************************
Function `cnvt_z_sample`:
  Sample redshifts uniformly in the redshift range of a catalog.
Arguments:
  * `data`:     structure for storing the input catalog;
  * `ndata`:    number of elements of the input catalog;
  * `num`:      number of sample points for the redshifts.
Return:
  Address of the sampled redshift array.
******************************************************************************/
static double *cnvt_z_sample(const DATA *data, const size_t ndata,
    const int num) {
  if (num < 1) return NULL;
  double zmax = -DBL_MAX;
  double zmin = DBL_MAX;

#ifdef OMP
  const int nomp = omp_get_max_threads();
  double *pmax;                 /* thread-private maximum redshift */
  double *pmin;                 /* thread-private minimum redshift */
  if (!(pmax = malloc(nomp * sizeof(double)))) {
    P_ERR("failed to allocate memory for thread-private redshift\n");
    return NULL;
  }
  if (!(pmin = malloc(nomp * sizeof(double)))) {
    P_ERR("failed to allocate memory for thread-private redshift\n");
    free(pmax); return NULL;
  }
  for (int i = 0; i < nomp; i++) {
    pmax[i] = -DBL_MAX;
    pmin[i] = DBL_MAX;
  }

  /* Determine the minimum and maximum redshift of the samples. */
#pragma omp parallel num_threads(nomp)
  {
    const int tid = omp_get_thread_num();
#pragma omp for
#endif
    for (size_t n = 0; n < ndata; n++) {
      double z = data[n].x[2];
      if (z < 0) {
        P_ERR("invalid negative redshift in the data catalog:\n"
            "(%g, %g, %g)\n", data[n].x[0], data[n].x[1], data[n].x[2]);
#ifdef OMP
        exit(FCFC_ERR_DATA);
#else
        return NULL;
#endif
      }

#ifdef OMP
      if (pmax[tid] < z) pmax[tid] = z;
      if (pmin[tid] > z) pmin[tid] = z;
#else
      if (zmax < z) zmax = z;
      if (zmin > z) zmin = z;
#endif
    }
#ifdef OMP
  }
  /* Gather the largest redshift from all threads. */
  for (int i = 0; i < nomp; i++) {
    if (zmax < pmax[i]) zmax = pmax[i];
    if (zmin > pmin[i]) zmin = pmin[i];
  }
  free(pmax);
  free(pmin);
#endif

  if (zmin > zmax) {
    P_ERR("invalid redshift value in the catalogs\n");
    return NULL;
  }

  /* Generate sample points. */
  double *zsp = malloc(num * sizeof(double));
  if (!zsp) {
    P_ERR("failed to allocate memory for sample points of redshifts\n");
    return NULL;
  }
  for (int i = 0; i < num; i++)
    zsp[i] = zmin + i * (zmax - zmin) / (num - 1);

  return zsp;
}

/******************************************************************************
Function `cnvt_integrand`:
  Integrand for the redshift to radial comoving distance conversion.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `z`:        the redshift to be converted to comoving distance;
Return:
  The integrand for comoving distance integration.
******************************************************************************/
static inline double cnvt_integrand(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const double z) {
  double z1 = z + 1;
  double z2 = z1 * z1;
  double d = Omega_m * z2 * z1;
  if (Omega_k) d += Omega_k * z2;
  if (widx) d += Omega_L * pow(z1, widx);
  else d += Omega_L;
  d = SPEED_OF_LIGHT * 0.01 / sqrt(d);
  return d;
}

/******************************************************************************
Function `cnvt_legauss`:
  Convert redshift to comoving distance using the Legendre-Gauss integration.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `order`:    order of the integration;
  * `z`:        redshift to be converted to comoving distance;
Return:
  The radial comoving distance on success; HUGE_VAL on error.
******************************************************************************/
static inline double cnvt_legauss(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const int order, const double z) {
  /* Variable transformation for integration from 0 to z. */
  double zp = z * 0.5;
  double sum = 0;
  int i;
  for (i = LEGAUSS_IDX(order);
      i < LEGAUSS_IDX(order) + LEGAUSS_LEN_NONZERO(order); i++) {
    /* Look up the abscissas and weights. */
    double x = legauss_x[i];
    double w = legauss_w[i];

    /* Integrate for both positive and negative abscissas. */
    double z1 = zp * (1 + x);
    double z2 = zp * (1 - x);
    sum += w * (cnvt_integrand(Omega_m, Omega_L, Omega_k, widx, z1) +
        cnvt_integrand(Omega_m, Omega_L, Omega_k, widx, z2));
  }
  /* For odd orders, there is also the abscissas x = 0. */
  if (order & 1)
    sum += legauss_w[i] * cnvt_integrand(Omega_m, Omega_L, Omega_k, widx, zp);

  return sum * zp;
}

/******************************************************************************
Function `cnvt_legauss_order`:
  Check the order of Legendre-Gauss integration given the desired precision.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `err`:      maximum allowed error for the integration;
  * `data`:     structure for storing the input catalog;
  * `ndata`:    number of elements of the input catalog;
  * `num`:      the length of the redshift array.
Return:
  The order on success; INT_MAX on error.
******************************************************************************/
static int cnvt_legauss_order(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const double err,
    const DATA *data, const size_t ndata, const int num) {
  /* Uniformly sample redshifts from the data. */
  double *ztest = cnvt_z_sample(data, ndata, num);
  if (!ztest) return INT_MAX;

  int order = 0;
#ifdef OMP
#pragma omp parallel
  {
    int priv = 0;
#pragma omp for
#endif
    for (int i = 0; i < num; i++) {
      double oint, nint = 0;
      int n = LEGAUSS_MIN_ORDER - 1;
      do {
        if (n > LEGAUSS_MAX_ORDER) {
          n = INT_MAX;
          break;
        }
        oint = nint;
        nint = cnvt_legauss(Omega_m, Omega_L, Omega_k, widx, ++n, ztest[i]);
      }
      while (fabs(nint - oint) > nint * err);
#ifdef OMP
      if (priv < n) priv = n;
#else
      if (order < n) order = n;
#endif
    }
#ifdef OMP
#pragma omp critical
    if (order < priv) order = priv;
  }
#endif

  free(ztest);
  return order;
}

/******************************************************************************
Function `cnvt_cata_integr`:
  Apply coordinate conversion for a catalog with integration.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `order`:    order of the integration;
  * `data`:     data strucutre for the coordinates;
  * `ndata`:    number of objects.
******************************************************************************/
static void cnvt_coord_integr(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const int order, DATA *data,
    const size_t ndata) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    double ra = data[i].x[0] * DEGREE_2_RAD;
    double dec = data[i].x[1] * DEGREE_2_RAD;
    double dist =
      cnvt_legauss(Omega_m, Omega_L, Omega_k, widx, order, data[i].x[2]);

    data[i].x[0] = dist * cos(dec) * cos(ra);
    data[i].x[1] = dist * cos(dec) * sin(ra);
    data[i].x[2] = dist * sin(dec);
  }
}


/*============================================================================*\
                      Interfaces for coordinate conversion
\*============================================================================*/

/******************************************************************************
Function `cnvt_init`:
  Initialise the structure for coordinate conversion.
Return:
  Address of the structure.
******************************************************************************/
COORD_CNVT *cnvt_init(void) {
  COORD_CNVT *cnvt = malloc(sizeof *cnvt);
  if (!cnvt) return NULL;
  cnvt->nsp = 0;
  cnvt->z = cnvt->d = cnvt->ypp = NULL;
  return cnvt;
}

/******************************************************************************
Function `cnvt_destroy`:
  Deconstruct the structure for coordinate conversion.
Arguments:
  * `cnvt`:     the structure to be deconstrcuted.
******************************************************************************/
void cnvt_destroy(COORD_CNVT *cnvt) {
  if (!cnvt) return;
  if (cnvt->z) free(cnvt->z);
  if (cnvt->d) free(cnvt->d);
  if (cnvt->ypp) free(cnvt->ypp);
  free(cnvt);
}

/******************************************************************************
Function `cnvt_coord`:
  Interface for applying coordinate conversion.
Arguments:
  * `conf`:     the structure for configurations;
  * `data`:     structure for storing the input catalog;
  * `ndata`:    number of elements of the input catalog;
  * `coord`:    structure for redshift to comoving distance interpolation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cnvt_coord(const CONF *conf, DATA *data, const size_t ndata,
    COORD_CNVT *coord) {
  if (!conf) {
    P_ERR("configuration parameters not loaded\n");
    return FCFC_ERR_CONF;
  }
  if (!data) {
    P_ERR("catalog not read before applying coordinate concersion\n");
    return FCFC_ERR_CNVT;
  }

  /* Apply coordinate conversion. */
  if (conf->fcnvt) {    /* binary search and cubic spline interpolation */
    if (!coord) {
      P_ERR("coordinate interpolation not initialised\n");
      return FCFC_ERR_CNVT;
    }
    /* Setup the interpolation sample. */
    if (!coord->nsp) {
      /* Read the sample from file. */
      if (read_sample(conf->fcnvt, coord)) return FCFC_ERR_FILE;
      /* Compute the second derivative for cubic spline interpolation. */
      if (!(coord->ypp = cspline_ypp(coord->z, coord->d, coord->nsp))) {
        P_ERR("failed to interpolate the sample coordinates for conversion\n");
        return FCFC_ERR_CNVT;
      }
    }
    /* Apply cubic spline interpolation. */
    cnvt_coord_interp(coord, data, ndata);
    if (conf->verbose)
      printf("  Coordinates converted using cubic spline interpolation\n");
  }
  else {                /* Legendre-Gauss integration */
    /* Pre-compute the power index for Lambda. */
    double widx = (conf->dew == -1) ? 0 : 3 * (1 + conf->dew);
    /* Choose the integration order via convergency test. */
    int order = cnvt_legauss_order(conf->omega_m, conf->omega_l, conf->omega_k,
        widx, conf->ecnvt, data, ndata, FCFC_INT_NUM_ZSP);
    if (order == INT_MAX) {
      P_ERR("failed to perform the convergency test for integrations.\n");
      return FCFC_ERR_CNVT;
    }

    /* Apply Legendre-Gauss integration. */
    cnvt_coord_integr(conf->omega_m, conf->omega_l, conf->omega_k, widx, order,
        data, ndata);
    if (conf->verbose)
      printf("  Coordinates converted using Legendre-Gauss integration "
          "with order %d\n", order);
  }

  return 0;
}

