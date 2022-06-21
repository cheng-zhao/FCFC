/*******************************************************************************
* 2pt/eval_cf.h: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __EVAL_CF_H__
#define __EVAL_CF_H__

#include "define.h"
#include "load_conf.h"
#include "cnvt_coord.h"
#include "libast.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

/*============================================================================*\
              Data structure for storing information of the 2PCFs
\*============================================================================*/

/* Structure for the input catalog. */
typedef struct {
  size_t n;             /* number of objects          */
  real *x[FCFC_XDIM];   /* coordinates                */
  real *w;              /* weights                    */
  double wt;            /* weighted number of objects */
} DATA;

typedef struct {
  int bintype;          /* binning scheme: iso, smu, or spi               */
  real *s2bin;          /* edges of squared separation (or s_perp) bins   */
  real *p2bin;          /* edges of squared pi bins                       */
  void *stab;           /* lookup table for separation (or s_perp) bins   */
  void *ptab;           /* lookup table for squared pi bins               */
  uint8_t *mutab;       /* lookup table for mu bins                       */
  int swidth;           /* separation (or s_perp) lookup table entry size */
  int pwidth;           /* pi lookup table entry size                     */
  int tabtype;          /* type of the lookup tables (integer vs. hybrid) */
  int ns;               /* number of separation (or s_perp) bins          */
  int np;               /* number of pi bins                              */
  int nmu;              /* number of mu bins                              */
#if FCFC_SIMD  <  FCFC_SIMD_AVX512
  COUNT *pcnt;          /* thread-private array for counting in parallel  */
#else
  void *pcnt;           /* vector-private array for counting in parallel  */
#endif
  int treetype;         /* type of the tree structure                     */

  int ncat;             /* number of catalogues to be read                */
  const char *label;    /* labels of the input catalogues                 */
  DATA *data;           /* structures for saving the input catalogues     */
  real rescale;         /* rescaling factor for the input coordinates     */
  real *sbin_raw;       /* unrescaled separation (or s_perp) bin edges    */
  real *pbin_raw;       /* unrescaled pi bin edges                        */

  int nthread;          /* number of threads to be used                   */
  size_t ntot;          /* total number of bins                           */
  real *sbin;           /* edges of separation (or s_perp) bins           */
  real *pbin;           /* edges of pi bins                               */
  int verbose;          /* indicate whether to show detailed outputs      */

  int npc;              /* number of pair counts to be evaluated          */
  int *pc_idx[2];       /* pairs to be counted, defined as input indices  */
#ifdef MPI
  bool *comp_pc;        /* indicate whether to evaluate the pair counts   */
#else
  const bool *comp_pc;  /* indicate whether to evaluate the pair counts   */
#endif
  bool *wt;             /* indicate whether using weights for pair counts */
  bool *cat_wt;         /* indicate whether using weights in catalogues   */
  const bool *cnvt;     /* indicate whether to run coordinate conversion  */
  COORD_CNVT *coord;    /* structure for coordinate interpolation         */

  COUNT **cnt;          /* array for storing evaluated pair counts        */
  double *norm;         /* normalisation factors for pair counts          */
  double **ncnt;        /* array for normalised pair counts               */
  int ncf;              /* number of correlation functions to be computed */
  char **cf_exp;        /* expression for correlation function estimators */
  ast_t **ast_cf;       /* abstract syntax trees for 2PCF estimators      */
  double **cf;          /* array for storing correlation functions        */

  int nl;               /* number of multipoles to be evaluated           */
  const int *poles;     /* orders of Legendre polynomials to be evaluated */
  double **mp;          /* array for storing 2PCF multipoles              */
  bool comp_wp;         /* indicate whether to compute projected 2PCF     */
  double **wp;          /* array for storing projected 2PCFs              */
} CF;


/*============================================================================*\
                Interfaces for correlation function calculations
\*============================================================================*/

/******************************************************************************
Function `cf_setup`:
  Initialise the structure for correlation function evaluations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `nthread`:  number of OpenMP threads.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
CF *cf_setup(const CONF *conf
#ifdef OMP
    , const int nthread
#endif
    );

#ifdef MPI
/******************************************************************************
Function `cf_setup_worker`:
  Setup the correlation function configurations for MPI workers.
Arguments:
  * `cf`:       structure for correlation function configurations;
  * `rank`:     ID of MPI task.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
void cf_setup_worker(CF **cf, const int rank);
#endif

/******************************************************************************
Function `cf_destroy`:
  Release memory allocated for the correlation function calculator.
Arguments:
  * `cf`:       structure for correlation function evaluations.
******************************************************************************/
void cf_destroy(CF *cf);

/******************************************************************************
Function `eval_cf`:
  Evaluate correlation functions.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations;
  * `ntask`:    number of MPI tasks;
  * `rank`:     ID of MPI task.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int eval_cf(const CONF *conf, CF *cf
#ifdef MPI
    , const int ntask, const int rank
#endif
    );

#endif
