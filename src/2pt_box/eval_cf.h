/*******************************************************************************
* 2pt_box/eval_cf.h: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>
 
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
#include "libast.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

/*============================================================================*\
              Data structure for storing information of the 2PCFs
\*============================================================================*/

typedef struct {
  /* frequently accessed variables at the beginning */
  real bsize;           /* side length of the periodic box                */
  real s2min;           /* minimum squared separation (or s_perp)         */
  real s2max;           /* maximum squared separation (or s_perp)         */
  real pmin;            /* minimum pi of interest                         */
  real pmax;            /* maximum pi of interest                         */
  real prec;            /* precision for truncating (squared) distances   */

  real *s2bin;          /* edges of squared separation (or s_perp) bins   */
  int ns;               /* number of separation (or s_perp) bins          */
  real *pbin;           /* edges of pi bins                               */
  int np;               /* number of pi bins                              */
  int nmu;              /* number of mu bins                              */
  size_t nmu2;          /* squared number of mu bins                      */
  size_t *stab;         /* lookup table for separation (or s_perp) bins   */
  size_t *ptab;         /* lookup table for pi bins                       */
  size_t *mutab;        /* lookup table for mu bins                       */


  int nthread;          /* number of threads to be used                   */
  int bintype;          /* binning scheme: iso, smu, or spi               */
  size_t ntot;          /* total number of bins                           */
  real *sbin;           /* edges of separation (or s_perp) bins           */

  int ncat;             /* number of catalogues to be read                */
  const char *label;    /* labels of the input catalogues                 */
  DATA **data;          /* structures for saving the input catalogues     */
  size_t *ndata;        /* number of objects in the input catalogues      */

  int npc;              /* number of pair counts to be evaluated          */
  int *pc_idx[2];       /* pairs to be counted, defined as input indices  */
  const bool *comp_pc;  /* indicate whether to evaluate the pair counts   */
  size_t **cnt;         /* array for storing evaluated pair counts        */
  double *norm;         /* normalisation factors for pair counts          */
  double **ncnt;        /* array for normalised pair counts               */
  double *rr;           /* array for analytical RR counts                 */
  size_t *pcnt;         /* thread-private array for counting in parallel  */
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
Function `cf_init`:
  Initialise the structure for correlation function evaluations.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
CF *cf_init(const CONF *conf);

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
  * `cf`:       structure for correlation function evaluations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int eval_cf(const CONF *conf, CF *cf);

#endif
