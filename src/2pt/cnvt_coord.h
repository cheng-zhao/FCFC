/*******************************************************************************
* 2pt/conv_coord.h: this file is part of the FCFC program.

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

#ifndef _CONV_COORD_H_
#define _CONV_COORD_H_

#include "define.h"
#include "load_conf.h"
#include <stdio.h>

/*============================================================================*\
                 Data structure for cubic spline interpolation
\*============================================================================*/

typedef struct {
  size_t nsp;           /* number of sample points, excluding (0,0)      */
  double *z;            /* redshifts                                     */
  double *d;            /* radial comoving distances                     */
  double *ypp;          /* second derivative for spline interpolation    */
} COORD_CNVT;


/*============================================================================*\
                      Interface for coordinate conversion
\*============================================================================*/

/******************************************************************************
Function `cnvt_init`:
  Initialise the structure for coordinate conversion.
Return:
  Address of the structure.
******************************************************************************/
COORD_CNVT *cnvt_init(void);

/******************************************************************************
Function `cnvt_destroy`:
  Deconstruct the structure for coordinate conversion.
Arguments:
  * `cnvt`:     the structure to be deconstrcuted.
******************************************************************************/
void cnvt_destroy(COORD_CNVT *cnvt);

/******************************************************************************
Function `cnvt_coord`:
  Interface for applying coordinate conversion.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for storing the input catalog;
  * `ndata`:    number of elements of the input catalog;
  * `coord`:    structure for redshift to comoving distance interpolation.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cnvt_coord(const CONF *conf, DATA *data, const size_t ndata,
    COORD_CNVT *coord);

#endif
