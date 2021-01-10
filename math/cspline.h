/*******************************************************************************
* cspline.h: this file is part of the FCFC program.

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

#ifndef __CSPLINE_H__
#define __CSPLINE_H__

#include <stdio.h>

/*******************************************************************************
  Implementation of the "natural" cubic spline interpolation algorithm.
  ref: https://doi.org/10.5281/zenodo.3611922
  see also: https://arxiv.org/abs/2001.09253

  The original source codes are released under a CC0 license by Haysn Hornbeck.
*******************************************************************************/

/*============================================================================*\
                    Interface for cubic spline interpolation
\*============================================================================*/

/******************************************************************************
Function `cspline_ypp`:
  Compute the second derivative of sample points.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `y`:        y coordinates of the sample points;
  * `n`:        number of sample points.
Return:
  The second derivative of `y`.
******************************************************************************/
double *cspline_ypp(const double *x, const double *y, const size_t n);

/******************************************************************************
Function `cspline_eval`:
  Evaluate the cubic spline interpolation.
  It can be further optimised for uniformly spaced sample points.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `y`:        y coordinates of the sample points;
  * `ypp`:      second derivative of y;
  * `xv`:       x coordinate of the value to be evaluated;
  * `i`:        index of the value to be evaluated.
Return:
  The value for the given x coordinate.
******************************************************************************/
double cspline_eval(const double *x, const double *y, const double *ypp,
    const double xv, const size_t i);

#endif

