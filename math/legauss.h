/*******************************************************************************
* legauss.h: this file is part of the FCFC program.

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

#ifndef _LEGAUSS_H_
#define _LEGAUSS_H_

/*******************************************************************************
  Pre-computed abscissas and weights for Legendre-Gauss integration.
  See e.g. https://mathworld.wolfram.com/Legendre-GaussQuadrature.html
*******************************************************************************/

/* Minimum and maximum order for the pre-computed orders. */
#define LEGAUSS_MIN_ORDER       4
#define LEGAUSS_MAX_ORDER       32

/* Caution: the table is represented as 1-D array.
 * Since the root of the Legendre polynomials are symmetric about 0, only
 *   (N + 1) / 2 values are stored for each order N.
 */

/* Starting index of order N. */
#define LEGAUSS_IDX(N)                  \
  ((((N) * (N)) >> 2) - ((LEGAUSS_MIN_ORDER * LEGAUSS_MIN_ORDER) >> 2))

/* Length of order N. */
#define LEGAUSS_LEN(N)                  (((N) + 1) >> 1)

/* Length for non-zero abscissas of order N. */
#define LEGAUSS_LEN_NONZERO(N)          ((N) >> 1)

extern const double legauss_x[LEGAUSS_IDX(LEGAUSS_MAX_ORDER + 1)];
extern const double legauss_w[LEGAUSS_IDX(LEGAUSS_MAX_ORDER + 1)];

#endif
