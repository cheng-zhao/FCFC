/*******************************************************************************
* legpoly.h: this file is part of the FCFC program.

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

#ifndef __LEGPOLY_H__
#define __LEGPOLY_H__

/*============================================================================*\
                    Inline function for Legendre polynomials
\*============================================================================*/

/******************************************************************************
Functions `legpoly`:
  Compute the Legengre polynomial.
Arguments:
  * `ell`:      the order;
  * `x`:        the variable.
Return:
  P_ell (x)
******************************************************************************/
inline double legpoly(const int ell, const double x) {
  if (ell == 0) return 1;
  const double x2 = x * x;
  switch (ell) {
    /* quadrupole and hexadecapole are used more frequently than the others */
    case 2: return 1.5 * x2 - 0.5;
    case 4: return 4.375 * x2 * x2 - 3.75 * x2 + 0.375;
    case 6: return 14.4375 * x2 * x2 * x2 - 19.6875 * x2 * x2
              + 6.5625 * x2 - 0.3125;
    case 1: return x;
    case 3: return 2.5 * x * (x2 - 0.6);
    case 5: return 8.75 * x * (0.9 * x2 * x2 - x2 + 0x1.b6db6db6db6dbp-3);
    default: return 0;  /* should not happen */
  }
}

#endif

