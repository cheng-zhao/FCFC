/*******************************************************************************
* benchmark/common/mt19937.h: this file is part of the FCFC program.
 
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

#ifndef __MT19937_H__
#define __MT19937_H__

#include <stdint.h>

/*============================================================================*\
                         Definition of the random state
\*============================================================================*/

#define MT19937_N 624

typedef struct {
  uint32_t mt[MT19937_N];
  int idx;
  double g;     /* for a pre-generated Gaussian random number */
} mt19937_t;


/*============================================================================*\
                    Interface of the random number generator
\*============================================================================*/

/******************************************************************************
Function `mt19937_seed`:
  Initialise the state with an integer.
Arguments:
  * `state`:    the state to be intialised;
  * `seed`:     a positive integer for the initialisation.
******************************************************************************/
void mt19937_seed(mt19937_t *state, uint64_t seed);

/******************************************************************************
Function `mt19937_get`:
  Generate an integer and update the state.
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random integer.
******************************************************************************/
uint64_t mt19937_get(mt19937_t *state);

/******************************************************************************
Function `mt19937_get_double`:
  Generate a double-precision floating-point number from the integer.
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
double mt19937_get_double(mt19937_t *state);

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
double mt19937_get_gauss(mt19937_t *state, const double cen, const double sig);

#endif
