/*******************************************************************************
* benchmark/struct/create_data.h: this file is part of the FCFC program.

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

#ifndef __CREATE_DATA_H__
#define __CREATE_DATA_H__

#include <stddef.h>
#include <stdbool.h>

/*============================================================================*\
                     Data structure for the input particles
\*============================================================================*/

typedef struct {
  real *x;
  real *y;
  real *z;
  size_t n;
  real bsize;
  bool isbox;
} DATA;


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
DATA *create_data(const double bsize, const size_t npar, const bool isbox);

/******************************************************************************
Function `destroy_data`:
  Deconstuct the data catalogue.
Arguments:
  * `data`:     pointer to the data catalogue.
******************************************************************************/
void destroy_data(DATA *data);

#endif
