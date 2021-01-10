/*******************************************************************************
* timsort.h: this file is part of the FCFC program.

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

#ifndef __TIMSORT_H__
#define __TIMSORT_H__

#include <stdio.h>

/*******************************************************************************
  Implementation of the timsort algorithm for sorting two arrays, based on
  https://github.com/swenson/sort

  The original source codes are released under the MIT license, by
  Copyright (c) 2010-2019 Christopher Swenson and others as listed in
  https://github.com/swenson/sort/blob/master/CONTRIBUTORS.md
*******************************************************************************/

/*============================================================================*\
                              Interface of timsort
\*============================================================================*/

/******************************************************************************
Function `tim_sort`:
  Sort two arrays with equal length, based on elements of the first one,
  using the timsort algorithm.
Arguments:
  * `x`:        pointer to the first array;
  * `y`:        pointer to the second array;
  * `size`:     size of the input arrays.
******************************************************************************/
void tim_sort(double *x, double *y, const size_t size);

#endif

