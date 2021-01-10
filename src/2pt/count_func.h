/*******************************************************************************
* 2pt/count_func.h: this file is part of the FCFC program.

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

#ifndef __COUNT_FUNC_H__
#define __COUNT_FUNC_H__

#include "eval_cf.h"

/*============================================================================*\
                          Interface for pair counting
\*============================================================================*/

/******************************************************************************
Function `kdtree_count_cross_pairs`:
  Count pairs based on the k-D tree data structure.
Arguments:
  * `tree1`:    pointer to the root of the first k-D tree;
  * `tree2`:    pointer to the root of the second k-D tree;
  * `cf`:       structure for congifurations of correlation functions;
  * `cnt`:      array for storing pair counts;
  * `isauto`:   true for counting auto pairs;
  * `usewt`:    apply weights for pair counts.
******************************************************************************/
void count_pairs(const void *tree1, const void *tree2, CF *cf,
    pair_count_t *cnt, bool isauto, bool usewt);

#endif
