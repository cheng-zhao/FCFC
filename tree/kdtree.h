/*******************************************************************************
* kdtree.h: this file is part of the FCFC program.

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

#ifndef __KDTREE_H__
#define __KDTREE_H__

#include "define.h"
#include <stdio.h>

/*============================================================================*\
                        Definition of the tree structure
\*============================================================================*/

typedef struct kdtree_struct {
  size_t n;                             /* number of objects       */
  DATA *data;                           /* pointer to the data     */
  DATA min;                             /* lower corner of the box */
  DATA max;                             /* upper corner of the box */
  struct kdtree_struct *left;           /* left child              */
  struct kdtree_struct *right;          /* right child             */
} KDT;

/*============================================================================*\
                        Interfaces of the tree structure
\*============================================================================*/

/******************************************************************************
Function `kdtree_build`:
  Construct the k-D tree from a data set.
Arguments:
  * `data`:     the input dataset;
  * `ndata`:    number of elements of the input dataset;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `err`:      error indicator.
Return:
  Root of the constructed tree.
******************************************************************************/
KDT* kdtree_build(DATA *data, const size_t ndata, DATA *buf, int *err);

/******************************************************************************
Function `kdtree_free`:
  Release the memory allocated for the k-D tree.
Arguments:
  * `node`:     pointer to a root of the k-D tree.
******************************************************************************/
void kdtree_free(KDT *node);

#endif
