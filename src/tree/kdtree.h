/*******************************************************************************
* kdtree.h: this file is part of the FCFC program.

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

#ifndef __KDTREE_H__
#define __KDTREE_H__

#include "define.h"
#include <stddef.h>
#include <stdbool.h>

/*============================================================================*\
                        Definition of the tree structure
\*============================================================================*/

typedef struct kdtree_struct {
  real *cen;                            /* centre of the box       */
  real *wid;                            /* half-width of the box   */
  size_t n;                             /* number of objects       */
  real *x[FCFC_XDIM];                   /* pointers to coordinates */
  real *w;                              /* pointers to weights     */
  struct kdtree_struct *left;           /* left child              */
  struct kdtree_struct *right;          /* right child             */
} KDT;

/*============================================================================*\
                        Interfaces of the tree structure
\*============================================================================*/

/******************************************************************************
Function `create_kdtree`:
  Construct the k-d tree for a input dataset.
Arguments:
  * `x`:        coordinates of the input dataset;
  * `w`:        weights of the input dataset;
  * `ndata`:    number of input data points;
  * `bsize`:    side lengths of the periodic box;
  * `nleaf`:    maximum number of data points on leaf nodes;
  * `nnode`:    number of tree nodes created successfully.
Return:
  Address of the root node on success; NULL on error.
******************************************************************************/
KDT *create_kdtree(real *x[static FCFC_XDIM], real *w, const size_t ndata,
#ifdef FCFC_METRIC_PERIODIC
    const real bsize[static 3],
#endif
    const size_t nleaf, size_t *nnode);

#ifdef MPI
/******************************************************************************
Function `kdtree_broadcast`:
  Broadcast the k-d tree to all MPI tasks.
Arguments:
  * `root`:     pointer to the root of the k-d tree;
  * `nnode`:    total number of tree nodes;
  * `wt`:       indicate whether to broadcast weights;
  * `src`:      the source for the broadcast;
  * `rank`:     ID of MPI task.
******************************************************************************/
void kdtree_broadcast(KDT **root, size_t *nnode, const bool wt, const int src,
    const int rank);
#endif

#if defined(MPI) || defined(OMP)
/******************************************************************************
Function `kdtree_get_nodes`:
  Get all tree nodes at the level with at least the specific number of nodes.
Arguments:
  * `root`:     pointer to the root of the k-d tree;
  * `nmin`:     the minimum number of nodes to be requested;
  * `nnode`:    the number of requested nodes.
Return:
  Address of the node list on success; NULL on error.
******************************************************************************/
KDT **kdtree_get_nodes(const KDT *root, const uint32_t nmin, size_t *nnode);
#endif

/******************************************************************************
Function `kdtree_free`:
  Release the memory allocated for the k-d tree.
Arguments:
  * `root`:     pointer to the root of the k-d tree.
******************************************************************************/
void kdtree_free(KDT *root);

#endif
