/*******************************************************************************
* kdtree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "kdtree.h"
#include <stdlib.h>

/*============================================================================*\
                         Functions for data processing
\*============================================================================*/

/******************************************************************************
Function `kth_compare`:
  Compare the k-th coordinate of two data points in composite space.
  (cf. section 5.5 of de Berg et al. 2008, Computing Geometry 3rd Edition)
Arguments:
  * `a`:        the first data point;
  * `b`:        the second data point;
  * `k`:        the direction of the coordinate to be compared.
Return:
  The order of the two data points.
******************************************************************************/
static inline int kth_compare(const DATA *a, const DATA *b, int k) {
  if (a->x[k] > b->x[k]) return 1;
  if (a->x[k] < b->x[k]) return -1;

  k = (k + 1) % 3;
  if (a->x[k] > b->x[k]) return 1;
  if (a->x[k] < b->x[k]) return -1;

  k = (k + 1) % 3;
  if (a->x[k] > b->x[k]) return 1;
  if (a->x[k] < b->x[k]) return -1;
  return 0;
}

/* Import the function for partition data by the median. */

#ifdef QSELECT_COMPARE
  #undef QSELECT_COMPARE
#endif
#ifdef QSELECT_DTYPE
  #undef QSELECT_DTYPE
#endif

#define QSELECT_COMPARE(a,b,k)  kth_compare(a,b,*((int *)(k)))
#define QSELECT_DTYPE           DATA

#include "quick_select.c"

#undef QSELECT_DTYPE
#undef QSELECT_COMPARE


/*============================================================================*\
                        Functions for tree construction
\*============================================================================*/

/******************************************************************************
Function `kdtree_init`:
  Initialise a node of the k-D tree.
Arguments:
  * `data`:     the dataset for this node;
  * `ndata`:    number of elements of the input dataset.
Return:
  Adress to the k-D tree node.
******************************************************************************/
static KDT *kdtree_init(DATA *data, const size_t ndata) {
  KDT *node = malloc(sizeof(KDT));
  if (!node) return NULL;

  node->n = ndata;
  node->data = data;
  node->left = node->right = NULL;

  return node;
}


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
KDT* kdtree_build(DATA *data, const size_t ndata, DATA *buf, int *err) {
  if (*err) return NULL;
  if (!data || !ndata) {
    *err = FCFC_ERR_ARG;
    return NULL;
  }
  KDT *node = kdtree_init(data, ndata);
  if (!node) {
    *err = FCFC_ERR_MEMORY;
    return NULL;
  }

  /* Compute only the corner of `data` if this is a leaf node. */
  if (ndata <= KDTREE_LEAF_SIZE) {
    for (int k = 0; k < 3; k++) {
      real min, max;
      min = max = data[0].x[k];
      for (size_t i = 0; i < ndata; i++) {
        if (min > data[i].x[k]) min = data[i].x[k];
        if (max < data[i].x[k]) max = data[i].x[k];
      }
      node->min.x[k] = min;
      node->max.x[k] = max;
    }
    return node;
  }

  /* Find the direction with the largest variance, and corners of `data`. */
  int dir = 0;          /* direction with the largest variance */
  real var_max = 0;
  for (int k = 0; k < 3; k++) {
    real mean = 0;
    real min, max;
    min = max = data[0].x[k];

    /* Mean and min/max. */
    for (size_t i = 0; i < ndata; i++) {
      real x = data[i].x[k];
      mean += x;
      if (min > x) min = x;
      if (max < x) max = x;
    }
    mean /= (real) ndata;
    node->min.x[k] = min;
    node->max.x[k] = max;

    real var = 0;
    for (size_t i = 0; i < ndata; i++) {
      real d = data[i].x[k] - mean;
      var += d * d;
    }

    if (var > var_max) {
      dir = k;
      var_max = var;
    }
  }

  size_t n = ndata >> 1;        /* index of the median point */
  /* Split the dataset by the median in the direction with largest variance. */
  qselect(data, n, ndata, buf, &dir);
  node->left = kdtree_build(data, n, buf, err);
  node->right = kdtree_build(data + n, ndata - n, buf, err);

  return node;
}

/******************************************************************************
Function `kdtree_free`:
  Release the memory allocated for the k-D tree.
Arguments:
  * `node`:     pointer to a root of the k-D tree.
******************************************************************************/
void kdtree_free(KDT *node) {
  if(!node) return;
  kdtree_free(node->left);
  kdtree_free(node->right);
  free(node);
}
