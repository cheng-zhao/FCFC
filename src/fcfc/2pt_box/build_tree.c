/*******************************************************************************
* 2pt_box/build_tree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "build_tree.h"
#include "read_file.h"
#include "kdtree.h"
#include "balltree.h"
#include <stdio.h>
#include <stdlib.h>

/*============================================================================*\
                 Functions for tree creation and deconstruction
\*============================================================================*/

/******************************************************************************
Function `tree_create`:
  Construct the tree from an input catalogue for pair counting.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the catalogue to be processed;
  * `para`:     structure for parallelisms.
Return:
  Address of the tree on success; NULL on error.
******************************************************************************/
void *tree_create(const CONF *conf, CF *cf, const int idx
#ifdef MPI
    , const PARA *para
#endif
    ) {
  void *tree = NULL;
  size_t nnode = 0;
  if (!cf) {
    P_ERR("correlation function evaluation has not been initialised\n");
    return NULL;
  }
  DATA *data = cf->data + idx;

#ifdef MPI
  if (para->rank == para->root) {
#endif
    if (!conf) {
      P_ERR("configuration parameters are not loaded\n");
      return NULL;
    }
    if (idx < 0 || idx > conf->ninput) {
      P_ERR("unexpected index of the catalog: %d\n", idx);
      return NULL;
    }

    printf("Construct the tree for catalog '%c' ...", cf->label[idx]);
    if (conf->verbose)  printf("\n");
    fflush(stdout);

    /* Read catalogue from file. */
    const size_t skip = (conf->skip) ? conf->skip[idx] : DEFAULT_ASCII_SKIP;
    const char cmt = (conf->comment) ? conf->comment[idx] :
        DEFAULT_ASCII_COMMENT;
    const char *sel = (conf->sel) ? conf->sel[idx] : NULL;
    if (sel && ((sel[0] == '\'' && sel[1] == '\'') ||
        (sel[0] == '"' && sel[1] == '"')) && sel[2] == '\0') sel = NULL;

    char *col_id[4];
    int ncol = 3;
    for (int i = 0; i < 3; i++) col_id[i] = conf->pos[idx * 3 + i];
    if (conf->has_wt[idx]) {
      ncol = 4;
      col_id[3] = conf->wt[idx];
    }
    else col_id[3] = NULL;
    real **res = NULL;

    int ftype = conf->ftype ? conf->ftype[idx] : DEFAULT_FILE_TYPE;
    switch (ftype) {
      case FCFC_FFMT_ASCII:
        if (read_ascii_data(conf->input[idx], skip, cmt, conf->fmtr[idx],
            col_id, ncol, sel, &res, &data->n, conf->verbose)) return NULL;
        break;
      case FCFC_FFMT_FITS:
#ifdef WITH_CFITSIO
        if (read_fits_data(conf->input[idx], col_id, ncol, sel, &res,
            &data->n, conf->verbose)) return NULL;
        break;
#else
        P_ERR("FITS format is not enabled\n"
            "Please re-compile the code with option -DWITH_CFITSIO\n");
        return NULL;
#endif
      case FCFC_FFMT_HDF5:
#ifdef WITH_HDF5
        if (read_hdf5_data(conf->input[idx], col_id, ncol, sel, &res,
            &data->n, conf->verbose)) return NULL;
        break;
#else
        P_ERR("HDF5 format is not enabled\n"
            "Please re-compile the code with option -DWITH_HDF5\n");
        return NULL;
#endif
      default:
        P_ERR("unexpected format (%d) of file: `%s'\n", conf->ftype[idx],
            conf->input[idx]);
        return NULL;
    }

    data->x[0] = res[0];
    data->x[1] = res[1];
    data->x[2] = res[2];
    data->w = (conf->has_wt[idx]) ? res[3] : NULL;
    free(res);

    /* Rescale the input coordinates if necessary. */
    if (cf->rescale != 1) {
  #ifdef OMP
    #pragma omp parallel for default(none) shared(data,cf)
  #endif
      for (size_t i = 0; i < data->n; i++) {
        data->x[0][i] *= cf->rescale;
        data->x[1][i] *= cf->rescale;
        data->x[2][i] *= cf->rescale;
      }
    }
    /* Process weights. */
    if (conf->has_wt[idx]) {
      double sumw = 0;
  #ifdef OMP
    #pragma omp parallel for reduction(+:sumw) default(none) shared(data)
  #endif
      for (size_t i = 0; i < data->n; i++) sumw += data->w[i];
      data->wt = sumw;
    }
    else if (cf->cat_wt[idx]) {
  #if FCFC_SIMD  ==  FCFC_SIMD_NONE
      if (!(data->w = malloc(data->n * sizeof(real))))
  #else
      if (!(data->w = malloc((data->n + FCFC_NUM_REAL) * sizeof(real))))
  #endif
      {
        P_ERR("failed to allocate memory for the weights\n");
        return NULL;
      }
  #ifdef OMP
    #pragma omp parallel for default(none) shared(data)
  #endif
      for (size_t i = 0; i < data->n; i++) data->w[i] = 1;
      data->wt = (double) data->n;
    }

    /* Construct the tree. */
    switch (cf->treetype) {
      case FCFC_STRUCT_KDTREE:
        tree = create_kdtree(data->x, data->w, data->n, cf->bsize,
            FCFC_KDTREE_LEAF_SIZE, &nnode);
        if (!tree || !nnode) return NULL;
        if (conf->verbose) printf("  k-d tree (%zu nodes) "
            "constructed for the catalog\n", nnode);
        break;
      case FCFC_STRUCT_BALLTREE:
        tree = create_balltree(data->x, data->w, data->n, cf->bsize,
            FCFC_BALLTREE_LEAF_SIZE, &nnode);
        if (!tree || !nnode) return NULL;
        if (conf->verbose) printf("  ball tree (%zu nodes) "
            "constructed for the catalog\n", nnode);
        break;
      default:
        P_ERR("unsupported tree type\n");
        return NULL;
    }
#ifdef MPI
  }

  /* Broadcast the data catalog and tree. */
  switch (cf->treetype) {
    case FCFC_STRUCT_KDTREE:
      kdtree_broadcast((KDT **) (&tree), &nnode, cf->wt[idx], para);
      if (para->rank != para->root) {
        KDT *root = (KDT *) tree;
        cf->data[idx].n = root->n;
        cf->data[idx].x[0] = root->x[0];
        cf->data[idx].x[1] = root->x[1];
        cf->data[idx].x[2] = root->x[2];
        cf->data[idx].w = root->w;
      }
      break;
    case FCFC_STRUCT_BALLTREE:
      balltree_broadcast((BLT **) (&tree), &nnode, cf->wt[idx], para);
      if (para->rank != para->root) {
        BLT *root = (BLT *) tree;
        cf->data[idx].n = root->n;
        cf->data[idx].x[0] = root->x[0];
        cf->data[idx].x[1] = root->x[1];
        cf->data[idx].x[2] = root->x[2];
        cf->data[idx].w = root->w;
      }
      break;
    default:
      P_ERR("unsupported tree type\n");
      return NULL;
  }

  if (para->rank == para->root) {
#endif
    printf(FMT_DONE);
#ifdef MPI
    fflush(stdout);
  }
#endif

  return tree;
}

/******************************************************************************
Function `tree_destroy`:
  Deconstruct a tree used for pair counting.
Arguments:
  * `tree`:     address of the tree;
  * `type`:     type of the tree.
******************************************************************************/
void tree_destroy(void *tree, const int type) {
  if (!tree) return;
  switch (type) {
    case FCFC_STRUCT_KDTREE:
      kdtree_free((KDT *) tree);
      break;
    case FCFC_STRUCT_BALLTREE:
      balltree_free((BLT *) tree);
      break;
    default:
      P_WRN("unsupported tree type: %d\n", type);
  }
}
