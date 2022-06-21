/*******************************************************************************
* 2pt_box/dual_tree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/*******************************************************************************
  Implementation of the dual tree algorithm.
  ref: http://dx.doi.org/10.1007/10849171_5
*******************************************************************************/

/* Macros for the template functions. */
#if defined(FCFC_TREE_TYPE) && defined(FCFC_CNT_TYPE) && \
  defined(FCFC_BIN_TYPE) && defined(FCFC_TAB_TYPE) && \
  defined(FCFC_STAB_WIDTH) && defined(FCFC_PTAB_WIDTH) && \
  defined(FCFC_BIN_SMIN) && defined(FCFC_BIN_PMIN) && defined(FCFC_CNT_WT)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef FCFC_TREE_NAME
  #undef FCFC_TREE_NAME
#endif
#ifdef FCFC_TREE_DTYPE
  #undef FCFC_TREE_DTYPE
#endif
#ifdef FCFC_CNT_NAME
  #undef FCFC_CNT_NAME
#endif
#ifdef FCFC_BIN_NAME
  #undef FCFC_BIN_NAME
#endif
#ifdef FCFC_TAB_NAME
  #undef FCFC_TAB_NAME
#endif
#ifdef FCFC_SWIDTH_NAME
  #undef FCFC_SWIDTH_NAME
#endif
#ifdef FCFC_PWIDTH_NAME
  #undef FCFC_PWIDTH_NAME
#endif
#ifdef FCFC_STAB_DTYPE
  #undef FCFC_STAB_DTYPE
#endif
#ifdef FCFC_PTAB_DTYPE
  #undef FCFC_PTAB_DTYPE
#endif
#ifdef FCFC_STAB_MASK
  #undef FCFC_STAB_MASK
#endif
#ifdef FCFC_PTAB_MASK
  #undef FCFC_PTAB_MASK
#endif
#ifdef FCFC_SMIN_NAME
  #undef FCFC_SMIN_NAME
#endif
#ifdef FCFC_PMIN_NAME
  #undef FCFC_PMIN_NAME
#endif
#ifdef FCFC_WT_NAME
  #undef FCFC_WT_NAME
#endif

#if     FCFC_TREE_TYPE  ==  FCFC_STRUCT_KDTREE
  #define       FCFC_TREE_NAME          kdtree
  #define       FCFC_TREE_DTYPE         KDT
#elif   FCFC_TREE_TYPE  ==  FCFC_STRUCT_BALLTREE
  #define       FCFC_TREE_NAME          balltree
  #define       FCFC_TREE_DTYPE         BLT
#else
  #error unexpected definition of `FCFC_TREE_TYPE`
#endif

#if     FCFC_CNT_TYPE  ==  FCFC_COUNT_AUTO
  #define       FCFC_CNT_NAME           auto
#elif   FCFC_CNT_TYPE  ==  FCFC_COUNT_CROSS
  #define       FCFC_CNT_NAME           cross
#else
  #error unexpected definition of `FCFC_CNT_TYPE`
#endif

#if     FCFC_BIN_TYPE  ==  FCFC_BIN_ISO
  #define       FCFC_BIN_NAME           iso
#elif   FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  #define       FCFC_BIN_NAME           smu
#elif   FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  #define       FCFC_BIN_NAME           spi
#else
  #error unexpected definition of `FCFC_BIN_TYPE`
#endif

#if     FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_INT
  #define       FCFC_TAB_NAME           itab
#elif   FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  #define       FCFC_TAB_NAME           htab
#else
  #error unexpected definition of `FCFC_TAB_TYPE`
#endif

#if     FCFC_STAB_WIDTH  ==  FCFC_LOOKUP_TABLE_W8
  #define       FCFC_SWIDTH_NAME        s8
  #define       FCFC_STAB_DTYPE         uint8_t
  #define       FCFC_STAB_MASK          0xFF
#elif   FCFC_STAB_WIDTH  ==  FCFC_LOOKUP_TABLE_W16
  #define       FCFC_SWIDTH_NAME        s16
  #define       FCFC_STAB_DTYPE         uint16_t
  #define       FCFC_STAB_MASK          0xFFFF
#else
  #error unexpected definition of `FCFC_STAB_WIDTH`
#endif

#if     FCFC_PTAB_WIDTH  ==  FCFC_LOOKUP_TABLE_W8
  #define       FCFC_PWIDTH_NAME        p8
  #define       FCFC_PTAB_DTYPE         uint8_t
  #define       FCFC_PTAB_MASK          0xFF
#elif   FCFC_PTAB_WIDTH  ==  FCFC_LOOKUP_TABLE_W16
  #define       FCFC_PWIDTH_NAME        p16
  #define       FCFC_PTAB_DTYPE         uint16_t
  #define       FCFC_PTAB_MASK          0xFFFF
#else
  #error unexpected definition of `FCFC_PTAB_WIDTH`
#endif

#if     FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_ZERO
  #define       FCFC_SMIN_NAME          _smin0
#elif   FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  #define       FCFC_SMIN_NAME
#else
  #error unexpected definition of `FCFC_BIN_SMIN`
#endif

#if     FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_ZERO
  #define       FCFC_PMIN_NAME          _pmin0
#elif   FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  #define       FCFC_PMIN_NAME
#else
  #error unexpected definition of `FCFC_BIN_PMIN`
#endif

#if     FCFC_CNT_WT  ==  FCFC_COUNT_NO_WT
  #define       FCFC_WT_NAME
#elif   FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
  #define       FCFC_WT_NAME            _wt
#else
  #error unexpected definition of `FCFC_CNT_WT`
#endif

/* Macros for generating function names. */
#ifndef CONCAT_STRING_10
  #define CONCAT_STRING_10(a,b,c,d,e,f,g,h,i,j)                         \
    a##_##b##_##c##_##d##_##e##_##f##_##g##h##i##j
#endif
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c,d,e,f,g,h,i,j)                             \
    CONCAT_STRING_10(a,b,c,d,e,f,g,h,i,j)
#endif

#ifdef PRIVATE_NAME
  #undef PRIVATE_NAME
#endif
#define PRIVATE_NAME(x)                                                 \
  CONCAT_FNAME(x, FCFC_TREE_NAME, FCFC_CNT_NAME, FCFC_BIN_NAME,         \
      FCFC_TAB_NAME, FCFC_SWIDTH_NAME, FCFC_PWIDTH_NAME,                \
      FCFC_SMIN_NAME, FCFC_PMIN_NAME, FCFC_WT_NAME)


/*============================================================================*\
             Data structure and functions for distance evaluations
\*============================================================================*/

/* Structure with the least configuration parameters for pair counting. */
typedef struct {
  real *bsize;          /* side lengths of the periodic box */
  real *hsize;          /* half side lengths of the periodic box */
#if     FCFC_TREE_TYPE  ==  FCFC_STRUCT_BALLTREE
  real smax;            /* maximum s (s_perp) of interest */
#endif
  real s2max;           /* maximum squared s (s_perp) of interest */
#if     FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO
  #if   FCFC_TREE_TYPE  ==  FCFC_STRUCT_BALLTREE
  real smin;            /* minimum s (s_perp) of interest */
  #endif
  real s2min;           /*minimum squared s (s_perp) of interest */
#endif
#if     FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  real pmax;            /* maximum pi of interest */
  #if   FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO
  real pmin;            /* minimum pi of interest */
  #endif
#endif
  FCFC_STAB_DTYPE *stab;        /* s (s_perp) lookup table */
#if     FCFC_BIN_TYPE  !=  FCFC_BIN_ISO  ||  \
        FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  int ns;               /* number of s (s_perp) bins */
#endif
#if     FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  real *s2bin;          /* edges of squared s (s_perp) bins */
#endif
#if     FCFC_BIN_TYPE  ==  FCFC_BIN_SMU
  uint8_t *mutab;       /* mu lookup table */
  int nmu2;             /* squared number of mu bins */
#endif
#if     FCFC_BIN_TYPE  ==  FCFC_BIN_SPI
  FCFC_PTAB_DTYPE *ptab;        /* pi lookup table */
  #if   FCFC_TAB_TYPE  ==  FCFC_LOOKUP_TYPE_HYBRID
  int np;
  real *pbin;           /* edges of pi bins */
  #endif
#endif
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
  int ntot;             /* total number of pair count bins */
#endif
} PRIVATE_NAME(conf);

/* Function for checking the distances between two nodes. */
#if     FCFC_TREE_TYPE  ==  FCFC_STRUCT_KDTREE
  #include "metric_kdtree.c"
#else
  #include "metric_balltree.c"
#endif

/* Functions for distance evaluation and pair count update. */
#include "metric_common.c"


/*============================================================================*\
                       Function for dual tree pair counts
\*============================================================================*/

/******************************************************************************
Function `paircnt_<IDENTIFIERS>`:
  Count pairs of objects given the tree structure.
Arguments:
  * `stack`:    stack for dual nodes;
  * `config`:   configurations for pair counting;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static void PRIVATE_NAME(paircnt) (STACK_DUAL_NODE *stack, const void *config,
    void *cnt) {
  /* Retrieve two nodes from the stack. */
  DUAL_NODE *nodes = stack_pop(stack);
  if (!nodes) return;
  FCFC_TREE_DTYPE *node1 = (FCFC_TREE_DTYPE *) ((*nodes)[0]);
  FCFC_TREE_DTYPE *node2 = (FCFC_TREE_DTYPE *) ((*nodes)[1]);

#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
  if (node1 == node2) {                 /* the same node */
    /* Check the separation range of points on the node. */
    const fcfc_node_dist_e flag =
      PRIVATE_NAME(compare_single_node) (node1, config);
  #if          (FCFC_BIN_SMIN  ==  FCFC_BIN_MIN_NONZERO) ||             \
               (FCFC_BIN_TYPE  ==  FCFC_BIN_SPI  &&                     \
                FCFC_BIN_PMIN  ==  FCFC_BIN_MIN_NONZERO)
    /* Terminate if the separation range is not of interest. */
    if (flag == FCFC_NODE_DIST_OUT) return;
    if (node1->left == NULL)
  #else
    /* Traverse the node if the separation range is fully inside the distance
     * range of interest, or if it is a leaf node. */
    if (flag == FCFC_NODE_DIST_IN || node1->left == NULL)
  #endif
    {
      PRIVATE_NAME(count_single_node) (node1->x,
  #if           FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          node1->w,
  #endif
          node1->n, config, cnt);
    }
    /* Visit children nodes. */
    else {
      /* Avoid double counting for auto pairs. */
      if (node1->left->x[0] <= node2->left->x[0])
        stack_push(stack, node1->left, node2->left);
      if (node1->left->x[0] <= node2->right->x[0])
        stack_push(stack, node1->left, node2->right);
      if (node1->right->x[0] <= node2->left->x[0])
        stack_push(stack, node1->right, node2->left);
      if (node1->right->x[0] <= node2->right->x[0])
        stack_push(stack, node1->right, node2->right);
    }
  }
  else {                                /* node1->x[0] < node2->x[0] */
#endif       /* FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO */

    /* Check the separation range between two nodes. */
#if             FCFC_SIMD  ==  FCFC_SIMD_NONE
    real shift[3] = {0,0,0};
#else        /* FCFC_SIMD  !=  FCFC_SIMD_NONE */
    real shift[3] __attribute__ ((aligned (FCFC_MEMALIGN_BYTE))) = {0,0,0};
#endif
    const fcfc_node_dist_e flag =
      PRIVATE_NAME(compare_dual_node) (node1, node2, config, shift);

    /* The separation range between the two nodes is not of interest. */
    if (flag == FCFC_NODE_DIST_OUT) return;
    /* The separation range is fully inside the distance range of interest,
     * or the nodes are both leaves. */
    else if (flag == FCFC_NODE_DIST_IN) {
      PRIVATE_NAME(count_dual_node) (node1->x, node2->x,
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
          node1->w, node2->w,
#endif
          node1->n, node2->n, shift, config, cnt);
      return;
    }
    if ((node1->left == NULL && node2->left == NULL)) {
      if (flag != FCFC_NODE_DIST_X_WRAP) {
        PRIVATE_NAME(count_dual_node) (node1->x, node2->x,
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
            node1->w, node2->w,
#endif
            node1->n, node2->n, shift, config, cnt);
      }
      else {
        PRIVATE_NAME(count_dual_node_wrap) (node1->x, node2->x,
#if             FCFC_CNT_WT  ==  FCFC_COUNT_WITH_WT
            node1->w, node2->w,
#endif
            node1->n, node2->n, config, cnt);
      }
      return;
    }
    /* Go deeper for non-leaf nodes. */
    else if (node1->left && node2->left) {
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->right->x[0] <= node2->right->x[0])
#endif
        stack_push(stack, node1->right, node2->right);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->right->x[0] <= node2->left->x[0])
#endif
        stack_push(stack, node1->right, node2->left);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->left->x[0] <= node2->right->x[0])
#endif
        stack_push(stack, node1->left, node2->right);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->left->x[0] <= node2->left->x[0])
#endif
        stack_push(stack, node1->left, node2->left);
    }
    else if (node1->left) {
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->right->x[0] <= node2->x[0])
#endif
        stack_push(stack, node1->right, node2);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->left->x[0] <= node2->x[0])
#endif
        stack_push(stack, node1->left, node2);
    }
    else {
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->x[0] <= node2->right->x[0])
#endif
        stack_push(stack, node1, node2->right);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->x[0] <= node2->left->x[0])
#endif
        stack_push(stack, node1, node2->left);
    }

#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
  }
#endif
}


/*============================================================================*\
                              Recycle definitions
\*============================================================================*/

#undef FCFC_TREE_NAME
#undef FCFC_TREE_DTYPE
#undef FCFC_CNT_NAME
#undef FCFC_BIN_NAME
#undef FCFC_TAB_NAME
#undef FCFC_SWIDTH_NAME
#undef FCFC_PWIDTH_NAME
#undef FCFC_STAB_DTYPE
#undef FCFC_PTAB_DTYPE
#undef FCFC_STAB_MASK
#undef FCFC_PTAB_MASK
#undef FCFC_SMIN_NAME
#undef FCFC_PMIN_NAME
#undef FCFC_WT_NAME

#undef FCFC_TREE_TYPE
#undef FCFC_CNT_TYPE
#undef FCFC_BIN_TYPE
#undef FCFC_TAB_TYPE
#undef FCFC_STAB_WIDTH
#undef FCFC_PTAB_WIDTH
#undef FCFC_BIN_SMIN
#undef FCFC_BIN_PMIN
#undef FCFC_CNT_WT

#endif

