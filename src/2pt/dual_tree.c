/*******************************************************************************
* 2pt/dual_tree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/*******************************************************************************
  Implementation of the dual tree algorithm.
  ref: http://dx.doi.org/10.1007/10849171_5
*******************************************************************************/

/* Macros for the template functions. */
#if defined(FCFC_TREE_TYPE) && defined(FCFC_CNT_TYPE) && \
  defined(FCFC_BIN_TYPE) && defined(FCFC_BIN_PREC) && \
  defined(FCFC_BIN_SMIN) && defined(FCFC_BIN_PMIN) && defined(FCFC_CNT_WT)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef FCFC_TREE_TYPE_NAME
  #undef FCFC_TREE_TYPE_NAME
#endif
#ifdef FCFC_TREE_DTYPE
  #undef FCFC_TREE_DTYPE
#endif
#ifdef FCFC_CNT_TYPE_NAME
  #undef FCFC_CNT_TYPE_NAME
#endif
#ifdef FCFC_BIN_TYPE_NAME
  #undef FCFC_BIN_TYPE_NAME
#endif
#ifdef FCFC_BIN_PREC_NAME
  #undef FCFC_BIN_PREC_NAME
#endif
#ifdef FCFC_BIN_SMIN_NAME
  #undef FCFC_BIN_SMIN_NAME
#endif
#ifdef FCFC_BIN_PMIN_NAME
  #undef FCFC_BIN_PMIN_NAME
#endif
#ifdef FCFC_CNT_WT_NAME
  #undef FCFC_CNT_WT_NAME
#endif

#if     FCFC_TREE_TYPE == FCFC_TREE_TYPE_KDTREE
  #define FCFC_TREE_TYPE_NAME   kdtree
  #define FCFC_TREE_DTYPE       KDT
#else
  #error "unexpected definition of `FCFC_TREE_TYPE`"
#endif

#if     FCFC_CNT_TYPE == FCFC_PAIR_COUNT_AUTO
  #define FCFC_CNT_TYPE_NAME    auto
#elif   FCFC_CNT_TYPE == FCFC_PAIR_COUNT_CROSS
  #define FCFC_CNT_TYPE_NAME    cross
#else
  #error "unexpected definition of `FCFC_CNT_TYPE`"
#endif

#if     FCFC_BIN_TYPE == FCFC_BIN_ISO
  #define FCFC_BIN_TYPE_NAME    iso
#elif   FCFC_BIN_TYPE == FCFC_BIN_SMU
  #define FCFC_BIN_TYPE_NAME    smu
#elif   FCFC_BIN_TYPE == FCFC_BIN_SPI
  #define FCFC_BIN_TYPE_NAME    spi
#else
  #error "unexpected definition of `FCFC_BIN_TYPE`"
#endif

#if     FCFC_BIN_PREC == FCFC_BIN_EXACT
  #define FCFC_BIN_PREC_NAME    exact
#elif   FCFC_BIN_PREC == FCFC_BIN_INTEG
  #define FCFC_BIN_PREC_NAME    intbin
#elif   FCFC_BIN_PREC == FCFC_BIN_TRUNC
  #define FCFC_BIN_PREC_NAME    trunc
#else
  #error "unexpected definition of `FCFC_BIN_PREC`"
#endif

#if     FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO
  #define FCFC_BIN_SMIN_NAME    _smin0
#elif   FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
  #define FCFC_BIN_SMIN_NAME
#else
  #error "unexpected definition of `FCFC_BIN_SMIN`"
#endif

#if     FCFC_BIN_TYPE == FCFC_BIN_SPI   /* pi only used for xi(s_perp,pi) */
  #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO
    #define FCFC_BIN_PMIN_NAME  _pmin0
  #elif FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
    #define FCFC_BIN_PMIN_NAME
  #else
    #error "unexpected definition of `FCFC_BIN_PMIN`"
  #endif
#else
  #define FCFC_BIN_PMIN_NAME
#endif

#if     FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
  #define FCFC_CNT_WT_NAME
#elif   FCFC_CNT_WT == FCFC_CNT_WITH_WEIGHT
  #define FCFC_CNT_WT_NAME      _wt
#else
  #error "unexpected definition of `FCFC_CNT_WT`"
#endif

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c,d,e,f,g)   a##_##b##_##c##_##d##e##f##g
#endif

#ifndef PAIR_COUNT_FUNC
  #define PAIR_COUNT_FUNC(a,b,c,d,e,f,g)        CONCAT_FNAME(a,b,c,d,e,f,g)
#endif


/*============================================================================*\
                       Function for dual tree pair counts
\*============================================================================*/

/******************************************************************************
Function `<FCFC_TREE_TYPE>_<FCFC_CNT_TYPE_NAME>_<FCFC_BIN_TYPE_NAME>_
  <FCFC_BIN_PREC_NAME><FCFC_BIN_SMIN_NAME><FCFC_BIN_PMIN_NAME>
  <FCFC_CNT_WT_NAME>`:
  Count pairs of objects given the tree structure.
Arguments:
  * `stack`:    stack of dual nodes to be visited;
  * `cf`:       structure for storing 2PCF configurations;
  * `cnt`:      array for storing pair counts.
******************************************************************************/
static void PAIR_COUNT_FUNC(FCFC_TREE_TYPE_NAME, FCFC_CNT_TYPE_NAME,
    FCFC_BIN_TYPE_NAME, FCFC_BIN_PREC_NAME, FCFC_BIN_SMIN_NAME,
    FCFC_BIN_PMIN_NAME, FCFC_CNT_WT_NAME)
    (STACK_DUAL_NODE *stack, const CF *cf, pair_count_t *cnt) {

  /* Retrieve two nodes from the stack. */
  DUAL_NODE *nodes = stack_pop(stack);
  if (!nodes) return;
  FCFC_TREE_DTYPE *node1 = (FCFC_TREE_DTYPE *) nodes->a;
  FCFC_TREE_DTYPE *node2 = (FCFC_TREE_DTYPE *) nodes->b;

#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO

  if (node1 == node2) {                 /* the same node */
    /* Return if the maximum distance of pairs in this node is outside the
       range of interest. */
#if     FCFC_BIN_TYPE == FCFC_BIN_SPI
    /* There is no good way to check separately the two projected distances. */
  #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO && \
          FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
    if (squared_distance_alter(&node1->min, &node1->max) < cf->s2min) return;
  #elif  FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
    if (squared_distance_alter(&node1->min, &node1->max) < cf->p2min) return;
  #elif  FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
    if (squared_distance_alter(&node1->min, &node1->max) < cf->sp2min) return;
  #endif
#else   /* FCFC_BIN_TYPE == FCFC_BIN_ISO or FCFC_BIN_SMU */
  #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
    if (squared_distance_alter(&node1->min, &node1->max) < cf->s2min) return;
  #endif
#endif
    /* Compare all the pairs if the node is a leaf. */
    if (node1->left == NULL) {
      for (size_t i = 0; i + 1 < node1->n; i++) {
        for (size_t j = i + 1; j < node1->n; j++) {
          DATA *a = node1->data + i;
          DATA *b = node1->data + j;

          /* Check if the distance is inside the range of interest. */
#if     FCFC_BIN_TYPE == FCFC_BIN_SPI
          real pi, sp;
          sp = squared_distance_par_with_pi(a, b, &pi);
  #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
          if (pi >= cf->p2max || pi < cf->p2min) continue;
  #else   /* FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO */
          if (pi >= cf->p2max) continue;
  #endif
  #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          if (sp >= cf->sp2max || sp < cf->sp2min) continue;
  #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          if (sp >= cf->sp2max) continue;
  #endif
#else   /* FCFC_BIN_TYPE == FCFC_BIN_ISO or FCFC_BIN_SMU */
  #if     FCFC_BIN_TYPE == FCFC_BIN_ISO
          real dist = squared_distance(a, b);
  #else   /* FCFC_BIN_TYPE == FCFC_BIN_SMU */
          real pi;
          real dist = squared_distance_with_pi(a, b, &pi);
  #endif
  #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          if (dist >= cf->s2max || dist < cf->s2min) continue;
  #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          if (dist >= cf->s2max) continue;
  #endif
#endif

          /* Find the desired distance bin. */
#if     FCFC_BIN_TYPE == FCFC_BIN_SPI
          size_t idx1, idx2;
  #if     FCFC_BIN_PREC == FCFC_BIN_EXACT
          if ((idx1 = find_dist_bin(pi, cf->p2bin, cf->np)) == SIZE_MAX)
            continue;
          if ((idx2 = find_dist_bin(sp, cf->s2bin, cf->ns)) == SIZE_MAX)
            continue;
  #else   /* FCFC_BIN_PREC == FCFC_BIN_INTEG or FCFC_BIN_TRUNC */
    #if     FCFC_BIN_PREC == FCFC_BIN_INTEG
      #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
          idx1 = (size_t) (pi - cf->p2min);
      #else   /* FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO */
          idx1 = (size_t) pi;
      #endif
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx2 = (size_t) (sp - cf->s2min);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx2 = (size_t) sp;
      #endif
    #else   /* FCFC_BIN_PREC == FCFC_BIN_TRUNC */
      #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
          idx1 = (size_t) ((pi - cf->p2min) * cf->prec);
      #else   /* FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO */
          idx1 = (size_t) (pi * cf->prec);
      #endif
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx2 = (size_t) ((sp - cf->s2min) * cf->prec);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx2 = (size_t) (sp * cf->prec);
      #endif
    #endif
          idx1 = cf->ptab[idx1];
          idx2 = cf->stab[idx2];
  #endif
          /* Increment the pair count. */
  #if   FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
          cnt[idx2 + idx1 * cf->ns].i += 1;
  #else
          cnt[idx2 + idx1 * cf->ns].d += a->w * b->w;
  #endif
#else   /* FCFC_BIN_TYPE == FCFC_BIN_ISO or FCFC_BIN_SMU */
          size_t idx;
  #if     FCFC_BIN_PREC == FCFC_BIN_EXACT
          if ((idx = find_dist_bin(dist, cf->s2bin, cf->ns)) == SIZE_MAX)
            continue;
  #else   /* FCFC_BIN_PREC == FCFC_BIN_INTEG or FCFC_BIN_TRUNC */
    #if     FCFC_BIN_PREC == FCFC_BIN_INTEG
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx = (size_t) (dist - cf->s2min);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx = (size_t) dist;
      #endif
    #else   /* FCFC_BIN_PREC == FCFC_BIN_TRUNC */
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx = (size_t) ((dist - cf->s2min) * cf->prec);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx = (size_t) (dist * cf->prec);
      #endif
    #endif
          idx = cf->stab[idx];
  #endif

  #if     FCFC_BIN_TYPE == FCFC_BIN_ISO
          /* Increment the pair count. */
    #if     FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
          cnt[idx].i += 1;
    #else
          cnt[idx].d += a->w * b->w;
    #endif
  #else   /* FCFC_BIN_TYPE == FCFC_BIN_SMU */
          /* Compute mu and find the corresponding bin. */
          if (dist < REAL_EPS) continue;
          size_t muidx = (size_t) (pi * cf->nmu2 / dist);
          if (muidx >= cf->nmu2) continue;
          muidx = cf->mutab[muidx];
          /* Increment the pair count. */
    #if     FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
          cnt[idx + muidx * cf->ns].i += 1;
    #else
          cnt[idx + muidx * cf->ns].d += a->w * b->w;
    #endif
  #endif
#endif
        }
      }
    }
    /* Check the children if the node is not a leaf. */
    else {
      /* Avoid double counting for auto pairs. */
      if (node1->left->data <= node2->left->data)
        stack_push(stack, node1->left, node2->left);
      if (node1->left->data <= node2->right->data)
        stack_push(stack, node1->left, node2->right);
      if (node1->right->data <= node2->left->data)
        stack_push(stack, node1->right, node2->left);
      if (node1->right->data <= node2->right->data)
        stack_push(stack, node1->right, node2->right);
    }
  }
  else {                                /* node1->data < node2->data */

#endif    /*    FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO    */

    /* Return if the minimum/maximum squared distance of the two nodes is
       outside the range of interest. */
#if     FCFC_BIN_TYPE == FCFC_BIN_SPI
    if (min_squared_dist_between_box(&node1->min, &node1->max,
        &node2->min, &node2->max) >= cf->s2max) return;
  #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO && \
          FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
    if (max_squared_dist_between_box(&node1->min, &node1->max,
        &node2->min, &node2->max) < cf->s2min) return;
  #elif   FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
    if (max_squared_dist_between_box(&node1->min, &node1->max,
        &node2->min, &node2->max) < cf->p2min) return;
  #elif   FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
    if (max_squared_dist_between_box(&node1->min, &node1->max,
        &node2->min, &node2->max) < cf->sp2min) return;
  #endif
#else   /* FCFC_BIN_TYPE == FCFC_BIN_ISO or FCFC_BIN_SMU */
    if (min_squared_dist_between_box(&node1->min, &node1->max,
        &node2->min, &node2->max) >= cf->s2max) return;
  #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
    if (max_squared_dist_between_box(&node1->min, &node1->max,
        &node2->min, &node2->max) < cf->s2min) return;
  #endif
#endif

    /* Compare all the pairs if the nodes are both leaves. */
    if (node1->left == NULL && node2->left == NULL) {
      for (size_t i = 0; i < node1->n; i++) {
        for (size_t j = 0; j < node2->n; j++) {
          DATA *a = node1->data + i;
          DATA *b = node2->data + j;

          /* Check if the distance is inside the range of interest. */
#if     FCFC_BIN_TYPE == FCFC_BIN_SPI
          real pi, sp;
          sp = squared_distance_par_with_pi(a, b, &pi);
  #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
          if (pi >= cf->p2max || pi < cf->p2min) continue;
  #else   /* FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO */
          if (pi >= cf->p2max) continue;
  #endif
  #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          if (sp >= cf->sp2max || sp < cf->sp2min) continue;
  #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          if (sp >= cf->sp2max) continue;
  #endif
#else   /* FCFC_BIN_TYPE == FCFC_BIN_ISO or FCFC_BIN_SMU */
  #if     FCFC_BIN_TYPE == FCFC_BIN_ISO
          real dist = squared_distance(a, b);
  #else   /* FCFC_BIN_TYPE == FCFC_BIN_SMU */
          real pi;
          real dist = squared_distance_with_pi(a, b, &pi);
  #endif
  #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          if (dist >= cf->s2max || dist < cf->s2min) continue;
  #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          if (dist >= cf->s2max) continue;
  #endif
#endif

          /* Find the desired distance bin. */
#if     FCFC_BIN_TYPE == FCFC_BIN_SPI
          size_t idx1, idx2;
  #if     FCFC_BIN_PREC == FCFC_BIN_EXACT
          if ((idx1 = find_dist_bin(pi, cf->p2bin, cf->np)) == SIZE_MAX)
            continue;
          if ((idx2 = find_dist_bin(sp, cf->s2bin, cf->ns)) == SIZE_MAX)
            continue;
  #else   /* FCFC_BIN_PREC == FCFC_BIN_INTEG or FCFC_BIN_TRUNC */
    #if     FCFC_BIN_PREC == FCFC_BIN_INTEG
      #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
          idx1 = (size_t) (pi - cf->p2min);
      #else   /* FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO */
          idx1 = (size_t) pi;
      #endif
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx2 = (size_t) (sp - cf->s2min);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx2 = (size_t) sp;
      #endif
    #else   /* FCFC_BIN_PREC == FCFC_BIN_TRUNC */
      #if     FCFC_BIN_PMIN == FCFC_BIN_MIN_NONZERO
          idx1 = (size_t) ((pi - cf->p2min) * cf->prec);
      #else   /* FCFC_BIN_PMIN == FCFC_BIN_MIN_ZERO */
          idx1 = (size_t) (pi * cf->prec);
      #endif
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx2 = (size_t) ((sp - cf->s2min) * cf->prec);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx2 = (size_t) (sp * cf->prec);
      #endif
    #endif
          idx1 = cf->ptab[idx1];
          idx2 = cf->stab[idx2];
  #endif
          /* Increment the pair count. */
  #if   FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
          cnt[idx2 + idx1 * cf->ns].i += 1;
  #else
          cnt[idx2 + idx1 * cf->ns].d += a->w * b->w;
  #endif
#else   /* FCFC_BIN_TYPE == FCFC_BIN_ISO or FCFC_BIN_SMU */
          size_t idx;
  #if     FCFC_BIN_PREC == FCFC_BIN_EXACT
          if ((idx = find_dist_bin(dist, cf->s2bin, cf->ns)) == SIZE_MAX)
            continue;
  #else   /* FCFC_BIN_PREC == FCFC_BIN_INTEG or FCFC_BIN_TRUNC */
    #if     FCFC_BIN_PREC == FCFC_BIN_INTEG
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx = (size_t) (dist - cf->s2min);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx = (size_t) dist;
      #endif
    #else   /* FCFC_BIN_PREC == FCFC_BIN_TRUNC */
      #if     FCFC_BIN_SMIN == FCFC_BIN_MIN_NONZERO
          idx = (size_t) ((dist - cf->s2min) * cf->prec);
      #else   /* FCFC_BIN_SMIN == FCFC_BIN_MIN_ZERO */
          idx = (size_t) (dist * cf->prec);
      #endif
    #endif
          idx = cf->stab[idx];
  #endif

  #if     FCFC_BIN_TYPE == FCFC_BIN_ISO
          /* Increment the pair count. */
    #if     FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
          cnt[idx].i += 1;
    #else
          cnt[idx].d += a->w * b->w;
    #endif
  #else   /* FCFC_BIN_TYPE == FCFC_BIN_SMU */
          /* Compute mu and find the corresponding bin. */
          if (dist < REAL_EPS) continue;
          size_t muidx = (size_t) (pi * cf->nmu2 / dist);
          if (muidx >= cf->nmu2) continue;
          muidx = cf->mutab[muidx];
          /* Increment the pair count. */
    #if     FCFC_CNT_WT == FCFC_CNT_NO_WEIGHT
          cnt[idx + muidx * cf->ns].i += 1;
    #else
          cnt[idx + muidx * cf->ns].d += a->w * b->w;
    #endif
  #endif
#endif
        }
      }
    }
    /* Check the children of non-leaf nodes. */
    else if (node1->left == NULL) {
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->data <= node2->left->data)
#endif
        stack_push(stack, node1, node2->left);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->data <= node2->right->data)
#endif
        stack_push(stack, node1, node2->right);
    }
    else if (node2->left == NULL) {
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->left->data <= node2->data)
#endif
        stack_push(stack, node1->left, node2);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->right->data <= node2->data)
#endif
        stack_push(stack, node1->right, node2);
    }
    /* Both nodes are not leaves. */
    else {
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->left->data <= node2->left->data)
#endif
        stack_push(stack, node1->left, node2->left);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->left->data <= node2->right->data)
#endif
        stack_push(stack, node1->left, node2->right);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->right->data <= node2->left->data)
#endif
        stack_push(stack, node1->right, node2->left);
#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
      if (node1->right->data <= node2->right->data)
#endif
        stack_push(stack, node1->right, node2->right);
    }

#if             FCFC_CNT_TYPE   ==   FCFC_PAIR_COUNT_AUTO
  }
#endif

}

#undef FCFC_TREE_TYPE_NAME
#undef FCFC_TREE_DTYPE
#undef FCFC_CNT_TYPE_NAME
#undef FCFC_BIN_TYPE_NAME
#undef FCFC_BIN_PREC_NAME
#undef FCFC_BIN_SMIN_NAME
#undef FCFC_BIN_PMIN_NAME
#undef FCFC_CNT_WT_NAME

#undef FCFC_CNT_TYPE
#undef FCFC_BIN_TYPE
#undef FCFC_BIN_PREC
#undef FCFC_BIN_SMIN
#undef FCFC_BIN_PMIN
#undef FCFC_CNT_WT

#endif

