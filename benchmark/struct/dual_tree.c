/*******************************************************************************
* benchmark/struct/dual_tree.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/*******************************************************************************
  Implementation of the dual tree algorithm.
  ref: http://dx.doi.org/10.1007/10849171_5

  Please #include "metric_<treetype>.c" before this file.
*******************************************************************************/

/* Macros for the template functions. */
#if defined(BENCHMARK_TREE_TYPE) && defined(BENCHMARK_PAIRCOUNT_TYPE) && \
  defined(BENCHMARK_BIN_SMIN)

#if !defined(BENCHMARK_TREE_NAME) || !defined(BENCHMARK_TREE_DTYPE) || \
  !defined(BENCHMARK_PAIRCOUNT_NAME) || !defined(BENCHMARK_BIN_SMIN_NAME) || \
  !defined(BENCHMARK_GET_NODES) || !defined(BENCHMARK_LCHILD) || \
  !defined(BENCHMARK_RCHILD) || !defined(BENCHMARK_SELF_NODE) || \
  !defined(BENCHMARK_NOT_LEAF)
  #error please include `metric_<treetype>.c` before `dual_tree.c`
#endif


/*============================================================================*\
                       Function for dual tree pair counts
\*============================================================================*/

/******************************************************************************
Function `pair_<BENCHMARK_TREE_NAME><BENCHMARK_PAIRCOUNT_NAME>
    <BENCHMARK_BIN_SMIN_NAME>`:
  Count pairs of objects given the tree structure.
Arguments:
  * `stack`:    stack for dual nodes;
  * `config`:   configurations for pair counting.
******************************************************************************/
static void PRIVATE_NAME(pair, BENCHMARK_TREE_NAME,
    BENCHMARK_PAIRCOUNT_NAME, BENCHMARK_BIN_SMIN_NAME)
    (STACK_DUAL_NODE *stack, void *config) {

  /* Retrieve two nodes from the stack. */
  DUAL_NODE *nodes = stack_pop(stack);
  if (!nodes) return;
  BENCHMARK_GET_NODES(node1, node2);

  /* Check the separation range between two nodes. */
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
  #if   BENCHMARK_SIMD  ==  BENCHMARK_SIMD_NONE
  real shift[3] = {0,0,0};
  #else
  real shift[3] __attribute__ ((aligned (BENCHMARK_MEMALIGN_BYTE)));
  #endif
  const int flag = PRIVATE_NAME(compare_node_dist, BENCHMARK_TREE_NAME,
      BENCHMARK_PAIRCOUNT_NAME, BENCHMARK_BIN_SMIN_NAME)
      (node1, node2, config, shift);
#else   /* BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_NOBOX */
  const int flag = PRIVATE_NAME(compare_node_dist, BENCHMARK_TREE_NAME,
      BENCHMARK_PAIRCOUNT_NAME, BENCHMARK_BIN_SMIN_NAME)
      (node1, node2, config);
#endif

  /* The separation range between the two nodes is not of interest. */
  if (flag == 1) return;
  /* The separation range is fully inside the distance range of interest,
   * or the nodes are both leaves. */
  if (flag == -1 || (!(BENCHMARK_NOT_LEAF(node1)) &&
      !(BENCHMARK_NOT_LEAF(node2)))) {
    PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
        BENCHMARK_BIN_SMIN_NAME) *conf =
        (PRIVATE_NAME(conf, BENCHMARK_TREE_NAME, BENCHMARK_PAIRCOUNT_NAME,
        BENCHMARK_BIN_SMIN_NAME) *) config;

    PRIVATE_NAME(count_dual_node, BENCHMARK_PAIRCOUNT_NAME,
        BENCHMARK_BIN_SMIN_NAME, )
        (node1->n, node1->x, node2->n, node2->x,
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
        shift,
#endif
#if     BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
        conf->r2min,
#endif
        conf->r2max, &conf->npair);

#ifndef BENCHMARK_TIMING
    conf->ndist += node1->n * node2->n;
#endif
    return;
  }
  /* Go deeper for non-leaf nodes. */
  if (BENCHMARK_NOT_LEAF(node1) && BENCHMARK_NOT_LEAF(node2)) {
    stack_push(stack, BENCHMARK_RCHILD(node1), BENCHMARK_RCHILD(node2));
    stack_push(stack, BENCHMARK_RCHILD(node1), BENCHMARK_LCHILD(node2));
    stack_push(stack, BENCHMARK_LCHILD(node1), BENCHMARK_RCHILD(node2));
    stack_push(stack, BENCHMARK_LCHILD(node1), BENCHMARK_LCHILD(node2));
  }
  else if (BENCHMARK_NOT_LEAF(node1)) {
    stack_push(stack, BENCHMARK_RCHILD(node1), BENCHMARK_SELF_NODE(node2));
    stack_push(stack, BENCHMARK_LCHILD(node1), BENCHMARK_SELF_NODE(node2));
  }
  else {
    stack_push(stack, BENCHMARK_SELF_NODE(node1), BENCHMARK_RCHILD(node2));
    stack_push(stack, BENCHMARK_SELF_NODE(node1), BENCHMARK_LCHILD(node2));
  }
}

#undef BENCHMARK_TREE_NAME
#undef BENCHMARK_TREE_DTYPE
#undef BENCHMARK_PAIRCOUNT_NAME
#undef BENCHMARK_BIN_SMIN_NAME

#undef BENCHMARK_PAIRCOUNT_TYPE
#undef BENCHMARK_BIN_SMIN

#endif

