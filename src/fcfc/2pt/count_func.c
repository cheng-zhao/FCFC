/*******************************************************************************
* 2pt/count_func.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "count_func.h"
#include "kdtree.h"
#include "balltree.h"
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                            Macro for error handling
\*============================================================================*/

#define FCFC_ALLOC_CNT_CONF(conf, dtype)                                \
  dtype *conf = calloc(1, sizeof(dtype));                               \
  if (!conf) {                                                          \
    P_ERR("failed to allocate memory for pair count settings\n");       \
    return FCFC_ERR_MEMORY;                                             \
  }

/*============================================================================*\
                   Data structures for processing dual nodes
\*============================================================================*/

/* Structure for pointers to the two nodes. */
typedef const void * DUAL_NODE[2];

/* Stack for dual nodes to be processed. */
typedef struct {
  DUAL_NODE *nodes;
  size_t size;
  size_t capacity;
} STACK_DUAL_NODE;

/* Enumeration for the spatial relationship between two nodes. */
typedef enum {
  FCFC_NODE_DIST_IN,            /* node separation inside distance range */
  FCFC_NODE_DIST_OUT,           /* node separation outside distance range */
  FCFC_NODE_DIST_X,             /* node separation cross distance range */
} fcfc_node_dist_e;


/*============================================================================*\
                        Functions for stack manipulation
\*============================================================================*/

/******************************************************************************
Function `stack_push`:
  Push an element to the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack;
  * `a`:        the first node to be pushed to the stack;
  * `b`:        the second node to be pushed to the stack.
******************************************************************************/
static void stack_push(STACK_DUAL_NODE *s, const void *a, const void *b) {
  if (s->size >= s->capacity) {
    /* Enlarge the memory allocated for the stack. */
    if (s->capacity) {
      if ((FCFC_STACK_MAX_SIZE >> 1) <= s->capacity) {
        P_EXT("too many elements to be pushed to the dual-tree stack\n");
        exit(FCFC_ERR_MEMORY);
      }
      s->capacity <<= 1;
    }
    else {      /* initialise the stack */
      s->nodes = NULL;
      s->capacity = FCFC_STACK_INIT_SIZE;
    }

    if (s->capacity <= s->size) {
      P_EXT("unable to expand the size of dual-tree stack\n");
      exit(FCFC_ERR_UNKNOWN);
    }

    DUAL_NODE *tmp = realloc(s->nodes, s->capacity * sizeof *tmp);
    if (!tmp) {
      P_EXT("failed to allocate memory for the dual-tree stack\n");
      exit(FCFC_ERR_MEMORY);
    }
    s->nodes = tmp;
  }

  s->nodes[s->size][0] = a;
  s->nodes[s->size++][1] = b;
}

/******************************************************************************
Function `stack_pop`:
  Pop one pair of nodes from the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack.
Return:
  Address of the dual nodes on the top of the stack on success; NULL on error.
******************************************************************************/
static inline DUAL_NODE *stack_pop(STACK_DUAL_NODE *s) {
  if (!s->size) return NULL;
  s->size -= 1;
  return s->nodes + s->size;
}

#if defined(MPI) && defined(OMP)
/******************************************************************************
Function `stack_clear`:
  Clear the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack.
******************************************************************************/
static inline void stack_clear(STACK_DUAL_NODE *s) {
  s->size = 0;
}
#endif

/******************************************************************************
Function `stack_destroy`:
  Deconstruct the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack.
******************************************************************************/
static void stack_destroy(STACK_DUAL_NODE *s) {
  if (!s->capacity) return;
  s->size = s->capacity = 0;
  if (s->nodes) free(s->nodes);
}


/*============================================================================*\
                     Pair counting functions from templates
\*============================================================================*/

/* Clean all the relevant macros first */
#ifdef FCFC_TREE_TYPE
  #undef FCFC_TREE_TYPE
#endif

#ifdef FCFC_CNT_TYPE
  #undef FCFC_CNT_TYPE
#endif
#ifdef FCFC_BIN_TYPE
  #undef FCFC_BIN_TYPE
#endif
#ifdef FCFC_TAB_TYPE
  #undef FCFC_TAB_TYPE
#endif
#ifdef FCFC_STAB_WIDTH
  #undef FCFC_STAB_WIDTH
#endif
#ifdef FCFC_PTAB_WIDTH
  #undef FCFC_PTAB_WIDTH
#endif
#ifdef FCFC_BIN_SMIN
  #undef FCFC_BIN_SMIN
#endif
#ifdef FCFC_BIN_PMIN
  #undef FCFC_BIN_PMIN
#endif
#ifdef FCFC_CNT_WT
  #undef FCFC_CNT_WT
#endif

/*******************************************************************************
                                k-d tree -- auto
*******************************************************************************/

/************************************ ISO ************************************/

/* kdtree_auto_iso_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_iso_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SMU ************************************/

/* kdtree_auto_smu_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_smu_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SPI ************************************/

/* kdtree_auto_spi_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_itab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_auto_spi_htab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"


/*******************************************************************************
                               k-d tree -- cross
*******************************************************************************/

/************************************ ISO ************************************/

/* kdtree_cross_iso_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_iso_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SMU ************************************/

/* kdtree_cross_smu_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_smu_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SPI ************************************/

/* kdtree_cross_spi_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_itab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* kdtree_cross_spi_htab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_KDTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"


/*******************************************************************************
                                k-d tree -- auto
*******************************************************************************/

/************************************ ISO ************************************/

/* balltree_auto_iso_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_iso_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SMU ************************************/

/* balltree_auto_smu_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_smu_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SPI ************************************/

/* balltree_auto_spi_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_itab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_auto_spi_htab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_AUTO
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"


/*******************************************************************************
                               k-d tree -- cross
*******************************************************************************/

/************************************ ISO ************************************/

/* balltree_cross_iso_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_iso_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_ISO
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SMU ************************************/

/* balltree_cross_smu_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_smu_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SMU
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/************************************ SPI ************************************/

/* balltree_cross_spi_itab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_itab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_INT
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s8_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p8_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W8
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_smin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_smin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_NONZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_smin0_pmin0 */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_NO_WT
#include "dual_tree.c"

/* balltree_cross_spi_htab_s16_p16_smin0_pmin0_wt */
#define FCFC_TREE_TYPE          FCFC_STRUCT_BALLTREE
#define FCFC_CNT_TYPE           FCFC_COUNT_CROSS
#define FCFC_BIN_TYPE           FCFC_BIN_SPI
#define FCFC_TAB_TYPE           FCFC_LOOKUP_TYPE_HYBRID
#define FCFC_STAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_PTAB_WIDTH         FCFC_LOOKUP_TABLE_W16
#define FCFC_BIN_SMIN           FCFC_BIN_MIN_ZERO
#define FCFC_BIN_PMIN           FCFC_BIN_MIN_ZERO
#define FCFC_CNT_WT             FCFC_COUNT_WITH_WT
#include "dual_tree.c"


/*============================================================================*\
                          Interface for pair counting
\*============================================================================*/

/******************************************************************************
Function `count_pairs`:
  Count cross pairs based on the k-D tree data structure.
Arguments:
  * `tree1`:    pointer to the root of the first k-D tree;
  * `tree2`:    pointer to the root of the second k-D tree;
  * `cf`:       structure for congifurations of correlation functions;
  * `cnt`:      array for storing pair counts;
  * `isauto`:   true for counting auto pairs;
  * `withwt`:   true for enabling weights;
  * `para`:     structure for parallelisms.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int count_pairs(const void *tree1, const void *tree2, CF *cf, COUNT *cnt,
    const bool isauto, const bool withwt
#ifdef MPI
    , const PARA *para
#endif
    ) {
  void (*paircnt_func) (STACK_DUAL_NODE *, const void *, void *) = NULL;
  void *conf = NULL;

#if     defined(OMP)  ||  FCFC_SIMD  >=  FCFC_SIMD_AVX512
  /* Clean the array for storing private pair counts. */
  #if   FCFC_SIMD  >=  FCFC_SIMD_AVX512
  memset(cf->pcnt, 0, cf->ntot * cf->nthread * FCFC_SIMD_BYTES);
  #else
  memset(cf->pcnt, 0, sizeof(COUNT) * cf->ntot * cf->nthread);
  #endif
#endif
#ifdef MPI
  memset(cnt, 0, sizeof(COUNT) * cf->ntot);
#endif

  /* Choose the optimal pair counting function. */
  bool smin0 = (cf->s2bin[0] == 0);
  if (cf->treetype == FCFC_STRUCT_KDTREE) {
    if (cf->bintype == FCFC_BIN_ISO) {
      if (cf->tabtype == FCFC_LOOKUP_TYPE_INT) {
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_iso_itab_s8_p8_smin0);
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_itab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_itab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_iso_itab_s8_p8);
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
      else {                    /* cf->tabtype == FCFC_LOOKUP_TYPE_HYBRID */
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_iso_htab_s8_p8_smin0);
          cfg->ns = cf->ns;
          cfg->s2bin = cf->s2bin;
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_iso_htab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_iso_htab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_iso_htab_s8_p8);
          cfg->ns = cf->ns;
          cfg->s2bin = cf->s2bin;
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
    }
    else if (cf->bintype == FCFC_BIN_SMU) {
      if (cf->tabtype == FCFC_LOOKUP_TYPE_INT) {
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_smu_itab_s8_p8_smin0);
          cfg->ns = cf->ns;
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_itab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_itab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_smu_itab_s8_p8);
          cfg->ns = cf->ns;
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
      else {                    /* cf->tabtype == FCFC_LOOKUP_TYPE_HYBRID */
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_smu_htab_s8_p8_smin0);
          cfg->ns = cf->ns;
          cfg->s2bin = cf->s2bin;
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_kdtree_auto_smu_htab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_kdtree_cross_smu_htab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_smu_htab_s8_p8);
          cfg->ns = cf->ns;
          cfg->s2bin = cf->s2bin;
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
    }
    else {                      /* cf->bintype == FCFC_BIN_SPI */
      bool pmin0 = (cf->p2bin[0] == 0);
      if (cf->tabtype == FCFC_LOOKUP_TYPE_INT) {
        if (smin0) {
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_smin0_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_smin0_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_kdtree_auto_spi_itab_s8_p8_smin0_pmin0);
            cfg->ns = cf->ns;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_smin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_smin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_spi_itab_s8_p8_smin0);
            cfg->ns = cf->ns;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
        else {                  /* smin != 0 */
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_kdtree_auto_spi_itab_s8_p8_pmin0);
            cfg->ns = cf->ns;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s8_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s8_p16;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_itab_s16_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_itab_s16_p16;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_spi_itab_s8_p8);
            cfg->ns = cf->ns;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2min = cfg->s2min + cfg->p2min;
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
      }
      else {                    /* cf->tabtype == FCFC_LOOKUP_TYPE_HYBRID */
        if (smin0) {
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_smin0_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_smin0_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_kdtree_auto_spi_htab_s8_p8_smin0_pmin0);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_smin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_smin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_spi_htab_s8_p8_smin0);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
        else {                  /* smin != 0 */
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_kdtree_auto_spi_htab_s8_p8_pmin0);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s8_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s8_p16;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_kdtree_auto_spi_htab_s16_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_kdtree_cross_spi_htab_s16_p16;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_kdtree_auto_spi_htab_s8_p8);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->ms2max = FCFC_MIN(cfg->s2max, cfg->p2max);
            cfg->ts2min = cfg->s2min + cfg->p2min;
            cfg->ts2max = cfg->s2max + cfg->p2max;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
      }
    }
  }
  else {                        /* cf->treetype == FCFC_STRUCT_BALLTREE */
    if (cf->bintype == FCFC_BIN_ISO) {
      if (cf->tabtype == FCFC_LOOKUP_TYPE_INT) {
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_iso_itab_s8_p8_smin0);
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_itab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_itab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_iso_itab_s8_p8);
          cfg->smin = cf->sbin[0];
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
      else {                    /* cf->tabtype == FCFC_LOOKUP_TYPE_HYBRID */
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_iso_htab_s8_p8_smin0);
          cfg->ns = cf->ns;
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2bin = cf->s2bin;
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_iso_htab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_iso_htab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_iso_htab_s8_p8);
          cfg->ns = cf->ns;
          cfg->smin = cf->sbin[0];
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2bin = cf->s2bin;
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
    }
    else if (cf->bintype == FCFC_BIN_SMU) {
      if (cf->tabtype == FCFC_LOOKUP_TYPE_INT) {
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_smu_itab_s8_p8_smin0);
          cfg->ns = cf->ns;
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_itab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_itab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_smu_itab_s8_p8);
          cfg->ns = cf->ns;
          cfg->smin = cf->sbin[0];
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
      else {                    /* cf->tabtype == FCFC_LOOKUP_TYPE_HYBRID */
        if (smin0) {
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s8_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s8_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s8_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s8_p8_smin0;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s16_p8_smin0_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s16_p8_smin0_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s16_p8_smin0;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s16_p8_smin0;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_smu_htab_s8_p8_smin0);
          cfg->ns = cf->ns;
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2bin = cf->s2bin;
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
        else {                  /* smin != 0 */
          if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s8_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s8_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s8_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s8_p8;
            }
          }
          else {                /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
            if (withwt) {
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s16_p8_wt;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s16_p8_wt;
            }
            else {              /* no weight */
              if (isauto)       /* auto pair count */
                paircnt_func = paircnt_balltree_auto_smu_htab_s16_p8;
              else              /* cross pair count */
                paircnt_func = paircnt_balltree_cross_smu_htab_s16_p8;
            }
          }

          FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_smu_htab_s8_p8);
          cfg->ns = cf->ns;
          cfg->smin = cf->sbin[0];
          cfg->smax = cf->sbin[cf->ns];
          cfg->s2bin = cf->s2bin;
          cfg->s2min = cf->s2bin[0];
          cfg->s2max = cf->s2bin[cf->ns];
          cfg->stab = cf->stab;
          cfg->nmu2 = cf->nmu * cf->nmu;
          cfg->mutab = cf->mutab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
          cfg->ntot = cf->ntot;
#endif
          conf = cfg;
        }
      }
    }
    else {                      /* cf->bintype == FCFC_BIN_SPI */
      bool pmin0 = (cf->p2bin[0] == 0);
      if (cf->tabtype == FCFC_LOOKUP_TYPE_INT) {
        if (smin0) {
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_smin0_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_smin0_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_balltree_auto_spi_itab_s8_p8_smin0_pmin0);
            cfg->ns = cf->ns;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_smin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_smin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_spi_itab_s8_p8_smin0);
            cfg->ns = cf->ns;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->pmin = cf->pbin[0];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
        else {                  /* smin != 0 */
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_balltree_auto_spi_itab_s8_p8_pmin0);
            cfg->ns = cf->ns;
            cfg->smin = cf->sbin[0];
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s8_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s8_p16;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_itab_s16_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_itab_s16_p16;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_spi_itab_s8_p8);
            cfg->ns = cf->ns;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmin = REAL_SQRT(cfg->s2min + cfg->p2min) - REAL_EPS;
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
      }
      else {                    /* cf->tabtype == FCFC_LOOKUP_TYPE_HYBRID */
        if (smin0) {
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_smin0_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_smin0_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_smin0_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_smin0_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_smin0_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_smin0_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_balltree_auto_spi_htab_s8_p8_smin0_pmin0);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_smin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_smin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_smin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_smin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_smin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_smin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_spi_htab_s8_p8_smin0);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->pmin = cf->pbin[0];
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
        else {                  /* smin != 0 */
          if (pmin0) {
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_pmin0;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_pmin0;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_pmin0_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_pmin0_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_pmin0;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_pmin0;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg,
                conf_balltree_auto_spi_htab_s8_p8_pmin0);
            cfg->ns = cf->ns;
            cfg->smin = cf->sbin[0];
            cfg->s2bin = cf->s2bin;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
          else {                /* pmin != 0 */
            if (cf->swidth == FCFC_LOOKUP_TABLE_W8) {
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s8_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s8_p16;
                }
              }
            }
            else {              /* cf->swidth == FCFC_LOOKUP_TABLE_W16 */
              if (cf->pwidth == FCFC_LOOKUP_TABLE_W8) {
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p8;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p8;
                }
              }
              else {            /* cf->pwidth == FCFC_LOOKUP_TABLE_W16 */
                if (withwt) {
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16_wt;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16_wt;
                }
                else {          /* no weight */
                  if (isauto)   /* auto pair count */
                    paircnt_func =
                        paircnt_balltree_auto_spi_htab_s16_p16;
                  else          /* cross pair count */
                    paircnt_func =
                        paircnt_balltree_cross_spi_htab_s16_p16;
                }
              }
            }

            FCFC_ALLOC_CNT_CONF(cfg, conf_balltree_auto_spi_htab_s8_p8);
            cfg->ns = cf->ns;
            cfg->s2bin = cf->s2bin;
            cfg->s2min = cf->s2bin[0];
            cfg->s2max = cf->s2bin[cf->ns];
            cfg->stab = cf->stab;
            cfg->np = cf->np;
            cfg->p2bin = cf->p2bin;
            cfg->p2min = cf->p2bin[0];
            cfg->p2max = cf->p2bin[cf->np];
            cfg->msmax = FCFC_MIN(cf->sbin[cf->ns], cf->pbin[cf->np]);
            cfg->tsmin = REAL_SQRT(cfg->s2min + cfg->p2min) - REAL_EPS;
            cfg->tsmax = REAL_SQRT(cfg->s2max + cfg->p2max) + REAL_EPS;
            cfg->ptab = cf->ptab;
#if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
            cfg->ntot = cf->ntot;
#endif
            conf = cfg;
          }
        }
      }
    }
  }

  /* Initialise the stack for dual nodes. */
  STACK_DUAL_NODE stack;
  stack.size = stack.capacity = 0;

  /* Setup the initial stack for dual nodes. */
#ifdef MPI
  /* Prepare for dynamic MPI scheduler. */
  uint32_t *addr = NULL;
  MPI_Aint win_size = 0;
  if (para->rank == para->root) {
    win_size = sizeof(uint32_t);
    if (MPI_Alloc_mem(win_size, MPI_INFO_NULL, &addr)) {
      P_ERR("failed to allocate memory for MPI window\n");
      free(conf); FCFC_QUIT(FCFC_ERR_MPI);
    }
    *addr = 0;
  }
  MPI_Win win;
  if (MPI_Win_create(addr, win_size, sizeof(uint32_t), MPI_INFO_NULL,
      para->comm, &win)) {
    P_ERR("failed to create the MPI window\n");
    free(conf); FCFC_QUIT(FCFC_ERR_MPI);
  }
  uint32_t start, step = cf->nthread * FCFC_STACK_SIZE_PER_THREAD;

  /* Compute the tree level for dual-node distribution. */
  size_t nnodes = para->ntask * FCFC_STACK_SIZE_PER_TASK;
  #ifdef OMP    /* both MPI and OMP */
  nnodes *= cf->nthread * FCFC_STACK_SIZE_PER_THREAD;
  #endif
  if (nnodes > (UINT32_MAX >> 1)) nnodes = UINT32_MAX >> 1;

  if (cf->treetype == FCFC_STRUCT_KDTREE) {     /* k-d tree */
    if (isauto) {       /* auto pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      KDT **nodes = kdtree_get_nodes(tree1, nnodes, &nnodes);
      if (!nodes) {
        free(conf); FCFC_QUIT(FCFC_ERR_MEMORY);
      }

      /* Distribute dual nodes to tasks dynamically. */
      do {
        if (MPI_Win_lock(MPI_LOCK_EXCLUSIVE, para->root, 0, win) ||
            MPI_Fetch_and_op(&step, &start, MPI_UINT32_T, para->root, 0,
            MPI_SUM, win) || MPI_Win_unlock(para->root, win)) {
          P_ERR("failed to fetch dual nodes with MPI tasks\n");
          free(conf); free(nodes);
          FCFC_QUIT(FCFC_ERR_MPI);
        }

        uint32_t end = start + step;
        if (end > (uint32_t) nnodes) end = nnodes;
        for (uint32_t i = start; i < end; i++) {
          for (size_t j = 0; j < nnodes; j++) {
            if (nodes[i]->x[0] <= nodes[j]->x[0])
              stack_push(&stack, nodes[i], nodes[j]);
          }
        }

        /* Perform pair counts. */
  #ifdef OMP
    #pragma omp parallel num_threads(cf->nthread)
        {
          STACK_DUAL_NODE ps;           /* thread-private stack */
          ps.size = ps.capacity = 0;
          const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
          for (size_t i = 0; i < stack.size; i++) {
            stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
            while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
              paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
              paircnt_func(&ps, conf,
                  ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
          }
          stack_destroy(&ps);
        }
        stack_clear(&stack);
  #else
        while (stack.size)
    #if FCFC_SIMD >= FCFC_SIMD_AVX512
          paircnt_func(&stack, conf, cf->pcnt);
    #else
          paircnt_func(&stack, conf, cnt);
    #endif
  #endif
      } while (start + step < (uint32_t) nnodes);

      free(nodes);
    }
    else {              /* cross pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      size_t nnodes2 = nnodes;
      KDT **nodes1, **nodes2;
      if (!(nodes1 = kdtree_get_nodes(tree1, nnodes, &nnodes)) ||
          !(nodes2 = kdtree_get_nodes(tree2, nnodes2, &nnodes2))) {
        free(conf);
        if (nodes1) free(nodes1);
        FCFC_QUIT(FCFC_ERR_MEMORY);
      }

      /* Distribute dual nodes to tasks dynamically. */
      do {
        if (MPI_Win_lock(MPI_LOCK_EXCLUSIVE, para->root, 0, win) ||
            MPI_Fetch_and_op(&step, &start, MPI_UINT32_T, para->root, 0,
            MPI_SUM, win) || MPI_Win_unlock(para->root, win)) {
          P_ERR("failed to fetch dual nodes with MPI tasks\n");
          free(conf); free(nodes1); free(nodes2);
          FCFC_QUIT(FCFC_ERR_MPI);
        }

        uint32_t end = start + step;
        if (end > (uint32_t) nnodes) end = nnodes;
        for (uint32_t i = start; i < end; i++) {
          for (size_t j = 0; j < nnodes2; j++)
            stack_push(&stack, nodes1[i], nodes2[j]);
        }

        /* Perform pair counts. */
  #ifdef OMP
    #pragma omp parallel num_threads(cf->nthread)
        {
          STACK_DUAL_NODE ps;           /* thread-private stack */
          ps.size = ps.capacity = 0;
          const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
          for (size_t i = 0; i < stack.size; i++) {
            stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
            while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
              paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
              paircnt_func(&ps, conf,
                  ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
          }
          stack_destroy(&ps);
        }
        stack_clear(&stack);
  #else
        while (stack.size)
    #if FCFC_SIMD >= FCFC_SIMD_AVX512
          paircnt_func(&stack, conf, cf->pcnt);
    #else
          paircnt_func(&stack, conf, cnt);
    #endif
  #endif
      } while (start + step < (uint32_t) nnodes);

      free(nodes1);
      free(nodes2);
    }
  }
  else {                                        /* ball tree */
    if (isauto) {       /* auto pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      BLT **nodes = balltree_get_nodes(tree1, nnodes, &nnodes);
      if (!nodes) {
        free(conf); FCFC_QUIT(FCFC_ERR_MEMORY);
      }

      /* Distribute dual nodes to tasks dynamically. */
      do {
        if (MPI_Win_lock(MPI_LOCK_EXCLUSIVE, para->root, 0, win) ||
            MPI_Fetch_and_op(&step, &start, MPI_UINT32_T, para->root, 0,
            MPI_SUM, win) || MPI_Win_unlock(para->root, win)) {
          P_ERR("failed to fetch dual nodes with MPI tasks\n");
          free(conf); free(nodes);
          FCFC_QUIT(FCFC_ERR_MPI);
        }

        uint32_t end = start + step;
        if (end > (uint32_t) nnodes) end = nnodes;
        for (uint32_t i = start; i < end; i++) {
          for (size_t j = 0; j < nnodes; j++) {
            if (nodes[i]->x[0] <= nodes[j]->x[0])
              stack_push(&stack, nodes[i], nodes[j]);
          }
        }

        /* Perform pair counts. */
  #ifdef OMP
    #pragma omp parallel num_threads(cf->nthread)
        {
          STACK_DUAL_NODE ps;           /* thread-private stack */
          ps.size = ps.capacity = 0;
          const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
          for (size_t i = 0; i < stack.size; i++) {
            stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
            while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
              paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
              paircnt_func(&ps, conf,
                  ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
          }
          stack_destroy(&ps);
        }
        stack_clear(&stack);
  #else
        while (stack.size)
    #if FCFC_SIMD  >=  FCFC_SIMD_AVX512
          paircnt_func(&stack, conf, cf->pcnt);
    #else
          paircnt_func(&stack, conf, cnt);
    #endif
  #endif
      } while (start + step < (uint32_t) nnodes);

      free(nodes);
    }
    else {              /* cross pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      size_t nnodes2 = nnodes;
      BLT **nodes1, **nodes2;
      if (!(nodes1 = balltree_get_nodes(tree1, nnodes, &nnodes)) ||
          !(nodes2 = balltree_get_nodes(tree2, nnodes2, &nnodes2))) {
        free(conf);
        if (nodes1) free(nodes1);
        FCFC_QUIT(FCFC_ERR_MEMORY);
      }

      /* Distribute dual nodes to tasks dynamically. */
      do {
        if (MPI_Win_lock(MPI_LOCK_EXCLUSIVE, para->root, 0, win) ||
            MPI_Fetch_and_op(&step, &start, MPI_UINT32_T, para->root, 0,
            MPI_SUM, win) || MPI_Win_unlock(para->root, win)) {
          P_ERR("failed to fetch dual nodes with MPI tasks\n");
          free(conf); free(nodes1); free(nodes2);
          FCFC_QUIT(FCFC_ERR_MPI);
        }

        uint32_t end = start + step;
        if (end > (uint32_t) nnodes) end = nnodes;
        for (uint32_t i = start; i < end; i++) {
          for (size_t j = 0; j < nnodes2; j++)
            stack_push(&stack, nodes1[i], nodes2[j]);
        }

        /* Perform pair counts. */
  #ifdef OMP
    #pragma omp parallel num_threads(cf->nthread)
        {
          STACK_DUAL_NODE ps;           /* thread-private stack */
          ps.size = ps.capacity = 0;
          const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
          for (size_t i = 0; i < stack.size; i++) {
            stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
            while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
              paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
              paircnt_func(&ps, conf,
                  ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
          }
          stack_destroy(&ps);
        }
        stack_clear(&stack);
  #else
        while (stack.size)
    #if FCFC_SIMD  >=  FCFC_SIMD_AVX512
          paircnt_func(&stack, conf, cf->pcnt);
    #else
          paircnt_func(&stack, conf, cnt);
    #endif
  #endif
      } while (start + step < (uint32_t) nnodes);

      free(nodes1);
      free(nodes2);
    }
  }

  /* Cleanup MPI scheduler. */
  if (MPI_Win_free(&win) ||
      (para->rank == para->root && MPI_Free_mem(addr))) {
    P_ERR("failed to release memory for MPI scheduler\n");
    free(conf); stack_destroy(&stack);
    FCFC_QUIT(FCFC_ERR_MPI);
  }
#else
  #ifdef OMP    /* OMP only */
  size_t nnodes = cf->nthread * FCFC_STACK_SIZE_PER_THREAD;
  if (nnodes > (UINT32_MAX >> 1)) nnodes = UINT32_MAX >> 1;

  /* Expand the tree and push the node pairs at the same level. */
  if (cf->treetype == FCFC_STRUCT_KDTREE) {     /* k-d tree */
    if (isauto) {       /* auto pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      KDT **nodes = kdtree_get_nodes(tree1, nnodes, &nnodes);
      if (!nodes) {
        free(conf);
        return FCFC_ERR_MEMORY;
      }

      /* Distribute dual nodes to threads. */
      for (size_t i = 0; i < nnodes; i++) {
        for (size_t j = 0; j < nnodes; j++) {
          if (nodes[i]->x[0] <= nodes[j]->x[0])
            stack_push(&stack, nodes[i], nodes[j]);
        }
      }

      /* Perform pair counts. */
    #pragma omp parallel num_threads(cf->nthread)
      {
        STACK_DUAL_NODE ps;             /* thread-private stack */
        ps.size = ps.capacity = 0;
        const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < stack.size; i++) {
          stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
          while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
            paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
            paircnt_func(&ps, conf,
                ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
        }
        stack_destroy(&ps);
      }
      free(nodes);
    }
    else {              /* cross pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      size_t nnodes2 = nnodes;
      KDT **nodes1, **nodes2;
      if (!(nodes1 = kdtree_get_nodes(tree1, nnodes, &nnodes)) ||
          !(nodes2 = kdtree_get_nodes(tree2, nnodes2, &nnodes2))) {
        free(conf);
        if (nodes1) free(nodes1);
        return FCFC_ERR_MEMORY;
      }

      /* Distribute dual nodes to threads. */
      for (size_t i = 0; i < nnodes; i++) {
        for (size_t j = 0; j < nnodes2; j++)
          stack_push(&stack, nodes1[i], nodes2[j]);
      }

      /* Perform pair counts. */
    #pragma omp parallel num_threads(cf->nthread)
      {
        STACK_DUAL_NODE ps;             /* thread-private stack */
        ps.size = ps.capacity = 0;
        const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < stack.size; i++) {
          stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
          while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
            paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
            paircnt_func(&ps, conf,
                ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
        }
        stack_destroy(&ps);
      }
      free(nodes1);
      free(nodes2);
    }
  }
  else {                                        /* ball tree */
    if (isauto) {       /* auto pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      BLT **nodes = balltree_get_nodes(tree1, nnodes, &nnodes);
      if (!nodes) {
        free(conf);
        return FCFC_ERR_MEMORY;
      }

      /* Distribute dual nodes to threads. */
      for (size_t i = 0; i < nnodes; i++) {
        for (size_t j = 0; j < nnodes; j++) {
          if (nodes[i]->x[0] <= nodes[j]->x[0])
            stack_push(&stack, nodes[i], nodes[j]);
        }
      }

      /* Perform pair counts. */
    #pragma omp parallel num_threads(cf->nthread)
      {
        STACK_DUAL_NODE ps;             /* thread-private stack */
        ps.size = ps.capacity = 0;
        const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < stack.size; i++) {
          stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
          while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
            paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
            paircnt_func(&ps, conf,
                ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
        }
        stack_destroy(&ps);
      }
      free(nodes);
    }
    else {              /* cross pair counts */
      /* Expand the tree and push the node pairs at the same level. */
      size_t nnodes2 = nnodes;
      BLT **nodes1, **nodes2;
      if (!(nodes1 = balltree_get_nodes(tree1, nnodes, &nnodes)) ||
          !(nodes2 = balltree_get_nodes(tree2, nnodes2, &nnodes2))) {
        free(conf);
        if (nodes1) free(nodes1);
        return FCFC_ERR_MEMORY;
      }

      /* Distribute dual nodes to threads. */
      for (size_t i = 0; i < nnodes; i++) {
        for (size_t j = 0; j < nnodes2; j++)
          stack_push(&stack, nodes1[i], nodes2[j]);
      }

      /* Perform pair counts. */
    #pragma omp parallel num_threads(cf->nthread)
      {
        STACK_DUAL_NODE ps;             /* thread-private stack */
        ps.size = ps.capacity = 0;
        const int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < stack.size; i++) {
          stack_push(&ps, stack.nodes[i][0], stack.nodes[i][1]);
          while (ps.size)
    #if FCFC_SIMD  <  FCFC_SIMD_AVX512
            paircnt_func(&ps, conf, cf->pcnt + tid * cf->ntot);
    #else
            paircnt_func(&ps, conf,
                ((char *) cf->pcnt) + tid * cf->ntot * FCFC_SIMD_BYTES);
    #endif
        }
        stack_destroy(&ps);
      }
      free(nodes1);
      free(nodes2);
    }
  }
  #else         /* serial code */
  stack_push(&stack, tree1, tree2);
  while (stack.size)
    #if FCFC_SIMD  >=  FCFC_SIMD_AVX512
    paircnt_func(&stack, conf, cf->pcnt);
    #else
    paircnt_func(&stack, conf, cnt);
    #endif
  #endif
#endif
  free(conf);
  stack_destroy(&stack);

#ifdef OMP
  /* Gather pair counts from threads. */
  if (withwt) {
  #pragma omp parallel for num_threads(cf->nthread) default(none) shared(cf,cnt)
    for (size_t i = 0; i < cf->ntot; i++) {
      for (int j = 0; j < cf->nthread; j++) {
#if             FCFC_SIMD  <  FCFC_SIMD_AVX512
        cnt[i].d += cf->pcnt[i + j * cf->ntot].d;
#else        /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
        for (int k = 0; k < FCFC_SIMD_BYTES / (int) sizeof(double); k++) {
          cnt[i].d += ((double *) cf->pcnt)[i +
              (j * FCFC_SIMD_BYTES / (int) sizeof(double) + k) * cf->ntot];
        }
#endif
      }
    }
  }
  else {
  #pragma omp parallel for num_threads(cf->nthread) default(none) shared(cf,cnt)
    for (size_t i = 0; i < cf->ntot; i++) {
      for (int j = 0; j < cf->nthread; j++) {
#if             FCFC_SIMD  <  FCFC_SIMD_AVX512
        cnt[i].i += cf->pcnt[i + j * cf->ntot].i;
#else        /* FCFC_SIMD  >=  FCFC_SIMD_AVX512 */
        for (int k = 0; k < FCFC_SIMD_BYTES / (int) sizeof(int64_t); k++) {
          cnt[i].i += ((int64_t *) cf->pcnt)[i +
              (j * FCFC_SIMD_BYTES / (int) sizeof(int64_t) + k) * cf->ntot];
        }
#endif
      }
    }
  }
#else
  #if   FCFC_SIMD  >=  FCFC_SIMD_AVX512
  /* Gather pair counts */
  if (withwt) {
    for (size_t i = 0; i < cf->ntot; i++) {
      for (int k = 0; k < FCFC_SIMD_BYTES / (int) sizeof(double); k++)
        cnt[i].d += ((double *) cf->pcnt)[i + k * cf->ntot];
    }
  }
  else {
    for (size_t i = 0; i < cf->ntot; i++) {
      for (int k = 0; k < FCFC_SIMD_BYTES / (int) sizeof(double); k++)
        cnt[i].i += ((int64_t *) cf->pcnt)[i + k * cf->ntot];
    }
  }
  #endif
#endif

#ifdef MPI
  /* Gather pair counts from different MPI processes. */
  MPI_Request req;
  if (withwt) {
    double *pcnt;
    if (sizeof(double) == sizeof(COUNT)) pcnt = (double *) cnt;
    else {
      if (!(pcnt = malloc(cf->ntot * sizeof(double)))) {
        P_ERR("failed to allocate memory for gathering pair counts\n");
        FCFC_QUIT(FCFC_ERR_MEMORY);
      }
      for (size_t i = 0; i < cf->ntot; i++) pcnt[i] = cnt[i].d;
    }

    if (para->rank == para->root) {
      if (MPI_Ireduce(MPI_IN_PLACE, pcnt, cf->ntot, MPI_DOUBLE, MPI_SUM,
          para->root, para->comm, &req)) {
        P_ERR("failed to gather pair counts from MPI tasks\n");
        FCFC_QUIT(FCFC_ERR_MPI);
      }
    }
    else {
      if (MPI_Ireduce(pcnt, NULL, cf->ntot, MPI_DOUBLE, MPI_SUM,
                      para->root, para->comm, &req)) {
        P_ERR("failed to gather pair counts from MPI tasks\n");
        FCFC_QUIT(FCFC_ERR_MPI);
      }
    }

    if (MPI_Wait(&req, MPI_STATUS_IGNORE)) {
      P_ERR("failed to gather pair counts from MPI tasks\n");
      FCFC_QUIT(FCFC_ERR_MPI);
    }
    if (sizeof(double) != sizeof(COUNT)) free(pcnt);
  }
  else {        /* no weight */
    int64_t *pcnt;
    if (sizeof(int64_t) == sizeof(COUNT)) pcnt = (int64_t *) cnt;
    else {
      if (!(pcnt = malloc(cf->ntot * sizeof(int64_t)))) {
        P_ERR("failed to allocate memory for gathering pair counts\n");
        FCFC_QUIT(FCFC_ERR_MEMORY);
      }
      for (size_t i = 0; i < cf->ntot; i++) pcnt[i] = cnt[i].i;
    }

    if (para->rank == para->root) {
      if (MPI_Ireduce(MPI_IN_PLACE, pcnt, cf->ntot, MPI_INT64_T, MPI_SUM,
          para->root, para->comm, &req)) {
        P_ERR("failed to gather pair counts from MPI tasks\n");
        FCFC_QUIT(FCFC_ERR_MPI);
      }
    }
    else {
      if (MPI_Ireduce(pcnt, NULL, cf->ntot, MPI_INT64_T, MPI_SUM,
          para->root, para->comm, &req)) {
        P_ERR("failed to gather pair counts from MPI tasks\n");
        FCFC_QUIT(FCFC_ERR_MPI);
      }
    }

    if (MPI_Wait(&req, MPI_STATUS_IGNORE)) {
      P_ERR("failed to gather pair counts from MPI tasks\n");
      FCFC_QUIT(FCFC_ERR_MPI);
    }
    if (sizeof(int64_t) != sizeof(COUNT)) free(pcnt);
  }
#endif

  return 0;
}

