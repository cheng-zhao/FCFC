/*******************************************************************************
* benchmark/struct/data_struct.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "data_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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
#ifdef BENCHMARK_TREE_AS_ARRAY
void stack_push(STACK_DUAL_NODE *s, const size_t a, const size_t b)
#else
void stack_push(STACK_DUAL_NODE *s, const void *a, const void *b)
#endif
{
  if (s->size >= s->capacity) {
    /* Enlarge the memory allocated for the stack. */
    if (s->capacity) {
      if ((BENCHMARK_STACK_MAX_SIZE >> 1) <= s->capacity) {
        P_ERR("too many elements to be pushed to the stack of dual nodes\n");
        exit(EXIT_FAILURE);
      }
      s->capacity <<= 1;
    }
    else {      /* initialise the stack */
      s->nodes = NULL;
      s->capacity = BENCHMARK_STACK_INIT_SIZE;
    }

    if (s->capacity <= s->size) {
      P_ERR("unable to expand the size of the stack of dual nodes\n");
      exit(EXIT_FAILURE);
    }

    DUAL_NODE *tmp = realloc(s->nodes, s->capacity * sizeof *tmp);
    if (!tmp) {
      P_ERR("failed to allocate memory for the stack of dual nodes\n");
      exit(EXIT_FAILURE);
    }
    s->nodes = tmp;
  }

  s->nodes[s->size].a = a;
  s->nodes[s->size++].b = b;
}

/******************************************************************************
Function `stack_destroy`:
  Deconstruct the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack.
******************************************************************************/
void stack_destroy(STACK_DUAL_NODE *s) {
  if (!s->capacity) return;
  s->size = s->capacity = 0;
  if (s->nodes) free(s->nodes);
}
