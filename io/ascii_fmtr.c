/*******************************************************************************
* ascii_fmtr.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define_comm.h"
#include "ascii_fmtr.h"
#include <limits.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>

/******************************************************************************
Function `print_fmtr_error`:
  Print the error message for the formatter.
Arguments:
  * `fmtr`:     the formatter string;
  * `ptr`:      pointer to the position on error.
******************************************************************************/
static void print_fmtr_error(const char *fmtr, const char *ptr) {
  fprintf(stderr, "%s\n", fmtr);
  if (ptr > fmtr) {
    for (ptrdiff_t i = 0; i < ptr - fmtr; i++) fprintf(stderr, " ");
    fprintf(stderr, "^\n");
  }
}

/******************************************************************************
Function `arg_insert`:
  Insert one parsed argument to the argument array.
Arguments:
  * `arg`:      address of the argument array;
  * `num`:      current number of arguments;
  * `dtype`:    data type of the argument;
  * `fmtr`:     the string for the argument specifier;
  * `len`:      length of the specifier string.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int arg_insert(asc_arg_t **arg, int *num, const asc_dtype_t dtype,
    const char *fmtr, const int len) {
  if (*num == INT_MAX) {        /* there is no more space for the insertion */
    P_ERR("the formatter string is too long\n");
    return FCFC_ERR_ASCII;
  }

  asc_arg_t *tmp;
  /* Check if the allocated space is enough. */
  if ((*num & (*num - 1)) == 0) {       /* num is 0 or power of 2 */
    int size = 1;
    if (INT_MAX / 2 < *num) *num = INT_MAX;
    else if (*num) size = *num << 1;    /* double the size */

    tmp = realloc(*arg, size * sizeof(asc_arg_t));
    if (!tmp) {
      P_ERR("failed to allocate memory for parsing the formatter\n");
      return FCFC_ERR_MEMORY;
    }
    *arg = tmp;
  }

  /* Append the argument to the end of the array. */
  tmp = *arg + *num;
  tmp->dtype = dtype;

  /* Copy the specifier, add "%n", and null terminate it. */
  tmp->fmtr = calloc(len + 3, sizeof(char));
  if (!tmp->fmtr) {
    P_ERR("failed to allocate memory for parsing the formatter\n");
    return FCFC_ERR_MEMORY;
  }
  memcpy(tmp->fmtr, fmtr, len);
  tmp->fmtr[len] = '%';
  tmp->fmtr[len + 1] = 'n';
  *num += 1;
  return 0;
}

/******************************************************************************
Function `ascii_arg_destroy`:
  Deconstruct the structure array for parsed arguments.
Arguments:
  * `arg`:      parsed arguments from a formatter string;
  * `num`:      number of parsed elements in the array.
******************************************************************************/
void ascii_arg_destroy(asc_arg_t *arg, const int num) {
  if (!arg) return;
  for (int i = 0; i < num; i++) if (arg[i].fmtr) free(arg[i].fmtr);
  free(arg);
}

/******************************************************************************
Function `parse_ascii_fmtr`:
  Parse the formatter string for reading ASCII files.
  The formatter is compliant with fscanf, see C99 standard: section 7.19.6.2.
Arguments:
  * `fmtr`:     the formatter string;
  * `num`:      number of parsed arguments;
  * `rnum`:     number of arguments that are not suppressed.
Return:
  Arguments parsed from the formatter on success; NULL on error.
******************************************************************************/
asc_arg_t *parse_ascii_fmtr(const char *fmtr, int *num, int *rnum) {
  if (!fmtr || !num || !rnum) return NULL;
  if (!(*fmtr)) {
    P_ERR("empty formatter string\n");
    return NULL;
  }
  int n, rn;
  const char *start = fmtr;
  const char *fmt = fmtr;
  asc_arg_t *args = NULL;
  n = rn = 0;

  while (*fmt) {
    char c = *fmt++;
    if (c != '%') continue;

    /* Omit the literal '%'. */
    c = *fmt++;
    if (c == '%') continue;

    /* Check the assignment-suppressing charater '*'. */
    bool skip = false;
    if (c == '*') {
      skip = true;
      c = *fmt++;
    }

    const char *ptr = fmt - 1;
    /* Check the maximum field width (a decimal integer greater than 0). */
    if (isdigit(c)) {
      int width = 0;
      /* Read the width. */
      while (isdigit(c)) {
        int digit = c - '0';
        if (INT_MAX / 10 < width) {
          P_ERR("field width overflow:\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        width *= 10;
        if (INT_MAX - digit < width) {
          P_ERR("field width overflow:\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        width += digit;
        c = *fmt++;
      }
      if (!width) {
        P_ERR("zero field width not allowed:\n");
        print_fmtr_error(fmtr, ptr);
        ascii_arg_destroy(args, n);
        return NULL;
      }
    }

    /* Check the length modifier. */
    char lmod = '\0';
    ptr = fmt - 1;
    switch (c) {
      case 'h':
        if (*fmt == 'h') {
          fmt++;
          lmod = 'H';
        }
        else lmod = 'h';
        break;
      case 'l':
        if (*fmt == 'l') {
          fmt++;
          lmod = 'L';
        }
        else lmod = 'l';
        break;
      case 'j':
      case 'z':
      case 't': lmod = c; break;
      case 'L': lmod = 'D'; break;
      default:          /* Not a C99 modifier. */
        fmt--;
        break;
    }

    /* Check the conversion specifier. */
    asc_dtype_t dtype = ASCII_DTYPE_SKIP;
    c = *fmt++;
    switch (c) {
      case 'n':
      case 'd':
      case 'i':
        switch (lmod) {
          case '\0': dtype = ASCII_DTYPE_INT; break;
          case 'H': dtype = ASCII_DTYPE_SCHAR; break;
          case 'h': dtype = ASCII_DTYPE_SHRT; break;
          case 'L': dtype = ASCII_DTYPE_LLNG; break;
          case 'l': dtype = ASCII_DTYPE_LONG; break;
          case 'j': dtype = ASCII_DTYPE_INTMX; break;
          case 'z': dtype = ASCII_DTYPE_SIZE; break;
          case 't': dtype = ASCII_DTYPE_PTR; break;
          default:
            P_ERR("invalid length modifier for specifier '%c':\n", c);
            print_fmtr_error(fmtr, ptr);
            ascii_arg_destroy(args, n);
            return NULL;
        }
        break;
      case 'o':
      case 'u':
      case 'x':
      case 'X':
        switch (lmod) {
          case '\0': dtype = ASCII_DTYPE_UINT; break;
          case 'H': dtype = ASCII_DTYPE_UCHAR; break;
          case 'h': dtype = ASCII_DTYPE_USHRT; break;
          case 'L': dtype = ASCII_DTYPE_ULLNG; break;
          case 'l': dtype = ASCII_DTYPE_ULONG; break;
          case 'j': dtype = ASCII_DTYPE_UINTMX; break;
          case 'z': dtype = ASCII_DTYPE_SIZE; break;
          case 't': dtype = ASCII_DTYPE_PTR; break;
          default:
            P_ERR("invalid length modifier for specifier '%c':\n", c);
            print_fmtr_error(fmtr, ptr);
            ascii_arg_destroy(args, n);
            return NULL;
        }
        break;
      case 'a':
      case 'A':
      case 'e':
      case 'E':
      case 'f':
      case 'F':
      case 'g':
      case 'G':
        switch (lmod) {
          case '\0': dtype = ASCII_DTYPE_FLT; break;
          case 'l': dtype = ASCII_DTYPE_DBL; break;
          case 'D': dtype = ASCII_DTYPE_LDBL; break;
          default:
            P_ERR("invalid length modifier for specifier '%c':\n", c);
            print_fmtr_error(fmtr, ptr);
            ascii_arg_destroy(args, n);
            return NULL;
        }
        break;
      case 'c':
        if (lmod == '\0') dtype = ASCII_DTYPE_CHAR;
        else if (lmod == 'l') dtype = ASCII_DTYPE_WCHAR;
        else {
          P_ERR("invalid length modifier for specifier 'c':\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        break;
      case 's':
        if (lmod == '\0') dtype = ASCII_DTYPE_STR;
        else if (lmod == 'l') dtype = ASCII_DTYPE_WSTR;
        else {
          P_ERR("invalid length modifier for specifier 's':\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        break;
      case '[':
        if (lmod == '\0') dtype = ASCII_DTYPE_STR;
        else if (lmod == 'l') dtype = ASCII_DTYPE_WSTR;
        else {
          P_ERR("invalid length modifier for specifier '[]':\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        ptr = fmt - 1;
        if (*fmt == '^') ++fmt;
        if (*fmt == ']') ++fmt;
        while ((c = *fmt++) != '\0' && c != ']');
        if (c == '\0') {
          P_ERR("unclosed conversion specifier:\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        break;
      case 'p':
        if (lmod == '\0') dtype = ASCII_DTYPE_VOID;
        else {
          P_ERR("unexpected length modifier for specifier 'p':\n");
          print_fmtr_error(fmtr, ptr);
          ascii_arg_destroy(args, n);
          return NULL;
        }
        break;
      case '\0':
      case ' ':
        P_ERR("conversion specifier not found:\n");
        print_fmtr_error(fmtr, fmt - 1);
        ascii_arg_destroy(args, n);
        return NULL;
      default:
        P_ERR("unexpected conversion specifier:\n");
        print_fmtr_error(fmtr, fmt - 1);
        ascii_arg_destroy(args, n);
        return NULL;
    }

    if (skip) {
      if (arg_insert(&args, &n, ASCII_DTYPE_SKIP, start, fmt - start)) {
        ascii_arg_destroy(args, n);
        return NULL;
      }
    }
    else {
      if (arg_insert(&args, &n, dtype, start, fmt - start)) {
        ascii_arg_destroy(args, n);
        return NULL;
      }
      rn++;
    }
    start = fmt;
  }

  /* Check uninterpreted part of the formatter string. */
  if (start != fmt) {
    const char *tmp = start;
    while (isspace(*tmp)) tmp++;
    if (tmp != fmt) {
      P_WRN("uninterpreted part of the formatter:\n");
      print_fmtr_error(fmtr, start);
    }
  }

  /* Release memory that are not necessary any more. */
  asc_arg_t *tmp = realloc(args, sizeof(asc_arg_t) * n);
  if (!tmp) P_WRN("failed to reduce the allocated memory\n");
  else (args = tmp);

  *num = n;
  *rnum = rn;
  return args;
}

