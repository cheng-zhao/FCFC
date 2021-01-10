/*******************************************************************************
* quick_select.c: this file is part of the FCFC program.

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

#if defined(QSELECT_COMPARE) && defined(QSELECT_DTYPE)

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

/*******************************************************************************
  Implementation of the adaptive quick selection algorithm, based on
  the algorithm median-of-ninthers
  ref: http://dx.doi.org/10.4230/LIPIcs.SEA.2017.24

  The original implementation is publicly available:
      https://github.com/andralex/MedianOfNinthers
  with the following license statements:

  Copyright Andrei Alexandrescu, 2016-.
  Distributed under the Boost Software License, Version 1.0.
  (See copy at https://www.boost.org/LICENSE_1_0.txt)

  This variant has a similar interface as `qsort`.
  But the comparison function is passed as a macro, for efficiency.
*******************************************************************************/

#define QSELECT_SWAP(a,b) {                                \
  memcpy(buf, (a), sizeof(QSELECT_DTYPE));                 \
  memcpy((a), (b), sizeof(QSELECT_DTYPE));                 \
  memcpy((b), buf, sizeof(QSELECT_DTYPE));                 \
}

/* Declaration of the quick selection function. */
static void qselect(QSELECT_DTYPE *, size_t, size_t, void *, void *);

/******************************************************************************
Functions `expand_partition_l`, `expand_partition_r`, `expand_partition`:
  Expand the partition in the range a <= p <= b, to all elements of the array.
Arguments:
  * `A`:        pointer to the array;
  * `len`:      length of the array `A`;
  * `a`:        lower limit of the current partition;
  * `p`:        the pivot;
  * `b`:        upper limit of the current partition;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
Return:
  Index of the partition point (pivot).
******************************************************************************/
static size_t expand_partition_l(QSELECT_DTYPE *A, size_t a, size_t p,
    void *buf, void *arg) {
  size_t l = 0;
  const size_t oldp = p;
  bool done = false;
  for (; a < p; l++) {
    if (l == a) {
      done = true;
      break;
    }
    if (QSELECT_COMPARE(A + oldp, A + l, arg) >= 0) continue;
    p -= 1;
    QSELECT_SWAP(A + l, A + p);
  }
  if (!done) {
    for (;; l++) {
      if (l == p) break;
      if (QSELECT_COMPARE(A + oldp, A + l, arg) >= 0) continue;
      for (;;) {
        if (l == p) {
          done = true;
          break;
        }
        p -= 1;
        if (QSELECT_COMPARE(A + p, A + oldp, arg) < 0) {
          QSELECT_SWAP(A + p, A + l);
          break;
        }
      }
      if (done) break;
    }
  }

  QSELECT_SWAP(A + p, A + oldp);
  return p;
}

static size_t expand_partition_r(QSELECT_DTYPE *A, size_t b, size_t r,
    void *buf, void *arg) {
  size_t p = 0;
  bool done = false;
  for (; p < b; r--) {
    if (r == b) {
      done = true;
      break;
    }
    if (QSELECT_COMPARE(A, A + r, arg) <= 0) continue;
    p += 1;
    QSELECT_SWAP(A + r, A + p);
  }
  if (!done) {
    for (; r > p; r--) {
      if (QSELECT_COMPARE(A, A + r, arg) <= 0) continue;
      while (r > p) {
        p += 1;
        if (QSELECT_COMPARE(A, A + p, arg) < 0) {
          QSELECT_SWAP(A + r, A + p);
          break;
        }
      }
    }
  }

  QSELECT_SWAP(A, A + p);
  return p;
}

static size_t expand_partition(QSELECT_DTYPE *A, size_t len, size_t a,
    size_t p, size_t b, void *buf, void *arg) {
  size_t l = 0;
  len -= 1;
  b -= 1;
  for (;; l++, len--) {
    for (;; l++) {
      if (l == a)
        return p + expand_partition_r(A + p, b - p, len - p, buf, arg);
      if (QSELECT_COMPARE(A + p, A + l, arg) < 0) break;
    }
    for (;; len--) {
      if (len == b)
        return l + expand_partition_l(A + l, a - l, p - l, buf, arg);
      if (QSELECT_COMPARE(A + p, A + len, arg) >= 0) break;
    }
    QSELECT_SWAP(A + l, A + len);
  }
}

/******************************************************************************
Functions `Hoare_partition`, `median_of_minima`, `median_of_maxima`:
  Move the nth smallest elements of an array A, to the left of A[n].
Arguments:
  * `A`:        pointer to the array;
  * `n`:        the number of smallest elements to be moved;
  * `len`:      length of the array `A`;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
Return:
  Index of the partition point (pivot).
******************************************************************************/
static size_t Hoare_partition(QSELECT_DTYPE *A, size_t n, size_t len,
    void *buf, void *arg) {
  QSELECT_SWAP(A, A + n);
  size_t a = 1;
  size_t b = len - 1;
  bool quit_loop = false;
  for (;;) {
    for (;;) {
      if (a > b) {
        quit_loop = true;
        break;
      }
      if (QSELECT_COMPARE(A, A + a, arg) <= 0) break;
      a += 1;
    }
    if (quit_loop) break;
    while (QSELECT_COMPARE(A, A + b, arg) < 0) b -= 1;
    if (a >= b) break;
    QSELECT_SWAP(A + a, A + b);
    a += 1;
    b -= 1;
  }
  a -= 1;
  QSELECT_SWAP(A, A + a);
  return a;
}

static size_t median_of_minima(QSELECT_DTYPE *A, size_t n, size_t len,
    void *buf, void *arg) {
  const size_t twon = n << 1;
  const size_t gamma1 = (len - twon) / twon;    /* gamma - 1 */
  size_t j = twon;
  for (size_t i = 0; i < twon; i++) {
    const size_t limit = j + gamma1;
    size_t imin = j;
    while (++j < limit) {
      if (QSELECT_COMPARE(A + j, A + imin, arg) < 0) imin = j;
    }
    if (QSELECT_COMPARE(A + imin, A + i, arg) < 0)
      QSELECT_SWAP(A + i, A + imin);
  }
  qselect(A, n, twon, buf, arg);
  return expand_partition(A, len, 0, n, twon, buf, arg);
}

static size_t median_of_maxima(QSELECT_DTYPE *A, size_t n, size_t len,
    void *buf, void *arg) {
  const size_t twon = (len - n) << 1;
  const size_t istart = len - twon;
  const size_t gamma1 = istart / twon;
  size_t j = istart - twon * gamma1;
  for (size_t i = istart; i < len; i++) {
    const size_t limit = j + gamma1;
    size_t imax = j;
    while (++j < limit) {
      if (QSELECT_COMPARE(A + j, A + imax, arg) > 0) imax = j;
    }
    if (QSELECT_COMPARE(A + imax, A + i, arg) > 0)
      QSELECT_SWAP(A + i, A + imax);
  }
  qselect(A + istart, len - n, twon, buf, arg);
  return expand_partition(A, len, istart, n, len, buf, arg);
}

/******************************************************************************
Function `median3`:
  Compute the median of an array with 3 elements.
Arguments:
  * `A`:        pointer to the array;
  * `i1-i3`:    the three elements of the array;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
******************************************************************************/
static inline size_t median3(QSELECT_DTYPE *A, size_t i1, size_t i2, size_t i3,
    void *buf, void *arg) {
  if (QSELECT_COMPARE(A + i1, A + i3, arg) > 0) {
    size_t tmp = i1;
    i1 = i3;
    i3 = tmp;
  }
  if (QSELECT_COMPARE(A + i2, A + i3, arg) > 0) return i3;
  if (QSELECT_COMPARE(A + i2, A + i1, arg) < 0) return i1;
  return i2;
}

/******************************************************************************
Function `ninther`:
  Compute the median of an array with 9 elements.
Arguments:
  * `A`:        pointer to the array;
  * `i1-i9`:    the nine elements of the array;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
******************************************************************************/
static inline void ninther(QSELECT_DTYPE *A, size_t i1, size_t i2, size_t i3,
    size_t i4, size_t i5, size_t i6, size_t i7, size_t i8, size_t i9,
    void *buf, void *arg) {
  i2 = median3(A, i1, i2, i3, buf, arg);
  i8 = median3(A, i7, i8, i9, buf, arg);
  if (QSELECT_COMPARE(A + i2, A + i8, arg) > 0) {
    size_t tmp = i2;
    i2 = i8;
    i8 = tmp;
  }
  if (QSELECT_COMPARE(A + i4, A + i6, arg) > 0) {
    size_t tmp = i4;
    i4 = i6;
    i6 = tmp;
  }
  if (QSELECT_COMPARE(A + i4, A + i5, arg) <= 0) {
    if (QSELECT_COMPARE(A + i5, A + i6, arg) > 0) i4 = i6;
    else {
      if (QSELECT_COMPARE(A + i2, A + i5, arg) > 0) {
        QSELECT_SWAP(A + i2, A + i5);
        return;
      }
      if (QSELECT_COMPARE(A + i5, A + i8, arg) > 0) {
        QSELECT_SWAP(A + i5, A + i8);
        return;
      }
      return;
    }
  }

  if (QSELECT_COMPARE(A + i2, A + i4, arg) > 0) i4 = i2;
  else if (QSELECT_COMPARE(A + i4, A + i8, arg) > 0) i4 = i8;
  QSELECT_SWAP(A + i4, A + i5);
}

/******************************************************************************
Function `median_of_ninthers`:
  Estimate the median of an array, and move smaller elements to its left.
Arguments:
  * `A`:        pointer to the array;
  * `len`:      length of the array `A`;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
******************************************************************************/
static size_t median_of_ninthers(QSELECT_DTYPE *A, size_t len, void *buf,
    void *arg) {
  size_t phi;
  if (len <= 1024) phi = len / 12;
  else if (len <= 131072) phi = len / 64;
  else phi = len / 1024;
  size_t p = phi >> 1;
  const size_t l = (len >> 1) - p;
  const size_t r = l + phi;
  const size_t gap = (len - 9 * phi) >> 2;
  size_t a = l - (phi << 2) - gap;
  size_t b = r + gap;
  for (size_t i = l; i < r; i++) {
    ninther(A, a, i - phi, b, a + 1, i, b + 1, a + 2, i + phi, b + 2,
        buf, arg);
    a += 3;
    b += 3;
  }
  qselect(A + l, p, phi, buf, arg);
  return expand_partition(A, len, l, l + p, r, buf, arg);
}

/******************************************************************************
Function `qselect`:
  Move the nth smallest elements of an array A, to the left of A[n].
Arguments:
  * `A`:        pointer to the array;
  * `n`:        the number of smallest elements to be moved;
  * `len`:      length of the array `A`;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
******************************************************************************/
static void qselect(QSELECT_DTYPE *A, size_t n, size_t len, void *buf,
    void *arg) {
  if (n >= len) return;
  for (;;) {
    if (n == 0) {       /* move the minimum element to the beginning */
      QSELECT_DTYPE *tmp = A;
      for (size_t i = 1; i < len; i++)
        if (QSELECT_COMPARE(A + i, tmp, arg) < 0) tmp = A + i;
      if (tmp != A) QSELECT_SWAP(A, tmp);
      return;
    }
    if (n == len - 1) { /* move the maximum element to the end */
      QSELECT_DTYPE *tmp = A;
      for (size_t i = 1; i < len; i++)
        if (QSELECT_COMPARE(A + i, tmp, arg) > 0) tmp = A + i;
      QSELECT_SWAP(A + len - 1, tmp);
      return;
    }

    size_t p;   /* index of the pivot */
    if (len <= 16) p = Hoare_partition(A, n, len, buf, arg);
    else if (n <= len / 6) p = median_of_minima(A, n, len, buf, arg);
    else if (n / 5 >= len / 6) p = median_of_maxima(A, n, len, buf, arg);
    else p = median_of_ninthers(A, len, buf, arg);

    if (p == n) return;
    if (p > n) len = p;
    else {
      p += 1;
      A += p;
      n -= p;
      len -= p;
    }
  }
}

#endif

