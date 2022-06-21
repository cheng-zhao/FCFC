/*******************************************************************************
* define_simd.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_SIMD_H__
#define __DEFINE_SIMD_H__

#include <stdint.h>

#ifdef WITH_SIMD
  #if !defined(WITH_AVX512) && defined(__AVX512F__)
    #define WITH_AVX512
  #endif
  #if !defined(WITH_AVX2) && defined(__AVX2__)
    #define WITH_AVX2
  #endif
  #if !defined(WITH_AVX) && defined(__AVX__)
    #define WITH_AVX
  #endif
#endif

#ifdef FCFC_SIMD
  #undef FCFC_SIMD
#endif
#ifdef FCFC_SIMD_FMA
  #undef FCFC_SIMD_FMA
#endif

#define FCFC_SIMD_NONE          0
#define FCFC_SIMD_AVX           1
#define FCFC_SIMD_AVX2          2
#define FCFC_SIMD_AVX512        3


#if defined(WITH_AVX512)
  #ifndef __AVX512F__
    #error `WITH_AVX512` is set but `__AVX512F__` is unavailable
  #endif

  #define FCFC_SIMD             FCFC_SIMD_AVX512
  #define FCFC_SIMD_FMA

  #ifdef __AVX512DQ__
    #define FCFC_SIMD_AVX512DQ
  #endif
#elif defined(WITH_AVX2)
  #ifndef __AVX2__
    #error `WITH_AVX2` is set but `__AVX2__` is unavailable
  #endif

  #define FCFC_SIMD             FCFC_SIMD_AVX2
  #ifdef __FMA__
    #define FCFC_SIMD_FMA
  #endif
#elif defined(WITH_AVX)
  #ifndef __AVX__
    #error `WITH_AVX` is set but `__AVX__` is unavailable
  #endif

  #define FCFC_SIMD             FCFC_SIMD_AVX
  #ifdef __FMA__
    #define FCFC_SIMD_FMA
  #endif
#else
  #define FCFC_SIMD             FCFC_SIMD_NONE
#endif

#if     FCFC_SIMD != FCFC_SIMD_NONE
#include <immintrin.h>
#endif


/* Constants and function aliases. */
#if     FCFC_SIMD == FCFC_SIMD_AVX  ||  FCFC_SIMD == FCFC_SIMD_AVX2

  #define FCFC_MEMALIGN_BYTE            32
  #define FCFC_SIMD_BYTES               32

  typedef __m256i                       vec_i;
  #define SIMD_ZEROS_INT                _mm256_setzero_si256
  #define SIMD_AND_INT(a,b)             _mm256_and_si256(a,b)
  #define SIMD_ANDNOT_INT(a,b)          _mm256_andnot_si256(a,b)
  #if   FCFC_SIMD  ==  FCFC_SIMD_AVX2
    #define SIMD_MVMSK_INT(a)           _mm256_movemask_epi8(a)
  #endif

  #ifdef        SINGLE_PREC
    #define FCFC_NUM_REAL               8
    #define FCFC_NUM_MASK               0xFFFFFFFFFFFFFFF8ULL
    #define FCFC_REM_MASK               0xFF

    typedef __m256                      vec_r;
    typedef __m256i                     vec_r2i;
    typedef __m256                      vec_msk;
    #define SIMD_ZEROS_REAL             _mm256_setzero_ps
    #define SIMD_SET1_REAL(a)           _mm256_set1_ps(a)
    #define SIMD_LOAD_REAL(a)           _mm256_load_ps(a)
    #define SIMD_LOADU_REAL(a)          _mm256_loadu_ps(a)
    #define SIMD_STORE_REAL(m,a)        _mm256_store_ps(m,a)
    #define SIMD_STOREU_REAL(m,a)       _mm256_storeu_ps(m,a)
    #define SIMD_MSKLOAD_REAL(a,m)      _mm256_maskload_ps(a,m)
    #define SIMD_MSKSTORE_REAL(a,m,v)   _mm256_maskstore_ps(a,m,v)
    #define SIMD_FLOOR_REAL(a)          _mm256_floor_ps(a)
    #define SIMD_R2I(a)                 _mm256_cvttps_epi32(a)
    #define SIMD_CMP_REAL(a,b,i)        _mm256_cmp_ps(a,b,i)
    #define SIMD_AND_REAL(a,b)          _mm256_and_ps(a,b)
    #define SIMD_OR_REAL(a,b)           _mm256_or_ps(a,b)
    #define SIMD_ANDNOT_REAL(a,b)       _mm256_andnot_ps(a,b)
    #define SIMD_XOR_REAL(a,b)          _mm256_xor_ps(a,b)
    #define SIMD_ADD_REAL(a,b)          _mm256_add_ps(a,b)
    #define SIMD_SUB_REAL(a,b)          _mm256_sub_ps(a,b)
    #define SIMD_ABS_REAL(a)            SIMD_ANDNOT_REAL(SIMD_SET1_REAL(-0.0),a)
    #define SIMD_NEG_REAL(a)            SIMD_XOR_REAL(SIMD_SET1_REAL(-0.0),a)
    #define SIMD_MUL_REAL(a,b)          _mm256_mul_ps(a,b)
    #define SIMD_DIV_REAL(a,b)          _mm256_div_ps(a,b)
    #define SIMD_MVMSK_REAL(a)          _mm256_movemask_ps(a)
    #define SIMD_FMADD_REAL(a,b,c)      _mm256_fmadd_ps(a,b,c)
    #define SIMD_CAST_R2I(a)            _mm256_castps_si256(a)
    #define SIMD_CAST_I2R(a)            _mm256_castsi256_ps(a)
    #define SIMD_SET1_INT(a)            _mm256_set1_epi32(a)
    #define SIMD_SETR_INT               _mm256_setr_epi32
    #if FCFC_SIMD  ==  FCFC_SIMD_AVX2
      #define SIMD_SUB_R2I(a,b)         _mm256_sub_epi32(a,b)
      #define SIMD_MSK_GATHER_REAL(a,b,c,d,e)                           \
                                        _mm256_mask_i32gather_ps(a,b,c,d,e)
      #define SIMD_MSK_R2I_GATHER_INT(a,b,c,d,e)                        \
                                        _mm256_mask_i32gather_epi32(a,b,c,d,e)
      #define SIMD_CMPGT_INT(a,b)       _mm256_cmpgt_epi32(a,b)
      #define SIMD_SUB_INT(a,b)         _mm256_sub_epi32(a,b)
    #endif
  #else      /* DOUBLE_PREC */
    #define FCFC_NUM_REAL               4
    #define FCFC_NUM_MASK               0xFFFFFFFFFFFFFFFCULL
    #define FCFC_REM_MASK               0xF

    typedef __m256d                     vec_r;
    typedef __m128i                     vec_r2i;
    typedef __m256d                     vec_msk;
    #define SIMD_ZEROS_REAL             _mm256_setzero_pd
    #define SIMD_SET1_REAL(a)           _mm256_set1_pd(a)
    #define SIMD_LOAD_REAL(a)           _mm256_load_pd(a)
    #define SIMD_LOADU_REAL(a)          _mm256_loadu_pd(a)
    #define SIMD_STORE_REAL(m,a)        _mm256_store_pd(m,a)
    #define SIMD_STOREU_REAL(m,a)       _mm256_storeu_pd(m,a)
    #define SIMD_MSKLOAD_REAL(a,m)      _mm256_maskload_pd(a,m)
    #define SIMD_MSKSTORE_REAL(a,m,v)   _mm256_maskstore_pd(a,m,v)
    #define SIMD_FLOOR_REAL(a)          _mm256_floor_pd(a)
    #define SIMD_R2I(a)                 _mm256_cvttpd_epi32(a)
    #define SIMD_CMP_REAL(a,b,i)        _mm256_cmp_pd(a,b,i)
    #define SIMD_AND_REAL(a,b)          _mm256_and_pd(a,b)
    #define SIMD_OR_REAL(a,b)           _mm256_or_pd(a,b)
    #define SIMD_ANDNOT_REAL(a,b)       _mm256_andnot_pd(a,b)
    #define SIMD_XOR_REAL(a,b)          _mm256_xor_pd(a,b)
    #define SIMD_ADD_REAL(a,b)          _mm256_add_pd(a,b)
    #define SIMD_SUB_REAL(a,b)          _mm256_sub_pd(a,b)
    #define SIMD_ABS_REAL(a)            SIMD_ANDNOT_REAL(SIMD_SET1_REAL(-0.0),a)
    #define SIMD_NEG_REAL(a)            SIMD_XOR_REAL(SIMD_SET1_REAL(-0.0),a)
    #define SIMD_MUL_REAL(a,b)          _mm256_mul_pd(a,b)
    #define SIMD_DIV_REAL(a,b)          _mm256_div_pd(a,b)
    #define SIMD_MVMSK_REAL(a)          _mm256_movemask_pd(a)
    #define SIMD_FMADD_REAL(a,b,c)      _mm256_fmadd_pd(a,b,c)
    #define SIMD_CAST_R2I(a)            _mm256_castpd_si256(a)
    #define SIMD_CAST_I2R(a)            _mm256_castsi256_pd(a)
    #define SIMD_SET1_INT(a)            _mm256_set1_epi64x(a)
    #define SIMD_SETR_INT               _mm256_setr_epi64x
    #define SIMD_SUB_R2I(a,b)           _mm_sub_epi32(a,b)
    #if FCFC_SIMD  ==  FCFC_SIMD_AVX2
      #define SIMD_MSK_GATHER_REAL(a,b,c,d,e)                           \
                                        _mm256_mask_i64gather_pd(a,b,c,d,e)
      #define SIMD_MSK_R2I_GATHER_INT(a,b,c,d,e)                        \
                                        _mm256_mask_i32gather_epi64(a,b,c,d,e)
      #define SIMD_CMPGT_INT(a,b)       _mm256_cmpgt_epi64(a,b)
      #define SIMD_SUB_INT(a,b)         _mm256_sub_epi64(a,b)
    #endif
  #endif

#elif   FCFC_SIMD == FCFC_SIMD_AVX512

  #define FCFC_MEMALIGN_BYTE            64
  #define FCFC_SIMD_BYTES               64

  typedef __m512i                       vec_i;
  #define SIMD_ZEROS_INT                _mm512_setzero_si512
  #define SIMD_ANDNOT_INT(a,b)          _mm512_andnot_si512(a,b)

  #ifdef        SINGLE_PREC
    #define FCFC_NUM_REAL               16
    #define FCFC_NUM_MASK               0xFFFFFFFFFFFFFFF0ULL
    #define FCFC_REM_MASK               0xFFFF

    typedef __m512                      vec_r;
    typedef __m512i                     vec_r2i;
    typedef __mmask16                   vec_msk;
    typedef __m512d                     vec_r64;
    #define SIMD_ZEROS_REAL             _mm512_setzero_ps
    #define SIMD_SET1_REAL(a)           _mm512_set1_ps(a)
    #define SIMD_LOAD_REAL(a)           _mm512_load_ps(a)
    #define SIMD_LOADU_REAL(a)          _mm512_loadu_ps(a)
    #define SIMD_STORE_REAL(m,a)        _mm512_store_ps(m,a)
    #define SIMD_STOREU_REAL(m,a)       _mm512_storeu_ps(m,a)
    #define SIMD_MSK0_LOADU_REAL(k,a)   _mm512_maskz_loadu_ps(k,a)
    #define SIMD_R2I(a)                 _mm512_cvttps_epi32(a)
    #define SIMD_CMP_REAL(a,b,i)        _mm512_cmp_ps_mask(a,b,i)
    #define SIMD_MSK_CMP_REAL(k,a,b,i)  _mm512_mask_cmp_ps_mask(k,a,b,i)
    #define SIMD_MSK_MOV_REAL(s,k,a)    _mm512_mask_mov_ps(s,k,a)
    #define SIMD_MSK0_MOV_REAL(k,a)     _mm512_maskz_mov_ps(k,a)
    #define SIMD_MSKSTORE_REAL(m,k,a)   _mm512_mask_store_ps(m,k,a)
    #define SIMD_ADD_REAL(a,b)          _mm512_add_ps(a,b)
    #define SIMD_SUB_REAL(a,b)          _mm512_sub_ps(a,b)
    #define SIMD_XOR_REAL(a,b)          _mm512_xor_ps(a,b)
    #define SIMD_ABS_REAL(a)            _mm512_abs_ps(a)
    #define SIMD_MUL_REAL(a,b)          _mm512_mul_ps(a,b)
    #define SIMD_MSK_SUB_REAL(s,k,a,b)  _mm512_mask_sub_ps(s,k,a,b)
    #define SIMD_MSK_ADD_REAL(s,k,a,b)  _mm512_mask_add_ps(s,k,a,b)
    #define SIMD_DIV_REAL(a,b)          _mm512_div_ps(a,b)
    #define SIMD_MSK0_MUL_REAL(k,a,b)   _mm512_maskz_mul_ps(k,a,b)
    #define SIMD_FMADD_REAL(a,b,c)      _mm512_fmadd_ps(a,b,c)
    #define SIMD_AND_INT(a,b)           _mm512_and_epi32(a,b)
    #define SIMD_MSK_SUB_INT(s,k,a,b)   _mm512_mask_sub_epi32(s,k,a,b)
    #define SIMD_SET1_INT(a)            _mm512_set1_epi32(a)
    #define SIMD_SET1_I64(a)            _mm512_set1_epi64(a)
    #define SIMD_ADD_INT(a,b)           _mm512_add_epi32(a,b)
    #define SIMD_ADD_I64(a,b)           _mm512_add_epi64(a,b)
    #define SIMD_ADD_R64(a,b)           _mm512_add_pd(a,b)
    #define SIMD_SUB_R2I(a,b)           _mm512_sub_epi32(a,b)
    #define SIMD_CAST_R2I(a)            _mm512_castps_si512(a)
    #define SIMD_CAST_I2R(a)            _mm512_castsi512_ps(a)
    #define SIMD_MSK_CMPGE_INT(k,a,b)   _mm512_mask_cmpge_epi32_mask(k,a,b)
    #define SIMD_MSK_CMPNEQ_INT(k,a,b)  _mm512_mask_cmpneq_epi32_mask(k,a,b)
    #define SIMD_MUL_INT(a,b)           _mm512_mullo_epi32(a,b)
    #ifdef FCFC_SIMD_AVX512DQ
      #define SIMD_NEG_REAL(a)          _mm512_xor_ps(SIMD_SET1_REAL(-0.0),a)
    #else
      #define SIMD_NEG_REAL(a)          SIMD_CAST_I2R(_mm512_xor_si512( \
          SIMD_CAST_R2I(SIMD_SET1_REAL(-0.0)),SIMD_CAST_R2I(a)))
    #endif
    #define SIMD_MSK0_DIV_ROUND_REAL(k,a,b,r)                           \
                                        _mm512_maskz_div_round_ps(k,a,b,r)
    #define SIMD_MSK_GATHER_REAL(a,b,c,d,e)                             \
                                        _mm512_mask_i32gather_ps(a,b,c,d,e)
    #define SIMD_MSK_R2I_GATHER_INT(a,b,c,d,e)                          \
                                        _mm512_mask_i32gather_epi32(a,b,c,d,e)
    #define SIMD_MSK_GATHER_I64(a,b,c,d,e)                              \
                                        _mm512_mask_i32gather_epi64(a,b,c,d,e)
    #define SIMD_MSK_SCATTER_I64(a,b,c,d,e)                             \
                                        _mm512_mask_i32scatter_epi64(a,b,c,d,e)
    #define SIMD_MSK_GATHER_R64(a,b,c,d,e)                              \
                                        _mm512_mask_i32gather_pd(a,b,c,d,e)
    #define SIMD_MSK_SCATTER_R64(a,b,c,d,e)                             \
                                        _mm512_mask_i32scatter_pd(a,b,c,d,e)
  #else      /* DOUBLE_PREC */
    #define FCFC_NUM_REAL               8
    #define FCFC_NUM_MASK               0xFFFFFFFFFFFFFFF8ULL
    #define FCFC_REM_MASK               0xFF

    typedef __m512d                     vec_r;
    typedef __m256i                     vec_r2i;
    typedef __mmask8                    vec_msk;
    #define SIMD_ZEROS_REAL             _mm512_setzero_pd
    #define SIMD_SET1_REAL(a)           _mm512_set1_pd(a)
    #define SIMD_LOAD_REAL(a)           _mm512_load_pd(a)
    #define SIMD_LOADU_REAL(a)          _mm512_loadu_pd(a)
    #define SIMD_STORE_REAL(m,a)        _mm512_store_pd(m,a)
    #define SIMD_STOREU_REAL(m,a)       _mm512_storeu_pd(m,a)
    #define SIMD_MSK0_LOADU_REAL(k,a)   _mm512_maskz_loadu_pd(k,a)
    #define SIMD_R2I(a)                 _mm512_cvttpd_epi32(a)
    #define SIMD_CMP_REAL(a,b,i)        _mm512_cmp_pd_mask(a,b,i)
    #define SIMD_MSK_CMP_REAL(k,a,b,i)  _mm512_mask_cmp_pd_mask(k,a,b,i)
    #define SIMD_MSK_MOV_REAL(s,k,a)    _mm512_mask_mov_pd(s,k,a)
    #define SIMD_MSK0_MOV_REAL(k,a)     _mm512_maskz_mov_pd(k,a)
    #define SIMD_MSKSTORE_REAL(m,k,a)   _mm512_mask_store_pd(m,k,a)
    #define SIMD_ADD_REAL(a,b)          _mm512_add_pd(a,b)
    #define SIMD_SUB_REAL(a,b)          _mm512_sub_pd(a,b)
    #define SIMD_MSK_SUB_REAL(s,k,a,b)  _mm512_mask_sub_pd(s,k,a,b)
    #define SIMD_MSK_ADD_REAL(s,k,a,b)  _mm512_mask_add_pd(s,k,a,b)
    #define SIMD_ABS_REAL(a)            _mm512_abs_pd(a)
    #define SIMD_MUL_REAL(a,b)          _mm512_mul_pd(a,b)
    #define SIMD_DIV_REAL(a,b)          _mm512_div_pd(a,b)
    #define SIMD_MSK0_MUL_REAL(k,a,b)   _mm512_maskz_mul_pd(k,a,b)
    #define SIMD_FMADD_REAL(a,b,c)      _mm512_fmadd_pd(a,b,c)
    #define SIMD_AND_INT(a,b)           _mm512_and_epi64(a,b)
    #define SIMD_SET1_INT(a)            _mm512_set1_epi64(a)
    #define SIMD_SET1_I64(a)            _mm512_set1_epi64(a)
    #define SIMD_ADD_INT(a,b)           _mm512_add_epi64(a,b)
    #define SIMD_ADD_I64(a,b)           _mm512_add_epi64(a,b)
    #define SIMD_MSK_SUB_INT(s,k,a,b)   _mm512_mask_sub_epi64(s,k,a,b)
    #define SIMD_SUB_R2I(a,b)           _mm256_sub_epi32(a,b)
    #define SIMD_CAST_R2I(a)            _mm512_castpd_si512(a)
    #define SIMD_CAST_I2R(a)            _mm512_castsi512_pd(a)
    #define SIMD_MSK_CMPGE_INT(k,a,b)   _mm512_mask_cmpge_epi64_mask(k,a,b)
    #define SIMD_MSK_CMPNEQ_INT(k,a,b)  _mm512_mask_cmpneq_epi64_mask(k,a,b)
    #ifdef FCFC_SIMD_AVX512DQ
      #define SIMD_NEG_REAL(a)          _mm512_xor_pd(SIMD_SET1_REAL(-0.0),a)
      #define SIMD_MUL_INT(a,b)         _mm512_mullo_epi64(a,b)
    #else
      #define SIMD_NEG_REAL(a)          SIMD_CAST_I2R(_mm512_xor_si512( \
          SIMD_CAST_R2I(SIMD_SET1_REAL(-0.0)),SIMD_CAST_R2I(a)))
      #define SIMD_MUL_INT(a,b)         _mm512_mullo_epi32(a,b)
    #endif
    #define SIMD_MSK0_DIV_ROUND_REAL(k,a,b,r)                           \
                                        _mm512_maskz_div_round_pd(k,a,b,r)
    #define SIMD_MSK_GATHER_REAL(a,b,c,d,e)                             \
                                        _mm512_mask_i64gather_pd(a,b,c,d,e)
    #define SIMD_MSK_R2I_GATHER_INT(a,b,c,d,e)                          \
                                        _mm512_mask_i32gather_epi64(a,b,c,d,e)
    #define SIMD_MSK_GATHER_I64(a,b,c,d,e)                              \
                                        _mm512_mask_i64gather_epi64(a,b,c,d,e)
    #define SIMD_MSK_SCATTER_I64(a,b,c,d,e)                             \
                                        _mm512_mask_i64scatter_epi64(a,b,c,d,e)
    #define SIMD_MSK_GATHER_R64(a,b,c,d,e)                              \
                                        _mm512_mask_i64gather_pd(a,b,c,d,e)
    #define SIMD_MSK_SCATTER_R64(a,b,c,d,e)                             \
                                        _mm512_mask_i64scatter_pd(a,b,c,d,e)
  #endif

#endif


/* Integer types with the same width as floating-point numbers. */
#ifdef          SINGLE_PREC
  typedef int32_t                       int_wr;
  typedef uint32_t                      uint_wr;
#else        /* DOUBLE_PREC */
  typedef int64_t                       int_wr;
  typedef uint64_t                      uint_wr;
#endif

#endif
