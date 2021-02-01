#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CUSTOM_PARAMS_H
#define PBRT_CUSTOM_PARAMS_H

#include "pbrt.h"

namespace pbrt {
    
#if defined(CHUNKSIZE64)
  const int chunksize = 64;
#elif defined(CHUNKSIZE128)
  const int chunksize = 128;
#endif

#if defined(BFSIZE8)
  typedef uint8_t bftype;
#elif defined(BFSIZE16)
  typedef uint16_t bftype;
#elif defined(BFSIZE32)
  typedef uint32_t bftype;
#elif defined(BFSIZE64)
  typedef uint64_t bftype;
#endif

#if defined(_MSC_VER)
  #include "intrin.h"
  #if defined(BFSIZE8) || defined(BFSIZE16)
    #define popcnt __popcnt16
  #elif defined(BFSIZE32)
    #define popcnt __popcnt
  #elif defined(BFSIZE64)
    #define popcnt __popcnt64
  #endif
#else
  #if defined(__clang__)
    #include "popcntintrin.h"
  #elif defined (__linux__)
    #include <nmmintrin.h>
  #endif
  #if defined(BFSIZE8) || defined(BFSIZE16) || defined(BFSIZE32)
    #define popcnt _mm_popcnt_u32
  #elif defined(BFSIZE64)
    #define popcnt _mm_popcnt_u64
  #endif
#endif

const int bfsize = 8 * sizeof(bftype);
const bftype bftzero = (bftype)0;
const bftype bftone = (bftype)1;
}

#endif