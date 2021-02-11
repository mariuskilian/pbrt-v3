#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CUSTOM_PARAMS_H
#define PBRT_CUSTOM_PARAMS_H

#include "pbrt.h"

namespace pbrt {

enum class metric { PRIMITIVES, NODES, LEAFNODES, TIME };

#define COUNT_STATS //debug

#if defined (COUNT_STATS)
  static const bool count_stats = true;
  static std::vector<int> *prim_isects_per_ray;
  static std::vector<int> *node_isects_per_ray;
  static std::vector<int> *leafNode_isects_per_ray;
  static std::vector<int> *chunk_isects_per_ray;
#else
  static const bool count_stats = false;
#endif

#if defined (REL_KEYS)
  static const bool relative_keys = true;
#else
  static const bool relative_keys = false;
#endif
    
#if defined (CHUNKSIZE64)
  static const int chunksize = 64;
#elif defined (CHUNKSIZE128)
  static const int chunksize = 128;
#elif defined (CHUNKSIZE256)
  static const int chunksize = 256;
#endif

#if defined (BFSIZE8)
  typedef uint8_t bftype;
#elif defined (BFSIZE16)
  typedef uint16_t bftype;
#elif defined (BFSIZE32)
  typedef uint32_t bftype;
#elif defined (BFSIZE64)
  typedef uint64_t bftype;
#endif

#if defined (_MSC_VER)
  #include "intrin.h"
  #if defined (BFSIZE8) || defined (BFSIZE16)
    #define popcnt __popcnt16
  #elif defined (BFSIZE32)
    #define popcnt __popcnt
  #elif defined (BFSIZE64)
    #define popcnt __popcnt64
  #endif
#elif defined (__clang__)
  #include "popcntintrin.h"
  #if defined (BFSIZE8) || defined (BFSIZE16) || defined (BFSIZE32)
    #define popcnt _mm_popcnt_u32
  #elif defined (BFSIZE64)
    #define popcnt _mm_popcnt_u64
  #endif
#elif defined (__linux__)
  #if defined (BFSIZE8) || defined (BFSIZE16) || defined (BFSIZE32)
    #define popcnt __builtin_popcount
  #elif defined (BFSIZE64)
    #define popcnt __builtin_popcountll
  #endif
#endif

static const int bfsize = 8 * sizeof(bftype);
static const bftype bftzero = (bftype)0;
static const bftype bftone = (bftype)1;
}

#endif