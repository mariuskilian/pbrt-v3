
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_OCTREEACCEL_H
#define PBRT_ACCELERATORS_OCTREEACCEL_H

// accelerators/octree.h*
#include "pbrt.h"
#include "primitive.h"
#include "octree-basic.h"
#include <array>

namespace pbrt {

#if defined(_MSC_VER)
    #include "intrin.h"
    #define POPCNT __popcnt64
#elif defined(__clang__)
    #include "popcntintrin.h"
    #define POPCNT _mm_popcnt_u64
#endif

#define BITFIELD_64 // Change this to change bitfield sizes
  #if defined BITFIELD_8
    #define BF_SIZE 8
    #define BITFIELD_TYPE uint8_t
  #elif defined BITFIELD_16
    #define BF_SIZE 16
    #define BITFIELD_TYPE uint16_t
  #elif defined BITFIELD_32
    #define BF_SIZE 32
    #define BITFIELD_TYPE uint32_t
  #elif defined BITFIELD_64
    #define BF_SIZE 64
    #define BITFIELD_TYPE uint64_t
  #endif
const int CHUNK_DEPTH = 5; // Size of chunk array of type BITFIELD_X below
const int BITFIELD_SIZE = BF_SIZE;
const int NUM_SETS_PER_BITFIELD = BITFIELD_SIZE / 8;
const int NUM_SETS_PER_CHUNK = NUM_SETS_PER_BITFIELD * CHUNK_DEPTH;

// OcteeAccel Declarations

class OctreeAccel : public Aggregate {
  public:
    // KdTreeAccel Public Methods
    OctreeAccel(std::vector<std::shared_ptr<Primitive>> p);
    Bounds3f WorldBound() const { return wb; }
    ~OctreeAccel();
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;

  private:

    // OctreeAccel Private Data
    std::vector<std::shared_ptr<Primitive>> primitives;
    Bounds3f wb; // World Bounds
    OctreeBasicAccel oba;

    struct chunk {
      uint32_t child_chunk_offset;
      uint32_t leaf_offset; // with size
      std::array<BITFIELD_TYPE, CHUNK_DEPTH> node_type; // 1 inner node, 0 leaf node
    };
    std::vector<chunk> octree;
    
    Bounds3f octreeDivide(Bounds3f bounds, int idx) const;
    void Recurse(uint32_t root_node_offset, int chunk_idx);  
    void RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t chunk_offset, Bounds3f parent_bounds, bool &hit) const;
    void lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds);
    void lh_dump_dfs(const char *path);
    void lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds);
    void lh_dump(const char *path);

};

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

} // namespace pbrt

#endif  // PBRT_ACCELERATORS_KDTREEACCEL_H
