
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

const int chunk_depth = 5;

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
      uint32_t leaf_offset;
      std::array<uint8_t, chunk_depth> node_type; // 0 inner node, 1 leaf node
      std::array<uint8_t, chunk_depth> leaf_type; // 0 primitive leaf, 1 chunk leaf
    };
    std::vector<chunk> octree;
    
    Bounds3f octreeDivide(Bounds3f bounds, int idx) const;
    void Recurse(uint32_t root_node_offset, int chunk_idx);  
    void RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t offset, Bounds3f bounds, bool &hit) const;
    void lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds);
    void lh_dump(const char *path);

};

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

} // namespace pbrt

#endif  // PBRT_ACCELERATORS_KDTREEACCEL_H
