
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

namespace pbrt {

// OcteeAccel Declarations
struct OctreeNode;

class OctreeAccel : public Aggregate {
  public:
    // KdTreeAccel Public Methods
    OctreeAccel(std::vector<std::shared_ptr<Primitive>> p);
    Bounds3f WorldBound() const { return bounds; }
    ~OctreeAccel();
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;

  private:
    // OctreeAccel Private Methods
    void buildTree();

    // OctreeAccel Private Data
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::vector<int> primitiveIndices;
    OctreeNode *nodes;
    int nAllocedNodes, nextFreeNode;
    Bounds3f bounds;
};

struct OcToDo {
    const OctreeNode *node;
};

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

} // namespace pbrt

#endif  // PBRT_ACCELERATORS_KDTREEACCEL_H
