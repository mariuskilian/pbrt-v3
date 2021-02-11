
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

#ifndef PBRT_ACCELERATORS_OCTREEBASICACCEL_H
#define PBRT_ACCELERATORS_OCTREEBASICACCEL_H

// accelerators/octree.h*
#include "pbrt.h"
#include "primitive.h"
#include "custom_params.h"
#include <array>

namespace pbrt {

static int MAX_PRIMS;
static Float PRM_THRESH;
static Float VOL_THRESH;

// Size of chunk array of type BITFIELD_X below

#if defined (CHUNKSIZE64)
  const int bytes_free = 56;
#elif defined (CHUNKSIZE128)
  const int bytes_free = 120;
#endif
#if defined (BFSIZE8)
  const int bytes_per_bf = 1;
#elif defined (BFSIZE16)
  const int bytes_per_bf = 2;
#elif defined (BFSIZE32)
  const int bytes_per_bf = 4;
#elif defined (BFSIZE64)
  const int bytes_per_bf = 8;
#endif
const int BFS_CHUNK_DEPTH = bytes_free / bytes_per_bf;

// DON'T change these!
const int DFS_CHUNK_DEPTH = BFS_CHUNK_DEPTH / 2;
const int NUM_SETS_PER_BITFIELD = sizeof(bftype);
const int BFS_NUM_SETS_PER_CHUNK = NUM_SETS_PER_BITFIELD * BFS_CHUNK_DEPTH;
const int DFS_NUM_SETS_PER_CHUNK = NUM_SETS_PER_BITFIELD * DFS_CHUNK_DEPTH;

struct ChildHit { uint32_t idx; float tMin; };
struct ChildTraversal { std::array<ChildHit, 4> nodes; uint32_t size; };

Vector3f BoundsHalf(Bounds3f &b);
Bounds3f DivideBounds(Bounds3f b, int idx, Vector3f b_half);
ChildTraversal FindTraversalOrder(const Ray &ray, Bounds3f &b, Float tMin, Vector3f &invDir);
int Rank(bftype bits, int n = bfsize);
bool BoundsContainPrim(Bounds3f &b, std::shared_ptr<Primitive> p);
bool BoundsContainPoint(Bounds3f &b, Point3f &p);
bool MakeLeafNode(Bounds3f &b, std::vector<std::shared_ptr<Primitive>>);

// OcteeAccel Declarations

class OctreeBasicAccel : public Aggregate {
  public:
    // KdTreeAccel Public Methods
    OctreeBasicAccel();
    OctreeBasicAccel(std::vector<std::shared_ptr<Primitive>> p, int maxPrims = 32, float prm_thresh = 0.9, float vol_thresh = 0.9);
    Bounds3f WorldBound() const { return wb; }
    std::vector<uint32_t> &Nodes() { return nodes; }
    std::vector<uint32_t> &Sizes() { return sizes; }
    std::vector<std::shared_ptr<Primitive>> &Leaves() { return leaves; }
    ~OctreeBasicAccel();
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    float IntersectMetric(const Ray &ray, metric m) const;
    bool IntersectP(const Ray &ray) const;

  private:

    // OctreeAccel Private Data
    std::vector<std::shared_ptr<Primitive>> primitives;
    Bounds3f wb; // World Bounds
    
    void Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth);
    void RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t offset, Bounds3f bounds, Float tMin, bool &hit) const;
    void lh_dump(const char *path);
    void lh_dump_rec(FILE *f, uint32_t *vcnt_, int offset, Bounds3f bounds);

    std::vector<uint32_t> nodes, sizes; 
    std::vector<std::shared_ptr<Primitive>> leaves;
};

std::shared_ptr<OctreeBasicAccel> CreateOctreeBasicAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);
    
} // namespace pbrt

#endif  // PBRT_ACCELERATORS_KDTREEACCEL_H
