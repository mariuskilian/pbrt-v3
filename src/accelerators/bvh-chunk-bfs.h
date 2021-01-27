
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

#ifndef PBRT_ACCELERATORS_BVH_CHUNK_BFS_H
#define PBRT_ACCELERATORS_BVH_CHUNK_BFS_H

// accelerators/bvh-chunk-bfs.h*
#include "pbrt.h"
#include "primitive.h"
#include "bvh.h"
#include <atomic>

namespace pbrt {

const bool relative_keys = false;

typedef BVHAccel::SplitMethod SplitMethod;
struct Bounds3k { uint8_t min[3]; uint8_t max[3]; };
struct BVHBFSNodeInfo { Bounds3k bk; uint8_t axis; }; //IDEA: take 2 bits of non-axis min to save axis. 
struct ChildrenBuildInfo { uint32_t idx[2]; Bounds3f b_comp[2]; Bounds3k b_key[2]; };

// BVHChunkBFSAccel Forward Declarations

// BVHAccel Declarations
class BVHChunkBFSAccel : public Aggregate {
  public:
    // BVHAccel Public Types

    // BVHAccel Public Methods
    BVHChunkBFSAccel(std::vector<std::shared_ptr<Primitive>> p,
             int maxPrimsInNode = 1,
             SplitMethod splitMethod = SplitMethod::SAH);
    Bounds3f WorldBound() const;
    ~BVHChunkBFSAccel();
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;

  private:
    struct BVHChunkBFS;

    typedef uint64_t bf_type;
    static const int bf_size = sizeof(bf_type) * 8;
    static const int chunk_depth = 3;
    static const int node_pairs_per_chunk = bf_size * chunk_depth / 2;
    static const int node_pairs_per_bitfield = node_pairs_per_chunk / chunk_depth;
    int num_chunk_layers = 0;

    // BVHAccel Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<std::shared_ptr<Primitive>> primitives;

    std::vector<BVHChunkBFS> bvh_chunks;
    std::vector<BVHBFSNodeInfo> node_info;
    std::vector<uint16_t> sizes;

    BVHAccel *bvh;

    uint32_t Rank(bf_type bitfield[chunk_depth], int n) const;

    Bounds3k FindBoundsKey(Bounds3f b_root, Bounds3f b, Vector3f factor) const;
    Bounds3f FindCompressedBounds(Bounds3f b_root, Bounds3k b_k, Vector3f factor) const;
    ChildrenBuildInfo GetChildrenBuildInfo(Bounds3f b_root, uint32_t node_offset) const;

    void Recurse(uint32_t chunk_offset, uint32_t root_node_idx, Bounds3f b_root, int chunk_layer);

    void lh_dump(const char *path, bool dfs = false);
    void lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds);
    void lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds);
};

std::shared_ptr<BVHChunkBFSAccel> CreateBVHChunkBFSAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_BVH_CHUNK_BFS_H
