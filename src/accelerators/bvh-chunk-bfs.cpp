
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


// accelerators/bvh-chunk-bfs.cpp*
#include "accelerators/bvh-chunk-bfs.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>
#include <queue>
#include <cmath>

namespace pbrt {

struct alignas(64) BVHChunkBFSAccel::BVHChunkBFS {
    Bounds3f b_root;
    uint32_t multipliers_offset;
    uint32_t primitive_offset;
    uint32_t child_chunk_offset; // TODO now have 32 bits more than 64B :(
    uint64_t bitfield[chunk_depth];
};
std::vector<Bounds3k> multipliers;
std::vector<uint32_t> sizes;

Bounds3k BVHChunkBFSAccel::FindBoundsKey(Bounds3f b_root, Bounds3f b) {
    Bounds3k b_k = Bounds3k{};
    for (int i = 0; i < 3; i++) {
        Float dist = b_root.pMax[i] - b_root.pMin[i];
        b_k.min[i] = (int)floor(255 * (b.pMin[i] - b_root.pMin[i]) / dist);
        b_k.max[i] = (int)ceil(255 * (b.pMax[i] - b_root.pMin[i]) / dist);
    }
    return b_k;
}

Bounds3f BVHChunkBFSAccel::FindCompressedBounds(Bounds3f b_root, Bounds3k b_k) {
    Bounds3f b;
    for (int i = 0; i < 3; i++) {
        Float dist = b_root.pMax[i] - b_root.pMin[i];
        b.pMin[i] = (dist * b_k.min[i] / 255) + b.pMin[i];
        b.pMax[i] = (dist * b_k.max[i] / 255) + b.pMin[i];
    }
    return b;

}

ChildrenBuildInfo BVHChunkBFSAccel::GetChildrenBuildInfo(
        Bounds3f b_root, uint32_t node_offset) {
    ChildrenBuildInfo cbi = ChildrenBuildInfo{};
    cbi.idx[0] = node_offset + 1;
    cbi.idx[1] = bvh->GetNodes()[node_offset].secondChildOffset;
    for (int child = 0; child < 2; child++) {
        Bounds3f b_child = bvh->GetNodes()[cbi.idx[child]].bounds;
        cbi.b_key[child] = FindBoundsKey(b_root, b_child);
        cbi.b_comp[child] = FindCompressedBounds(b_root, cbi.b_key[child]);
    }
    return cbi;
}

// BVHAccel Method Definitions
BVHChunkBFSAccel::BVHChunkBFSAccel(std::vector<std::shared_ptr<Primitive>> p,
                   int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      splitMethod(splitMethod),
      primitives(std::move(p)) {
    ProfilePhase _(Prof::AccelConstruction);
    // Let BVHAccel build BVH Structure
    bvh = new BVHAccel(p, maxPrimsInNode, splitMethod);
    if (!bvh->GetNodes()) return;
    // Chunk struct for building the chunks
    struct Chunk { uint32_t root_node_idx; Bounds3f b_root; };
    uint32_t chunk_idx = 0;
    std::queue<Chunk> chunk_q;
    bvh_chunks.push_back(BVHChunkBFS{});
    std::vector<std::shared_ptr<Primitive>> orderedPrims;
    // Build a full chunk within each iteration of the loop
    while (!chunk_q.empty()) {
        Chunk c = chunk_q.front();
        // === Building Chunk ===
        // Set meta information of chunk
        bvh_chunks[chunk_idx].b_root = c.b_root;
        bvh_chunks[chunk_idx].multipliers_offset = multipliers.size();
        bvh_chunks[chunk_idx].primitive_offset = orderedPrims.size();
        bvh_chunks[chunk_idx].child_chunk_offset = bvh_chunks.size();
        // Make a node queue and push the root nodes 2 children
        struct BVHChunkBFSNodeBuild { uint32_t node_idx; Bounds3f b_comp; Bounds3k b_key; };
        std::queue<BVHChunkBFSNodeBuild> node_q;
        ChildrenBuildInfo cbi = GetChildrenBuildInfo(c.b_root, c.root_node_idx);
        for (int a = 0; a < 2; a++) 
            node_q.push(BVHChunkBFSNodeBuild{cbi.idx[a], cbi.b_comp[a], cbi.b_key[a]});
        // Initialize variables needed to fill bitfield
        int chunk_fill = 1; // number of pairs of nodes processed and/or in queue
        int nodes_processed = 0; // number of nodes already processed
        bool firstLeafInChunk = true;
        // Fill the bitfield of this chunk, processing one node per iteration
        while (!node_q.empty()) {
            // Variables for easier access to information
            BVHChunkBFSNodeBuild build_node = node_q.front();
            LinearBVHNode *node = &bvh->GetNodes()[build_node.node_idx];
            bool is_inner_node = node->nPrimitives == 0;
            // Find correct idx and bit position in bitfield
            int idx = nodes_processed / bf_size;
            int bit_pos = nodes_processed % bf_size;
            // When starting a new index, make sure the initial value of the bitcode is 0
            if (bit_pos == 0) bvh_chunks[chunk_idx].bitfield[idx] = (bf_type)0;
            // Determine node type, and process accordingly
            if (is_inner_node) {
                bvh_chunks[chunk_idx].bitfield[idx] |= (bf_type)1 << bit_pos;
                if (chunk_fill >= node_pairs_per_chunk) {
                    // Chunk Ptr Inner Node
                    chunk_q.push(Chunk{build_node.node_idx, node->bounds});
                    bvh_chunks.push_back(BVHChunkBFS{});
                } else {
                    // Real Inner Node
                    ChildrenBuildInfo cbi = GetChildrenBuildInfo(build_node.b_comp, build_node.node_idx);
                    for (int a = 0; a < 2; a++) 
                        node_q.push(BVHChunkBFSNodeBuild{cbi.idx[a], cbi.b_comp[a], cbi.b_key[a]});
                    chunk_fill++;
                }
            } else {
                // Leaf Node
                if (firstLeafInChunk) {
                    sizes.push_back(0);
                    firstLeafInChunk = false;
                }
                sizes.push_back(sizes.back() + node->nPrimitives);
                for (int i = 0; i < node->nPrimitives; i++)
                    orderedPrims.push_back(bvh->GetPrimitives()[node->primitivesOffset + i]);
            }
            multipliers.push_back(build_node.b_key);
            // Update node counter and -queue
            nodes_processed++;
            node_q.pop();
        }
        // Update chunk counter and -queue
        chunk_idx++;
        chunk_q.pop();
    }
    primitives.swap(orderedPrims);
}

Bounds3f BVHChunkBFSAccel::WorldBound() const {
    return bvh->WorldBound();
}

BVHChunkBFSAccel::~BVHChunkBFSAccel() {
    bvh->~BVHAccel();
    FreeAligned(bvh);
}

std::shared_ptr<BVHChunkBFSAccel> CreateBVHChunkBFSAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    std::string splitMethodName = ps.FindOneString("splitmethod", "sah");
    SplitMethod splitMethod;
    if (splitMethodName == "sah")
        splitMethod = SplitMethod::SAH;
    else if (splitMethodName == "hlbvh")
        splitMethod = SplitMethod::HLBVH;
    else if (splitMethodName == "middle")
        splitMethod = SplitMethod::Middle;
    else if (splitMethodName == "equal")
        splitMethod = SplitMethod::EqualCounts;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                splitMethodName.c_str());
        splitMethod = SplitMethod::SAH;
    }

    int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return std::make_shared<BVHChunkBFSAccel>(std::move(prims), maxPrimsInNode, splitMethod);
}

bool BVHChunkBFSAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    if (true) return false;
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
    return hit;
}

bool BVHChunkBFSAccel::IntersectP(const Ray &ray) const {
    if (true) return false;
    ProfilePhase p(Prof::AccelIntersectP);
    return false;
}

// VISUALISATION
void BVHChunkBFSAccel::lh_dump(const char *path, bool dfs) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    if (dfs) lh_dump_rec_dfs(f, &vcnt, 0, WorldBound());
    else lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void BVHChunkBFSAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    BVHChunkBFS c = bvh_chunks[chunk_offset];

    std::queue<Bounds3f> bounds_q;
    bounds_q.push(bounds);

    int child_chunk_cnt = 0;
    int current_pair_idx = 0;

    int num_node_pairs = 1;
    for (int idx = 0; idx < chunk_depth; idx++) {
        for (bf_type pair = 0; pair < num_node_pairs; pair++) {
            if (current_pair_idx++ >= num_node_pairs) break;

            
        }
    }
}

void BVHChunkBFSAccel::lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {

}

}  // namespace pbrt
