
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

#include <algorithm>
#include <array>
#include <cmath>
#include <queue>

#include "interaction.h"
#include "parallel.h"
#include "paramset.h"
#include "stats.h"

#if defined(_MSC_VER)
#include "intrin.h"
#define POPCNT __popcnt64
#elif defined(__clang__)
#include "popcntintrin.h"
#define POPCNT _mm_popcnt_u64
#endif

namespace pbrt {

STAT_COUNTER("BVH-BFS/Total nodes", total_nodes);

struct alignas(64) BVHChunkBFSAccel::BVHChunkBFS {
    Bounds3f b_root;
    uint32_t primitive_offset;
    uint32_t sizes_offset;
    uint32_t node_info_offset;
    uint32_t child_chunk_offset;  // TODO now have 32 bits more than 64B :(
    bf_type bitfield[chunk_depth];
};

uint32_t BVHChunkBFSAccel::Rank(bf_type bitfield[chunk_depth], int n) const {
    int idx = n;
    int count = 0;
    bf_type bits;
    for (int i = 0; i < chunk_depth; i++) {
        bits = bitfield[i];
        if ((idx -= bf_size) < 0) break;
        count += POPCNT(bits);
    }
    return count + POPCNT(bits & (((bf_type)1 << n) - (bf_type)1));
}

Bounds3k BVHChunkBFSAccel::FindBoundsKey(Bounds3f b_root, Bounds3f b,
                                         Vector3f b2k) const {
    Bounds3k b_k = Bounds3k{};
    for (int i = 0; i < 3; i++) {
        // b2k = 255 / dist
        b_k.min[i] = (uint8_t)floor(b2k[i] * (b.pMin[i] - b_root.pMin[i]));
        b_k.max[i] = (uint8_t)ceil(b2k[i] * (b.pMax[i] - b_root.pMin[i]));
    }
    return b_k;
}

Bounds3f BVHChunkBFSAccel::FindCompressedBounds(Bounds3f b_root, Bounds3k b_k,
                                                Vector3f k2b) const {
    Bounds3f b;
    for (int i = 0; i < 3; i++) {
        // k2b = dist / 255
        b.pMin[i] = (k2b[i] * b_k.min[i]) + b_root.pMin[i];
        b.pMax[i] = (k2b[i] * b_k.max[i]) + b_root.pMin[i];
    }
    return b;
}

struct BoundConversionFactors {
    Vector3f b2k;
    Vector3f k2b;
};
BoundConversionFactors CalcBoundConversionFactors(Bounds3f b) {
    Vector3f b2k, k2b;
    for (int i = 0; i < 3; i++) {
        b2k[i] = 255 / (b.pMax[i] - b.pMin[i]);
        k2b[i] = 1 / b2k[i];
    }
    return BoundConversionFactors{b2k, k2b};
}

ChildrenBuildInfo BVHChunkBFSAccel::GetChildrenBuildInfo(
    Bounds3f b_root, uint32_t node_offset) const {
    ChildrenBuildInfo cbi = ChildrenBuildInfo{};
    cbi.idx[0] = node_offset + 1;
    cbi.idx[1] = bvh->GetNodes()[node_offset].secondChildOffset;
    BoundConversionFactors bcf = CalcBoundConversionFactors(b_root);
    for (int child = 0; child < 2; child++) {
        Bounds3f b_child = bvh->GetNodes()[cbi.idx[child]].bounds;
        cbi.b_key[child] = FindBoundsKey(b_root, b_child, bcf.b2k);
        cbi.b_comp[child] =
            FindCompressedBounds(b_root, cbi.b_key[child], bcf.k2b);
    }
    return cbi;
}

// BVHAccel Method Definitions
BVHChunkBFSAccel::BVHChunkBFSAccel(std::vector<std::shared_ptr<Primitive>> p,
                                   int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      splitMethod(splitMethod) {
    ProfilePhase _(Prof::AccelConstruction);
    // Let BVHAccel build BVH Structure
    printf("BVHBFS: Starting original BVH build!\n");
    bvh = new BVHAccel(p, maxPrimsInNode, splitMethod);
    if (!bvh->GetNodes()) return;
    printf("BVHBFS: Original BVH build done!\n");
    // Chunk struct for building the chunks
    bvh_chunks.push_back(BVHChunkBFS{});
    printf("BVHBFS: Starting build!\n");
    total_nodes += 1;
    Recurse(0, 0, WorldBound(), 1);
    printf("BVHBFS: Build done!\n");
    // printf("Starting vis!\n");
    // Visualisation:
    // lh_dump("bvh_vis_bfs.obj");
}

void BVHChunkBFSAccel::Recurse(uint32_t chunk_offset, uint32_t root_node_idx,
                               Bounds3f b_root, int chunk_layer) {
    if (chunk_layer > num_chunk_layers) num_chunk_layers = chunk_layer;
    struct Chunk {
        uint32_t root_node_idx;
        Bounds3f b_root;
    };
    std::vector<Chunk> chunk_ptr_nodes;
    // === Building Chunk ===
    // Set meta information of chunk
    bvh_chunks[chunk_offset].b_root = b_root;
    bvh_chunks[chunk_offset].node_info_offset = node_info.size();
    bvh_chunks[chunk_offset].sizes_offset = sizes.size();
    bvh_chunks[chunk_offset].primitive_offset = primitives.size();
    bvh_chunks[chunk_offset].child_chunk_offset = bvh_chunks.size();
    // Make a node queue and push the root nodes 2 children
    struct BVHChunkBFSNodeBuild {
        uint32_t node_idx;
        Bounds3f b_comp;
        Bounds3k b_key;
    };
    std::queue<BVHChunkBFSNodeBuild> node_q;
    ChildrenBuildInfo cbi = GetChildrenBuildInfo(b_root, root_node_idx);
    for (int i = 0; i < 2; i++)
        node_q.push(
            BVHChunkBFSNodeBuild{cbi.idx[i], cbi.b_comp[i], cbi.b_key[i]});
    // Initialize variables needed to fill bitfield
    int chunk_fill = 1;  // number of pairs of nodes processed and/or in queue
    int nodes_processed = 0;  // number of nodes already processed
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
        // When starting a new index, make sure the initial value of the bitcode
        // is 0
        if (bit_pos == 0) bvh_chunks[chunk_offset].bitfield[idx] = (bf_type)0;
        // Determine node type, and process accordingly
        if (is_inner_node) {
            bvh_chunks[chunk_offset].bitfield[idx] |= ((bf_type)1 << bit_pos);
            if (chunk_fill < node_pairs_per_chunk) {
                // Real Inner Node
                Bounds3f root_bounds =
                    (relative_keys) ? build_node.b_comp : b_root;
                ChildrenBuildInfo cbi =
                    GetChildrenBuildInfo(root_bounds, build_node.node_idx);
                for (int i = 0; i < 2; i++)
                    node_q.push(BVHChunkBFSNodeBuild{cbi.idx[i], cbi.b_comp[i],
                                                     cbi.b_key[i]});
                chunk_fill++;
            } else {
                // Chunk Ptr Inner Node
                chunk_ptr_nodes.push_back(
                    Chunk{build_node.node_idx, node->bounds});
                bvh_chunks.push_back(BVHChunkBFS{});
            }
        } else {
            // Leaf Node
            if (firstLeafInChunk) {
                sizes.push_back(0);
                firstLeafInChunk = false;
            }
            sizes.push_back(sizes.back() + node->nPrimitives);
            for (int i = 0; i < node->nPrimitives; i++)
                primitives.push_back(
                    bvh->GetPrimitives()[node->primitivesOffset + i]);
        }
        node_info.push_back(BVHBFSNodeInfo{build_node.b_key, node->axis});
        // Update node counter and -queue
        nodes_processed++;
        node_q.pop();
    }  // nodes_q
    // Update chunk counter and -queue
    total_nodes += nodes_processed;

    for (int i = 0; i < chunk_ptr_nodes.size(); i++) {
        Recurse(bvh_chunks[chunk_offset].child_chunk_offset + i,
                chunk_ptr_nodes[i].root_node_idx, chunk_ptr_nodes[i].b_root,
                chunk_layer + 1);
    }
}

Bounds3f BVHChunkBFSAccel::WorldBound() const { return bvh->WorldBound(); }

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
    return std::make_shared<BVHChunkBFSAccel>(std::move(prims), maxPrimsInNode,
                                              splitMethod);
}

bool BVHChunkBFSAccel::Intersect(const Ray &ray,
                                 SurfaceInteraction *isect) const {
    if (!bvh->GetNodes()) return false;
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    // Follow ray through BVH nodes to find primitive intersections
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
    };
    // Variables that update whenever a new chunk is entered
    uint32_t chunk_offset = 99;  // Set to non-0 value to trigger variable updates in first iteration
    BVHChunkBFS current_chunk;
    BVHChunkBFSNode current_node;
    Vector3f k2b;
    // Initialize node stack
    BVHChunkBFSNode node_stack[64];
    uint32_t node_stack_offset = 0;
    if (dirIsNeg[bvh->GetNodes()[0].axis]) {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 0};
        current_node = BVHChunkBFSNode{0, 1};
    } else {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 1};
        current_node = BVHChunkBFSNode{0, 0};
    }
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_idx != chunk_offset) {
            chunk_offset = current_node.chunk_idx;
            current_chunk = bvh_chunks[chunk_offset];
            k2b = CalcBoundConversionFactors(current_chunk.b_root).k2b;
        }

        uint32_t rank = Rank(current_chunk.bitfield, current_node.node_idx);

        int idx = current_node.node_idx / bf_size;
        int bit_pos = current_node.node_idx % bf_size;
        bool is_inner_node = ((current_chunk.bitfield[idx] >> bit_pos) & (bf_type)1) == (bf_type)1;

        const BVHBFSNodeInfo *ni = &node_info[current_chunk.node_info_offset + current_node.node_idx];
        Bounds3f b_node = FindCompressedBounds(current_chunk.b_root, ni->bk, k2b);
        if (b_node.IntersectP(ray, invDir, dirIsNeg)) {
            if (is_inner_node) {
                // Update rank to include current node
                rank += 1;
                // Determine next chunk offset and next nodes offset
                uint32_t next_chunk_offset;
                uint32_t next_nodes_offset;
                if (rank < node_pairs_per_chunk) {
                    next_chunk_offset = chunk_offset;
                    next_nodes_offset = 2 * rank;
                } else {
                    next_chunk_offset = current_chunk.child_chunk_offset + rank - node_pairs_per_chunk;
                    next_nodes_offset = 0;
                }
                // Add nodes to stack depending on which axis was used to split BVH
                if (dirIsNeg[ni->axis]) {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset};
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1};
                } else {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1};
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset};
                }
            } else {
                uint32_t sizes_idx = current_chunk.sizes_offset + current_node.node_idx - rank;
                uint32_t prim_start = current_chunk.primitive_offset + sizes[sizes_idx];
                uint32_t prim_end = current_chunk.primitive_offset + sizes[sizes_idx + 1];
                for (uint32_t i = prim_start; i < prim_end; i++)
                    if (primitives[i].get()->Intersect(ray, isect)) hit = true;
                if (node_stack_offset == 0) break;
                current_node = node_stack[--node_stack_offset];
            }
        } else {
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
        }
    }
    return hit;
}

bool BVHChunkBFSAccel::IntersectP(const Ray &ray) const {
    if (!bvh->GetNodes()) return false;
    ProfilePhase p(Prof::AccelIntersect);
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    // Follow ray through BVH nodes to find primitive intersections
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
    };
    // Variables that update whenever a new chunk is entered
    uint32_t chunk_offset = 99;  // Set to non-0 value to trigger variable updates in first iteration
    BVHChunkBFS current_chunk;
    BVHChunkBFSNode current_node;
    Vector3f k2b;
    // Initialize node stack
    BVHChunkBFSNode node_stack[64];
    uint32_t node_stack_offset = 0;
    if (dirIsNeg[bvh->GetNodes()[0].axis]) {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 0};
        current_node = BVHChunkBFSNode{0, 1};
    } else {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 1};
        current_node = BVHChunkBFSNode{0, 0};
    }
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_idx != chunk_offset) {
            chunk_offset = current_node.chunk_idx;
            current_chunk = bvh_chunks[chunk_offset];
            k2b = CalcBoundConversionFactors(current_chunk.b_root).k2b;
        }

        uint32_t rank = Rank(current_chunk.bitfield, current_node.node_idx);

        int idx = current_node.node_idx / bf_size;
        int bit_pos = current_node.node_idx % bf_size;
        bool is_inner_node = ((current_chunk.bitfield[idx] >> bit_pos) & (bf_type)1) == (bf_type)1;

        const BVHBFSNodeInfo *ni = &node_info[current_chunk.node_info_offset + current_node.node_idx];
        Bounds3f b_node = FindCompressedBounds(current_chunk.b_root, ni->bk, k2b);
        if (b_node.IntersectP(ray, invDir, dirIsNeg)) {
            if (is_inner_node) {
                // Update rank to include current node
                rank += 1;
                // Determine next chunk offset and next nodes offset
                uint32_t next_chunk_offset;
                uint32_t next_nodes_offset;
                if (rank < node_pairs_per_chunk) {
                    next_chunk_offset = chunk_offset;
                    next_nodes_offset = 2 * rank;
                } else {
                    next_chunk_offset = current_chunk.child_chunk_offset + rank - node_pairs_per_chunk;
                    next_nodes_offset = 0;
                }
                // Add nodes to stack depending on which axis was used to split BVH
                if (dirIsNeg[ni->axis]) {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset};
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1};
                } else {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1};
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset};
                }
            } else {
                uint32_t sizes_idx = current_chunk.sizes_offset + current_node.node_idx - rank;
                uint32_t prim_start = current_chunk.primitive_offset + sizes[sizes_idx];
                uint32_t prim_end = current_chunk.primitive_offset + sizes[sizes_idx + 1];
                for (uint32_t i = prim_start; i < prim_end; i++)
                    if (primitives[i].get()->IntersectP(ray)) return true;
                if (node_stack_offset == 0) break;
                current_node = node_stack[--node_stack_offset];
            }
        } else {
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
        }
    }
    return false;
}

// VISUALISATION
void BVHChunkBFSAccel::lh_dump(const char *path, bool dfs) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    if (dfs)
        lh_dump_rec_dfs(f, &vcnt, 0, WorldBound());
    else
        lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void BVHChunkBFSAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_,
                                   uint32_t chunk_offset, Bounds3f bounds) {
    BVHChunkBFS c = bvh_chunks[chunk_offset];

    std::queue<Bounds3f> bounds_q;
    bounds_q.push(bounds);
    BoundConversionFactors bcf = CalcBoundConversionFactors(bounds);

    int child_chunk_cnt = 0;

    int num_chunk_pairs =
        1;  // number of pairs of nodes in this chunk (at least 1)
    for (int pair_id = 0; pair_id < node_pairs_per_chunk; pair_id++) {
        // If the chunk doesn't contain more nodes, stop
        if (pair_id >= num_chunk_pairs) break;
        // Basic information about where node is in bitfield
        int bf_idx = (pair_id * 2) / bf_size;  // which idx in bitfield array
        int pair_idx =
            pair_id % (bf_size / 2);  // how many-th pair inside single bitfield
        for (int bit = 0; bit < 2; bit++) {
            int node_id = 2 * pair_id + bit;  // how many-th node overall
            int node_idx =
                node_id % bf_size;  // how many-th node inside single bitfield
            bf_type one = (bf_type)1;  // macro
            // Bound calculations
            BVHBFSNodeInfo ni = node_info[c.node_info_offset + node_id];
            Bounds3f node_b_comp =
                FindCompressedBounds(bounds_q.front(), ni.bk, bcf.k2b);
            // Process note according to node type
            if (((c.bitfield[bf_idx] >> node_idx) & one) == one) {
                // Inner Node
                if (num_chunk_pairs < node_pairs_per_chunk) {
                    // Real Inner Node
                    if (relative_keys) bounds_q.push(node_b_comp);
                    num_chunk_pairs++;
                } else {
                    // Chunk Ptr Inner Node
                    lh_dump_rec(f, vcnt_,
                                c.child_chunk_offset + child_chunk_cnt++,
                                node_b_comp);
                }
            }
            // Vertices ausgeben
            for (uint32_t i = 0; i < 8; i++) {
                Float x =
                    ((i & 1) == 0) ? node_b_comp.pMin.x : node_b_comp.pMax.x;
                Float y =
                    ((i & 2) == 0) ? node_b_comp.pMin.y : node_b_comp.pMax.y;
                Float z =
                    ((i & 4) == 0) ? node_b_comp.pMin.z : node_b_comp.pMax.z;
                fprintf(f, "v %f %f %f\n", x, y, z);
            }
            // Vertex indices ausgeben
            uint32_t vcnt = *vcnt_;
            fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 1, vcnt + 5,
                    vcnt + 4);  // bottom
            fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7,
                    vcnt + 6);  // top
            fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 1, vcnt + 3,
                    vcnt + 2);  // front
            fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7,
                    vcnt + 6);  // back
            fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 4, vcnt + 6,
                    vcnt + 2);  // left
            fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7,
                    vcnt + 3);  // right
            *vcnt_ += 8;
        }
        if (relative_keys) bounds_q.pop();
    }
}

void BVHChunkBFSAccel::lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_,
                                       uint32_t chunk_offset, Bounds3f bounds) {

}

}  // namespace pbrt
