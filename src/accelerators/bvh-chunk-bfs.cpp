
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

#if defined (COUNT_STATS)
#include <set>
#endif

#include "interaction.h"
#include "parallel.h"
#include "paramset.h"
#include "stats.h"

namespace pbrt {

// Stats counted when building structure
STAT_MEMORY_COUNTER("Memory/BVH-BFS tree", bvhbfs_mem); // DONE
STAT_MEMORY_COUNTER("Memory/BVH-BFS topology", bvhbfs_mem_top);

STAT_COUNTER("BVH-BFS/Chunks - # Total", bvhbfs_stat_num_chunks); // DONE
STAT_COUNTER("BVH-BFS/Chunks - # Layers (excl. root chunk)", bvhbfs_stat_num_chunkLayers); // DONE
STAT_FLOAT_DISTRIBUTION("BVH-BFS/Chunks - Fill %", bvhbfs_dist_chunkFill); // DONE

STAT_COUNTER("BVH-BFS/Nodes - # Total (incl. implicit root node)", bvhbfs_stat_num_nodes); // DONE
STAT_COUNTER("BVH-BFS/Nodes - # Leaf", bvhbfs_stat_num_leafNode); // DONE
STAT_COUNTER("BVH-BFS/Primitives - # Total", bvhbfs_stat_num_prims);

// Stats counted when intersecting
#if defined (COUNT_STATS)
STAT_COUNTER("BVH-BFS - Intersects - Primitive/Total #", bvhbfs_stat_primIntersectsTotal); // DONE
STAT_INT_DISTRIBUTION("BVH-BFS - Intersects - Primitive/Distribution", bvhbfs_dist_primIntersects); // DONE

STAT_COUNTER("BVH-BFS - Intersects - Node - Leaf/Total #", bvhbfs_stat_leafNodeIntersectsTotal); // DONE
STAT_INT_DISTRIBUTION("BVH-BFS - Intersects - Node - Leaf/Distribution", bvhbfs_dist_leafNodeIntersects); // DONE

STAT_COUNTER("BVH-BFS - Intersects - Node - Total/Total #", bvhbfs_stat_nodeIntersectsTotal); // DONE
STAT_INT_DISTRIBUTION("BVH-BFS - Intersects - Node - Total/Distribution", bvhbfs_dist_nodeIntersects); // DONE

STAT_COUNTER("BVH-BFS - Intersects - Chunk/Total #", bvhbfs_stat_chunkIntersectsTotal); // DONE
STAT_INT_DISTRIBUTION("BVH-BFS - Intersects - Chunk/Distribution", bvhbfs_dist_chunkIntersects); // DONE

// STAT_COUNTER("BVH-BFS - Intersects - Primitive/Double Q1", bvhbfs_stat_primIntersectsDoubleQ1); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Primitive/Double Median", bvhbfs_stat_primIntersectsDoubleMedian); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Primitive/Double Q3", bvhbfs_stat_primIntersectsDoubleQ3); // DONE

// STAT_COUNTER("BVH-BFS - Intersects - Node - Leaf/Double Q1", bvhbfs_stat_leafNodeIntersectsDoubleQ1); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Node - Leaf/Double Median", bvhbfs_stat_leafNodeIntersectsDoubleMedian); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Node - Leaf/Double Q3", bvhbfs_stat_leafNodeIntersectsDoubleQ3); // DONE

// STAT_COUNTER("BVH-BFS - Intersects - Node - Total/Double Q1", bvhbfs_stat_nodeIntersectsDoubleQ1); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Node - Total/Double Median", bvhbfs_stat_nodeIntersectsDoubleMedian); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Node - Total/Double Q3", bvhbfs_stat_nodeIntersectsDoubleQ3); // DONE

// STAT_COUNTER("BVH-BFS - Intersects - Chunk/Double Q1", bvhbfs_stat_chunkIntersectsDoubleQ1); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Chunk/Double Median", bvhbfs_stat_chunkIntersectsDoubleMedian); // DONE
// STAT_COUNTER("BVH-BFS - Intersects - Chunk/Double Q3", bvhbfs_stat_chunkIntersectsDoubleQ3); // DONE
#endif

struct alignas(chunk_alignment) BVHChunkBFSAccel::BVHChunkBFS {
    Bounds3f b_root;
    uint32_t primitive_offset;
    uint32_t sizes_offset;
    uint32_t node_info_offset;
    uint32_t child_chunk_offset;
    bftype bitfield[chunk_depth];
};

uint32_t BVHChunkBFSAccel::Rank(const bftype bitfield[chunk_depth], int n) const {
    uint32_t count = 0;
    bftype bits;
    for (int i = 0; i < chunk_depth; i++) {
        bits = bitfield[i];
        if (n < bfsize) break;
        count += popcnt(bits);
        n -= bfsize;
    }
    return count + popcnt(bits & ((bftone << n) - bftone));
}

Bounds3k BVHChunkBFSAccel::FindBoundsKey(Bounds3f b_root, Bounds3f b,
                                         Vector3f b2k, Vector3f k2b) const {
    Bounds3k b_k = Bounds3k{};
    for (int i = 0; i < 3; i++) {
        // b2k = 255 / dist
        Float min = b2k[i] * (b.pMin[i] - b_root.pMin[i]);
        int floor_min = floor(min);
        floor_min = (floor_min < 0) ? 0 : floor_min > 255 ? 255 : floor_min;

        Float max = b2k[i] * (b.pMax[i] - b_root.pMin[i]);
        int ceil_max = ceil(max);
        ceil_max = (ceil_max < 0) ? 0 : ceil_max > 255 ? 255 : ceil_max;

        b_k.min[i] = (uint8_t)floor_min;
        b_k.max[i] = (uint8_t)ceil_max;
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
        b2k[i] = 255 / ((b.pMax[i] != b.pMin[i]) ? (b.pMax[i] - b.pMin[i]) : FLT_MIN);
        k2b[i] = 1 / b2k[i];
    }
    return BoundConversionFactors{b2k, k2b};
}

Vector3f CalcK2B(Bounds3f b) {
    Vector3f k2b;
    for (int i = 0; i < 3; i++) {
        k2b[i] = (b.pMax[i] - b.pMin[i]) / 255.0;
    }
    return k2b;
}

ChildrenBuildInfo BVHChunkBFSAccel::GetChildrenBuildInfo(
    Bounds3f b_root, uint32_t node_offset) const {
    ChildrenBuildInfo cbi = ChildrenBuildInfo{};
    cbi.idx[0] = node_offset + 1;
    cbi.idx[1] = bvh->GetNodes()[node_offset].secondChildOffset;
    BoundConversionFactors bcf = CalcBoundConversionFactors(b_root);
    for (int child = 0; child < 2; child++) {
        Bounds3f b_child = bvh->GetNodes()[cbi.idx[child]].bounds;
        cbi.b_key[child] = FindBoundsKey(b_root, b_child, bcf.b2k, bcf.k2b);
        cbi.b_comp[child] =
            FindCompressedBounds(b_root, cbi.b_key[child], bcf.k2b);
    }
    return cbi;
}

#if defined (COUNT_STATS)
// enum class bvhbfs_stat_name { PRIMS, NODES, LEAVES, CHUNKS }
// void Calc_Q1_Median_Q3(std::vector<int> data, std::string bvhbfs_stat_name) {
//     int bvhbfs_stat_q1, bvhbfs_stat_median, bvhbfs_stat_q3;
//     int num_rays = data.size();
//     std::sort(data.begin(), data.end());
//     // Median
//     int j = num_rays / 2;
//     int i = ((num_rays % 2) == 0) ? j-1 : j;
//     int double_median = data[i] + data[j];
//     // Q1
//     int max_q1 = (i == j) ? i-1 : i;
//     int j_q1 = max_q1 / 2;
//     int i_q1 = ((max_q1 % 2) == 0) ? j-1 : j;
//     int double_q1 = data[i_q1] + data[j_q1];
//     // Q3
//     int min_q3 = (i == j) ? j+1 : j;
//     int j_q3 = min_q3 + j_q1;
//     int i_q3 = min_q3 + i_q1;
//     int double_q3 = data[i_q3] + data[j_q3];

//     if (bvhbfs_stat_name == PRIMS) {
//         bvhbfs_stat_primIntersectsDoubleMedian = double_median;
//         bvhbfs_stat_primIntersectsDoubleQ1 = double_q1;
//         bvhbfs_stat_primIntersectsDoubleQ3 = double_q3;
//     } else if (bvhbfs_stat_name == LEAVES) {
//         bvhbfs_stat_leafNodeIntersectsDoubleMedian = double_median;
//         bvhbfs_stat_leafNodeIntersectsDoubleQ1 = double_q1;
//         bvhbfs_stat_leafNodeIntersectsDoubleQ3 = double_q3;
//     } else if (bvhbfs_stat_name == NODES) {
//         bvhbfs_stat_nodeIntersectsDoubleMedian = double_median;
//         bvhbfs_stat_nodeIntersectsDoubleQ1 = double_q1;
//         bvhbfs_stat_nodeIntersectsDoubleQ3 = double_q3;
//     } else if (bvhbfs_stat_name == CHUNKS) {
//         bvhbfs_stat_chunkIntersectsDoubleMedian = double_median;
//         bvhbfs_stat_chunkIntersectsDoubleQ1 = double_q1;
//         bvhbfs_stat_chunkIntersectsDoubleQ3 = double_q3;
//     }
// }
#endif

// BVHAccel Method Definitions
BVHChunkBFSAccel::BVHChunkBFSAccel(std::vector<std::shared_ptr<Primitive>> p,
                                   int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      splitMethod(splitMethod) {
    printf("Chosen Accelerator: BVH w/ BFS Chunks\n");

    #if defined (COUNT_STATS)
    // prim_isects_per_ray = new std::vector<int>();
    #endif

    ProfilePhase _(Prof::AccelConstruction);
    // Let BVHAccel build BVH Structure
    bvh = new BVHAccel(p, maxPrimsInNode, splitMethod);
    if (!bvh->GetNodes()) return;
    // Chunk struct for building the chunks
    bvh_chunks.push_back(BVHChunkBFS{});
    // We count the initial root node as a chunk pointer inner node
    bvhbfs_stat_num_nodes++;
    Recurse(0, 0, WorldBound(), 0);
    // Visualisation:
    // lh_dump("bvh_vis_bfs.obj");
    // lh_dump("bvh_vis_dfs.obj", true);

    bvhbfs_mem_top +=
            sizeof(bvh_chunks) + bvh_chunks.size() * sizeof(bvh_chunks[0]) + 
            sizeof(node_info) + node_info.size() * sizeof(node_info[0]); 
    
    bvhbfs_mem += sizeof(*this) - sizeof(bvh) + //Removing bvh since its not relevant
            bvh_chunks.size() * sizeof(bvh_chunks[0]) +
            node_info.size() * sizeof(node_info[0]) +
            primitives.size() * sizeof(primitives[0]) +
            sizes.size() * sizeof(sizes[0]);

}

void BVHChunkBFSAccel::Recurse(uint32_t chunk_offset, uint32_t root_node_idx,
                               Bounds3f b_root, int chunk_layer) {
    if (chunk_layer > bvhbfs_stat_num_chunkLayers) bvhbfs_stat_num_chunkLayers = chunk_layer;
    struct Chunk {
        uint32_t root_node_idx;
        Bounds3f b_root;
    };
    bvhbfs_stat_num_chunks++;
    std::vector<Chunk> chunk_ptr_nodes;
    // === Building Chunk ===
    // Set meta information of chunk
    bvh_chunks[chunk_offset].b_root = b_root;
    bvh_chunks[chunk_offset].node_info_offset = node_info.size();
    bvh_chunks[chunk_offset].sizes_offset = sizes.size();
    bvh_chunks[chunk_offset].primitive_offset = primitives.size();
    bvh_chunks[chunk_offset].child_chunk_offset = bvh_chunks.size();
    for (int i = 0; i < chunk_depth; i++)
        bvh_chunks[chunk_offset].bitfield[i] = bftzero;
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
        int idx = nodes_processed / bfsize;
        int bit_pos = nodes_processed % bfsize;
        // Determine node type, and process accordingly
        if (is_inner_node) {
            bvh_chunks[chunk_offset].bitfield[idx] |= (bftone << bit_pos);
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
                sizes.push_back(node->nPrimitives);
                firstLeafInChunk = false;
            } else sizes.push_back(sizes.back() + node->nPrimitives);
            for (int i = 0; i < node->nPrimitives; i++)
                primitives.push_back(bvh->GetPrimitives()[node->primitivesOffset + i]);
            bvhbfs_stat_num_prims += node->nPrimitives;
            bvhbfs_stat_num_leafNode++;
        }
        node_info.push_back(BVHBFSNodeInfo{build_node.b_key, node->axis});
        // Update node counter and -queue
        nodes_processed++;
        node_q.pop();
        bvhbfs_stat_num_nodes++;
    }  // nodes_q

    ReportValue(bvhbfs_dist_chunkFill, 100 * (float)nodes_processed / (float)(2 * node_pairs_per_chunk));

    for (int i = 0; i < chunk_ptr_nodes.size(); i++) {
        Recurse(bvh_chunks[chunk_offset].child_chunk_offset + i,
                chunk_ptr_nodes[i].root_node_idx, chunk_ptr_nodes[i].b_root,
                chunk_layer + 1);
    }
}

Bounds3f BVHChunkBFSAccel::WorldBound() const { return bvh->WorldBound(); }

BVHChunkBFSAccel::~BVHChunkBFSAccel() {
    #if defined (COUNT_STATS)
    // Calc_Q1_Median_Q3(*prim_isects_per_ray, PRIMITIVES);
    // Calc_Q1_Median_Q3(*node_isects_per_ray, ALLNODES);
    // Calc_Q1_Median_Q3(*leafNode_isects_per_ray, LEAFNODES);
    // Calc_Q1_Median_Q3(*chunk_isects_per_ray, CHUNKS);
    FreeAligned(prim_isects_per_ray);
    #endif
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
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    if (!WorldBound().IntersectP(ray, invDir, dirIsNeg)) return false;
    // === COUNT STATS INITIALIZE BEGIN ===
    #if defined (COUNT_STATS)
    int num_prims_intersected = 0;
    int num_nodes_intersected = 1;
    int num_leafNodes_intersected = 0;
    std::set<int> chunks_intersected;
    #endif
    // === COUNT STATS INITIALIZE END ===
    bool hit = false;
    // Follow ray through BVH nodes to find primitive intersections
    #if defined(REL_KEYS)
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
        Bounds3f root_bounds;
    };
    #else
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
    };
    #endif
    // Variables that update whenever a new chunk is entered
    uint32_t chunk_offset = 99;  // Set to non-0 value to trigger variable updates in first iteration
    const BVHChunkBFS *current_chunk;
    // More variables
    BVHChunkBFSNode current_node;
    Vector3f k2b;
    Bounds3f root_b;
    // Initialize node stack
    BVHChunkBFSNode node_stack[64];
    uint32_t node_stack_offset = 0;
    if (dirIsNeg[bvh->GetNodes()[0].axis]) {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 0
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
        current_node = BVHChunkBFSNode{0, 1
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
    } else {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 1
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
        current_node = BVHChunkBFSNode{0, 0
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
    }
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_idx != chunk_offset) {
            chunk_offset = current_node.chunk_idx;
            current_chunk = &(bvh_chunks[chunk_offset]);
            #if !defined(REL_KEYS)
            k2b = CalcK2B(current_chunk->b_root);
            root_b = current_chunk->b_root;
            #endif
            // === COUNT STATS BEGIN ===
            #if defined (COUNT_STATS)
            chunks_intersected.emplace(chunk_offset);
            #endif
            // === COUNT STATS END ===
        }

        uint32_t rank = Rank(current_chunk->bitfield, current_node.node_idx);

        bool is_inner_node = ((current_chunk->bitfield[current_node.node_idx / bfsize]
                >> (current_node.node_idx % bfsize)) & bftone) == bftone;

        const BVHBFSNodeInfo *ni = &node_info[current_chunk->node_info_offset + current_node.node_idx];
        #if defined(REL_KEYS)
        k2b = CalcK2B(current_node.root_bounds);
        root_b = current_node.root_bounds;
        #endif
        Bounds3f b_node = FindCompressedBounds(root_b, ni->bk, k2b);
        if (b_node.IntersectP(ray, invDir, dirIsNeg)) {
            // === COUNT STATS FOR ALL NODES BEGIN ===
            #if defined (COUNT_STATS)
            num_nodes_intersected++;
            #endif
            // === COUNT STATS FOR ALL NODES END ===
            if (is_inner_node) {
                // Update rank to include current node
                rank += 1;
                // Determine next chunk offset and next nodes offset
                uint32_t next_chunk_offset;
                uint32_t next_nodes_offset;
                #if defined(REL_KEYS)
                Bounds3f next_bounds;
                #endif
                if (rank < node_pairs_per_chunk) {
                    // Real Inner Node
                    next_chunk_offset = chunk_offset;
                    next_nodes_offset = 2 * rank;
                    #if defined(REL_KEYS)
                    next_bounds = b_node;
                    #endif
                } else {
                    // Chunk Ptr Inner Node
                    next_chunk_offset = current_chunk->child_chunk_offset + rank - node_pairs_per_chunk;
                    next_nodes_offset = 0;
                    #if defined(REL_KEYS)
                    next_bounds = bvh_chunks[next_chunk_offset].b_root;
                    #endif
                }
                // Add nodes to stack depending on which axis was used to split BVH
                if (dirIsNeg[ni->axis]) {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                } else {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                }
            } else {
                // Leaf Node
                uint32_t sizes_idx = current_chunk->sizes_offset + current_node.node_idx - rank;
                uint32_t prim_start = current_chunk->primitive_offset +
                        ((current_node.node_idx == rank) ? 0 : sizes[sizes_idx - 1]);
                uint32_t prim_end = current_chunk->primitive_offset + sizes[sizes_idx];
                for (uint32_t i = prim_start; i < prim_end; i++)
                    if (primitives[i].get()->Intersect(ray, isect)) hit = true;
                // === COUNT STATS FOR LEAF NODES BEGIN ===
                #if defined (COUNT_STATS)
                num_prims_intersected += prim_end - prim_start;
                num_leafNodes_intersected++;
                #endif
                // === COUNT STATS FOR LEAF NODES END ===
                if (node_stack_offset == 0) break;
                current_node = node_stack[--node_stack_offset];
            }
        } else {
            // Didn't intersect the nodes bounds
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
        }
    }
    // === COUNT STATS FOR RAY AND TOTAL STATS BEGIN ===
    #if defined (COUNT_STATS)
    bvhbfs_stat_primIntersectsTotal += num_prims_intersected;
    ReportValue(bvhbfs_dist_primIntersects, num_prims_intersected);
    // prim_isects_per_ray->push_back(num_prims_intersected);

    bvhbfs_stat_nodeIntersectsTotal += num_nodes_intersected;
    ReportValue(bvhbfs_dist_nodeIntersects, num_nodes_intersected);
    // node_isects_per_ray->push_back(num_nodes_intersected);

    bvhbfs_stat_leafNodeIntersectsTotal += num_leafNodes_intersected;
    ReportValue(bvhbfs_dist_leafNodeIntersects, num_leafNodes_intersected);
    // leafNode_isects_per_ray->push_back(num_leafNodes_intersected);

    bvhbfs_stat_chunkIntersectsTotal += chunks_intersected.size();
    ReportValue(bvhbfs_dist_chunkIntersects, chunks_intersected.size());
    // chunk_isects_per_ray->push_back(chunks_intersected.size());
    #endif
    // === COUNT STATS FOR RAY AND TOTAL STATS END ===
    return hit;
}

float BVHChunkBFSAccel::IntersectMetric(const Ray &ray, metric m) const {
    if (!bvh->GetNodes()) return 0;
    float metric_cnt = 0;
    SurfaceInteraction _isect;
    SurfaceInteraction *isect = &_isect;
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    if (!WorldBound().IntersectP(ray, invDir, dirIsNeg)) return false;
    bool hit = false;
    // Follow ray through BVH nodes to find primitive intersections
    #if defined(REL_KEYS)
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
        Bounds3f root_bounds;
    };
    #else
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
    };
    #endif
    // Variables that update whenever a new chunk is entered
    uint32_t chunk_offset = 99;  // Set to non-0 value to trigger variable updates in first iteration
    const BVHChunkBFS *current_chunk;
    // More variables
    BVHChunkBFSNode current_node;
    Vector3f k2b;
    Bounds3f root_b;
    // Initialize node stack
    BVHChunkBFSNode node_stack[64];
    uint32_t node_stack_offset = 0;
    if (dirIsNeg[bvh->GetNodes()[0].axis]) {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 0
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
        current_node = BVHChunkBFSNode{0, 1
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
    } else {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 1
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
        current_node = BVHChunkBFSNode{0, 0
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
    }
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_idx != chunk_offset) {
            chunk_offset = current_node.chunk_idx;
            current_chunk = &(bvh_chunks[chunk_offset]);
            #if !defined(REL_KEYS)
            k2b = CalcK2B(current_chunk->b_root);
            root_b = current_chunk->b_root;
            #endif
        }

        uint32_t rank = Rank(current_chunk->bitfield, current_node.node_idx);

        bool is_inner_node = ((current_chunk->bitfield[current_node.node_idx / bfsize]
                >> (current_node.node_idx % bfsize)) & bftone) == bftone;

        const BVHBFSNodeInfo *ni = &node_info[current_chunk->node_info_offset + current_node.node_idx];
        #if defined(REL_KEYS)
        k2b = CalcK2B(current_node.root_bounds);
        root_b = current_node.root_bounds;
        #endif
        Bounds3f b_node = FindCompressedBounds(root_b, ni->bk, k2b);
        if (b_node.IntersectP(ray, invDir, dirIsNeg)) {
            if (m == metric::NODES) metric_cnt++;
            if (is_inner_node) {
                // Update rank to include current node
                rank += 1;
                // Determine next chunk offset and next nodes offset
                uint32_t next_chunk_offset;
                uint32_t next_nodes_offset;
                #if defined(REL_KEYS)
                Bounds3f next_bounds;
                #endif
                if (rank < node_pairs_per_chunk) {
                    // Real Inner Node
                    next_chunk_offset = chunk_offset;
                    next_nodes_offset = 2 * rank;
                    #if defined(REL_KEYS)
                    next_bounds = b_node;
                    #endif
                } else {
                    // Chunk Ptr Inner Node
                    next_chunk_offset = current_chunk->child_chunk_offset + rank - node_pairs_per_chunk;
                    next_nodes_offset = 0;
                    #if defined(REL_KEYS)
                    next_bounds = bvh_chunks[next_chunk_offset].b_root;
                    #endif
                }
                // Add nodes to stack depending on which axis was used to split BVH
                if (dirIsNeg[ni->axis]) {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                } else {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                }
            } else {
                // Leaf Node
                if (m == metric::LEAFNODES) metric_cnt++;
                uint32_t sizes_idx = current_chunk->sizes_offset + current_node.node_idx - rank;
                uint32_t prim_start = current_chunk->primitive_offset +
                        ((current_node.node_idx == rank) ? 0 : sizes[sizes_idx - 1]);
                uint32_t prim_end = current_chunk->primitive_offset + sizes[sizes_idx];
                for (uint32_t i = prim_start; i < prim_end; i++)
                    if (primitives[i].get()->Intersect(ray, isect)) hit = true;
                if (m == metric::PRIMITIVES) metric_cnt += prim_end = prim_start;
                if (node_stack_offset == 0) break;
                current_node = node_stack[--node_stack_offset];
            }
        } else {
            // Didn't intersect the nodes bounds
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
        }
    }
    return metric_cnt;
}

bool BVHChunkBFSAccel::IntersectP(const Ray &ray) const {
    if (!bvh->GetNodes()) return false;
    ProfilePhase p(Prof::AccelIntersectP);
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    if (!WorldBound().IntersectP(ray, invDir, dirIsNeg)) return false;
    // Follow ray through BVH nodes to find primitive intersections
    #if defined(REL_KEYS)
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
        Bounds3f root_bounds;
    };
    #else
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
    };
    #endif
    // Variables that update whenever a new chunk is entered
    uint32_t chunk_offset = 99;  // Set to non-0 value to trigger variable updates in first iteration
    const BVHChunkBFS *current_chunk;
    // More variables
    BVHChunkBFSNode current_node;
    Vector3f k2b;
    Bounds3f root_b;
    // Initialize node stack
    BVHChunkBFSNode node_stack[64];
    uint32_t node_stack_offset = 0;
    if (dirIsNeg[bvh->GetNodes()[0].axis]) {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 0
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
        current_node = BVHChunkBFSNode{0, 1
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
    } else {
        node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 1
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
        current_node = BVHChunkBFSNode{0, 0
            #if defined(REL_KEYS)
            , WorldBound()
            #endif
            };
    }
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_idx != chunk_offset) {
            chunk_offset = current_node.chunk_idx;
            current_chunk = &(bvh_chunks[chunk_offset]);
            #if !defined(REL_KEYS)
            k2b = CalcK2B(current_chunk->b_root);
            root_b = current_chunk->b_root;
            #endif
        }

        uint32_t rank = Rank(current_chunk->bitfield, current_node.node_idx);

        bool is_inner_node = ((current_chunk->bitfield[current_node.node_idx / bfsize]
                >> (current_node.node_idx % bfsize)) & bftone) == bftone;

        const BVHBFSNodeInfo *ni = &node_info[current_chunk->node_info_offset + current_node.node_idx];
        #if defined(REL_KEYS)
        k2b = CalcK2B(current_node.root_bounds);
        root_b = current_node.root_bounds;
        #endif
        Bounds3f b_node = FindCompressedBounds(root_b, ni->bk, k2b);
        if (b_node.IntersectP(ray, invDir, dirIsNeg)) {
            if (is_inner_node) {
                // Update rank to include current node
                rank += 1;
                // Determine next chunk offset and next nodes offset
                uint32_t next_chunk_offset;
                uint32_t next_nodes_offset;
                #if defined(REL_KEYS)
                Bounds3f next_bounds;
                #endif
                if (rank < node_pairs_per_chunk) {
                    // Real Inner Node
                    next_chunk_offset = chunk_offset;
                    next_nodes_offset = 2 * rank;
                    #if defined(REL_KEYS)
                    next_bounds = b_node;
                    #endif
                } else {
                    // Chunk Ptr Inner Node
                    next_chunk_offset = current_chunk->child_chunk_offset + rank - node_pairs_per_chunk;
                    next_nodes_offset = 0;
                    #if defined(REL_KEYS)
                    next_bounds = bvh_chunks[next_chunk_offset].b_root;
                    #endif
                }
                // Add nodes to stack depending on which axis was used to split BVH
                if (dirIsNeg[ni->axis]) {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                } else {
                    node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                    current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset
                            #if defined(REL_KEYS)
                            , next_bounds
                            #endif
                            };
                }
            } else {
                // Leaf Node
                uint32_t sizes_idx = current_chunk->sizes_offset + current_node.node_idx - rank;
                uint32_t prim_start = current_chunk->primitive_offset +
                        ((current_node.node_idx == rank) ? 0 : sizes[sizes_idx - 1]);
                uint32_t prim_end = current_chunk->primitive_offset + sizes[sizes_idx];
                for (uint32_t i = prim_start; i < prim_end; i++)
                    if (primitives[i].get()->IntersectP(ray)) return true;
                if (node_stack_offset == 0) break;
                current_node = node_stack[--node_stack_offset];
            }
        } else {
            // Didn't intersect the nodes bounds
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

    if (chunk_offset == 0) {
        // Ver+= 8;tices ausgeben
        for (uint32_t i = 0; i < 8; i++) {
            Float x = ((i & 1) == 0) ? WorldBound().pMin.x : WorldBound().pMax.x;
            Float y = ((i & 2) == 0) ? WorldBound().pMin.y : WorldBound().pMax.y;
            Float z = ((i & 4) == 0) ? WorldBound().pMin.z : WorldBound().pMax.z;
            fprintf(f, "v %f %f %f\n", x, y, z);
        }
        // Vertex indices ausgeben
        uint32_t vcnt = *vcnt_;
        fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 1, vcnt + 5, vcnt + 4);  // bottom
        fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);  // top
        fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 1, vcnt + 3, vcnt + 2);  // front
        fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);  // back
        fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 4, vcnt + 6, vcnt + 2);  // left
        fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);  // right
        *vcnt_ += 8;
    }

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
        int bf_idx = (pair_id * 2) / bfsize;  // which idx in bitfield array
        int pair_idx =
            pair_id % (bfsize / 2);  // how many-th pair inside single bitfield
        for (int bit = 0; bit < 2; bit++) {
            int node_id = 2 * pair_id + bit;  // how many-th node overall
            int node_idx =
                node_id % bfsize;  // how many-th node inside single bitfield
            const bftype one = bftone;  // macro
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
                    node_b_comp = bvh_chunks[c.child_chunk_offset + child_chunk_cnt].b_root;
                }
            } else {
                uint32_t rank = Rank(c.bitfield, node_idx);
                uint32_t sizes_idx = c.sizes_offset + node_idx - rank;
                uint32_t prim_start = (node_idx == rank) ? 0 : sizes[sizes_idx - 1];
                uint32_t prim_end = sizes[sizes_idx];
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
    // Follow ray through BVH nodes to find primitive intersections
    struct BVHChunkBFSNode {
        uint32_t chunk_idx;
        uint32_t node_idx;
    };
    // Ver+= 8;tices ausgeben
    for (uint32_t i = 0; i < 8; i++) {
        Float x = ((i & 1) == 0) ? WorldBound().pMin.x : WorldBound().pMax.x;
        Float y = ((i & 2) == 0) ? WorldBound().pMin.y : WorldBound().pMax.y;
        Float z = ((i & 4) == 0) ? WorldBound().pMin.z : WorldBound().pMax.z;
        fprintf(f, "v %f %f %f\n", x, y, z);
    }
    // Vertex indices ausgeben
    uint32_t vcnt = *vcnt_;
    fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 1, vcnt + 5, vcnt + 4);  // bottom
    fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);  // top
    fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 1, vcnt + 3, vcnt + 2);  // front
    fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);  // back
    fprintf(f, "f %d %d %d %d\n", vcnt, vcnt + 4, vcnt + 6, vcnt + 2);  // left
    fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);  // right
    *vcnt_ += 8;
    // Variables that update whenever a new chunk is entered
    uint32_t current_chunk_idx = 99;  // Set to non-0 value to trigger variable updates in first iteration
    BVHChunkBFS current_chunk;
    BVHChunkBFSNode current_node;
    Vector3f k2b;
    // Initialize node stack
    BVHChunkBFSNode node_stack[64];
    uint32_t node_stack_offset = 0;
    node_stack[node_stack_offset++] = BVHChunkBFSNode{0, 0};
    current_node = BVHChunkBFSNode{0, 1};
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_idx != current_chunk_idx) {
            current_chunk_idx = current_node.chunk_idx;
            current_chunk = bvh_chunks[current_chunk_idx];
            k2b = CalcBoundConversionFactors(current_chunk.b_root).k2b;
        }

        uint32_t rank = Rank(current_chunk.bitfield, current_node.node_idx);

        int idx = current_node.node_idx / bfsize;
        int bit_pos = current_node.node_idx % bfsize;
        bool is_inner_node = ((current_chunk.bitfield[idx] >> bit_pos) & bftone) == bftone;

        const BVHBFSNodeInfo *ni = &node_info[current_chunk.node_info_offset + current_node.node_idx];
        Bounds3f b_node = FindCompressedBounds(current_chunk.b_root, ni->bk, k2b);
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
                b_node = bvh_chunks[next_chunk_offset].b_root;
            }
            // Add nodes to stack depending on which axis was used to split BVH
            node_stack[node_stack_offset++] = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset};
            current_node = BVHChunkBFSNode{next_chunk_offset, next_nodes_offset + 1};
        }
        // Vertices ausgeben
        for (uint32_t i = 0; i < 8; i++) {
            Float x =
                ((i & 1) == 0) ? b_node.pMin.x : b_node.pMax.x;
            Float y =
                ((i & 2) == 0) ? b_node.pMin.y : b_node.pMax.y;
            Float z =
                ((i & 4) == 0) ? b_node.pMin.z : b_node.pMax.z;
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
        if (!is_inner_node) {
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
        }
    }
}

}  // namespace pbrt
