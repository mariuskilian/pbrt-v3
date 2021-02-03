
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


// accelerators/kdtreeaccel.cpp*
#include "accelerators/octree-chunk-bfs.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>
#include <array>
#include <queue>

namespace pbrt {

// === HELPERS ===

int BitfieldRankOffset(std::array<bftype, BFS_CHUNK_DEPTH> bitfield, int n) {
    int count = 0;
    bftype bits;
    for (int i = 0; i < BFS_CHUNK_DEPTH; i++) {
        if ((n -= bfsize) < 0) break;
        count += Rank(bitfield[i]);
    }
    return count;
}

int BitfieldRankUnique(std::array<bftype, BFS_CHUNK_DEPTH> bitfield, int n) {
    return Rank(bitfield[n / bfsize], n % bfsize);
}

int RankAll(std::array<bftype, BFS_CHUNK_DEPTH> bitfield, int n) {
    uint32_t count = 0;
    bftype bits;
    for (int i = 0; i < BFS_CHUNK_DEPTH; i++) {
        bits = bitfield[i];
        if (n < bfsize) break;
        count += popcnt(bits);
        n -= bfsize;
    }
    return count + popcnt(bits & ((bftone << n) - bftone));
}

bool IsInnerNode(std::array<bftype, BFS_CHUNK_DEPTH> bitfield, int n) {
    return ((bitfield[n / bfsize] >> (n % bfsize)) & 1) == 1;
}

// === OCTREE STRUCT CREATION ==
OcChunkBFSAccel::OcChunkBFSAccel(std::vector<std::shared_ptr<Primitive>> p,
        int max_prims, float prm_thresh, float vol_thresh) 
        : primitives(std::move(p)) {
    printf("Chosen Accelerator: Octree w/ BFS Chunks\n");

    oba = OctreeBasicAccel(primitives, max_prims, prm_thresh, vol_thresh);
    wb = oba.WorldBound();
    
    if (oba.Nodes().size() > 1) {
        octree.push_back(Chunk{});
        sizes.push_back(0);
        Recurse(0, 0);
    }
    //lh_dump("visualize_bfs.obj");
    //lh_dump_dfs("visualize_dfs.obj");
}

void OcChunkBFSAccel::Recurse(uint32_t root_node_offset, int chunk_idx) {
    octree[chunk_idx].child_chunk_offset = octree.size();
    octree[chunk_idx].sizes_offset = sizes.size();

    uint32_t root_child_offset = oba.Nodes()[root_node_offset] >> 1;
    std::queue<uint32_t> bfs_nodes_q;
    for (int i = 0; i < 8; i++) bfs_nodes_q.push(root_child_offset + i);

    int chunk_fill_size = 1; // number of [set of 8 nodes] processed and/or in queue
    int num_nodes = 0; // number of individual nodes already processed
    std::vector<uint32_t> chunk_ptr_nodes; // nodes that point to other chunks

    uint32_t node_offset;
    bool is_inner_node;
    while (!bfs_nodes_q.empty()) {
        node_offset = bfs_nodes_q.front();
        is_inner_node = (oba.Nodes()[node_offset] & 1) == 0;
            
        int idx = num_nodes / bfsize;
        int bit_pos = num_nodes % bfsize;
        // When starting a new index, make sure the initial value of the bitcode is 0
        if (bit_pos == 0) octree[chunk_idx].nodes[idx] = bftzero;

        // If it's an inner node, properly reserve its children (new chunk or in current chunk)
        if (is_inner_node) {
            // Inner node
            if (chunk_fill_size == BFS_NUM_SETS_PER_CHUNK) {
                // Chunk is full, need to make new chunk
                chunk_ptr_nodes.push_back(node_offset);
                octree.push_back(Chunk{}); // reserve chunk slot
            } else {
                // Chunk has space for children
                uint32_t child_offset = oba.Nodes()[node_offset] >> 1;
                for (int i = 0; i < 8; i++) { bfs_nodes_q.push(child_offset + i); }
                chunk_fill_size++;
            }
            octree[chunk_idx].nodes[idx] |= bftone << bit_pos;
        } else {
            // Leaf Node
            uint32_t prim_start = oba.Nodes()[node_offset] >> 1;
            uint32_t prim_end = prim_start + oba.Sizes()[node_offset];
            sizes.push_back(sizes.back() + oba.Sizes()[node_offset]);
            leaves.insert(leaves.end(), oba.Leaves().begin() + prim_start, oba.Leaves().begin() + prim_end);
        }

        num_nodes++;
        bfs_nodes_q.pop();
    } 
    
    // Create additional chunks as needed
    for (int i = 0; i < chunk_ptr_nodes.size(); i++)
        Recurse(chunk_ptr_nodes[i], octree[chunk_idx].child_chunk_offset + i);

}

std::shared_ptr<OcChunkBFSAccel> CreateOcChunkBFSAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    int max_prims = ps.FindOneInt("maxprims", 32);
    float prm_thresh = ps.FindOneFloat("prmthresh", 0.9);
    float vol_thresh = ps.FindOneFloat("volthresh", 0.9);

    return std::make_shared<OcChunkBFSAccel>(std::move(prims), max_prims, prm_thresh, vol_thresh);
}

OcChunkBFSAccel::~OcChunkBFSAccel() { //FreeAligned(nodes2);
}

// === OCTREE RAY TRAVERSAL ===

// TODO Rekursion in Schleife umwandeln (schneller)
bool OcChunkBFSAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return false;
    bool hit = false;
    Vector3f invDir = {1/ray.d[0], 1/ray.d[1], 1/ray.d[2]};
    struct OctreeBFSTraversalNode { uint32_t chunk_offset; uint32_t node_idx; uint8_t bounds_stack_offset; Float tMin; };
    struct OctreeBFSTraversalBounds { Bounds3f bounds; Vector3f b_h; };
    // Variables that update whenever a new chunk is entered
    uint chunk_offset = 99; // Set to non-0 value to trigger variable updates in first iteration
    const Chunk *current_chunk;
    // More variables
    OctreeBFSTraversalNode current_node;
    // Initialize node stack and bounds stack
    OctreeBFSTraversalNode node_stack[128];
    OctreeBFSTraversalBounds bounds_stack[64];
    uint32_t node_stack_offset = 0;
    uint32_t bounds_stack_offset = 0;
    Bounds3f current_bounds = wb;
    // Add root nodes children to stack
    ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, tMin, invDir);
    for (int i = 3; i > 0; i--)
        if (i < traversal.size && traversal.nodes[i].tMin <= ray.tMax)
            node_stack[node_stack_offset++] = OctreeBFSTraversalNode{0, traversal.nodes[i].idx, 0, traversal.nodes[i].tMin};
    current_node = OctreeBFSTraversalNode{0, traversal.nodes[0].idx, 0, traversal.nodes[0].tMin};
    // Add root bounds to bounds stack
    bounds_stack[bounds_stack_offset] = OctreeBFSTraversalBounds{current_bounds, BoundsHalf(current_bounds)};
    uint16_t current_rank_offset = 0;
    uint16_t current_bitfield_idx = 0;
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_offset != chunk_offset) {
            chunk_offset = current_node.chunk_offset;
            current_chunk = &(octree[chunk_offset]);
        }
        
        current_bounds = DivideBounds(bounds_stack[bounds_stack_offset].bounds,
                                      current_node.node_idx & 7,
                                      bounds_stack[bounds_stack_offset].b_h);

        uint16_t bitfield_idx = current_node.node_idx / bfsize;
        bool is_inner_node = ((current_chunk->nodes[bitfield_idx]
                >> (current_node.node_idx % bfsize)) & bftone) == bftone;
        
        uint32_t rank = RankAll(current_chunk->nodes, current_node.node_idx);

        if (is_inner_node) {
            // Inner Node...
            rank++;
            ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, current_node.tMin, invDir);
            // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
            Vector3f b_h = BoundsHalf(current_bounds);
            bounds_stack[++bounds_stack_offset] = OctreeBFSTraversalBounds{current_bounds, b_h};
            uint32_t child_idx_offset;
            uint32_t child_chunk_offset;
            if (rank < BFS_NUM_SETS_PER_CHUNK) {
                // ... with children in same chunk
                child_idx_offset = 8 * rank;
                child_chunk_offset = current_node.chunk_offset;
            } else {
                // ... with children in different chunk
                child_idx_offset = 0;
                child_chunk_offset = current_chunk->child_chunk_offset + rank - BFS_NUM_SETS_PER_CHUNK;
            }
            // Add relevant child nodes to node stack
            for (int i = 3; i > 0; i--)
                if (i < traversal.size && traversal.nodes[i].tMin <= ray.tMax)
                    node_stack[node_stack_offset++] = OctreeBFSTraversalNode{
                            child_chunk_offset,
                            child_idx_offset + traversal.nodes[i].idx,
                            (uint8_t)bounds_stack_offset,
                            traversal.nodes[i].tMin};
            current_node = OctreeBFSTraversalNode{
                    child_chunk_offset,
                    child_idx_offset + traversal.nodes[0].idx,
                    (uint8_t)bounds_stack_offset,
                    traversal.nodes[0].tMin};
        } else {
            // Leaf Node
            int leaf_idx = current_chunk->sizes_offset + current_node.node_idx - rank;
            uint32_t prim_start = sizes[leaf_idx - 1];
            uint32_t prim_end = sizes[leaf_idx];
            for (uint32_t i = prim_start; i < prim_end; i++)
                // TODO LRU-Cache/Mailboxing, damit man nicht mehrmals dasselbe primitiv testen muss
                if (leaves[i].get()->Intersect(ray, isect)) hit = true;
            if (hit && BoundsContainPoint(current_bounds, isect->p)) return true;
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
            bounds_stack_offset = current_node.bounds_stack_offset;
        }
    }
    return hit;
}

bool OcChunkBFSAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return false;
    Vector3f invDir = {1/ray.d[0], 1/ray.d[1], 1/ray.d[2]};
    struct OctreeBFSTraversalNode { uint32_t chunk_offset; uint32_t node_idx; uint8_t bounds_stack_offset; Float tMin; };
    struct OctreeBFSTraversalBounds { Bounds3f bounds; Vector3f b_h; };
    // Variables that update whenever a new chunk is entered
    uint chunk_offset = 99; // Set to non-0 value to trigger variable updates in first iteration
    const Chunk *current_chunk;
    // More variables
    OctreeBFSTraversalNode current_node;
    // Initialize node stack and bounds stack
    OctreeBFSTraversalNode node_stack[128];
    OctreeBFSTraversalBounds bounds_stack[64];
    uint32_t node_stack_offset = 0;
    uint32_t bounds_stack_offset = 0;
    Bounds3f current_bounds = wb;
    // Add root nodes children to stack
    ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, tMin, invDir);
    for (int i = 3; i > 0; i--)
        if (i < traversal.size && traversal.nodes[i].tMin <= ray.tMax)
            node_stack[node_stack_offset++] = OctreeBFSTraversalNode{0, traversal.nodes[i].idx, 0, traversal.nodes[i].tMin};
    current_node = OctreeBFSTraversalNode{0, traversal.nodes[0].idx, 0, traversal.nodes[0].tMin};
    // Add root bounds to bounds stack
    bounds_stack[bounds_stack_offset] = OctreeBFSTraversalBounds{current_bounds, BoundsHalf(current_bounds)};
    uint16_t current_rank_offset = 0;
    uint16_t current_bitfield_idx = 0;
    while (true) {
        // If the current node is in a different chunk than the previous node,
        // update chunk
        if (current_node.chunk_offset != chunk_offset) {
            chunk_offset = current_node.chunk_offset;
            current_chunk = &(octree[chunk_offset]);
        }
        
        current_bounds = DivideBounds(bounds_stack[bounds_stack_offset].bounds,
                                      current_node.node_idx & 7,
                                      bounds_stack[bounds_stack_offset].b_h);

        uint16_t bitfield_idx = current_node.node_idx / bfsize;
        bool is_inner_node = ((current_chunk->nodes[bitfield_idx]
                >> (current_node.node_idx % bfsize)) & bftone) == bftone;
        
        uint32_t rank = RankAll(current_chunk->nodes, current_node.node_idx);

        if (is_inner_node) {
            // Inner Node...
            rank++;
            ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, current_node.tMin, invDir);
            // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
            Vector3f b_h = BoundsHalf(current_bounds);
            bounds_stack[++bounds_stack_offset] = OctreeBFSTraversalBounds{current_bounds, b_h};
            uint32_t child_idx_offset;
            uint32_t child_chunk_offset;
            if (rank < BFS_NUM_SETS_PER_CHUNK) {
                // ... with children in same chunk
                child_idx_offset = 8 * rank;
                child_chunk_offset = current_node.chunk_offset;
            } else {
                // ... with children in different chunk
                child_idx_offset = 0;
                child_chunk_offset = current_chunk->child_chunk_offset + rank - BFS_NUM_SETS_PER_CHUNK;
            }
            // Add relevant child nodes to node stack
            for (int i = 3; i > 0; i--)
                if (i < traversal.size && traversal.nodes[i].tMin <= ray.tMax)
                    node_stack[node_stack_offset++] = OctreeBFSTraversalNode{
                            child_chunk_offset,
                            child_idx_offset + traversal.nodes[i].idx,
                            (uint8_t)bounds_stack_offset,
                            traversal.nodes[i].tMin};
            current_node = OctreeBFSTraversalNode{
                    child_chunk_offset,
                    child_idx_offset + traversal.nodes[0].idx,
                    (uint8_t)bounds_stack_offset,
                    traversal.nodes[0].tMin};
        } else {
            // Leaf Node
            int leaf_idx = current_chunk->sizes_offset + current_node.node_idx - rank;
            uint32_t prim_start = sizes[leaf_idx - 1];
            uint32_t prim_end = sizes[leaf_idx];
            for (uint32_t i = prim_start; i < prim_end; i++)
                // TODO LRU-Cache/Mailboxing, damit man nicht mehrmals dasselbe primitiv testen muss
                if (leaves[i].get()->IntersectP(ray)) return true;
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
            bounds_stack_offset = current_node.bounds_stack_offset;
        }
    }
    return false;
}

// === VISUALIZATION ===

// Code to visualize octree
void OcChunkBFSAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void OcChunkBFSAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    Chunk c = octree[chunk_offset];

    std::queue<Bounds3f> bounds_q;
    bounds_q.push(bounds);

    int child_chunk_cnt = 0;
    int current_set_idx = 0;

    int num_node_sets = 1; // We have at least 1 node 'set' of 8 nodes in the chunk
    for (int idx = 0; idx < BFS_CHUNK_DEPTH; idx++) {
        for (bftype set = 0; set < NUM_SETS_PER_BITFIELD; set++) {
            if (current_set_idx++ >= num_node_sets) break; // This chunk wasn't completely filled
            
            // Go through current nodes children
            Vector3f b_h = BoundsHalf(bounds_q.front());
            for (bftype bit = 0; bit < 8; bit++) {
                Bounds3f b = DivideBounds(bounds_q.front(), bit, b_h);
                if (((c.nodes[idx] >> (set*8 + bit)) & bftone) == bftone) {
                    if (num_node_sets < BFS_NUM_SETS_PER_CHUNK) {
                        // Inner Node
                        bounds_q.push(b);
                        num_node_sets++;
                    } else {
                        // Inner Node pointing to another chunk
                        lh_dump_rec(f, vcnt_, c.child_chunk_offset + child_chunk_cnt++, b);
                    }
                    // Inner Node
                } else {
                    // Leaf Node
                    // Vertices ausgeben
                    for(uint32_t i = 0; i < 8; i++) {
                        Float x = ((i & 1) == 0) ? b.pMin.x : b.pMax.x;
                        Float y = ((i & 2) == 0) ? b.pMin.y : b.pMax.y;
                        Float z = ((i & 4) == 0) ? b.pMin.z : b.pMax.z;
                        fprintf(f, "v %f %f %f\n", x, y, z);
                    }
                    // Vertex indices ausgeben
                    uint32_t vcnt = *vcnt_;
                    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 5, vcnt + 4);//bottom
                    fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);//top
                    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 3, vcnt + 2);//front
                    fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);//back
                    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 4, vcnt + 6, vcnt + 2);//left
                    fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);//right
                    *vcnt_ += 8;
                }
            }
            bounds_q.pop();
        }
    }
}

void OcChunkBFSAccel::lh_dump_dfs(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec_dfs(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void OcChunkBFSAccel::lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    Chunk c = octree[chunk_offset];

    std::array<Node, 1 + bfsize * BFS_CHUNK_DEPTH> traversal; // Node stack
    traversal[0] = Node{0, bounds};
    int traversal_idx = 0;

    // TODO For schleife mit max möglichen knoten besser?
    while (traversal_idx >= 0) {
        Node node = traversal[traversal_idx--];

        if (node.bitcnt >= 0) {
            // Inner Node ... 
            if (node.bitcnt < BFS_NUM_SETS_PER_CHUNK) {
                // ... with children in same chunk
                Vector3f b_h = BoundsHalf(node.bounds);
                int child_base_rank = BitfieldRankOffset(c.nodes, 8 * node.bitcnt);
                for (int i = 7; i >= 0; i--) {
                    // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
                    // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
                    int idx = 8 * node.bitcnt + i;
                    int bitcnt = child_base_rank + BitfieldRankUnique(c.nodes, idx);
                    bitcnt += IsInnerNode(c.nodes, idx) ? 1 : -(idx + 1);
                    traversal[++traversal_idx] = Node{bitcnt, DivideBounds(node.bounds, i, b_h)};
                }
            } else {
                // ... with chilren in different chunk
                uint32_t child_chunk = c.child_chunk_offset + node.bitcnt - BFS_NUM_SETS_PER_CHUNK;
                lh_dump_rec_dfs(f, vcnt_, child_chunk, node.bounds);
            }
        } else {
            // Leaf Node
            // Vertices ausgeben
            for (uint32_t i = 0; i < 8; i++) {
                Float x = ((i & 1) == 0) ? node.bounds.pMin.x : node.bounds.pMax.x;
                Float y = ((i & 2) == 0) ? node.bounds.pMin.y : node.bounds.pMax.y;
                Float z = ((i & 4) == 0) ? node.bounds.pMin.z : node.bounds.pMax.z;
                fprintf(f, "v %f %f %f\n", x, y, z);
            }
            // Vertex indices ausgeben
            uint32_t vcnt = *vcnt_;
            fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 5, vcnt + 4);//bottom
            fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);//top
            fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 3, vcnt + 2);//front
            fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);//back
            fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 4, vcnt + 6, vcnt + 2);//left
            fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);//right
            *vcnt_ += 8;
        }
    }
}

}  // namespace pbrt
