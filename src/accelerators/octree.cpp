
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
#include "accelerators/octree.h"
#include "accelerators/octree-basic.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>
#include <array>
#include <queue>

#if defined(_MSC_VER)
    #include "intrin.h"
    #define POPCNT __popcnt64
#elif defined(__clang__)
    #include "popcntintrin.h"
    #define POPCNT _mm_popcnt_u64
#endif

namespace pbrt {

// === HELPERS ===

struct Node { int bitcnt; Bounds3f bounds; Float tMin; };
struct ChildHit { int idx; float tMin; };
struct ChildTraversal{ std::array<ChildHit, 4> nodes; int size; };


// TODO extract helper functions in octree-basic, then refer to those from here
// TODO give axisHalf as parameter; Inline
Vector3f BoundsHalf(Bounds3f b) {
    Vector3f h;
    for (int i = 0; i < 3; i++) h[i] = (b.pMin[i] + b.pMax[i]) / 2;
    return h;
}

Bounds3f DivideBounds(Bounds3f b, int idx, Vector3f b_half) {
    for (int i = 0; i < 3; i++) {
        if ((idx & (1<<i)) == 0) b.pMax[i] = b_half[i];
        else b.pMin[i] = b_half[i];
    }
    return b;
}

bool IsInnerNode(std::array<BITFIELD_TYPE, CHUNK_DEPTH> bitfield, int n) {
    int i = n / BITFIELD_SIZE;
    int offset = n % BITFIELD_SIZE;
    return ((bitfield[i] >> offset) & 1) == 1;
}

int Rank(std::array<BITFIELD_TYPE, CHUNK_DEPTH> bitfield, int n) {
    int count = 0;
    BITFIELD_TYPE bits;
    for (int i = 0; i < CHUNK_DEPTH; i++) {
        bits = bitfield[i];
        if (n < BITFIELD_SIZE) break;
        count += POPCNT(bits);
        n -= BITFIELD_SIZE;
    }
    count += POPCNT(bits & ((ONE << n) - ONE));
    return count;
}

ChildTraversal FindTraversalOrder(const Ray &ray, Bounds3f b, Float tMin) {
    int size = 1;
    std::array<ChildHit, 4> traversal;
    Vector3f b_h = BoundsHalf(b);

    // First child hit
    int idx = 0;
    Point3f init_point = ray.o + tMin * ray.d;
    for (int i = 0; i < 3; i++) if (init_point[i] > b_h[i]) idx |= (1<<i);
    traversal[0] = ChildHit{ idx, tMin };

    // Cut all bound-half-planes, and if intersection is within bounds, add to list
    for (int axis = 0; axis < 3; axis++) {
        Float t = (b_h[axis] - ray.o[axis]) / ray.d[axis];
        Point3f p = ray.o + t * ray.d;
        // Discard point if it's outside of the bounds
        for (int i = 0; i < 3; i++) if (p[i] < b.pMin[i] || p[i] > b.pMax[i]) continue;
        traversal[size++] = ChildHit{ 1<<axis, t };
    }

    // Sort list based on smallest tMin
    std::sort(traversal.begin() + 1, traversal.begin() + size, 
        [](const ChildHit &c1, const ChildHit &c2) {return c1.tMin < c2.tMin;});
    
    // Finally, determine idx for each child hit
    for (int i = 1; i < size; i++) traversal[i].idx ^= traversal[i-1].idx;

    return ChildTraversal{traversal, size};
}

// inline ChildTraversal FindTraversalOrder(const Ray &ray, Bounds3f b, Float tMin) {
//         int size = 0;
//         std::array<ChildHit, 8> traversal; // It can happen that more than 4 nodes are intersected when using intersect
//         // Pre calculate bounds half, since they are needed for every child
//         Vector3f b_h = BoundsHalf(b);
//         // 1st Step: Intersect all child bounding boxes and determine t parameter
//         for (int i = 0; i < 8; i++) {
//             Bounds3f child_bounds = DivideBounds(b, i, b_h);
//             ChildHit child_hit = { i, 0, child_bounds };
//             float tMax;
//             // TODO Optimierte Variante implementieren (3 Ebenentests)
//             if (child_bounds.IntersectP(ray, &child_hit.tMin, &tMax)) traversal[size++] = child_hit; 
//         }
//         // 2nd Step: Sort all children by smallest tMin parameter
//         // TODO eigene sortierung (bei 8 elementen ist ein naives insertionsort whr. besser)
//         std::sort(traversal.begin(), traversal.begin() + size,
//             [](const ChildHit &a, const ChildHit &b) {return a.tMin < b.tMin;});
//         return ChildTraversal{traversal, size};
// }

// === OCTREE STRUCT CREATION ==
OctreeAccel::OctreeAccel(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {
    oba = OctreeBasicAccel(primitives);
    wb = oba.WorldBound();
    
    if (oba.Nodes().size() > 1) {
        octree.push_back(Chunk{});
        sizes.push_back(0);
        Recurse(0, 0);
    }
    lh_dump("visualize_bfs.obj");
    lh_dump_dfs("visualize_dfs.obj");
}

void OctreeAccel::Recurse(uint32_t root_node_offset, int chunk_idx) {
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
            
        int idx = num_nodes / BITFIELD_SIZE;
        int bit_pos = num_nodes % BITFIELD_SIZE;
        // When starting a new index, make sure the initial value of the bitcode is 0
        if (bit_pos == 0) octree[chunk_idx].nodes[idx] = ZERO;

        // If it's an inner node, properly reserve its children (new chunk or in current chunk)
        if (is_inner_node) {
            // Inner node
            if (chunk_fill_size == NUM_SETS_PER_CHUNK) {
                // Chunk is full, need to make new chunk
                chunk_ptr_nodes.push_back(node_offset);
                octree.push_back(Chunk{}); // reserve chunk slot
            } else {
                // Chunk has space for children
                uint32_t child_offset = oba.Nodes()[node_offset] >> 1;
                for (int i = 0; i < 8; i++) { bfs_nodes_q.push(child_offset + i); }
                chunk_fill_size++;
            }
            octree[chunk_idx].nodes[idx] |= ONE << bit_pos;
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

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeAccel>(std::move(prims));
}

OctreeAccel::~OctreeAccel() { //FreeAligned(nodes2);
}

// === OCTREE RAY TRAVERSAL ===

// TODO Rekursion in Schleife umwandeln (schneller)
bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return false;
    RecurseIntersect(ray, isect, 0, wb, tMin, hit);
    return hit;
}

void OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t chunk_offset, Bounds3f parent_bounds, Float tMin, bool &hit) const {
    Chunk c = octree[chunk_offset];

    std::array<Node, 1 + BITFIELD_SIZE * CHUNK_DEPTH> traversal; // Node stack
    traversal[0] = Node{0, parent_bounds, tMin};
    int traversal_idx = 0;

    // TODO For schleife mit max möglichen knoten besser?
    while (traversal_idx >= 0) {
        Node node = traversal[traversal_idx--];

        if (node.bitcnt >= 0) {
            // Inner Node ... 
            if (node.bitcnt < NUM_SETS_PER_CHUNK) {
                // ... with children in same chunk
                ChildTraversal child_traversal = FindTraversalOrder(ray, node.bounds, node.tMin);
                Vector3f b_h = BoundsHalf(node.bounds);
                for (int i = child_traversal.size - 1; i >= 0; i--) {
                    // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
                    // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
                    if (child_traversal.nodes[i].tMin > ray.tMax) continue;

                    int idx = 8 * node.bitcnt + child_traversal.nodes[i].idx;
                    int bitcnt = Rank(c.nodes, idx) + (IsInnerNode(c.nodes, idx) ? 1 : -(idx + 1));
                    Bounds3f child_bounds = DivideBounds(node.bounds, child_traversal.nodes[i].idx, b_h);
                    traversal[++traversal_idx] = Node{bitcnt, child_bounds, child_traversal.nodes[i].tMin};
                }
            } else {
                // ... with children in different chunk
                uint32_t child_chunk = c.child_chunk_offset + node.bitcnt - NUM_SETS_PER_CHUNK;
                RecurseIntersect(ray, isect, child_chunk, node.bounds, node.tMin, hit);
            }
        } else {
            // Leaf Node
            int leaf_idx = c.sizes_offset - node.bitcnt - 1;
            uint32_t prim_start = sizes[leaf_idx - 1];
            uint32_t prim_end = sizes[leaf_idx];
            for (uint32_t i = prim_start; i < prim_end; i++)
                // TODO LRU-Cache/Mailboxing, damit man nicht mehrmals dasselbe primitiv testen muss
                if (leaves[i].get()->Intersect(ray, isect)) hit = true;
        }
    }
}

bool OctreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction isect;
    return Intersect(ray, &isect);
}

// === VISUALIZATION ===

// Code to visualize octree
void OctreeAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void OctreeAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    Chunk c = octree[chunk_offset];

    std::queue<Bounds3f> bounds_q;
    bounds_q.push(bounds);

    int child_chunk_cnt = 0;
    int current_set_idx = 0;

    int num_node_sets = 1; // We have at least 1 node 'set' of 8 nodes in the chunk
    for (int idx = 0; idx < CHUNK_DEPTH; idx++) {
        for (BITFIELD_TYPE set = 0; set < NUM_SETS_PER_BITFIELD; set++) {
            if (current_set_idx++ >= num_node_sets) break; // This chunk wasn't completely filled
            
            // Go through current nodes children
            Vector3f b_h = BoundsHalf(bounds_q.front());
            for (BITFIELD_TYPE bit = 0; bit < 8; bit++) {
                Bounds3f b = DivideBounds(bounds_q.front(), bit, b_h);
                if (((c.nodes[idx] >> (set*8 + bit)) & ONE) == ONE) {
                    if (num_node_sets < NUM_SETS_PER_CHUNK) {
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

void OctreeAccel::lh_dump_dfs(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec_dfs(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void OctreeAccel::lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    Chunk c = octree[chunk_offset];

    std::array<Node, 1 + BITFIELD_SIZE * CHUNK_DEPTH> traversal; // Node stack
    traversal[0] = Node{0, bounds};
    int traversal_idx = 0;

    // TODO For schleife mit max möglichen knoten besser?
    while (traversal_idx >= 0) {
        Node node = traversal[traversal_idx--];

        if (node.bitcnt >= 0) {
            // Inner Node ... 
            if (node.bitcnt < NUM_SETS_PER_CHUNK) {
                // ... with children in same chunk
                Vector3f b_h = BoundsHalf(node.bounds);
                for (int i = 7; i >= 0; i--) {
                    // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
                    // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
                    int idx = 8 * node.bitcnt + i;
                    int bitcnt = Rank(c.nodes, idx) + (IsInnerNode(c.nodes, idx) ? 1 : -(idx + 1));
                    traversal[++traversal_idx] = Node{bitcnt, DivideBounds(node.bounds, i, b_h)};
                }
            } else {
                // ... with chilren in different chunk
                uint32_t child_chunk = c.child_chunk_offset + node.bitcnt - NUM_SETS_PER_CHUNK;
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
