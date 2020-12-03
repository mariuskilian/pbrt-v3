
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
#elif defined(__clang__)
    #include "popcntintrin.h"
#endif

namespace pbrt {

const int MAX_DEPTH = 15; //15
const int MAX_PRIMS = 32; //30

Bounds3f OctreeAccel::octreeDivide(Bounds3f b, int idx) const {
    for (int i = 0; i < 3; i++) {
        Float axisHalf = (b.pMin[i] + b.pMax[i]) / 2;
        if ((idx & (1<<i)) == 0) b.pMax[i] = axisHalf;
        else b.pMin[i] = axisHalf;
    }
    return b;
}

// KdTreeAccel Method Definitions
OctreeAccel::OctreeAccel(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {
    oba = OctreeBasicAccel(primitives);
    wb = oba.WorldBound();
    
    if (oba.Nodes().size() > 1) {
        octree.push_back(chunk{});
        Recurse(0, 0);
    }

    lh_dump("visualize.obj");
}

struct InnerNodeHit { int children_offset; Bounds3f bounds; };
struct ChildHit { int idx; float tMin; };

// TODO Rekursion in Schleife umwandeln (schneller)
void OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t chunk_offset, Bounds3f parent_bounds, bool &hit) const {
    chunk c = octree[chunk_offset];

    std::deque<InnerNodeHit> traversal_dq; // Only Inner nodes with children in the same chunk will be in this dq
    traversal_dq.push_back(InnerNodeHit{0, parent_bounds}); // TODO: maybe change init idx and tMin values?

    while (!traversal_dq.empty()) {
        InnerNodeHit node = traversal_dq.front();
        traversal_dq.pop_front();

        int traversal_size = 0;
        std::array<ChildHit, 8> traversal; // It can happen that more than 4 nodes are intersected when using intersectP
        // 1. Step: Intersect all child bounding boxes and determine t parameter
        for (int i = 0; i < 8; i++) {
            Bounds3f child_bounds = octreeDivide(node.bounds, i);
            ChildHit child_hit = {i};
            float tMax;
            // TODO Optimierte Variante implementieren (3 Ebenentests)
            if (child_bounds.IntersectP(ray, &child_hit.tMin, &tMax)) {
                //if (traversal_size == 8) abort(); // sanity check
                traversal[traversal_size++] = child_hit; 
            }
            // Invert tMin of real inner nodes, so they're in the inverted order, and pushing them to the front \
                of traversal_dq later is easier. tMin is otherwise irrelevant for these nodes
            if (((c.node_type[node.children_offset] >> i) & 1) == 0 && node.children_offset < c.node_type.size())
                traversal[i].tMin = 1 / ++traversal[i].tMin; // ++ to avoid div by 0
        }
        // 2. Step: Sort all children by smallest tMin parameter
        // TODO eigene sortierung (bei 8 elementen ist ein naives insertionsort whr. besser)
        std::sort(traversal.begin(), traversal.begin() + traversal_size, [](const ChildHit &a, const ChildHit &b) {return a.tMin < b.tMin;});
        // 3. Step: Traverse child nodes in order
        for (int i = 0; i < 8; i++) {
            // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
            // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
            if (i >= traversal_size || traversal[i].tMin > ray.tMax) continue;
            // Otherwise handle node depending on type: \
                   - Inner node: insert into front of traversal_dq \
                   - Inner node pointing to chunk: Call recursion \
                   - Leaf node: check intersections with primitives
            if (((c.node_type[node.children_offset] >> i) & 1) == 0) {
                if (node.children_offset < c.node_type.size()) {
                    // Inner Node
                    InnerNodeHit inh = InnerNodeHit{0, octreeDivide(node.bounds, i)}; //TODO: replace 0 with rank(node.children_offset) \
                                                                                        and try to only do octreeDivide once
                    traversal_dq.push_front(inh);
                } else {
                    // Inner Node pointing to another chunk
                    
                }
            } else {
                // Leaf Node
            }
        }
    }
}

bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
    if (!wb.IntersectP(ray)) return false;
    RecurseIntersect(ray, isect, 0, wb, hit);
    return hit;
}

void OctreeAccel::Recurse(uint32_t root_node_offset, int chunk_idx) {
    octree[chunk_idx].child_chunk_offset = octree.size();
    // c.leaf_offset = ...

    uint32_t root_child_offset = oba.Nodes()[root_node_offset] >> 1;
    std::queue<uint32_t> bfs_nodes_q;
    for (int i = 0; i < 8; i++) { bfs_nodes_q.push(root_child_offset + i); }

    int chunk_fill_size = 1; // number of [set of 8 nodes] reserved in chunk
    int num_nodes = 0; // number of individual nodes already processed
    std::vector<uint32_t> chunk_ptr_nodes; // nodes that point to other chunks

    uint32_t node_offset;
    bool is_inner_node;

    while (!bfs_nodes_q.empty()) {
        node_offset = bfs_nodes_q.front();
        is_inner_node = (oba.Nodes()[node_offset] & 1) == 0;

        // If it's an inner node, properly reserve its children (new chunk or in current chunk)
        if (is_inner_node) {
            // Inner node
            if (chunk_fill_size == chunk_depth) {
                // Chunk is full, need to make new chunk
                chunk_ptr_nodes.push_back(node_offset);
                octree.push_back(chunk{}); // reserve chunk slot
            } else {
                // Chunk has space for children
                uint32_t child_offset = oba.Nodes()[node_offset] >> 1;
                for (int i = 0; i < 8; i++) { bfs_nodes_q.push(child_offset + i); }
                chunk_fill_size++;
            }
        } else {
            // Leaf node
            int arr_idx = num_nodes / 8;
            int bit_pos = num_nodes % 8;
            // When starting a new index, make sure the initial value of the bitcode is 0
            if (bit_pos == 0) octree[chunk_idx].node_type[arr_idx] = 0;
            octree[chunk_idx].node_type[arr_idx] |= (int)pow(2,bit_pos);
        }

        num_nodes++;
        bfs_nodes_q.pop();
    } 
    
    // Create additional chunks as needed
    for (int i = 0; i < chunk_ptr_nodes.size(); i++)
        Recurse(chunk_ptr_nodes[i], octree[chunk_idx].child_chunk_offset + i);

}

bool OctreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction isect;
    return Intersect(ray, &isect);
}

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeAccel>(std::move(prims));
}

OctreeAccel::~OctreeAccel() { //FreeAligned(nodes2);
}

// === VISUALIZATION ===

// Code to visualize octree
void OctreeAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    chunk c = octree[chunk_offset];

    std::queue<Bounds3f> bounds_q;
    bounds_q.push(bounds);

    int child_chunk_cnt = 0;

    int num_node_sets = 1; // We have at least 1 node 'set' of 8 nodes in the chunk
    for (int idx = 0; idx < c.node_type.size(); idx++) {
        if (idx == num_node_sets) break; // This chunk wasn't completely filled
        for (int bit = 0; bit < 8; bit++) {
            Bounds3f b = octreeDivide(bounds_q.front(), bit);
            if (((c.node_type[idx] >> bit) & 1) == 0) {
                if (num_node_sets < c.node_type.size()) {
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

void OctreeAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

}  // namespace pbrt
