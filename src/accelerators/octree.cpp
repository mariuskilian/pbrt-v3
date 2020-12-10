
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

    lh_dump_dfs("visualize_dfs.obj");
}

int rank(std::array<BITFIELD_TYPE, CHUNK_DEPTH> nodes, int node_idx) {
    BITFIELD_TYPE bits;
    BITFIELD_TYPE count = 0;
    for (int i = 0; i < CHUNK_DEPTH; i++) {
        bits = nodes[i];
        if (node_idx < 64) break;
        count += POPCNT(bits);
        node_idx -= 64;
    }
    return count + POPCNT(bits & ((1 << node_idx) - 1));
}

// TODO rename
struct InnerNodeHit { int children_offset; Bounds3f bounds; };
struct ChildHit { int idx; float tMin; };

void shift_array_at(std::array<InnerNodeHit, NUM_SETS_PER_CHUNK + 1> traversal, int fill_size, int idx, int amount) {
    for (int i = fill_size + amount - 1; i > idx; i--) traversal[i] = traversal[i-amount];
}

// TODO Rekursion in Schleife umwandeln (schneller)
void OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t chunk_offset, Bounds3f parent_bounds, bool &hit) const {
    chunk c = octree[chunk_offset];

    // TODO: vielleicht immer nur 1 child pro set auf einmal in array, und platz sparen
    std::array<InnerNodeHit, NUM_SETS_PER_CHUNK + 1> traversal;
    int traversal_idx = 0;
    traversal[traversal_idx] = InnerNodeHit{0, parent_bounds};
    int traversal_size = 1;

    // Handle all nodes in this chunk in the proper order (Depth First Search)
    while (traversal_idx < traversal_size) {
        InnerNodeHit node = traversal[traversal_idx++];

        // Handle note according to node type
        if (node.children_offset == -1) {
            // Leaf Node
            // Intersect primitives
            // if intersection found hit = true
            continue;
        } else if (node.children_offset >= NUM_SETS_PER_CHUNK) {
            // Chunk Pointer Inner Node
            uint32_t cco = c.child_chunk_offset + node.children_offset - c.node_type.size();
            RecurseIntersect(ray, isect, cco, node.bounds, hit);
            continue;
        }

        int child_traversal_size = 0;
        std::array<ChildHit, 8> child_traversal; // It can happen that more than 4 nodes are intersected when using intersectP
        // 1st Step: Intersect all child bounding boxes and determine t parameter
        for (int i = 0; i < 8; i++) {
            //Maybe TODO: try to only octreeDivide once (see above)
            Bounds3f child_bounds = octreeDivide(node.bounds, i);
            ChildHit child_hit = {i};
            float tMax;
            // TODO Optimierte Variante implementieren (3 Ebenentests)
            if (child_bounds.IntersectP(ray, &child_hit.tMin, &tMax)) child_traversal[child_traversal_size++] = child_hit; 
        }
        // 2nd Step: Sort all children by smallest tMin parameter
        // TODO eigene sortierung (bei 8 elementen ist ein naives insertionsort whr. besser)
        std::sort(child_traversal.begin(), child_traversal.begin() + child_traversal_size, [](const ChildHit &a, const ChildHit &b) {return a.tMin < b.tMin;});
        // 3rd Step: Traverse child nodes in order
        // Prepare array by shifting it to the right enough to insert children
        shift_array_at(traversal, traversal_size++, traversal_idx, child_traversal_size);
        int children_set_idx = node.children_offset / NUM_SETS_PER_BITFIELD;
        int children_bitfield_offset = 8 * (node.children_offset % NUM_SETS_PER_BITFIELD);
        for (int i = 0; i < child_traversal_size; i++) {
            // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
            // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
            if (child_traversal[i].tMin > ray.tMax) continue;

            int child_children_offset = -1;
            if (((c.node_type[children_set_idx] >> (children_bitfield_offset + child_traversal[i].idx)) & 1) == 1)
                child_children_offset = rank(c.node_type, node.children_offset * 8 + child_traversal[i].idx) + 1;
            Bounds3f child_bounds = octreeDivide(node.bounds, child_traversal[i].idx);
            traversal[traversal_idx + i] = InnerNodeHit{child_children_offset, child_bounds};
        }
    }
}

bool isInnerNode(BITFIELD_TYPE node_group, int idx) { return ((node_group >> idx) & 1) == 1; }

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
            
        int set_idx = num_nodes / BITFIELD_SIZE;
        int bit_pos = num_nodes % BITFIELD_SIZE;
        // When starting a new index, make sure the initial value of the bitcode is 0
        // TODO: koennte man vorher machen
        if (bit_pos == 0) octree[chunk_idx].node_type[set_idx] = 0;

        // If it's an inner node, properly reserve its children (new chunk or in current chunk)
        if (is_inner_node) {
            // Inner node
            if (chunk_fill_size == NUM_SETS_PER_CHUNK) {
                // Chunk is full, need to make new chunk
                chunk_ptr_nodes.push_back(node_offset);
                octree.push_back(chunk{}); // reserve chunk slot
            } else {
                // Chunk has space for children
                uint32_t child_offset = oba.Nodes()[node_offset] >> 1;
                for (int i = 0; i < 8; i++) { bfs_nodes_q.push(child_offset + i); }
                chunk_fill_size++;
            }
            octree[chunk_idx].node_type[set_idx] |= (int)pow(2,bit_pos);
        } else {
            // handle leaves offset etc.
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
void OctreeAccel::lh_dump_rec_dfs(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    chunk c = octree[chunk_offset];

    // TODO: vielleicht immer nur 1 child pro set auf einmal in array, und platz sparen
    std::array<InnerNodeHit, NUM_SETS_PER_CHUNK + 1> traversal;
    int traversal_idx = 0;
    traversal[traversal_idx] = InnerNodeHit{0, bounds};
    int traversal_size = 1;

    // Handle all nodes in this chunk in the proper order (Depth First Search)
    while (traversal_idx < traversal_size) {
        InnerNodeHit node = traversal[traversal_idx++];

        // Handle note according to node type
        if (node.children_offset == -1) {
            // Leaf Node
            // Vertices ausgeben
            for(uint32_t i = 0; i < 8; i++) {
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
            continue;
        } else if (node.children_offset >= NUM_SETS_PER_CHUNK) {
            // Chunk Pointer Inner Node
            uint32_t cco = c.child_chunk_offset + node.children_offset - NUM_SETS_PER_CHUNK;
            lh_dump_rec_dfs(f, vcnt_, cco, node.bounds);
            continue;
        }

        // Prepare array by shifting it to the right enough to insert children
        shift_array_at(traversal, traversal_size++, traversal_idx, 8);
        int children_set_idx = node.children_offset / NUM_SETS_PER_BITFIELD;
        int children_bitfield_offset = 8 * (node.children_offset % NUM_SETS_PER_BITFIELD);

        for (int i = 0; i < 8; i++) {
            int child_children_offset = -1;
            // Check children for node type and find rank if inner node
            if (((c.node_type[children_set_idx] >> (children_bitfield_offset + i)) & 1) == 1)
                child_children_offset = rank(c.node_type, node.children_offset * 8 + i) + 1;
            Bounds3f child_bounds = octreeDivide(node.bounds, i);
            traversal[traversal_idx + i] = InnerNodeHit{child_children_offset, child_bounds};
        }
    }
}

void OctreeAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, uint32_t chunk_offset, Bounds3f bounds) {
    chunk c = octree[chunk_offset];

    std::queue<Bounds3f> bounds_q;
    bounds_q.push(bounds);

    int child_chunk_cnt = 0;
    int current_set_idx = 0;

    int num_node_sets = 1; // We have at least 1 node 'set' of 8 nodes in the chunk
    for (int idx = 0; idx < CHUNK_DEPTH; idx++) {
        for (int set = 0; set < NUM_SETS_PER_BITFIELD; set++) {
            if (current_set_idx++ == num_node_sets) break; // This chunk wasn't completely filled

            for (int bit = 0; bit < 8; bit++) {
                Bounds3f b = octreeDivide(bounds_q.front(), bit);
                if (((c.node_type[idx] >> (set*8 + bit)) & 1) == 1) {
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
}

void OctreeAccel::lh_dump_dfs(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec_dfs(f, &vcnt, 0, WorldBound());
    fclose(f);
}

void OctreeAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

}  // namespace pbrt
