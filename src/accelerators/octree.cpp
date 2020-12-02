
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
#include "popcntintrin.h"

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
        octree.push_back({});
        Recurse(0, 0);
    }
}

// TODO Rekursion in Schleife umwandeln (schneller)
void OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t offset, Bounds3f bounds, bool &hit) const {
}

bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
}

void OctreeAccel::Recurse(uint32_t root_node_offset, int chunk_idx) { 
    chunk c = octree[chunk_idx];
    c.child_chunk_offset = octree.size();
    // c.leaf_offset = ...

    uint32_t root_child_offset = oba.Nodes()[root_node_offset >> 1];
    std::queue<uint32_t> bfs_nodes_q;
    for (int i = 0; i < 8; i++) { bfs_nodes_q.push(root_child_offset + i); }

    int chunk_fill_size = 1; // number of [set of 8 nodes] reserved in chunk
    int num_nodes = 0; // number of individual nodes already processed
    std::vector<uint32_t> inner_leaf_nodes; // nodes that point to other chunks

    uint32_t node_offset;
    bool is_leaf_node;
    bool is_inner_leaf; // "leaf" pointing to another chunk

    while (!bfs_nodes_q.empty()) {
        node_offset = bfs_nodes_q.front();
        is_leaf_node = false;
        is_inner_leaf = false;

        if (oba.Nodes()[node_offset] & 1 == 0) {
            // Inner node
            if (chunk_fill_size == chunk_depth) {
                // Chunk is full, need to make new chunk
                inner_leaf_nodes.push_back(node_offset);
                octree.push_back({}); // reserve chunk slot
                is_leaf_node = is_inner_leaf = true;
            } else {
                // Chunk has space for children
                uint32_t child_offset = oba.Nodes()[node_offset] >> 1;
                for (int i = 0; i < 8; i++) { bfs_nodes_q.push(child_offset + i); }
                chunk_fill_size++;
            }
        } else is_leaf_node = true;

        int arr_idx = num_nodes / 8;
        int bit_pos = num_nodes % 8;
        // When starting a new index, make sure the initial value of the bitcode is 0
        if (bit_pos == 0) { c.node_type[arr_idx] = 0; c.leaf_type[arr_idx] = 0; }
        // Set the correct bits accordingly
        if (is_leaf_node) c.node_type[arr_idx] |= (int)pow(2,bit_pos);
        if (is_inner_leaf) c.leaf_type[arr_idx] |= (int)pow(2,bit_pos);

        // Create additional chunks as needed
        for (int i = 0; i < inner_leaf_nodes.size(); i++)
            Recurse(inner_leaf_nodes[i], c.child_chunk_offset + i);

        num_nodes++;
        bfs_nodes_q.pop();
    } 

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

    std::queue<Bounds3f> bounds_queue;
    bounds_queue.push(bounds);
    for (int i = 0; i < 8; i++) {
        bounds_queue.push(bounds
        );
    }

    int num_nodes = 0;
    while (!bounds_queue.empty()) {

    }

    for (int idx = 0; idx < c.node_type.size(); idx++) {
        for (int bit_pos = 0; bit_pos < 8; bit_pos++) {

            if (((c.node_type[idx] >> bit_pos) & 1) == 0) {
                // Inner node

            } else if (((c.leaf_type[idx] >> bit_pos) & 1) == 0) {
                // Chunk pointer node
            } else {
                // Leaf node
            }

        }
    }

    // Vertices ausgeben
    for(uint32_t i = 0; i < 8; i++)
    {
        Float x = ((i & 1) == 0) ? bounds.pMin.x : bounds.pMax.x;
        Float y = ((i & 2) == 0) ? bounds.pMin.y : bounds.pMax.y;
        Float z = ((i & 4) == 0) ? bounds.pMin.z : bounds.pMax.z;
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

void OctreeAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

}  // namespace pbrt
