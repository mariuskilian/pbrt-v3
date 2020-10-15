
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
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>

namespace pbrt {

// - aabb bauen indem ich durch alle primitive iteriere - quadratisch machen - kann auch mal ohne testen

const int MAX_DEPTH = 4;

bool isOverlapping(Float min1, Float max1, Float min2, Float max2) {
    return max1 >= min2 && max2 >= min1;
}

Bounds3f OctreeAccel::octreeDivide(Bounds3f b, int idx) {
    Float xHalf = (b.pMin.x + b.pMax.x) / 2;
    Float yHalf = (b.pMin.y + b.pMax.y) / 2;
    Float zHalf = (b.pMin.z + b.pMax.z) / 2;

    if ((idx & 1) == 0) b.pMin.x = xHalf;
    else b.pMax.x = xHalf;

    if ((idx & 2) == 0) b.pMin.y = yHalf;
    else b.pMax.y = yHalf;

    if ((idx & 4) == 0) b.pMin.z = zHalf;
    else b.pMax.z = zHalf;

    return b;
}

void OctreeAccel::Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth) {
    if (depth > MAX_DEPTH) return;

    std::vector<std::shared_ptr<Primitive>> prims;
    for (int i = 0; i < primitives.size(); i++) {
        std::shared_ptr<Primitive> p = primitives.at(i);
        Bounds3f primBounds = p->WorldBound();
        if (isOverlapping(primBounds.pMin.x, primBounds.pMax.x, bounds.pMin.x, bounds.pMax.x) &&
            isOverlapping(primBounds.pMin.y, primBounds.pMax.y, bounds.pMin.y, bounds.pMax.y) &&
            isOverlapping(primBounds.pMin.z, primBounds.pMax.z, bounds.pMin.z, bounds.pMax.z)) {
                prims.push_back(p);
        }
    }

    if (prims.size() > 1 && depth <= MAX_DEPTH) { // Inner node
        uint32_t offset_children = nodes.size();

        std::vector<uint32_t> nodes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<uint32_t> sizes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        nodes.insert(nodes.end(), nodes_children.begin(), nodes_children.end());
        sizes.insert(sizes.end(), sizes_children.begin(), sizes_children.end());

        nodes[offset] = offset_children << 1 | 0;
        for (uint32_t i = 0; i < 8; i++) {
            Recurse(offset_children + i, prims, octreeDivide(bounds, i), depth + 1);
        }
    } else { // Leaf node
        uint32_t offset_leaves = leaves.size();

        leaves.insert(leaves.end(), prims.begin(), prims.end());

        nodes[offset] = offset_leaves << 1 | 1;
        sizes[offset] = prims.size();
    }
}

// Code to visualize octree
void OctreeAccel::lh_dump_rec(FILE *f, uint *vcnt_, int offset, Bounds3f bounds) {

    // Vertices ausgeben
    for(uint i = 0; i < 8; i++)
    {
        Float x = ((i & 1) == 0) ? bounds.pMin.x : bounds.pMax.x;
        Float y = ((i & 2) == 0) ? bounds.pMin.x : bounds.pMax.x;
        Float z = ((i & 4) == 0) ? bounds.pMin.x : bounds.pMax.x;
        fprintf(f, "v %f %f %f\n", x, y, z);
    }

    // Vertex indices ausgeben
    uint vcnt = *vcnt_;
    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 5, vcnt + 4);//bottom
    fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);//top
    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 3, vcnt + 2);//front
    fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);//back
    fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 4, vcnt + 6, vcnt + 2);//left
    fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);//right
    //fprintf(f, "l %d %d\n", vcnt + 8, vcnt + 10);
    *vcnt_ += 8;

    // Rekursion
    if ((offset & 1) == 0) { // Inner node
        for (uint32_t i = 0; i < 8; i++) {
            lh_dump_rec(f, vcnt_, nodes[offset + i] >> 1, octreeDivide(bounds, i));
        }
    } else { // Leaf node
        Float x = (bounds.pMin.x + bounds.pMax.x) / 2.0;
        Float y = (bounds.pMin.y + bounds.pMax.y) / 2.0;
        Float z = (bounds.pMin.z + bounds.pMax.z) / 2.0;
        fprintf(f, "v, %f %f %f\n", x, y, z);
        *vcnt_ += 1;
        return;
    }
}

void OctreeAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

// KdTreeAccel Method Definitions
OctreeAccel::OctreeAccel(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {
    primitives = p;

    wb.pMin.x = wb.pMax.x = 0;
    wb.pMin.y = wb.pMax.y = 0;
    wb.pMin.z = wb.pMax.z = 0;

    for (int i = 0; i < primitives.size(); i++) {
        Bounds3f b = primitives.at(i)->WorldBound();

        if (b.pMin.x < wb.pMin.x) wb.pMin.x = b.pMin.x;
        if (b.pMin.y < wb.pMin.y) wb.pMin.y = b.pMin.y;
        if (b.pMin.z < wb.pMin.z) wb.pMin.z = b.pMin.z;

        if (b.pMax.x > wb.pMax.x) wb.pMax.x = b.pMax.x;
        if (b.pMax.y > wb.pMax.y) wb.pMax.y = b.pMax.y;
        if (b.pMax.z > wb.pMax.z) wb.pMax.z = b.pMax.z;
    }

    Recurse(0, p, wb, 0);
    lh_dump("/Users/marius/visualization.obj");
}

OctreeAccel::~OctreeAccel() { //FreeAligned(nodes2);
}

bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent

    // Prepare to traverse kd-tree for ray

    // Traverse kd-tree nodes in order for ray
    return false;
}

bool OctreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    // Compute initial parametric range of ray inside kd-tree extent

    // Prepare to traverse kd-tree for ray
    return false;
}

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeAccel>(std::move(prims));
}

}  // namespace pbrt
