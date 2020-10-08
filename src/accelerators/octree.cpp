
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

// FRAGEN
// - ist bei sizes der default wert 0?
// - kann man sizes kleiner machen (viele 0en)?
// - wie zeigt leaves auf die primitive? was wenn mehrere primitive im knoten?
// - wie speicher ich bounding boxes? noch ein vector?
// - wie greife ich auf den raum zu, den ich aufteile? min/max xyz?
// - Kann ich if-Bedingung fuer bounds vereinfachen?
// - greif ich korrekt und auf die richtigen bounds zu? aka ist prim.WorldBound das prim aabb?

const int MAX_DEPTH = 10;
std::vector<uint32_t> nodes, sizes; 
std::vector<std::shared_ptr<Primitive>> leaves;

bool isOverlapping(Float min1, Float max1, Float min2, Float max2) {
    return max1 >= min2 && max2 >= min1;
}

Bounds3f octreeDivide(Bounds3f bounds, int idx) {
    Float xHalf = (bounds.pMin.x + bounds.pMax.x) / 2;
    Float yHalf = (bounds.pMin.y + bounds.pMax.y) / 2;
    Float zHalf = (bounds.pMin.z + bounds.pMax.z) / 2;

    if (idx % 2 == 0) bounds.pMin.x = xHalf;
    else bounds.pMax.x = xHalf;

    if ((idx >> 1) % 2 == 0) bounds.pMin.y = yHalf;
    else bounds.pMax.y = yHalf;

    if ((idx >> 2) % 2 == 0) bounds.pMin.z = zHalf;
    else bounds.pMax.z = zHalf;

    return bounds;
}

void Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth) {
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
        std::vector<uint32_t> nodes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<uint32_t> sizes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        nodes.insert(nodes.end(), nodes_children.begin(), nodes_children.end());
        sizes.insert(sizes.end(), sizes_children.begin(), sizes_children.end());

        uint32_t offset_children = nodes.size();
        nodes[offset] = offset_children << 1 | 0;
        for (int i = 0; i < 8; i++) {
            Recurse(offset_children, prims, octreeDivide(bounds, i), depth + 1);
        }
    } else { // Leaf node
        leaves.insert(leaves.end(), prims.begin(), prims.end());

        uint32_t offset_leaves = leaves.size();
        nodes[offset] = offset_leaves << 1 | 1;
        sizes[offset] = prims.size();
    }
}

// KdTreeAccel Method Definitions
OctreeAccel::OctreeAccel(std::vector<std::shared_ptr<Primitive>> p)
    : primitives(std::move(p)) {

    // Start recursive construction of kd-tree
    // buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
    //           maxDepth, edges, prims0.get(), prims1.get());
}

OctreeAccel::~OctreeAccel() { FreeAligned(nodes); }

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
