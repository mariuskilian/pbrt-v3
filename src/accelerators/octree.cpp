
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

const int MAX_DEPTH = 15;
const int MAX_PRIMS = 30;

bool isOverlapping(Float min1, Float max1, Float min2, Float max2) {
    return max1 >= min2 && max2 >= min1;
}

Bounds3f OctreeAccel::octreeDivide(Bounds3f b, int idx) {
    for (int i = 0; i < 3; i++) {
        Float axisHalf = (b.pMin[i] + b.pMax[i]) / 2;
        if ((idx & (int)pow(2,i)) == 0) b.pMax[i] = axisHalf;
        else b.pMin[i] = axisHalf;
    }
    return b;
}

void OctreeAccel::Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth) {        

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

    if (prims.size() > MAX_PRIMS && depth <= MAX_DEPTH) { // Inner node
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

// KdTreeAccel Method Definitions
OctreeAccel::OctreeAccel(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {

    wb.pMin.x = wb.pMax.x = 0;
    wb.pMin.y = wb.pMax.y = 0;
    wb.pMin.z = wb.pMax.z = 0;

    // Determine world bounds
    for (int i = 0; i < primitives.size(); i++) {
        Bounds3f b = primitives.at(i)->WorldBound();

        if (b.pMin.x < wb.pMin.x) wb.pMin.x = b.pMin.x;
        if (b.pMin.y < wb.pMin.y) wb.pMin.y = b.pMin.y;
        if (b.pMin.z < wb.pMin.z) wb.pMin.z = b.pMin.z;

        if (b.pMax.x > wb.pMax.x) wb.pMax.x = b.pMax.x;
        if (b.pMax.y > wb.pMax.y) wb.pMax.y = b.pMax.y;
        if (b.pMax.z > wb.pMax.z) wb.pMax.z = b.pMax.z;
    }

    // Set all dimension sizes of World Bounds to same size
    Float x = wb.pMax.x - wb.pMin.x;
    Float y = wb.pMax.y - wb.pMin.y;
    Float z = wb.pMax.z - wb.pMin.z;

    Float dMax = (x > y) ? ((x > z) ? x : z) : (y > z) ? y : z;
    Float dX = (dMax - x) / 2;
    Float dY = (dMax - y) / 2;
    Float dZ = (dMax - z) / 2;

    wb.pMin.x -= dX;
    wb.pMin.y -= dY;
    wb.pMin.z -= dZ;
    
    wb.pMax.x += dX;
    wb.pMax.y += dY;
    wb.pMax.z += dZ;

    nodes.push_back(0);
    sizes.push_back(0);
    Recurse(0, primitives, wb, 0);
    lh_dump("visualization.obj");
}

OctreeAccel::~OctreeAccel() { //FreeAligned(nodes2);
}

bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);

    Point3f pnt;
    Float length;
    
    // Check intersection with outer planes of WorldBound
    for (int i = 0; i < 6; i++) { // X=0,1; Y=2,3; Z=4,5; Even=Min; Odd=max
        int axis = i / 2; //axis: 0=x; 1=y; 2=z
        bool min = i % 2 == 0; //min: true=min; false=max

        if (ray.d[axis] == 0) continue;

        // Determine factor for formula: point = origin + factor * direction
        Float factor = (min) ? wb.pMin[axis] : wb.pMax[axis];
        factor = (factor - ray.o[axis]) / ray.d[axis];

        if (factor < 0) continue; // Intersection is behind ray

        Point3f pnt_tmp = ray.o + factor * ray.d; // The newly found point of intersection

        // Check if its within the World Bound box
        int otherAxis1 = (axis + 1) % 2;
        int otherAxis2 = (axis + 2) % 2;
        if (pnt_tmp[otherAxis1] > wb.pMax[otherAxis1] ||
            pnt_tmp[otherAxis2] > wb.pMax[otherAxis2] ||
            pnt_tmp[otherAxis1] < wb.pMin[otherAxis1] ||
            pnt_tmp[otherAxis2] < wb.pMin[otherAxis2]) {
                continue; 
        }

        // If its closer than any previously determined point then make this the saved point of intersection
        if (length == 0 || factor * ray.d.Length() < length) {
            pnt = ray.o + factor * ray.d;
            length = factor * ray.d.Length();
        }
    }

    if (length == 0) return false;

    //return RecurseIntersection(ray, isect, wb, pnt, 0);

    return false;
}

bool OctreeAccel::RecurseIntersection(const Ray &ray, SurfaceInteraction *isect,
        Bounds3f bounds, Point3f point, int offset) {

    // Figure out which child (index) the point is touching
    int nodeIdx;
    for (int i = 0; i < 3; i++) {
        Float axisHalf = (bounds.pMin[i] + bounds.pMax[i]) / 2;
        if (point[i] > axisHalf) nodeIdx = nodeIdx | (int)pow(2,i);
    }

    if (isInnerNode(nodes[offset])) { // inner node
        int child_offset = (nodes[offset] >> 1) + nodeIdx;
        if (RecurseIntersection(ray, isect, octreeDivide(bounds, nodeIdx), point, child_offset)) return true; // necessary?

    } else { // leaf node
        bool intersectedPrimitive = false;
        for (uint32_t i = 0; i < sizes[offset]; i++) {
            int leaf_offset = nodes[offset] >> 1;
            if (leaves[leaf_offset + i].get()->Intersect(ray, isect)) intersectedPrimitive = true;
        }
        if (intersectedPrimitive) return true;
    }

    // if no intersection has been found, call the recursive function with the 'next' cube of bounds (neighouring nodeIdx)
    Point3f next_point;
    Float length;
    int flipAxis;
    
    // Check intersection with outer planes of WorldBound
    for (int axis = 0; axis < 3; axis++) { // X=0; Y=1; Z=2

        // Determine factor for formula: point = origin + factor * direction
        Float factor = (bounds.pMin[axis] + bounds.pMax[axis]) / 2;
        factor = (factor - point[axis]) / ray.d[axis];

        Point3f pnt_tmp = point + factor * ray.d; // The newly found point of intersection

        // Check if its within the World Bound box
        int otherAxis1 = (axis + 1) % 2;
        int otherAxis2 = (axis + 2) % 2;
        if (pnt_tmp[otherAxis1] > bounds.pMax[otherAxis1] ||
            pnt_tmp[otherAxis2] > bounds.pMax[otherAxis2] ||
            pnt_tmp[otherAxis1] < bounds.pMin[otherAxis1] ||
            pnt_tmp[otherAxis2] < bounds.pMin[otherAxis2]) {
                continue; 
        }

        // If its closer than any previously determined point then make this the saved point of intersection
        if (length == 0 || factor * ray.d.Length() < length) {
            next_point = point + factor * ray.d;
            length = factor * ray.d.Length();
            flipAxis = axis;
        }
    }

    if (length == 0) return false;

    int next_node_idx = nodeIdx | flipAxis;
    int next_offset = offset - nodeIdx + next_node_idx;

    return RecurseIntersection(ray, isect, bounds, next_point, next_offset);
}

bool OctreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction *isect;
    return Intersect(ray, isect);
}

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeAccel>(std::move(prims));
}

bool isInnerNode(int node) {
    return node & 1 == 0;
}


// === VISUALIZATION ===
// Code to visualize octree
void OctreeAccel::lh_dump_rec(FILE *f, uint *vcnt_, int offset, Bounds3f bounds) {

    // Vertices ausgeben
    for(uint i = 0; i < 8; i++)
    {
        Float x = ((i & 1) == 0) ? bounds.pMin.x : bounds.pMax.x;
        Float y = ((i & 2) == 0) ? bounds.pMin.y : bounds.pMax.y;
        Float z = ((i & 4) == 0) ? bounds.pMin.z : bounds.pMax.z;
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
    *vcnt_ += 8;

    // Rekursion
    if ((nodes[offset] & 1) == 0) { // Inner node
        for (uint32_t i = 0; i < 8; i++) {
            lh_dump_rec(f, vcnt_, (nodes[offset] >> 1) + i, octreeDivide(bounds, i));
        }
    }
}

void OctreeAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

}  // namespace pbrt
