
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

Bounds3f OctreeAccel::octreeDivide(Bounds3f b, int idx) const {
    for (int i = 0; i < 3; i++) {
        Float axisHalf = (b.pMin[i] + b.pMax[i]) / 2;
        if ((idx & (int)pow(2,i)) == 0) b.pMax[i] = axisHalf;
        else b.pMin[i] = axisHalf;
    }
    return b;
}

bool BoundsContainPoint(Bounds3f b, Point3f p) {
    for (int i = 0; i < 2; i++) if (p[i] < b.pMin[i] || p[i] > b.pMax[i]) return false;
    return true;
}

int DetermineNodeIdx(Bounds3f b, Point3f p) {
    int node_idx = 0;
    for (int i = 0; i < 3; i++) {
        Float axisHalf = (b.pMin[i] + b.pMax[i]) / 2;
        if (p[i] > axisHalf) node_idx = node_idx | (int)pow(2,i);
    }
    return node_idx;
}

bool OctreeAccel::IntersectLeafPrims(const Ray &ray, SurfaceInteraction *isect, Bounds3f bounds, uint32_t offset) const{
    bool intersectedPrimitive = false;
    for (uint32_t i = 0; i < sizes[offset]; i++) {
        int leaf_offset = nodes[offset] >> 1;
        if (leaves[leaf_offset + i].get()->Intersect(ray, isect)) intersectedPrimitive = true;
    }
    return intersectedPrimitive && BoundsContainPoint(bounds, isect->p);
}

struct node_isect { int flip_mask; Float factor; };

bool CloserNode(node_isect n1, node_isect n2) {
    return n1.factor < n2.factor;
}

// Returns a sorted list of all cut half-axis
std::vector<node_isect> FindChildrenTraverseOrder(const Ray &ray, Bounds3f bounds) {
    std::vector<node_isect> node_order = {node_isect{0, 0}}; // "base" node
    for (int axis = 0; axis < 3; axis++) { 
        if (ray.d[axis] == 0) continue;

        node_isect node = {(int)pow(2,axis), ((bounds.pMin[axis] + bounds.pMax[axis]) / 2 - ray.o[axis]) / ray.d[axis]};

        if (node.factor >= 0 && node.factor <= ray.tMax && BoundsContainPoint(bounds, ray.o + node.factor * ray.d)) 
            node_order.push_back(node);
    }
    std::sort(node_order.begin(), node_order.end(), CloserNode);
    return node_order;
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

    Point3f point = ray.o;

    if (!BoundsContainPoint(wb, ray.o)) {
        // Check intersection with outer planes of WorldBound and find closest point
        Float factor = ray.tMax;
        for (int i = 0; i < 6; i++) { // X=0,1; Y=2,3; Z=4,5; Even=Min; Odd=max
            int axis = i / 2; bool min = i % 2 == 0;
            if (ray.d[axis] == 0) continue; 

            // Determine factor for formula: point = origin + factor * direction; Then check if factor is legal and closer
            Float f = (((min) ? wb.pMin[axis] : wb.pMax[axis]) - ray.o[axis]) / ray.d[axis];
            if (f >= 0 && f < factor && BoundsContainPoint(wb, ray.o + f * ray.d)) factor = f;
        }

        if (factor == ray.tMax) return false;
        point = ray.o + factor * ray.d;
    }

    return RecurseIntersect(ray, isect, wb, 0, point);
}

bool OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect,
        Bounds3f bounds, uint32_t offset, Point3f point) const {
    if ((nodes[offset] & 1) == 0) { // Inner node
        std::vector<node_isect> children = FindChildrenTraverseOrder(ray, bounds);
        int child_node_idx = DetermineNodeIdx(bounds, point);
        for (int i = 0; i < children.size(); i++) {
            child_node_idx ^= children[i].flip_mask;
            if (RecurseIntersect(ray, isect,
                    octreeDivide(bounds, child_node_idx),
                    (nodes[offset] >> 1) + child_node_idx,
                    ray.o + children[i].factor * ray.d))
                return true;
        }
    } else return IntersectLeafPrims(ray, isect, bounds, offset); // Leaf node

    return false;
    //int node_idx = (offset - 1) % 8;
    //Bounds3f node_bounds = octreeDivide(bounds, node_idx);
    //
    //if ((nodes[offset] & 1) == 0) {// inner node
    //    uint32_t child_offset = (nodes[offset] >> 1) + DetermineNodeIdx(node_bounds, point);
    //    if (RecurseIntersection(ray, isect, node_bounds, child_offset, point)) return true;
    //} else if (IntersectLeafPrims(ray, isect, node_bounds, offset)) return true;
    //
    //
    //// === TRAVERSE NEIGHBOURING NODES ===
    //struct node_isect { int flip_axis; Float factor; };
    //std::vector<node_isect> node_order; // holds order of indices of neighbouring nodes
    //for (int axis = 0; axis < 3; axis++) { // Find order of neighbouring nodes
    //    if (ray.d[axis] == 0) continue;
    //
    //    node_isect node;
    //
    //    // Point = Origin + factor * Direction
    //    node.factor = ((bounds.pMin[axis] + bounds.pMax[axis]) / 2 - ray.o[axis]) / ray.d[axis];
    //    if (node.factor < 0 || node.factor > ray.tMax) continue;
    //    if (!IsPointWithinBounds(node_bounds, ray.o + node.factor * ray.d)) continue;
    //    node.flip_axis = axis;
    //
    //    if (node_order.size() == 0) { node_order.push_back(node); continue; }
    //
    //    // Insert found point in the correct place into the 2 vectors
    //    for (int i = 0; i < node_order.size(); i++) {
    //        if (node.factor >= node_order[i].factor) continue;
    //
    //        // Insert node into index i
    //        node_order.push_back(node); // Placeholder
    //        for (int j = node_order.size() - 1; j > i; j--) node_order[j] = node_order[j-1];
    //        node_order[i] = node;
    //        break;
    //    }
    //}
    //
    //int next_node_idx = node_idx;
    //for (int i = 0; i < node_order.size(); i++) {
    //    next_node_idx = next_node_idx | (int)pow(2,node_order[i].flip_axis);
    //    uint32_t next_offset = offset - node_idx + next_node_idx;
    //    if (RecurseIntersection(ray, isect, bounds, next_offset, ray.o + node_order[i].factor * ray.d)) return true;
    //}
    //// === TRAVERSE NEIGHBOURING NODES END ===
    //
    //return false;
}

bool OctreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction *isect;
    return Intersect(ray, isect);
}

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeAccel>(std::move(prims));
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
