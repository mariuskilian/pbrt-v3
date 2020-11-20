
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

const int MAX_DEPTH = 15; //15
const int MAX_PRIMS = 30; //30

int idx; //vis
std::vector<Bounds3f> ray_leaves; //vis

bool BoundsOverlap(Bounds3f b1, Bounds3f b2) {
    for (int i = 0; i < 3; i++) if (b1.pMax[i] < b2.pMin[i] || b2.pMax[i] < b1.pMin[i]) return false;
    return true;
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
    for (int i = 0; i < 3; i++) if (p[i] < b.pMin[i] || p[i] > b.pMax[i]) return false;
    return true;
}

int DetermineNodeIdx(Bounds3f b, Point3f p) {
    int node_idx = 0;
    for (int i = 0; i < 3; i++) if (p[i] > (b.pMin[i] + b.pMax[i]) / 2) node_idx |= (int)pow(2,i);
    return node_idx;
}

bool OctreeAccel::IntersectLeafPrims(const Ray &ray, SurfaceInteraction *isect, Bounds3f bounds, uint32_t offset) const{
    bool hit = false;
    for (uint32_t i = 0; i < sizes[offset]; i++) {
        int leaf_offset = nodes[offset] >> 1;
        if (leaves[leaf_offset + i].get()->Intersect(ray, isect)) hit = true;
    }
    return hit && BoundsContainPoint(bounds, isect->p);
}

struct node_isect { int flip_mask; Float tMin; };

bool CloserNode(node_isect n1, node_isect n2) { return n1.tMin < n2.tMin; }

// Returns a sorted list of all cut half-axis
std::vector<node_isect> FindChildrenTraverseOrder(const Ray &ray, Bounds3f bounds) {
    std::vector<node_isect> node_order = {node_isect{0, 0}}; // "base" node
    for (int axis = 0; axis < 3; axis++) { 
        if (ray.d[axis] == 0) continue;

        node_isect node = {(int)pow(2,axis), (((bounds.pMin[axis] + bounds.pMax[axis]) / 2) - ray.o[axis]) / ray.d[axis]};

        if (node.tMin >= 0 && node.tMin <= ray.tMax && BoundsContainPoint(bounds, ray.o + node.tMin * ray.d)) 
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
        if (BoundsOverlap(primBounds, bounds)) prims.push_back(p);
    }

    if (prims.size() > MAX_PRIMS && depth < MAX_DEPTH) { // Inner node
        uint32_t offset_children = nodes.size();

        std::vector<uint32_t> nodes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<uint32_t> sizes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        nodes.insert(nodes.end(), nodes_children.begin(), nodes_children.end());
        sizes.insert(sizes.end(), sizes_children.begin(), sizes_children.end());

        nodes[offset] = offset_children << 1 | 0;
        for (uint32_t i = 0; i < 8; i++) Recurse(offset_children + i, prims, octreeDivide(bounds, i), depth + 1);
    } else { // Leaf node
        uint32_t offset_leaves = leaves.size();

        leaves.insert(leaves.end(), prims.begin(), prims.end());

        nodes[offset] = offset_leaves << 1 | 1;
        sizes[offset] = prims.size();
    }
}

// KdTreeAccel Method Definitions
OctreeAccel::OctreeAccel(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {

    for (int i = 0; i < 3; i++) wb.pMin[i] = wb.pMax[i] = 0;

    // Determine world bounds
    for (int i = 0; i < primitives.size(); i++) {
        Bounds3f b = primitives.at(i)->WorldBound();
        for (int i = 0; i < 3; i++) {
            if (b.pMin[i] < wb.pMin[i]) wb.pMin[i] = b.pMin[i];
            if (b.pMax[i] > wb.pMax[i]) wb.pMax[i] = b.pMax[i];
        }
    }

    // Set all dimension sizes of World Bounds to same size
    Vector3f size;
    Float max_size = 0;
    for (int i = 0; i < 3; i++) {
        size[i] = wb.pMax[i] - wb.pMin[i];
        if (size[i] > max_size) max_size = size[i];
    }

    for (int i = 0; i < 3; i++) {
        Float diff = (max_size - size[i]) / 2;
        wb.pMin[i] -= diff;
        wb.pMax[i] += diff;
    }
    
    nodes.push_back(0);
    sizes.push_back(0);
    Recurse(0, primitives, wb, 0);

    // Ray ray = Ray();
    // ray.o = Point3f(400, 20, 30);
    // ray.d = Vector3f(105 - 400, 140 - 20, 13 - 30);
    // SurfaceInteraction isect;
    // Intersect(ray, &isect);

    // lh_dump("visualization.obj");

    // Point3f t_p = Point3f(100, 115, -140);
    // uint32_t t_o = 0;
    // Bounds3f t_b = wb;
    // while (true) {
    //     if ((nodes[t_o] & 1) == 1) break;
    //     int nIdx = DetermineNodeIdx(t_b, t_p);
    //     t_b = octreeDivide(t_b, nIdx);
    //     t_o = (nodes[t_o] >> 1) + nIdx;
    // }
    // std::vector<std::shared_ptr<Primitive>> prims;
    // for (uint32_t i = 0; i < sizes[t_o]; i++) {
    //     int leaf_offset = nodes[t_o] >> 1;
    //     prims.push_back(leaves[leaf_offset + i]);
    // }
    // int x = 3;
}

OctreeAccel::~OctreeAccel() { //FreeAligned(nodes2);
}

bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);

    Float tMin;
    if (!wb.IntersectP(ray, &tMin, &(ray.tMax))) return false;
    if (tMin < 0) tMin = 0;
    return RecurseIntersect(ray, isect, wb, 0, tMin);

    //Point3f point = ray.o;
//
    //if (!BoundsContainPoint(wb, ray.o)) {
    //    // Check intersection with outer planes of WorldBound and find closest point
    //    Float factor = ray.tMax;
    //    for (int i = 0; i < 6; i++) { // X=0,1; Y=2,3; Z=4,5; Even=Min; Odd=max
    //        int axis = i / 2; bool min = i % 2 == 0;
    //        if (ray.d[axis] == 0) continue; 
//
    //        // Determine factor for formula: point = origin + factor * direction; Then check if factor is legal and closer
    //        Float f = (((min) ? wb.pMin[axis] : wb.pMax[axis]) - ray.o[axis]) / ray.d[axis];
    //        if (f >= 0 && f <= factor && BoundsContainPoint(wb, ray.o + f * ray.d)) factor = f;
    //    }
//
    //    if (factor == ray.tMax) return false;
    //    point = ray.o + factor * ray.d;
    //}
//
    //if (RecurseIntersect(ray, isect, wb, 0, point)) {
    //    return true;
    //} else {
    //    return false;
    //}
}

bool OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect,
        Bounds3f bounds, uint32_t offset, Float tMin) const {
    if ((nodes[offset] & 1) == 0) { // Inner node
        std::vector<node_isect> children = FindChildrenTraverseOrder(ray, bounds);
        int child_node_idx = DetermineNodeIdx(bounds, ray.o + tMin * ray.d);
        for (int i = 0; i < children.size(); i++) {
            child_node_idx ^= children[i].flip_mask;
            if (RecurseIntersect(ray, isect,
                    octreeDivide(bounds, child_node_idx),
                    (nodes[offset] >> 1) + child_node_idx,
                    children[i].tMin))
                return true;
        }
    } else {
        //ray_leaves.push_back(bounds); //vis
        if (IntersectLeafPrims(ray, isect, bounds, offset)) { // Leaf node
            //if (idx++ < 10) visualizeRayTraversal(idx, ray_leaves); //vis
            //for (int i = 0; i < ray_leaves.size(); i++) ray_leaves.pop_back(); //vis
            return true;
        } else {
            return false;
        }
    }

    //ray_leaves = std::vector<Bounds3f>(); //vis
    return false;
}

bool OctreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction isect;
    return Intersect(ray, &isect);
}

std::shared_ptr<OctreeAccel> CreateOctreeAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeAccel>(std::move(prims));
}


// === VISUALIZATION ===
// visualize ray traversal
void OctreeAccel::visualizeRayTraversal(int index, std::vector<Bounds3f> bounds) const {
    std::string str = "ray_vis_" + std::to_string(index) + ".obj";
    const char *path = str.c_str();
    FILE *f = fopen(path, "wb");
    for (int j = 0; j < bounds.size(); j++) {
        uint vcnt = 1;

        for(uint i = 0; i < 8; i++)
        {
            Float x = ((i & 1) == 0) ? bounds[j].pMin.x : bounds[j].pMax.x;
            Float y = ((i & 2) == 0) ? bounds[j].pMin.y : bounds[j].pMax.y;
            Float z = ((i & 4) == 0) ? bounds[j].pMin.z : bounds[j].pMax.z;
            fprintf(f, "v %f %f %f\n", x, y, z);
        }

        fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 5, vcnt + 4);//bottom
        fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);//top
        fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 3, vcnt + 2);//front
        fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);//back
        fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 4, vcnt + 6, vcnt + 2);//left
        fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);//right
        vcnt += 8;
    }
    fclose(f);
}

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
