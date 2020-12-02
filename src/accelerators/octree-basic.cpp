
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
#include "accelerators/octree-basic.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>
#include <array>

namespace pbrt {

const int MAX_DEPTH = 15; //15
const int MAX_PRIMS = 32; //30

Bounds3f OctreeBasicAccel::octreeDivide(Bounds3f b, int idx) const {
    for (int i = 0; i < 3; i++) {
        Float axisHalf = (b.pMin[i] + b.pMax[i]) / 2;
        if ((idx & (1<<i)) == 0) b.pMax[i] = axisHalf;
        else b.pMin[i] = axisHalf;
    }
    return b;
}

OctreeBasicAccel::OctreeBasicAccel(){}

// KdTreeAccel Method Definitions
OctreeBasicAccel::OctreeBasicAccel(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {
    // Hier hast du die Bounding Box falsch initialisiert. Sonst ist immer der Ursprung enthalten.
    for (int i = 0; i < 3; i++) { wb.pMin[i] = FLT_MAX; wb.pMax[i] = -FLT_MAX; };

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
}

struct ChildHit { uint32_t idx; float tMin; };

// TODO Rekursion in Schleife umwandeln (schneller)
void OctreeBasicAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t offset, Bounds3f bounds, bool &hit) const {
    if ((nodes[offset] & 1) == 0) {
        // Innerer Knoten
        uint32_t traversal_size = 0;
        // Besser anstatt std::vector zu nehmen, da man sich die Speicherallokation spart (sehr teuer)
        std::array<ChildHit, 8> traversal; // It can happen that more than 4 nodes are intersected when using intersectP
        // 1. Schritt: Alle Kind Bounding Boxen Schneiden und t Parameter bestimmen
        for (uint32_t i = 0; i < 8; i++) {
            Bounds3f child_bounds = octreeDivide(bounds, i);
            ChildHit childHit = {i};
            float tMax;
            // TODO Optimierte Variante implementieren (3 Ebenentests)
            if (child_bounds.IntersectP(ray, &childHit.tMin, &tMax)) {
                //if (traversal_size == 8) abort(); // sanity check
                traversal[traversal_size++] = childHit; 
            }
        }
        // 2. Schritt: Sortierung der Kinder nach t Parameter
        // TODO eigene sortierung (bei 8 elementen ist ein naives insertionsort whr. besser)
        std::sort(traversal.begin(), traversal.begin() + traversal_size, [](const ChildHit &a, const ChildHit &b) {return a.tMin < b.tMin;});
        // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
        for (uint32_t i = 0; i < 8; i++) {
            // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
            // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
            if (i < traversal_size && traversal[i].tMin <= ray.tMax) {
                uint32_t child_offset = (nodes[offset] >> 1) + traversal[i].idx;
                RecurseIntersect(ray, isect, child_offset, octreeDivide(bounds, traversal[i].idx), hit);
            }
        }
    } else {
        // Leaf Knoten
        uint32_t prim_start = nodes[offset] >> 1;
        uint32_t prim_end = prim_start + sizes[offset];
        for (uint32_t i = prim_start; i < prim_end; i++)
            // TODO LRU-Cache/Mailboxing, damit man nicht mehrmals dasselbe primitiv testen muss
            if (leaves[i].get()->Intersect(ray, isect)) hit = true;
    }
}

bool OctreeBasicAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
    Float tMin, tMax;
    if (!wb.IntersectP(ray, &tMin, &tMax)) return false;
    RecurseIntersect(ray, isect, 0, wb, hit);
    return hit;
}

void OctreeBasicAccel::Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth) {        

    std::vector<std::shared_ptr<Primitive>> prims;
    for (int i = 0; i < primitives.size(); i++) {
        std::shared_ptr<Primitive> p = primitives.at(i);
        Bounds3f primBounds = p->WorldBound();
        bool bounds_overlap = true;
        for (int i = 0; i < 3; i++)
            if (primBounds.pMax[i] < bounds.pMin[i] || bounds.pMax[i] < primBounds.pMin[i]) {
                bounds_overlap = false;
                break;
            }
        if (bounds_overlap) prims.push_back(p);
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

bool OctreeBasicAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction isect;
    return Intersect(ray, &isect);
}

std::shared_ptr<OctreeBasicAccel> CreateOctreeBasicAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<OctreeBasicAccel>(std::move(prims));
}

OctreeBasicAccel::~OctreeBasicAccel() { //FreeAligned(nodes2);
}

// === VISUALIZATION ===

// Code to visualize octree
void OctreeBasicAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, int offset, Bounds3f bounds) {

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

    // Rekursion
    if ((nodes[offset] & 1) == 0) { // Inner node
        for (uint32_t i = 0; i < 8; i++) {
            lh_dump_rec(f, vcnt_, (nodes[offset] >> 1) + i, octreeDivide(bounds, i));
        }
    }
}

void OctreeBasicAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

}  // namespace pbrt
