
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

namespace pbrt {

// === HELPERS ===

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
        if (ray.d[axis] == 0) continue;
        Float t = (b_h[axis] - ray.o[axis]) / ray.d[axis];
        if (t < tMin) continue;
        Point3f p = ray.o + t * ray.d;
        bool inside = true;
        for (int i = 0; i < 3; i++) if (p[i] < b.pMin[i] || p[i] > b.pMax[i]) { inside = false; break; }
        if (!inside) continue;
        // Add point to traversal array. It's index is currently only the to-be-flipped axis
        traversal[size++] = ChildHit{ 1<<axis, t };
    }
    // Sort list based on smallest tMin
    std::sort(traversal.begin() + 1, traversal.begin() + size, 
        [](const ChildHit &c1, const ChildHit &c2) {return c1.tMin < c2.tMin;});
    // Finally, determine idx for each child hit
    for (int i = 1; i < size; i++) traversal[i].idx ^= traversal[i-1].idx;
    return ChildTraversal{traversal, size};
}

int Rank(bftype bits, int n) {
    if (n == bfsize) return popcnt(bits);
    return popcnt(bits & ((bftone << n) - bftone));
}

bool BoundsContainPrim(Bounds3f b, std::shared_ptr<Primitive> p) {
    Bounds3f b_p = p->WorldBound();
    for (int i = 0; i < 3; i++)
        if (b_p.pMax[i] < b.pMin[i] || b.pMax[i] < b_p.pMin[i])
            return false;
    return true;
}

bool MakeLeafNode(Bounds3f b, std::vector<std::shared_ptr<Primitive>> prims) {
    // Check how many children contain at least 80% of the parents prims,
    // we call these 'cluster children'
    Vector3f b_h = BoundsHalf(b);
    std::vector<int> cluster_children;
    for (int idx = 0; idx < 8; idx++) {
        Bounds3f b_child = DivideBounds(b, idx, b_h);
        int num_prims = 0;
        for (int p = 0; p < prims.size(); p++)
            if (BoundsContainPrim(b_child, prims[p])) num_prims++;
        if (num_prims >= PRM_THRESH * prims.size()) cluster_children.push_back(idx);
    }
    if (cluster_children.size() < 2) return false;
    // Since there's at least 2 children with >80% of the parents prims, make a list
    // of the prims that are in ALL of these 'cluster children'. This is the primitive cluster
    std::vector<std::shared_ptr<Primitive>> cluster_prims;
    for (int p = 0; p < prims.size(); p++) {
        bool prim_in_all_cluster_children = true;
        for (int idx = 0; idx < cluster_children.size(); idx++) {
            Bounds3f b_child = DivideBounds(b, cluster_children[idx], b_h);
            if (!BoundsContainPrim(b_child, prims[p])) {
                prim_in_all_cluster_children = false;
                break;
            }
        }
        if (prim_in_all_cluster_children) cluster_prims.push_back(prims[p]);
    }
    // Now create bounds surrounding this cluster of primitives
    Bounds3f b_cluster;
    for (int i = 0; i < 3; i++) { b_cluster.pMin[i] = FLT_MAX; b_cluster.pMax[i] = -FLT_MAX; }
    for (int i = 0; i < cluster_prims.size(); i++) {
        Bounds3f b = cluster_prims[i]->WorldBound();
        for (int i = 0; i < 3; i++) {
            if (b.pMin[i] < b_cluster.pMin[i]) b_cluster.pMin[i] = b.pMin[i];
            if (b.pMax[i] > b_cluster.pMax[i]) b_cluster.pMax[i] = b.pMax[i];
        }
    }
    // Limit the cluster bounds to the current node's bounds
    for (int i = 0; i < 3; i++) {
        b_cluster.pMin[i] = (b_cluster.pMin[i] < b.pMin[i]) ? b.pMin[i] : b_cluster.pMin[i];
        b_cluster.pMax[i] = (b_cluster.pMax[i] > b.pMax[i]) ? b.pMax[i] : b_cluster.pMax[i];
    }
    // Now check the bounds of this cluster of prims, and if it takes up at least 80% of
    // the current bounds space, mark it as a leaf node
    return (b_cluster.Volume() >= VOL_THRESH * b.Volume());
}

// === OCTREE STRUCTURE CREATION ===

OctreeBasicAccel::OctreeBasicAccel(std::vector<std::shared_ptr<Primitive>> p, int max_prims, float prm_thresh, float vol_thresh)
        : primitives(std::move(p)) {
    MAX_PRIMS = max_prims;
    PRM_THRESH = prm_thresh;
    VOL_THRESH = vol_thresh;

    for (int i = 0; i < 3; i++) { wb.pMin[i] = FLT_MAX; wb.pMax[i] = -FLT_MAX; }

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
    printf("Octree: Starting creation!\n");
    Recurse(0, primitives, wb, 0);
    printf("Octree: Creation done!\n");

    printf("TEST: %i\n", max_prims);
    // printf("Starting visualization\n");
    // lh_dump("visualize_basic.obj");
    // printf("Visualization done!\n");
}

OctreeBasicAccel::OctreeBasicAccel(){}

void OctreeBasicAccel::Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth) {        

    // TODO primitives partitionieren in brauche ich/brauche ich nicht: std::partition
    std::vector<std::shared_ptr<Primitive>> prims;
    for (int i = 0; i < primitives.size(); i++) {
        std::shared_ptr<Primitive> p = primitives.at(i);
        if (BoundsContainPrim(bounds, p)) prims.push_back(p);
    }

    if (prims.size() > MAX_PRIMS && !MakeLeafNode(bounds, prims)) { // Inner node
        uint32_t offset_children = nodes.size();

        std::vector<uint32_t> nodes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<uint32_t> sizes_children = {0, 0, 0, 0, 0, 0, 0, 0};
        nodes.insert(nodes.end(), nodes_children.begin(), nodes_children.end());
        sizes.insert(sizes.end(), sizes_children.begin(), sizes_children.end());

        nodes[offset] = offset_children << 1 | 0;
        Vector3f b_h = BoundsHalf(bounds);
        for (uint32_t i = 0; i < 8; i++) Recurse(offset_children + i, prims, DivideBounds(bounds, i, b_h), depth + 1);
    } else { // Leaf node
        uint32_t offset_leaves = leaves.size();

        leaves.insert(leaves.end(), prims.begin(), prims.end());

        nodes[offset] = offset_leaves << 1 | 1;
        sizes[offset] = prims.size();
    }
}

std::shared_ptr<OctreeBasicAccel> CreateOctreeBasicAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    int max_prims = ps.FindOneInt("maxprims", 32);
    float prm_thresh = ps.FindOneFloat("prmthresh", 0.9);
    float vol_thresh = ps.FindOneFloat("volthresh", 0.9);
    return std::make_shared<OctreeBasicAccel>(std::move(prims), max_prims, prm_thresh, vol_thresh);
}

OctreeBasicAccel::~OctreeBasicAccel() { //FreeAligned(nodes2);
}

// === OCTREE RAY TRAVERSAL ===

// TODO Rekursion in Schleife umwandeln (schneller)
bool OctreeBasicAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return false;
    RecurseIntersect(ray, isect, 0, wb, tMin, hit);
    return hit;
}

void OctreeBasicAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t offset, Bounds3f bounds, Float tMin, bool &hit) const {
    if ((nodes[offset] & 1) == 0) {
        // Innerer Knoten
        ChildTraversal traversal = FindTraversalOrder(ray, bounds, tMin);
        // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
        Vector3f b_h = BoundsHalf(bounds);
        for (int i = 0; i < traversal.size; i++) {
            // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
            // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
            if (traversal.nodes[i].tMin <= ray.tMax) {
                uint32_t child_offset = (nodes[offset] >> 1) + traversal.nodes[i].idx;
                Bounds3f child_bounds = DivideBounds(bounds, traversal.nodes[i].idx, b_h);
                RecurseIntersect(ray, isect, child_offset, child_bounds, traversal.nodes[i].tMin, hit);
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

bool OctreeBasicAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    SurfaceInteraction isect;
    return Intersect(ray, &isect);
}

// === VISUALIZATION ===

// Code to visualize octree
void OctreeBasicAccel::lh_dump(const char *path) {
    FILE *f = fopen(path, "wb");
    uint32_t vcnt = 1;
    lh_dump_rec(f, &vcnt, 0, WorldBound());
    fclose(f);
}

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
        Vector3f b_h = BoundsHalf(bounds);
        for (uint32_t i = 0; i < 8; i++) {
            lh_dump_rec(f, vcnt_, (nodes[offset] >> 1) + i, DivideBounds(bounds, i, b_h));
        }
    }
}

}  // namespace pbrt
