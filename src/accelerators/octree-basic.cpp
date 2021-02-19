
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

// Stats counted when building structure
STAT_MEMORY_COUNTER("Memory/Octree tree", octree_mem);
STAT_MEMORY_COUNTER("Memory/Octree topology", octree_mem_top);

STAT_COUNTER("Octree/Nodes - # Total", octree_stat_num_nodes);
STAT_COUNTER("Octree/Nodes - # Leaf", octree_stat_num_leafNodes);
STAT_COUNTER("Octree/Primitives - # Total", octree_stat_num_prims);
STAT_RATIO("Octree/Duplicate Primitives", octree_num_totalPrims, octree_num_uniquePrims);

// Stats counted when intersecting
#if defined (COUNT_STATS)
STAT_COUNTER("Octree - Intersects - Primitive/Total #", octree_stat_primIntersectsTotal);
STAT_INT_DISTRIBUTION("Octree - Intersects - Primitive/Distribution", octree_dist_primIntersects);

STAT_COUNTER("Octree - Intersects - Node - Leaf/Total #", octree_stat_leafNodeIntersectsTotal);
STAT_INT_DISTRIBUTION("Octree - Intersects - Node - Leaf/Distribution", octree_dist_leafNodeIntersects);

STAT_COUNTER("Octree - Intersects - Node - Total/Total #", octree_stat_nodeIntersectsTotal);
STAT_INT_DISTRIBUTION("Octree - Intersects - Node - Total/Distribution", octree_dist_nodeIntersects);

STAT_COUNTER("Octree - Intersects - Chunk/Total #", octree_stat_chunkIntersectsTotal);
STAT_INT_DISTRIBUTION("Octree - Intersects - Chunk/Distribution", octree_dist_chunkIntersects);
#endif

// === HELPERS ===

Vector3f BoundsHalf(Bounds3f &b) {
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

ChildTraversal FindTraversalOrder(const Ray &ray, Bounds3f &b, Float tMin, Vector3f &invDir) {
    uint32_t size = 1;
    std::array<ChildHit, 4> traversal;
    Vector3f b_h = BoundsHalf(b);
    // First child hit
    uint32_t idx = 0;
    Point3f init_point = ray.o + tMin * ray.d;
    for (int i = 0; i < 3; i++) if (init_point[i] > b_h[i]) idx |= (1<<i);
    traversal[0] = ChildHit{ idx, tMin };
    // Cut all bound-half-planes, and if intersection is within bounds, add to list
    for (int axis = 0; axis < 3; axis++) {
        if (ray.d[axis] == 0) continue;
        Float t = invDir[axis] * (b_h[axis] - ray.o[axis]);
        if (t < tMin) continue;
        Point3f p = ray.o + t * ray.d;
        bool inside = true;
        for (int i = 0; i < 3; i++) if (p[i] < b.pMin[i] || p[i] > b.pMax[i]) { inside = false; break; }
        if (!inside) continue;
        // Add point to traversal array. It's index is currently only the to-be-flipped axis
        // 1. Find index to insert
        int insertion_idx = size;
        for (int i = 1; i < 4; i++) {
            if (i == size) break;
            if (t < traversal[i].tMin) {
                insertion_idx = i;
                break;
            }
        }
        // 2. Move every item above that index up
        for (int i = 3; i > 1; i--) {
            if (i <= size) {
                if (i == insertion_idx) break;
                traversal[i] = traversal[i-1];
            }
        }
        // 3. Insert new childhit in correct place
        traversal[insertion_idx] = ChildHit{(uint32_t)1<<axis, t};
        size++;
    }
    // Finally, determine idx for each child hit
    for (int i = 1; i < size; i++) traversal[i].idx ^= traversal[i-1].idx;
    return ChildTraversal{traversal, size};
}

int Rank(bftype bits, int n) {
    if (n == bfsize) return popcnt(bits);
    return popcnt(bits & ((bftone << n) - bftone));
}

bool BoundsContainPrim(Bounds3f &b, std::shared_ptr<Primitive> p) {
    Bounds3f b_p = p->WorldBound();
    for (int i = 0; i < 3; i++)
        if (b_p.pMax[i] < b.pMin[i] || b.pMax[i] < b_p.pMin[i])
            return false;
    return true;
}

bool BoundsContainPoint(Bounds3f &b, Point3f &p) {
    for (int i = 0; i < 3; i++)
        if (p[i] < b.pMin[i] || p[i] > b.pMax[i]) return false;
    return true;
}

// bool OctreeBasicAccel::MakeLeafNode(Bounds3f &b, uint32_t prim_count) {
//     // Check how many children contain at least 80% of the parents prims,
//     // we call these 'cluster children'
//     Vector3f b_h = BoundsHalf(b);
//     std::array<bool, 8> cluster_children;
//     for (int i = 0; i < 8; i++) {
//         cluster_children[i] = false;
//         Bounds3f b_child = DivideBounds(b, i, b_h);
//         int num_prims = 0;
//         for (int p = 0; p < prim_count; p++)
//                 if (BoundsContainPrim(b_child, primitives[p])) num_prims++;
//         if (num_prims >= PRM_THRESH * prim_count) cluster_children[i] = true;
//     }
//     if (std::count(cluster_children.begin(), cluster_children.end(), true) < 2) return false;
//     // Since there's at least 2 children with >80% of the parents prims, make a list
//     // of the prims that are in ALL of these 'cluster children'. This is the primitive cluster
//     uint32_t num_cluster_prims = 0;
//     std::partition(primitives.begin(), primitives.begin() + prim_count,
//         [&b, &b_h, &num_cluster_prims, &cluster_children](std::shared_ptr<Primitive> p){
//             bool prim_in_all_cluster_children = true;
//             for (int i = 0; i < 8; i++) {
//                 if (!cluster_children[i]) continue;
//                 Bounds3f b_child = DivideBounds(b, i, b_h);
//                 if (!BoundsContainPrim(b_child, p)) {
//                     prim_in_all_cluster_children = false;
//                     break;
//                 }
//             }
//             if (prim_in_all_cluster_children) {
//                 num_cluster_prims++;
//                 return true;
//             }
//             return false;
//         });
//     // Now create bounds surrounding this cluster of primitives
//     Bounds3f b_cluster;
//     for (int i = 0; i < 3; i++) { b_cluster.pMin[i] = b.pMax[i]; b_cluster.pMax[i] = b.pMin[i]; }
//     for (int i = 0; i < num_cluster_prims; i++) {
//         Bounds3f b_prim = primitives[i]->WorldBound();
//         for (int i = 0; i < 3; i++) {
//             if (b_prim.pMin[i] < b_cluster.pMin[i]) b_cluster.pMin[i] = b_prim.pMin[i];
//             if (b_prim.pMax[i] > b_cluster.pMax[i]) b_cluster.pMax[i] = b_prim.pMax[i];
//         }
//     }
//     // Limit the cluster bounds to the current node's bounds
//     for (int i = 0; i < 3; i++) {
//         b_cluster.pMin[i] = (b_cluster.pMin[i] < b.pMin[i]) ? b.pMin[i] : b_cluster.pMin[i];
//         b_cluster.pMax[i] = (b_cluster.pMax[i] > b.pMax[i]) ? b.pMax[i] : b_cluster.pMax[i];
//     }
//     // Now check the bounds of this cluster of prims, and if it takes up at least 80% of
//     // the current bounds space, mark it as a leaf node
//     return (b_cluster.Volume() >= VOL_THRESH * b.Volume());
// }

// === OCTREE STRUCTURE CREATION ===

OctreeBasicAccel::OctreeBasicAccel(std::vector<std::shared_ptr<Primitive>> p)
        : primitives(std::move(p)) {
    printf("Chosen Accelerator: Octree\n");

    std::vector<int> p_idxs;
    for (int i = 0; i < primitives.size(); i++) {
        p_idxs.push_back(i);
    }

    for (int i = 0; i < 3; i++) { wb.pMin[i] = FLT_MAX; wb.pMax[i] = -FLT_MAX; }

    // Determine world bounds
    for (int i = 0; i < primitives.size(); i++) {
        Bounds3f b = primitives[i]->WorldBound();
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

    struct OctreeBuildNode { uint32_t offset; uint32_t parent_prim_count; Bounds3f bounds; };
    OctreeBuildNode node_stack[128];
    OctreeBuildNode current_node = OctreeBuildNode{0, (uint32_t)primitives.size(), wb};
    nodes.push_back(0);
    sizes.push_back(0);
    uint32_t stack_size = 0;
    while (true) {
        octree_stat_num_nodes++;

        uint32_t prim_count = 0;
        std::partition(p_idxs.begin(), p_idxs.begin() + current_node.parent_prim_count,
                [this, &current_node, &prim_count](int i){
                    if (BoundsContainPrim(current_node.bounds, primitives[i])) { prim_count++; return true; }
                    return false; });

        // Determine whether or not this node should be an inner node or a leaf node
        uint32_t children_prim_count = 0;
        Vector3f b_h = BoundsHalf(current_node.bounds);
        for (int i = 0; i < prim_count; i++) {
            Bounds3f b_prim = primitives[i]->WorldBound();
            int prim_mult = 1; // how many times this prim gets multiplied
            for (int j = 0; j < 3; j++) if (b_prim.pMin[j] <= b_h[j] && b_prim.pMax[j] >= b_h[j]) prim_mult *= 2;
            children_prim_count += prim_mult;
        }
        bool make_leaf_node = children_prim_count > MULT_THRESH * prim_count;

        if (prim_count > MAX_PRIMS && !make_leaf_node) { // Inner node
            uint32_t offset_children = nodes.size();

            std::vector<uint32_t> nodes_children = {0, 0, 0, 0, 0, 0, 0, 0};
            std::vector<uint32_t> sizes_children = {0, 0, 0, 0, 0, 0, 0, 0};
            nodes.insert(nodes.end(), nodes_children.begin(), nodes_children.end());
            sizes.insert(sizes.end(), sizes_children.begin(), sizes_children.end());

            nodes[current_node.offset] = offset_children << 1 | 0;
            for (uint32_t i = 0; i < 7; i++)
                node_stack[stack_size++] = OctreeBuildNode{offset_children + i, prim_count, DivideBounds(current_node.bounds, i, b_h)};
            current_node = OctreeBuildNode{offset_children + 7, prim_count, DivideBounds(current_node.bounds, 7, b_h)};
        } else { // Leaf node
            octree_stat_num_leafNodes++;
            octree_stat_num_prims += prim_count;
            uint32_t offset_leaves = leaves.size();

            leaves.insert(leaves.end(), p_idxs.begin(), p_idxs.begin() + prim_count);
            nodes[current_node.offset] = offset_leaves << 1 | 1;
            sizes[current_node.offset] = prim_count;

            if (stack_size == 0) break;
            current_node = node_stack[--stack_size];
        }
    }

    octree_mem += sizeof(*this) - sizeof(primitives) +
            nodes.size() * sizeof(nodes[0]) +
            sizes.size() * sizeof(sizes[0]) +
            leaves.size() * sizeof(leaves[0]);
    octree_mem_top += sizeof(nodes) + nodes.size() * sizeof(nodes[0]);

    octree_num_totalPrims = octree_stat_num_prims;
    octree_num_uniquePrims = primitives.size();

    // lh_dump("visualize_basic.obj");
}

OctreeBasicAccel::OctreeBasicAccel(){}

// void OctreeBasicAccel::Recurse(int offset, uint32_t parent_prim_count, Bounds3f bounds, int depth) {        

//     uint32_t prim_count = 0;
//     std::partition(primitives.begin(), primitives.begin() + parent_prim_count,
//             [&bounds, &prim_count](std::shared_ptr<Primitive> p){
//                 if (BoundsContainPrim(bounds, p)) {
//                     prim_count++;
//                     return true;
//                 }
//                 return false;
//             });

//     octree_stat_num_nodes++;

//     if (prim_count > MAX_PRIMS && !MakeLeafNode(bounds, prim_count)) { // Inner node
//         uint32_t offset_children = nodes.size();

//         std::vector<uint32_t> nodes_children = {0, 0, 0, 0, 0, 0, 0, 0};
//         std::vector<uint32_t> sizes_children = {0, 0, 0, 0, 0, 0, 0, 0};
//         nodes.insert(nodes.end(), nodes_children.begin(), nodes_children.end());
//         sizes.insert(sizes.end(), sizes_children.begin(), sizes_children.end());

//         nodes[offset] = offset_children << 1 | 0;
//         Vector3f b_h = BoundsHalf(bounds);
//         for (uint32_t i = 0; i < 8; i++) Recurse(offset_children + i, prim_count, DivideBounds(bounds, i, b_h), depth + 1);
//     } else { // Leaf node
//         octree_stat_num_leafNodes++;
//         octree_stat_num_prims += prim_count;
//         uint32_t offset_leaves = leaves.size();

//         leaves.insert(leaves.end(), primitives.begin(), primitives.begin() + prim_count);
//         nodes[offset] = offset_leaves << 1 | 1;
//         sizes[offset] = prim_count;
//     }
// }

std::shared_ptr<OctreeBasicAccel> CreateOctreeBasicAccelerator(std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    int max_prims = ps.FindOneInt("maxprims", 32);
    float mult_thresh = ps.FindOneFloat("multthresh", 4);
    MAX_PRIMS = max_prims;
    MULT_THRESH = mult_thresh;
    return std::make_shared<OctreeBasicAccel>(std::move(prims));
}

OctreeBasicAccel::~OctreeBasicAccel() { //FreeAligned(nodes2);
}

// === OCTREE RAY TRAVERSAL ===
static thread_local int mailbox_gen = 0;
static thread_local std::vector<int> mailbox;

// === OCTREE RAY TRAVERSAL ===
// TODO Rekursion in Schleife umwandeln (schneller)
bool OctreeBasicAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    mailbox_gen++;
    if (mailbox.size() == 0) {
        mailbox = std::vector<int>(primitives.size(), 0);
    }

    ProfilePhase p(Prof::AccelIntersect);
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return false;
    bool hit = false;
    // === COUNT STATS INITIALIZE BEGIN ===
    #if defined (COUNT_STATS)
    int num_prims_intersected = 0;
    int num_nodes_intersected = 0;
    int num_leafNodes_intersected = 0;
    #endif
    // === COUNT STATS INITIALIZE END ===
    Vector3f invDir = {1/ray.d[0], 1/ray.d[1], 1/ray.d[2]};
    struct OctreeTraversalNode { uint32_t offset; uint8_t bounds_stack_offset; Float tMin; };
    struct OctreeTraversalBounds { Bounds3f bounds; Vector3f b_h; };
    OctreeTraversalNode node_stack[128];
    OctreeTraversalBounds bounds_stack[64];
    int node_stack_offset = 0;
    int bounds_stack_offset = -1;
    OctreeTraversalNode current_node = OctreeTraversalNode{0, 0, tMin};
    Bounds3f current_bounds = wb;
    while (true) {
        if (bounds_stack_offset >= 0) {
            current_bounds = 
                    DivideBounds(bounds_stack[bounds_stack_offset].bounds,
                                 (current_node.offset - 1) & 7,
                                 bounds_stack[bounds_stack_offset].b_h);
        }
        // === COUNT STATS FOR ALL NODES BEGIN ===
        #if defined (COUNT_STATS)
        num_nodes_intersected++;
        #endif
        // === COUNT STATS FOR ALL NODES END ===
        if ((nodes[current_node.offset] & 1) == 0) {
            // Innerer Knoten
            ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, current_node.tMin, invDir);
            // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
            Vector3f b_h = BoundsHalf(current_bounds);
            bounds_stack[++bounds_stack_offset] = OctreeTraversalBounds{current_bounds, b_h};
            ChildHit *traversal_node;
            for (int i = 3; i > 0; i--) {
                // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
                // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
                traversal_node = &traversal.nodes[i];
                if (i < traversal.size && traversal_node->tMin <= ray.tMax) {
                    uint32_t child_offset = (nodes[current_node.offset] >> 1) + traversal_node->idx;
                    node_stack[node_stack_offset++] = OctreeTraversalNode{child_offset, (uint8_t)bounds_stack_offset, traversal_node->tMin};
                }
            }
            traversal_node = &traversal.nodes[0];
            uint32_t child_offset = (nodes[current_node.offset] >> 1) + traversal_node->idx;
            current_node = OctreeTraversalNode{child_offset, (uint8_t)bounds_stack_offset, traversal_node->tMin};
        } else {
            // Leaf Knoten
            uint32_t prim_start = nodes[current_node.offset] >> 1;
            uint32_t prim_end = prim_start + sizes[current_node.offset];
            for (uint32_t i = prim_start; i < prim_end; i++) {
                // TODO LRU-Cache/Mailboxing, damit man nicht mehrmals dasselbe primitiv testen muss
                int prim_idx = leaves[i];
                if (mailbox[prim_idx] != mailbox_gen) {
                    mailbox[prim_idx] = mailbox_gen;
                    if (primitives[prim_idx].get()->Intersect(ray, isect)) hit = true;
                }
            }
            // === COUNT STATS FOR LEAF NODES BEGIN ===
            #if defined (COUNT_STATS)
            num_prims_intersected += prim_end - prim_start;
            num_leafNodes_intersected++;
            #endif
            // === COUNT STATS FOR LEAF NODES END ===
            if (hit && BoundsContainPoint(current_bounds, isect->p)) break;
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
            bounds_stack_offset = current_node.bounds_stack_offset;
        }
    }
    // === COUNT STATS FOR RAY AND TOTAL STATS BEGIN ===
    #if defined (COUNT_STATS)
    octree_stat_primIntersectsTotal += num_prims_intersected;
    ReportValue(octree_dist_primIntersects, num_prims_intersected);
    // prim_isects_per_ray->push_back(num_prims_intersected);

    octree_stat_nodeIntersectsTotal += num_nodes_intersected;
    ReportValue(octree_dist_nodeIntersects, num_nodes_intersected);
    // node_isects_per_ray->push_back(num_nodes_intersected);

    octree_stat_leafNodeIntersectsTotal += num_leafNodes_intersected;
    ReportValue(octree_dist_leafNodeIntersects, num_leafNodes_intersected);
    // leafNode_isects_per_ray->push_back(num_leafNodes_intersected);
    #endif
    // === COUNT STATS FOR RAY AND TOTAL STATS END ===
    return hit;
}

float OctreeBasicAccel::IntersectMetric(const Ray &ray, metric m) const {
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return 0;
    bool hit = false;
    float metric_cnt = 0;
    SurfaceInteraction _isect;
    SurfaceInteraction *isect = &_isect;
    Vector3f invDir = {1/ray.d[0], 1/ray.d[1], 1/ray.d[2]};
    // boundskey = 29 bits - index in bounds_stack to this nodes parent bounds; 3 bits - node idx [0-7]
    struct OctreeTraversalNode { uint32_t offset; uint8_t bounds_stack_offset; Float tMin; };
    struct OctreeTraversalBounds { Bounds3f bounds; Vector3f b_h; };
    OctreeTraversalNode node_stack[128];
    OctreeTraversalBounds bounds_stack[64];
    int node_stack_offset = 0;
    int bounds_stack_offset = -1;
    OctreeTraversalNode current_node = OctreeTraversalNode{0, 0, tMin};
    Bounds3f current_bounds = wb;
    while (true) {
        if (bounds_stack_offset >= 0) {
            current_bounds = 
                    DivideBounds(bounds_stack[bounds_stack_offset].bounds,
                                 (current_node.offset - 1) & 7,
                                 bounds_stack[bounds_stack_offset].b_h);
        }
        if (m == metric::NODES) metric_cnt++;
        if ((nodes[current_node.offset] & 1) == 0) {
            // Innerer Knoten
            ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, current_node.tMin, invDir);
            // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
            Vector3f b_h = BoundsHalf(current_bounds);
            bounds_stack[++bounds_stack_offset] = OctreeTraversalBounds{current_bounds, b_h};
            ChildHit *traversal_node;
            for (int i = 3; i > 0; i--) {
                // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
                // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
                traversal_node = &traversal.nodes[i];
                if (i < traversal.size && traversal_node->tMin <= ray.tMax) {
                    uint32_t child_offset = (nodes[current_node.offset] >> 1) + traversal_node->idx;
                    node_stack[node_stack_offset++] = OctreeTraversalNode{child_offset, (uint8_t)bounds_stack_offset, traversal_node->tMin};
                }
            }
            traversal_node = &traversal.nodes[0];
            uint32_t child_offset = (nodes[current_node.offset] >> 1) + traversal_node->idx;
            current_node = OctreeTraversalNode{child_offset, (uint8_t)bounds_stack_offset, traversal_node->tMin};
        } else {
            if (m == metric::LEAFNODES) metric_cnt++;
            // Leaf Knoten
            uint32_t prim_start = nodes[current_node.offset] >> 1;
            uint32_t prim_end = prim_start + sizes[current_node.offset];
            for (uint32_t i = prim_start; i < prim_end; i++) {
                // TODO LRU-Cache/Mailboxing, damit man nicht mehrmals dasselbe primitiv testen muss
                int prim_idx = leaves[i];
                if (mailbox[prim_idx] != mailbox_gen) {
                    mailbox[prim_idx] = mailbox_gen;
                    if (primitives[prim_idx].get()->Intersect(ray, isect)) hit = true;
                }
            }
            if (m == metric::PRIMITIVES) metric_cnt += prim_end - prim_start;
            if (hit && BoundsContainPoint(current_bounds, isect->p)) break;
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
            bounds_stack_offset = current_node.bounds_stack_offset;
        }
    }
    return metric_cnt;
}

bool OctreeBasicAccel::IntersectP(const Ray &ray) const {
    mailbox_gen++;
    if (mailbox.size() == 0) {
        mailbox = std::vector<int>(primitives.size(), 0);
    }

    ProfilePhase p(Prof::AccelIntersectP);
    Float tMin, _;
    if (!wb.IntersectP(ray, &tMin, &_)) return false;
    Vector3f invDir = {1/ray.d[0], 1/ray.d[1], 1/ray.d[2]};
    // boundskey = 29 bits - index in bounds_stack to this nodes parent bounds; 3 bits - node idx [0-7]
    struct OctreeTraversalNode { uint32_t offset; uint8_t bounds_stack_offset; Float tMin; };
    struct OctreeTraversalBounds { Bounds3f bounds; Vector3f b_h; };
    OctreeTraversalNode node_stack[128];
    OctreeTraversalBounds bounds_stack[64];
    int node_stack_offset = 0;
    int bounds_stack_offset = -1;
    OctreeTraversalNode current_node = OctreeTraversalNode{0, 0, tMin};
    Bounds3f current_bounds = wb;
    while (true) {
        if (bounds_stack_offset >= 0) {
            current_bounds = 
                    DivideBounds(bounds_stack[bounds_stack_offset].bounds,
                                 (current_node.offset - 1) & 7,
                                 bounds_stack[bounds_stack_offset].b_h);
        }
        if ((nodes[current_node.offset] & 1) == 0) {
            // Innerer Knoten
            ChildTraversal traversal = FindTraversalOrder(ray, current_bounds, current_node.tMin, invDir);
            // 3. Schritt: Traversierung der Kindknoten in sortierter Reihenfolge
            Vector3f b_h = BoundsHalf(current_bounds);
            bounds_stack[++bounds_stack_offset] = OctreeTraversalBounds{current_bounds, b_h};
            ChildHit *traversal_node;
            for (int i = 3; i > 0; i--) {
                // Kindknoten werden dann nicht mehr traversiert, wenn bereits ein näherer Schnitt ermittelt wurde
                // Dadurch deckt man auch den Fall ab, dass zwar ein Schnitt gefunden wurde, dieser aber außerhalb der Knotens liegt
                traversal_node = &traversal.nodes[i];
                if (i < traversal.size && traversal_node->tMin <= ray.tMax) {
                    uint32_t child_offset = (nodes[current_node.offset] >> 1) + traversal_node->idx;
                    node_stack[node_stack_offset++] = OctreeTraversalNode{child_offset, (uint8_t)bounds_stack_offset, traversal_node->tMin};
                }
            }
            traversal_node = &traversal.nodes[0];
            uint32_t child_offset = (nodes[current_node.offset] >> 1) + traversal_node->idx;
            current_node = OctreeTraversalNode{child_offset, (uint8_t)bounds_stack_offset, traversal_node->tMin};
        } else {
            // Leaf Knoten
            uint32_t prim_start = nodes[current_node.offset] >> 1;
            uint32_t prim_end = prim_start + sizes[current_node.offset];
            for (uint32_t i = prim_start; i < prim_end; i++) {
                int prim_idx = leaves[i];
                if (mailbox[prim_idx] != mailbox_gen) {
                    mailbox[prim_idx] = mailbox_gen;
                    if (primitives[prim_idx].get()->IntersectP(ray)) return true;
                }
            }
            if (node_stack_offset == 0) break;
            current_node = node_stack[--node_stack_offset];
            bounds_stack_offset = current_node.bounds_stack_offset;
        }
    }
    return false;
}

// === VISUALIZATION ===

// // Code to visualize octree
// void OctreeBasicAccel::lh_dump(const char *path) {
//     FILE *f = fopen(path, "wb");
//     uint32_t vcnt = 1;
//     lh_dump_rec(f, &vcnt, 0, WorldBound());
//     fclose(f);
// }

// void OctreeBasicAccel::lh_dump_rec(FILE *f, uint32_t *vcnt_, int offset, Bounds3f bounds) {

//     // Vertices ausgeben
//     for(uint32_t i = 0; i < 8; i++)
//     {
//         Float x = ((i & 1) == 0) ? bounds.pMin.x : bounds.pMax.x;
//         Float y = ((i & 2) == 0) ? bounds.pMin.y : bounds.pMax.y;
//         Float z = ((i & 4) == 0) ? bounds.pMin.z : bounds.pMax.z;
//         fprintf(f, "v %f %f %f\n", x, y, z);
//     }

//     // Vertex indices ausgeben
//     uint32_t vcnt = *vcnt_;
//     fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 5, vcnt + 4);//bottom
//     fprintf(f, "f %d %d %d %d\n", vcnt + 2, vcnt + 3, vcnt + 7, vcnt + 6);//top
//     fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 1, vcnt + 3, vcnt + 2);//front
//     fprintf(f, "f %d %d %d %d\n", vcnt + 4, vcnt + 5, vcnt + 7, vcnt + 6);//back
//     fprintf(f, "f %d %d %d %d\n", vcnt    , vcnt + 4, vcnt + 6, vcnt + 2);//left
//     fprintf(f, "f %d %d %d %d\n", vcnt + 1, vcnt + 5, vcnt + 7, vcnt + 3);//right
//     *vcnt_ += 8;

//     // Rekursion
//     if ((nodes[offset] & 1) == 0) { // Inner node
//         Vector3f b_h = BoundsHalf(bounds);
//         for (uint32_t i = 0; i < 8; i++) {
//             lh_dump_rec(f, vcnt_, (nodes[offset] >> 1) + i, DivideBounds(bounds, i, b_h));
//         }
//     }
// }

}  // namespace pbrt
