
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
#include <array>
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
    
}

struct ChildHit { uint32_t idx; float tMin; };

// TODO Rekursion in Schleife umwandeln (schneller)
void OctreeAccel::RecurseIntersect(const Ray &ray, SurfaceInteraction *isect, uint32_t offset, Bounds3f bounds, bool &hit) const {
}

bool OctreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
}

void OctreeAccel::Recurse(int offset, std::vector<std::shared_ptr<Primitive>> primitives, Bounds3f bounds, int depth) {        
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

}  // namespace pbrt
