
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_EMBREE_H
#define PBRT_ACCELERATORS_EMBREE_H

// accelerators/bvh.h*
#include <embree3/rtcore.h>

#include <atomic>

#include "pbrt.h"
#include "primitive.h"

namespace pbrt {

typedef std::vector<std::shared_ptr<Primitive>> Primitives;

// BVHAccel Declarations
class EmbreeAccel : public Aggregate {
  public:
    static void BoundsCallback(const struct RTCBoundsFunctionArguments *args);

    // BVHAccel Public Methods
    EmbreeAccel(std::vector<std::shared_ptr<Primitive>> p);
    Bounds3f WorldBound() const;
    ~EmbreeAccel();
    static void IntersectCallback(
        const struct RTCIntersectFunctionNArguments *args);
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;

    static void OccludedCallback(
        const struct RTCOccludedFunctionNArguments *args);
    bool IntersectP(const Ray &ray) const;

  private:
    RTCDevice device;
    RTCScene scene;
    Primitives primitives;
};

std::shared_ptr<EmbreeAccel> CreateEmbreeAccelerator(Primitives prims,
                                                     const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_EMBREE_H
