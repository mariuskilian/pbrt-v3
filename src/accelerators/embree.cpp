
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

#include "accelerators/embree.h"

#include <algorithm>

#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "custom_params.h"

namespace pbrt {

void EmbreeAccel::BoundsCallback(
    const struct RTCBoundsFunctionArguments *args) {
    const Primitives *primitives = (Primitives *)args->geometryUserPtr;
    Bounds3f bounds = (*primitives)[args->primID]->WorldBound();
    args->bounds_o->lower_x = bounds.pMin.x;
    args->bounds_o->lower_y = bounds.pMin.y;
    args->bounds_o->lower_z = bounds.pMin.z;
    args->bounds_o->upper_x = bounds.pMax.x;
    args->bounds_o->upper_y = bounds.pMax.y;
    args->bounds_o->upper_z = bounds.pMax.z;
}

// BVHAccel Method Definitions
EmbreeAccel::EmbreeAccel(std::vector<std::shared_ptr<Primitive>> p)
    : primitives(std::move(p)) {
    ProfilePhase _(Prof::AccelConstruction);
    if (primitives.empty()) return;
    device = rtcNewDevice("");
    RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    rtcSetGeometryBoundsFunction(geometry, EmbreeAccel::BoundsCallback,
                                 (void *)this);
    rtcSetGeometryIntersectFunction(geometry, EmbreeAccel::IntersectCallback);
    rtcSetGeometryOccludedFunction(geometry, EmbreeAccel::OccludedCallback);
    rtcSetGeometryUserData(geometry, (void *)&primitives);
    rtcSetGeometryUserPrimitiveCount(geometry, primitives.size());
    rtcCommitGeometry(geometry);
    scene = rtcNewScene(device);
    rtcAttachGeometry(scene, geometry);
    rtcReleaseGeometry(geometry);
    rtcCommitScene(scene);
}

Bounds3f EmbreeAccel::WorldBound() const {
    RTCBounds bounds;
    rtcGetSceneBounds(scene, &bounds);
    return Bounds3f(Point3f(bounds.lower_x, bounds.lower_y, bounds.lower_z),
                    Point3f(bounds.upper_x, bounds.upper_y, bounds.upper_z));
}

EmbreeAccel::~EmbreeAccel() {}

struct IntersectContext {
    RTCIntersectContext rtc;
    bool hit;
    const Ray *ray;
    SurfaceInteraction *isect;
    float metric_cnt;
};

void EmbreeAccel::IntersectCallback(
    const struct RTCIntersectFunctionNArguments *args) {
    // assumes N = 1 and valid ray
    assert(args->N == 1 && args->valid[0]);
    const Primitives *primitives = (Primitives *)args->geometryUserPtr;
    IntersectContext *context = (IntersectContext *)args->context;
    if ((*primitives)[args->primID]->Intersect(*context->ray, context->isect)) {
        context->hit = true;
        RTCRayN *ray = RTCRayHitN_RayN(args->rayhit, args->N);
        RTCRayN_tfar(ray, args->N, 0) = context->ray->tMax;
    }
    context->metric_cnt++;
}

bool EmbreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);

    IntersectContext context{.hit = false, .ray = &ray, .isect = isect, .metric_cnt = 0.0};
    rtcInitIntersectContext(&context.rtc);
    RTCRayHit rayhit{.ray =
                         {
                             .org_x = ray.o.x,
                             .org_y = ray.o.y,
                             .org_z = ray.o.z,
                             .tnear = 0,
                             .dir_x = ray.d.x,
                             .dir_y = ray.d.y,
                             .dir_z = ray.d.z,
                             .time = 0,
                             .tfar = ray.tMax,
                             .mask = 0,
                             .id = 0,
                             .flags = 0,
                         },
                     .hit = {.geomID = RTC_INVALID_GEOMETRY_ID,
                             .instID = {RTC_INVALID_GEOMETRY_ID}}};
    rtcIntersect1(scene, (RTCIntersectContext *)&context, &rayhit);
    return context.hit;
}

float EmbreeAccel::IntersectMetric(const Ray &ray, metric m) const {
    ProfilePhase p(Prof::AccelIntersect);

    SurfaceInteraction _isect;
    SurfaceInteraction* isect = &_isect;

    IntersectContext context{.hit = false, .ray = &ray, .isect = isect, .metric_cnt = 0.0};
    rtcInitIntersectContext(&context.rtc);
    RTCRayHit rayhit{.ray =
                         {
                             .org_x = ray.o.x,
                             .org_y = ray.o.y,
                             .org_z = ray.o.z,
                             .tnear = 0,
                             .dir_x = ray.d.x,
                             .dir_y = ray.d.y,
                             .dir_z = ray.d.z,
                             .time = 0,
                             .tfar = ray.tMax,
                             .mask = 0,
                             .id = 0,
                             .flags = 0,
                         },
                     .hit = {.geomID = RTC_INVALID_GEOMETRY_ID,
                             .instID = {RTC_INVALID_GEOMETRY_ID}}};
    rtcIntersect1(scene, (RTCIntersectContext *)&context, &rayhit);
    return context.metric_cnt;
}

struct OccludedContext {
    RTCIntersectContext rtc;
    const Ray *ray;
};

void EmbreeAccel::OccludedCallback(
    const struct RTCOccludedFunctionNArguments *args) {
    // assumes N = 1 and valid ray
    assert(args->N == 1 && args->valid[0]);
    const Primitives *primitives = (Primitives *)args->geometryUserPtr;
    OccludedContext *context = (OccludedContext *)args->context;
    if ((*primitives)[args->primID]->IntersectP(*context->ray)) {
        RTCRayN_tfar(args->ray, args->N, 0) = -INFINITY;
    }
}

bool EmbreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);

    OccludedContext context{.ray = &ray};
    rtcInitIntersectContext(&context.rtc);
    RTCRay ray_ = {
        .org_x = ray.o.x,
        .org_y = ray.o.y,
        .org_z = ray.o.z,
        .tnear = 0,
        .dir_x = ray.d.x,
        .dir_y = ray.d.y,
        .dir_z = ray.d.z,
        .time = 0,
        .tfar = ray.tMax,
        .mask = 0,
        .id = 0,
        .flags = 0,
    };
    rtcOccluded1(scene, (RTCIntersectContext *)&context, &ray_);
    return ray_.tfar == -INFINITY;
}

std::shared_ptr<EmbreeAccel> CreateEmbreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    return std::make_shared<EmbreeAccel>(std::move(prims));
}

}  // namespace pbrt

