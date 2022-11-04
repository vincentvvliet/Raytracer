#pragma once
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "shading.h"



 std::list<PointLight> sampleSegmentLight(const SegmentLight& segmentLight,std::list<PointLight> points, glm::vec3& position, glm::vec3& color);

std::list<PointLight> sampleParallelogramLight(const ParallelogramLight& parallelogramLight, std::list<PointLight> parallelogramPoints, glm::vec3& position, glm::vec3& color);

glm::vec3 testVisibilityLightSample(const glm::vec3& samplePos, glm::vec3 debugColor,const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo);

