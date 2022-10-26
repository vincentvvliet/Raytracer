#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>


// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
    float random = rand() / RAND_MAX;
  
    position = segmentLight.endpoint0 + ((segmentLight.endpoint1 - segmentLight.endpoint0) * random);
    color = segmentLight.color0 + ((segmentLight.color1 - segmentLight.color0) * random);

}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
    float random = rand() / RAND_MAX;
    position = parallelogramLight.v0 + (parallelogramLight.edge01 * random) + parallelogramLight.v0 + (parallelogramLight.edge02 * random);
    color = parallelogramLight.color0 * (parallelogramLight.edge01 * random) 
        + parallelogramLight.color1 * (parallelogramLight.edge01 * random) 
        + parallelogramLight.color2 * (parallelogramLight.edge02 * random) 
        + parallelogramLight.color3 * (parallelogramLight.edge02 * random);

}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 lightVector = glm::normalize(samplePos - intersection);
    glm::vec3 camVector = glm::normalize(intersection - ray.origin);
    float distance = glm::distance(intersection, samplePos);
    HitInfo shadowInfo;
    Ray shadowRay = Ray { intersection, lightVector, distance };   
    
    if (glm::dot(-lightVector, hitInfo.normal) * glm::dot(camVector, hitInfo.normal) >= 0.0f) {
        if (bvh.intersect(shadowRay, shadowInfo, features) && shadowRay.t < distance) {
            // Shadow ray intersect, therefore show no colour (return 0.0)
            drawRay(shadowRay, glm::vec3 { 0.0f, 0.0f, 1.0f });
            return 0.0f;
        } else {
            // No shadow ray intersect, therefore show colour as per usual (return 1.0)
            drawRay(shadowRay, glm::vec3 { 1.0f, 1.0f, 1.0f });
            return 1.0f;
        }
    } else {
        return 0.0f;
    }
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 total = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 sample = glm::vec3(1.0f, 1.0f, 1.0f);
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Perform your calculations for a point light.
                glm::vec3 intersection = ray.origin + ray.t * ray.direction;

                float shadowFactor = 1.0f;
                if (features.enableHardShadow) {
                    shadowFactor = testVisibilityLightSample(pointLight.position, bvh, features, ray, hitInfo);
                }

                total += shadowFactor * computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                // Perform your calculations for a segment light.
              /*  for (int i = 0; i < 100; i++) {
                     sampleSegmentLight(segmentLight, sample, sample);
                }*/
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                // Perform your calculations for a parallelogram light.
             /*   for (int i = 0; i < 100; i++) {
                    sampleParallelogramLight(parallelogramLight, sample, sample);
                }*/
            }
        }

        drawRay(ray, total);
        return total;
    }
    // If shading is disabled, return the albedo of the material.
    return hitInfo.material.kd;
}
