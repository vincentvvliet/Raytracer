#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <tuple>

int segmentLightPoints = 10;
int parallelogramLightPoints = 100;

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
std::list<PointLight> interpolateLine(const SegmentLight& light, std::list<PointLight> points, int segmentLightPoints) {
    for (int i = 0; i < segmentLightPoints; i++) {
        
        float random = rand() / RAND_MAX;
        glm::vec3 point = light.endpoint0 + ((light.endpoint1 - light.endpoint0) * random);
        glm::vec3 colour = light.color0 + ((light.color1 - light.color0) * random);
        points.insert(points.begin(), PointLight(point, colour));
        
    }
    return points;
}
    std::list<PointLight> sampleSegmentLight(const SegmentLight& segmentLight, std::list<PointLight> points,glm::vec3& position, glm::vec3& color)
{
    
   
    return interpolateLine(segmentLight,points,segmentLightPoints);
}


std::list<PointLight> interpolate(const ParallelogramLight& light, std::list<PointLight> points, int sqrtSampleCount, glm::vec3 v1, glm::vec3 v2, float ratio_x, float ratio_y, float dist_x, float dist_y)
{
    // Uses bilinear interpolation for samples
    for (int i = 0; i <= sqrtSampleCount; i++) {
        float x1 = (light.v0.x - v1.x) * (ratio_x * i) + v1.x;
        float y1 = (light.v0.y - v1.y) * (ratio_x * i) + v1.y;
        float z1 = (light.v0.z - v1.z) * (ratio_x * i) + v1.z;
        glm::vec3 point1 = glm::vec3 { x1, y1, z1 };
        float alpha = (float)glm::distance(light.v0, point1) / dist_x;
        for (int j = 0; j <= sqrtSampleCount; j++) {
            glm::vec3 color1 = (1 - alpha) * light.color0 + (alpha)*light.color1;
            glm::vec3 color2 = (1 - alpha) * light.color2 + (alpha)*light.color3;
            float x2 = (light.v0.x - v2.x) * (ratio_y * j) + v2.x;
            float y2 = (light.v0.y - v2.y) * (ratio_y * j) + v2.y;
            float z2 = (light.v0.z - v2.z) * (ratio_y * j) + v2.z;
            glm::vec3 point2 = glm::vec3 { x2, y2, z2 };
            float beta = (float)glm::distance(light.v0, point2) / dist_y;
            glm::vec3 finalColor = (1 - beta) * color2 + (beta) * color1;
            float final_x = -(light.v0.x - v2.x) * (ratio_y * j) + point1.x;
            float final_y = -(light.v0.y - v2.y) * (ratio_y * j) + point1.y;
            float final_z = -(light.v0.z - v2.z) * (ratio_y * j) + point1.z;

            points.insert(points.begin(), PointLight(glm::vec3(final_x, final_y, final_z), finalColor));
        }
    }

    return points;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
std::list<PointLight> sampleParallelogramLight(const ParallelogramLight& parallelogramLight, std::list<PointLight> parallelogramPoints, glm::vec3& position, glm::vec3& color)
{
    glm::vec3 v1 = parallelogramLight.v0 + parallelogramLight.edge01;
    glm::vec3 v2 = parallelogramLight.v0 + parallelogramLight.edge02;

    int sqrtSampleCount = int(std::sqrt(parallelogramLightPoints));

    float dist_x = glm::distance(parallelogramLight.v0, v1);
    float dist_y = glm::distance(parallelogramLight.v0, v2);
    float normalizedDist_x = dist_x / sqrtSampleCount;
    float normalizedDist_y = dist_y / sqrtSampleCount;

    // Percentage in x and y directions
    float ratio_x = normalizedDist_x / dist_x;
    float ratio_y = normalizedDist_y / dist_y;

    return interpolate(parallelogramLight, parallelogramPoints, sqrtSampleCount, v1, v2, ratio_x, ratio_y, dist_x, dist_y);

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
        glm::vec3 intersection = ray.origin + ray.t * ray.direction;
        std::list<PointLight> segmentPoints; 
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Perform your calculations for a point light.
                float shadowFactor = 1.0f;
                float transparency = 1.0f;
               
                if (features.enableHardShadow) {
                    shadowFactor = testVisibilityLightSample(pointLight.position, bvh, features, ray, hitInfo);
                    
                }
                total += shadowFactor * computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);              
           
            } else if (std::holds_alternative<SegmentLight>(light) && features.enableSoftShadow) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                std::list<PointLight> linepoints;
                
                linepoints = sampleSegmentLight(segmentLight,linepoints,sample,sample);
                glm::vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);
                
                // interate over the sampled points and add all the colors together
                for (std::list<PointLight>::iterator it = linepoints.begin(); it != linepoints.end(); it++) {
                    // If in shadow, apply shadowFactor
                    PointLight light = *it;
                    float shadowFactor = testVisibilityLightSample(light.position, bvh, features, ray, hitInfo);
                    if (shadowFactor == 1.0f) {
                        color += shadowFactor * computeShading(light.position, light.color, features, ray, hitInfo);
                    }
                }

                
                total += color / glm::vec3 { segmentLightPoints, segmentLightPoints , segmentLightPoints  };
                

            } else if (std::holds_alternative<ParallelogramLight>(light) && features.enableSoftShadow) {
                // TODO: add check for enableSoftShadows -> what to return when softShadows is disabled?
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                std::list<PointLight> points;
                int sqrtSampleCount = int(std::sqrt(parallelogramLightPoints));

                points = sampleParallelogramLight(parallelogramLight, points, sample, sample);

                glm::vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);

                // interate over the sampled points and add all the colors together
                for (std::list<PointLight>::iterator it = points.begin(); it != points.end(); it++) {
                    // If in shadow, apply shadowFactor
                    PointLight light = *it;
                    float shadowFactor = testVisibilityLightSample(light.position, bvh, features, ray, hitInfo);
                    if (shadowFactor == 1.0f) {
                        color += shadowFactor * computeShading(light.position, light.color, features, ray, hitInfo);
                    }
                }

                float denom = (sqrtSampleCount + 1) * (sqrtSampleCount + 1);
                total += color / glm::vec3 { sqrtSampleCount * sqrtSampleCount, sqrtSampleCount * sqrtSampleCount, sqrtSampleCount * sqrtSampleCount };
            }
        }

        drawRay(ray, total);
        return total;
    }
    // If shading is disabled, return the albedo of the material.
    return hitInfo.material.kd;
}
