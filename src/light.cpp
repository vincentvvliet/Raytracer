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

// Global variables
int segmentLightPoints = 50;
int parallelogramLightPoints = 100;
float focalLength = 3.0f;
int depthOfFieldPoints = 100;
float aperture = 0.5f;

Ray depthOfField(Ray ray, int i, int points)
{
    glm::vec3 focalPoint = ray.origin + (float)focalLength * ray.direction;
    double pi = 2 * acos(0.0);
    float rho = sqrt((i + 0.001f) / (points));
    float theta = pi * (3.0f - sqrt(5.0f)) * (i - 1.0f);

    glm::vec3 secondaryPosition = ray.origin + ((float)aperture) * glm::vec3 { cos(theta) * rho, sin(theta) * rho, 0.0f };

    Ray secondaryRay { secondaryPosition, (focalPoint - secondaryPosition) };
    drawRay(secondaryRay, glm::vec3(0.0f, 1.0f, 0.0f));
    return secondaryRay;
}

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
std::list<PointLight> interpolateLine(const SegmentLight& light, std::list<PointLight> points, int segmentLightPoints) {

    float scale = 1 / sqrt(segmentLightPoints);
    for (int i = 0; i < sqrt(segmentLightPoints); i++) {
        glm::vec3 position = (light.endpoint0 - light.endpoint1) * (scale * i) + light.endpoint1;
        float alpha = (float)glm::distance(light.endpoint0, position) / glm::distance(light.endpoint0,light.endpoint1);
        glm::vec3 color = (1 - alpha) * light.color0 + alpha * light.color1;
        points.insert(points.begin(),PointLight(position,color));
    }

    return points;
}
  
std::list<PointLight> sampleSegmentLight(const SegmentLight& segmentLight, std::list<PointLight> points,glm::vec3& position, glm::vec3& color)
{   
    return interpolateLine(segmentLight,points,segmentLightPoints);
}


std::list<PointLight> interpolate(const ParallelogramLight& light, std::list<PointLight> points, int sqrtSampleCount, glm::vec3 v1, glm::vec3 v2, float ratio, float dist_x, float dist_y)
{
    // Uses bilinear interpolation for samples
    for (int i = 0; i <= sqrtSampleCount; i++) {
        glm::vec3 point1 = (light.v0 - v1) * (ratio * i) + v1;
        float alpha = (float)glm::distance(light.v0, point1) / dist_x;
        for (int j = 0; j <= sqrtSampleCount; j++) {
            glm::vec3 color1 = (1 - alpha) * light.color0 + (alpha)*light.color1;
            glm::vec3 color2 = (1 - alpha) * light.color2 + (alpha)*light.color3;
            glm::vec3 point2 = (light.v0 - v2) * (ratio * j) + v2;
            float beta = (float)glm::distance(light.v0, point2) / dist_y;
            glm::vec3 finalColor = (1 - beta) * color2 + (beta) * color1;
            glm::vec3 finalPoint = -(light.v0 - v2) * (ratio * j) + point1;

            points.insert(points.begin(), PointLight(finalPoint, finalColor));
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
    float ratio = 1.0f / sqrtSampleCount;

    return interpolate(parallelogramLight, parallelogramPoints, sqrtSampleCount, v1, v2, ratio, dist_x, dist_y);

}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
glm::vec3 testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3 debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 lightVector = glm::normalize(intersection - samplePos);
    glm::vec3 camVector = intersection - samplePos;
    float distance = glm::length(camVector) - 0.00001f;

    Ray shadowRay = Ray { samplePos, lightVector, distance };
    HitInfo shadowRayInfo;

    glm::vec3 colour = debugColor;
    float t = shadowRay.t;
    while (bvh.intersect(shadowRay, shadowRayInfo, features) && shadowRay.t < distance) {
        if (features.extra.enableTransparency && shadowRayInfo.material.transparency < 1.0f) {
            colour = colour * shadowRayInfo.material.transparency + (1 - shadowRayInfo.material.transparency) * shadowRayInfo.material.kd;
            shadowRay = { shadowRay.origin + shadowRay.t * shadowRay.direction,
                glm::normalize(ray.origin + ray.t * ray.direction - samplePos),
                glm::length((shadowRay.origin + shadowRay.t * shadowRay.direction - (ray.origin + ray.t * ray.direction))) - 0.00001f };
        } else return { 0, 0, 0 };
    }
    return colour;
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
        bool hardShadows = features.enableHardShadow;
        std::list<PointLight> segmentPoints; 
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Perform your calculations for a point light.
                      
                glm::vec3 lightColour = pointLight.color;
                if (features.enableHardShadow) {
                    lightColour = testVisibilityLightSample(pointLight.position, lightColour, bvh, features, ray, hitInfo);
                }
               
                total += computeShading(pointLight.position, lightColour, features, ray, hitInfo);              
           
            } else if (std::holds_alternative<SegmentLight>(light) && features.enableSoftShadow) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                std::list<PointLight> linepoints;
                
                linepoints = sampleSegmentLight(segmentLight,linepoints,sample,sample);
                glm::vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);
                
                // interate over the sampled points and add all the colors together
                for (std::list<PointLight>::iterator it = linepoints.begin(); it != linepoints.end(); it++) {
                    // If in shadow, apply shadowFactor
                    PointLight light = *it;
                    glm::vec3 lightColour = light.color;
                    color += computeShading(light.position, lightColour, features, ray, hitInfo);
                }

                total += color / glm::vec3 { sqrt(50), sqrt(50), sqrt(50) };
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
                    glm::vec3 lightColour = light.color;
                    if (features.enableSoftShadow) {
                        lightColour = testVisibilityLightSample(light.position, lightColour, bvh, features, ray, hitInfo);
                    }
                    color += computeShading(light.position, lightColour, features, ray, hitInfo);
                }
                
                total += color / (float) parallelogramLightPoints;
            }
        }

        if (features.extra.enableDepthOfField) {
            for (int i = 0; i < depthOfFieldPoints; i++) {
                depthOfField(ray, i, depthOfFieldPoints);
            }
        }

        drawRay(ray, total);
        return total;
    }
    // If shading is disabled, return the albedo of the material.
    return hitInfo.material.kd;
}
