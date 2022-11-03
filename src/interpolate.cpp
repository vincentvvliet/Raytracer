#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    float total = glm::length(glm::cross((v1 - v0), (v2 - v0))) / 2.0f;
    float area_v0 = glm::length(glm::cross((v1 - p), (v2 - p))) / 2.0f;
    float area_v1 = glm::length(glm::cross((v0 - p), (v2 - p))) / 2.0f;

    float alpha = area_v0 / total;
    float beta = area_v1 / total;
    float gamma = 1.0f - alpha - beta;

    return glm::vec3 { alpha, beta, gamma };
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
}

glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord.x * t0 + barycentricCoord.y * t1 + barycentricCoord.z * t2;
}
