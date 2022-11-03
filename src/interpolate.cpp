#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // Barycentric coordinates using Christer Ericson's Real-Time Collision Detection algorithm
    glm::vec3 a = v1 - v0, b = v2 - v0, c = p - v0;
    float d00 = glm::dot(a, a);
    float d01 = glm::dot(a, b);
    float d11 = glm::dot(b, b);
    float d20 = glm::dot(c, a);
    float d21 = glm::dot(c, b);
    float denom = d00 * d11 - d01 * d01;

    float alpha = (d11 * d20 - d01 * d21) / denom;
    float beta = (d00 * d21 - d01 * d20) / denom;
    float gamma = 1.0f - alpha - beta;

    return glm::normalize(glm::vec3 { alpha, beta, gamma });
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
}

glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord.x * t0 + barycentricCoord.y * t1 + barycentricCoord.z * t2;
}
