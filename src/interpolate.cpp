#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // Barycentric coordinates using Christer Ericson's Real-Time Collision Detection algorithm
    glm::vec3 a = v1 - v0, b = v2 - v0, c = p - v0;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;

    float alpha = (d11 * d20 - d01 * d21) / denom;
    float beta = (d00 * d21 - d01 * d20) / denom;
    float gamma = 1.0f - alpha - beta;

    return glm::vec3 { alpha, beta, gamma };
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord[0] * n0 + barycentricCoord[2] * n1 + barycentricCoord[1] * n2;
}

glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return barycentricCoord[0] * t0 + barycentricCoord[2] * t1 + barycentricCoord[1] * t2;
}
