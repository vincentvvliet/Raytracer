#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>
#include <iostream>
#include <interpolate.cpp>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    // Point in triangle test using barycentric coordinates.
    glm::vec3 bary = computeBarycentricCoord(v0, v1, v2, p);    

    if (bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0) {
        return true;
    }

    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    if (glm::dot((ray.origin + ray.direction * ray.t), plane.normal) - plane.D != 0) {
        float t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal);
        if (t >= 0.0f && t <= ray.t) {
            ray.t = t;
            return true;
        }
    }
    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 normal = glm::normalize(glm::cross(v0 - v2, v1 - v2) / glm::length(glm::cross(v0 - v2, v1 - v2)));
    plane.D = glm::dot(normal, v0);
    plane.normal = normal;

    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo, Features features)
{    
    Plane plane = trianglePlane(v0, v1, v2);
    float t = ray.t;
    if (intersectRayWithPlane(plane, ray)) {
        glm::vec3 p = ray.origin + ray.t * ray.direction;

        glm::vec3 normal = plane.normal;
        if (features.enableNormalInterp) {
            glm::vec3 n0 = v0 + plane.normal;
            glm::vec3 n1 = v1 + plane.normal;
            glm::vec3 n2 = v2 + plane.normal;
            normal = interpolateNormal(n0, n1, n2, computeBarycentricCoord(v0, v1, v2, p));
        }

        if (pointInTriangle(v0, v1, v2, plane.normal, p)) {
            // Point must be in triangle, therefore hit
            hitInfo.normal = normal;
            return true;
        } else {
            ray.t = t;
        }
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // from real-time collision detection by C. Ericson
    glm::vec3 center = sphere.center;
    float radius = sphere.radius;
    glm::vec3 p = ray.origin;
    glm::vec3 m = p - center;
    glm::vec3 d = glm::normalize(ray.direction);
    float b = glm::dot<3, float, glm::qualifier::highp>(m, d);
    float c = glm::dot<3, float, glm::qualifier::highp>(m, m) - radius * radius;
    if (c > 0.0f && b > 0.0f)
        return false;
    float discr = b * b - c;
    if (discr < 0.0f)
        return false;
    ray.t = -b - sqrt(discr);
    // If t is negative, ray started inside sphere so clamp t to zero
    return true;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{   
    glm::vec3 bmin = box.lower;
    glm::vec3 bmax = box.upper;

    float tx1 = (bmin.x - ray.origin.x) / ray.direction.x, tx2 = (bmax.x - ray.origin.x) / ray.direction.x;
    float tmin = fmin(tx1, tx2), tmax = fmax(tx1, tx2);
    float ty1 = (bmin.y - ray.origin.y) / ray.direction.y, ty2 = (bmax.y - ray.origin.y) / ray.direction.y;
    tmin = fmax(tmin, fmin(ty1, ty2)), tmax = fmin(tmax, fmax(ty1, ty2));
    float tz1 = (bmin.z - ray.origin.z) / ray.direction.z, tz2 = (bmax.z - ray.origin.z) / ray.direction.z;
    tmin = fmax(tmin, fmin(tz1, tz2)), tmax = fmin(tmax, fmax(tz1, tz2));

    return tmax >= tmin && tmin < ray.t && tmax > 0;
}
