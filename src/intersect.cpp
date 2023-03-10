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
    glm::vec3 side_0 = v1 - v0;
    glm::vec3 side_1 = v2 - v1;
    glm::vec3 side_2 = v0 - v2;
    glm::vec3 point_vec_0 = p - v0;
    glm::vec3 point_vec_1 = p - v1;
    glm::vec3 point_vec_2 = p - v2;

    float alpha = glm::dot(n, glm::cross(side_0, point_vec_0));
    float beta = glm::dot(n, glm::cross(side_1, point_vec_1));
    float gamma = glm::dot(n, glm::cross(side_2, point_vec_2));
    if (alpha >= 0 && beta >= 0 && gamma >= 0) {
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
bool intersectRayWithTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2, const Plane plane, Ray& ray, HitInfo& hitInfo, Features features)
{
    float t = ray.t;
    if (intersectRayWithPlane(plane, ray)) {
        glm::vec3 p = ray.origin + ray.t * ray.direction;

        if (pointInTriangle(v0.position, v1.position, v2.position, plane.normal, p)) {
            // Point must be in triangle, therefore hit
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
