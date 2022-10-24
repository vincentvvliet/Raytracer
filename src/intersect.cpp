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


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    // Point in triangle test using barycentric coordinates. Supposedly faster than using same side technique with crossproducts
   float dot0 =  glm::dot(v0, v0);
    float dot1 =glm::dot(v0, v1);
    float dot2 = glm::dot(v0, v2);
   float dot3 = glm::dot(v1, v1);
    float dot4 = glm::dot(v1, v2);
   float denominator = 1 / ((dot0 * dot3) - (dot1 * dot1) ;
   float u = (dot3 * dot2 - dot1 * dot4) * denominator;
    float v = (dot0 * dot4 - dot1 * dot2) * denominator;
    
    
    return (u >= 0) && (v >= 0) && (u + v < 1) ;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
 

    float d = plane.D;
    glm::vec3 n = plane.normal;
    glm::vec3 rayvec = ray.origin - ray.direction;
    ray.t = (d - Dot(n, ray.origin)) / Dot(n, rayvec);
    // If t in [0..1] compute and return intersection point
    if (t >= 0.0f && t <= 1.0f) {
       
        return true;
    }
    // Else no intersection
    return false;
    
    }
 


Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
 
    // TODO: implement this function.
    Plane plane;
    plane.normal = glm::cross((v0-v1)-(v2-v1));
    plane.D = v0;
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    Plane thisPlane = trianglePlane(v0, v1, v2);
    intersectRayWithPlane(thisPlane, ray);
    return intersectRayWithPlane(thisPlane, ray);
    
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // from real-time collision detection by C. Ericson
    glm::vec3 center = sphere.center;
    glm::vec3 radius = sphere.radius;
    glm::vec3 p = ray.origin;
    glm::vec3 m = p - center;
    glm::vec3 d = glm::normalize(ray.direction);
    float b = Dot(m, d);
    float c = Dot(m, m) - radius * radius;
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
    // TODO: implement this function.
    return false;
}
