#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <iostream>

// Helper to reflect vectors
glm::vec3 reflect(const glm::vec3& L, const glm::vec3& N)
{
    return L - 2 * glm::dot<3, float, glm::qualifier::highp>(L, N) * N;
}

const glm::vec3 computePhongSpecularity(const glm::vec3& lightPosition, const glm::vec3& lightColor, Ray ray, HitInfo hitInfo)
{
    // Phong specularity
    glm::vec3 vertexPosition = ray.origin + ray.direction * ray.t;
    glm::vec3 lightDirection = (lightPosition - vertexPosition);
    glm::vec3 R = reflect(-lightDirection, hitInfo.normal);
    glm::vec3 V = glm::normalize(ray.origin - vertexPosition);
    if (glm::dot(R, -V) < 1) {
        return { 0.0f, 0.0f, 0.0f };
    }
        
    return lightColor * hitInfo.material.ks * pow(fmax(dot(V, R), 0.0f), hitInfo.material.shininess);
}

const glm::vec3 computePhongDiffuse(const glm::vec3& lightPosition, const glm::vec3& lightColor, Ray ray, HitInfo hitInfo, Features features)
{
    // Phong diffuse
    glm::vec3 lightDirection = glm::normalize(lightPosition - (ray.origin + ray.direction * ray.t));

    /*if (hitInfo.material.kdTexture) {
        hitInfo.material.kdTexture->getTexel(acquireTexel(hitInfo.material.kdTexture.get(), hitInfo.texCoord, features));
    }*/

    return lightColor * hitInfo.material.kd * fmax(dot(hitInfo.normal, lightDirection), 0.0f);
}

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        return computePhongDiffuse(lightPosition, lightColor, ray, hitInfo, features) + computePhongSpecularity(lightPosition, lightColor, ray, hitInfo);
    }
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    //http://paulbourke.net/geometry/reflected/ 
    // ks check in getfinalcolour
 

    glm::vec3 vecray = glm::normalize(ray.direction);
    glm::vec3 normal = glm::normalize(hitInfo.normal);
  
    Ray reflectionRay
    {
 
    ray.origin + (ray.direction * ray.t) , 
    vecray + (normal * (2 * -glm::dot(normal,vecray))), 
    std::numeric_limits<float>::max() 
       
    };
   
    
    // TODO: implement the reflection ray computation

    return reflectionRay;
}