#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include <glm/glm.hpp>
#include <iostream>

std::vector<BVHNode> nodes;
int rootNodeId = 0, nodesUsed = 1;
std::vector<glm::uvec3> alltriangles;
std::vector<Mesh> allmeshes;
std::vector<int> meshpointer;

std::vector<int> triIdx;


BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    Scene ourscene = pScene[0];

    // nu alleen gedaan met triangels vgm moet het ook met spheres
    //  TODO: implement BVH construction
    int meshid = 0;
    int tri = 0;
    for (Mesh m : ourscene.meshes) {
        for (glm::uvec3 t : m.triangles) {
            alltriangles.push_back(t);
            triIdx.push_back(tri);
            meshpointer.push_back(meshid);
            tri++;
        }
        allmeshes.push_back(m);
        meshid++;
    }
    int N = tri;

    BVHNode& root = nodes[rootNodeId];
    root.leftChild = 0;
    root.rightChild = 0;
    root.firsttri = 0, root.triCount = N;
    UpdateNodeBounds(rootNodeId);
    subdivide(rootNodeId, 0);
}

void UpdateNodeBounds(int NodeId)
{
    BVHNode& node = nodes[NodeId];
    node.aabbMin[0] = std::numeric_limits<float>::max();
    node.aabbMin[1] = std::numeric_limits<float>::max();
    node.aabbMin[2] = std::numeric_limits<float>::max();

    node.aabbMax[0] = std::numeric_limits<float>::min();
    node.aabbMax[1] = std::numeric_limits<float>::min();
    node.aabbMax[2] = std::numeric_limits<float>::min();

    for (int first = node.firsttri, i = 0; i < node.triCount; i++) {
        int triId = triIdx[first + i];
        glm::uvec3& triangle = alltriangles[triId];
        int meshid = meshpointer[triId];
        Mesh mesh = allmeshes[meshid];


        node.aabbMin[0] = fminf(node.aabbMin[0], mesh.vertices[triangle.x].position[0]);
        node.aabbMin[1] = fminf(node.aabbMin[0], mesh.vertices[triangle.x].position[1]);
        node.aabbMin[2] = fminf(node.aabbMin[0], mesh.vertices[triangle.x].position[2]);

        node.aabbMin[0] = fminf(node.aabbMin[0], mesh.vertices[triangle.y].position[0]);
        node.aabbMin[1] = fminf(node.aabbMin[0], mesh.vertices[triangle.y].position[1]);
        node.aabbMin[2] = fminf(node.aabbMin[0], mesh.vertices[triangle.y].position[2]);

        node.aabbMin[0] = fminf(node.aabbMin[0], mesh.vertices[triangle.z].position[0]);
        node.aabbMin[1] = fminf(node.aabbMin[0], mesh.vertices[triangle.z].position[1]);
        node.aabbMin[2] = fminf(node.aabbMin[0], mesh.vertices[triangle.z].position[2]);

        
        
        node.aabbMax[0] = fmaxf(node.aabbMax[0], mesh.vertices[triangle.x].position[0]);
        node.aabbMax[1] = fmaxf(node.aabbMax[1], mesh.vertices[triangle.x].position[1]);
        node.aabbMax[2] = fmaxf(node.aabbMax[2], mesh.vertices[triangle.x].position[2]);

        node.aabbMax[0] = fmaxf(node.aabbMax[0], mesh.vertices[triangle.y].position[0]);
        node.aabbMax[1] = fmaxf(node.aabbMax[1], mesh.vertices[triangle.y].position[1]);
        node.aabbMax[2] = fmaxf(node.aabbMax[2], mesh.vertices[triangle.y].position[2]);

        node.aabbMax[0] = fmaxf(node.aabbMax[0], mesh.vertices[triangle.z].position[0]);
        node.aabbMax[1] = fmaxf(node.aabbMax[1], mesh.vertices[triangle.z].position[1]);
        node.aabbMax[2] = fmaxf(node.aabbMax[2], mesh.vertices[triangle.z].position[2]);

    }
}

void subdivide(int NodeId, int axis)
{
    BVHNode& node = nodes[NodeId];
    if (node.triCount <= 2)
        return;
    glm::vec3 extent = node.aabbMax - node.aabbMin;
    float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
    int i = node.firsttri;
    int j = i + node.triCount - 1;
    //--to be implemented--
    while (i <= j) {
        //calculate centroid triangle
        int triId = triIdx[i];
        glm::uvec3& triangle = alltriangles[triId];
        int meshid = meshpointer[triId];
        Mesh mesh = allmeshes[meshid];

        glm::vec3 centroid = (1.0f / 3.0f) * mesh.vertices[triangle.x].position + (1.0f / 3.0f) * mesh.vertices[triangle.y].position
        +(1.0f / 3.0f) * mesh.vertices[triangle.z].position; 

        if (centroid[axis] < splitPos) {
            i++;
        }else {
            std::swap(triIdx[i], triIdx[j--]);
        }
        int leftCount = i-node.firsttri;
        if (leftCount == 0 || leftCount == node.triCount)
            return;
        int leftChildid = nodesUsed++;
        int rightChildid = nodesUsed++;
        node.leftChild = leftChildid;
        nodes[leftChildid].firsttri = node.firsttri;
        nodes[leftChildid].triCount = leftCount;
        nodes[rightChildid].firsttri = i;
        nodes[rightChildid].triCount = node.triCount- leftCount;
        node.triCount = 0;
        UpdateNodeBounds(leftChildid);
        UpdateNodeBounds(rightChildid);
        if (axis == 0) {
            subdivide(leftChildid, 1);
            subdivide(rightChildid, 1);
        } else if (axis == 1) {
            subdivide(leftChildid, 2);
            subdivide(rightChildid, 2);
        } else if (axis == 2) {
            subdivide(leftChildid, 0);
            subdivide(rightChildid, 0);
        }
        

    }
    return;
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return 1;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    // drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    // drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}