#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include <glm/glm.hpp>
#include <iostream>
#include <stack>

#include <fmt/chrono.h>
#include <fmt/core.h>

#include <chrono>

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
    Scene ourscene = pScene[0];

    bool SAHbinning;

    // Checks if we have SAHbinning enabled
    if (features.extra.enableBvhSahBinning) {
        SAHbinning = true;
    } else {
        SAHbinning = false;
    }

    //  TODO: implement BVH construction
    // Fill the vectors which are gonna be needed to
    // acces the triangles
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
    // total amount of trianhgles
    int N = tri;

    // initialize the rootnode
    glm::vec3 aabbMin = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
    glm::vec3 aabbMax = { std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min() };
    BVHNode root = { aabbMin, aabbMax, 0, 0, N };

    // put too node in out tree vector
    nodes.push_back(root);
    // update the aabb bounds of the root node so it encompass the object
    UpdateNodeAABBBounds(0);
    // start dividing the scene
    divide(0, 0, SAHbinning);
    // calculate the height of the tree just created
    level = tree_height(0);

    const auto end = clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    fmt::print("Creating the BVH took {} ms\n", duration);
}

void BoundingVolumeHierarchy::UpdateNodeAABBBounds(int NodeId)
{

    BVHNode& node = nodes[NodeId];
    // for every triangle in the node we update the corresponding aabb to encompass all triangles
    // in the node
    for (int first = node.firsttri, i = 0; i < node.triCount; i++) {
        int triId = triIdx[first + i];
        glm::uvec3& triangle = alltriangles[triId];
        int meshid = meshpointer[triId];
        Mesh mesh = allmeshes[meshid];

        node.aabbMin[0] = fminf(node.aabbMin[0], mesh.vertices[triangle.x].position[0]);
        node.aabbMin[1] = fminf(node.aabbMin[1], mesh.vertices[triangle.x].position[1]);
        node.aabbMin[2] = fminf(node.aabbMin[2], mesh.vertices[triangle.x].position[2]);

        node.aabbMin[0] = fminf(node.aabbMin[0], mesh.vertices[triangle.y].position[0]);
        node.aabbMin[1] = fminf(node.aabbMin[1], mesh.vertices[triangle.y].position[1]);
        node.aabbMin[2] = fminf(node.aabbMin[2], mesh.vertices[triangle.y].position[2]);

        node.aabbMin[0] = fminf(node.aabbMin[0], mesh.vertices[triangle.z].position[0]);
        node.aabbMin[1] = fminf(node.aabbMin[1], mesh.vertices[triangle.z].position[1]);
        node.aabbMin[2] = fminf(node.aabbMin[2], mesh.vertices[triangle.z].position[2]);

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

void BoundingVolumeHierarchy::divide(int NodeId, int axis, bool SAHbinning)

{
    BVHNode& node = nodes[NodeId];

    if (node.triCount <= 2) {
        return;
    }

    glm::vec3 extent = node.aabbMax - node.aabbMin;
    // stansard position where we are gonna split the scene
    float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;

    // SAH
    if (SAHbinning) {
        // get the most coseffective way of splitting the scene on a certain axis
        // where cost is defined by the area of the whole aabb's times the triangles in the aabb's
        float bestCost = std::numeric_limits<float>::max();
        for (int i = 0; i < node.triCount; i++) {
            int triId = triIdx[node.firsttri + i];
            glm::uvec3& triangle = alltriangles[triId];
            int meshid = meshpointer[triId];
            Mesh mesh = allmeshes[meshid];

            glm::vec3 centroid = (1.0f / 3.0f) * mesh.vertices[triangle.x].position + (1.0f / 3.0f) * mesh.vertices[triangle.y].position
                + (1.0f / 3.0f) * mesh.vertices[triangle.z].position;

            // calculate the costo f having the split position at every triangle in the aabb
            float candidatePos = centroid[axis];
            float cost = costSAH(node, axis, candidatePos);

            if (cost < bestCost) {
                // get the split that optimizez this cost
                splitPos = candidatePos;
                bestCost = cost;
            }
        }
    }

    int i = node.firsttri;
    int j = i + node.triCount - 1;

    // keep swapping the triangles indices in the tri errar such that they are dvidied
    // between left of the split pos and right of the split positions
    while (i <= j) {
        // calculate centroid triangle
        int triId = triIdx[i];
        glm::uvec3& triangle = alltriangles[triId];
        int meshid = meshpointer[triId];
        Mesh mesh = allmeshes[meshid];

        glm::vec3 centroid = (1.0f / 3.0f) * mesh.vertices[triangle.x].position + (1.0f / 3.0f) * mesh.vertices[triangle.y].position
            + (1.0f / 3.0f) * mesh.vertices[triangle.z].position;

        if (centroid[axis] < splitPos) {
            i++;
        } else {
            std::swap(triIdx[i], triIdx[j--]);
        }
    }
    // totoal amount of triangles on the left aabb
    int leftCount = i - node.firsttri;
    if (leftCount == 0 || leftCount == node.triCount)
        return;

    // set the leftchild id to beign the nodesused+
    // same for the rightid
    // this will be used when traversing the scene tree
    int leftChildid = nodesUsed++;
    int rightChildid = nodesUsed++;

    node.leftChild = leftChildid;

    glm::vec3 aabbMin = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
    glm::vec3 aabbMax = { std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min() };

    BVHNode leftchild = { aabbMin, aabbMax, 0, node.firsttri, leftCount };
    BVHNode rightchild = { aabbMin, aabbMax, 0, i, node.triCount - leftCount };

    nodes.push_back(leftchild);
    nodes.push_back(rightchild);

    // the node value suddenly changed so i did this
    BVHNode& node1 = nodes[NodeId];
    node1.triCount = 0;

    // update the aabb bounds for the children just created
    UpdateNodeAABBBounds(leftChildid);
    UpdateNodeAABBBounds(rightChildid);

    // we go trhough the axis's one by one
    if (axis == 0) {
        divide(leftChildid, 1, SAHbinning);
        divide(rightChildid, 1, SAHbinning);
    } else if (axis == 1) {
        divide(leftChildid, 2, SAHbinning);
        divide(rightChildid, 2, SAHbinning);
    } else if (axis == 2) {
        divide(leftChildid, 0, SAHbinning);
        divide(rightChildid, 0, SAHbinning);
    }
    return;
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return level;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    // basicly counts the number of leaves in a binary tree
    int count = 0;
    for (BVHNode node : nodes) {
        if (node.isLeaf())
            count++;
    }
    return count;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    std::vector<BVHNode> nodesatlevel;

    LevelNodes(0, nodesatlevel, 0, level);

    for (BVHNode& node : nodesatlevel) {
        AxisAlignedBox aabb { node.aabbMin, node.aabbMax };
        drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    }
}

void BoundingVolumeHierarchy::LevelNodes(int NodeId, std::vector<BVHNode>& resultarray, int currentlevel, int level)
{
    // get the nodes corresponding to a certain level
    if (nodes[NodeId].isLeaf()) {
        resultarray.push_back(nodes[NodeId]);
        return;
    }
    if (currentlevel == level) {
        resultarray.push_back(nodes[NodeId]);
        return;
    }

    LevelNodes(nodes[NodeId].leftChild, resultarray, currentlevel + 1, level);
    LevelNodes(nodes[NodeId].leftChild + 1, resultarray, currentlevel + 1, level);
}

int BoundingVolumeHierarchy::tree_height(int NodeId)
{
    // calculated the height of the binary tree
    BVHNode& node = nodes[NodeId];
    if (node.isLeaf()) {
        return 1;
    } else {
        int index = nodes[NodeId].leftChild;
        int ldepth = tree_height(index);
        int rdepth = tree_height(index + 1);
        if (ldepth > rdepth) {
            return ldepth + 1;
        } else {

            return rdepth + 1;
        }
    }
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    int i = 0;
    AxisAlignedBox aabb;
    for (BVHNode node : nodes) {
        if (node.isLeaf() && i == leafIdx) {
            continue;
        } else if (node.isLeaf()) {
            AxisAlignedBox aabb { node.aabbMin, node.aabbMax };
            drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
            i++;
        }
    }
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
                Plane plane = trianglePlane(v0.position, v1.position, v2.position);
                if (intersectRayWithTriangle(v0, v1, v2, plane, ray, hitInfo, features)) {
                    hitInfo.material = mesh.material;
                    glm::vec3 normal = plane.normal;
                    glm::vec3 p = ray.origin + ray.t * ray.direction;
                    glm::vec3 bary = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

                    if (features.enableNormalInterp) {
                        if (bary.x <= 1.0f && bary.x >= 0.0f && bary.y <= 1.0f && bary.y >= 0.0f && bary.z <= 1.0f && bary.z >= 0.0f) {
                            normal = interpolateNormal(v0.normal, v1.normal, v2.normal, bary);
                        }
                    }

                    if (features.enableTextureMapping) {
                        if (hitInfo.material.kdTexture) {
                            hitInfo.material.kd = acquireTexel(*hitInfo.material.kdTexture, interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, bary));
                        }
                    }

                    hit = true;
                    hitInfo.normal = normal;
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
        return intersectBVH(ray, hitInfo, 0, features);
    }
}

bool BoundingVolumeHierarchy::intersectBVH(Ray& ray, HitInfo& hitInfo, int NodeId, Features features) const
{
    // BVH traversal code
    std::stack<BVHNode> stack;
    stack.push(nodes[NodeId]);
    float t = ray.t;
    while (!stack.empty()) {
        BVHNode curnode = stack.top();
        // For debugging
        // AxisAlignedBox aabb { curnode.aabbMin, curnode.aabbMax };
        // drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 0.1f);
        stack.pop();
        // if current node is  a lead we try to intersect with all the triangles
        if (curnode.isLeaf()) {
            for (int i = 0; i < curnode.triCount; i++) {
                int triId = triIdx[curnode.firsttri + i];
                glm::uvec3 triangle = alltriangles[triId];
                int meshid = meshpointer[triId];
                Mesh mesh = allmeshes[meshid];
                const auto v0 = mesh.vertices[triangle.x];
                const auto v1 = mesh.vertices[triangle.y];
                const auto v2 = mesh.vertices[triangle.z];
                Plane plane = trianglePlane(v0.position, v1.position, v2.position);
                if (intersectRayWithTriangle(v0, v1, v2, plane, ray, hitInfo, features)) {
                    hitInfo.material = mesh.material;
                }
            }
        } else {

            BVHNode leftnode = nodes[curnode.leftChild];
            BVHNode rightnode = nodes[curnode.leftChild + 1];

            AxisAlignedBox aabbleft { leftnode.aabbMin, leftnode.aabbMax };
            AxisAlignedBox aabbright { rightnode.aabbMin, rightnode.aabbMax };

            Ray leftray = ray;
            Ray rightray = ray;
            // If we interesect with both the left child and aabb and the right child  aabb then we push both
            if (intersectRayWithShape(aabbleft, leftray) && intersectRayWithShape(aabbright, rightray)) {
                // first push the one with the one with the heighest t
                // so that we end up rendering the one with the lowest t
                if (leftray.t > rightray.t) {
                    stack.push(rightnode);
                    stack.push(leftnode);
                } else {
                    stack.push(leftnode);
                    stack.push(rightnode);
                }
            }
            // if we only intersect with one child push that one
            else if (intersectRayWithShape(aabbleft, leftray) && !intersectRayWithShape(aabbright, rightray)) {
                stack.push(leftnode);
            } else if (!intersectRayWithShape(aabbleft, leftray) && intersectRayWithShape(aabbright, rightray)) {
                stack.push(rightnode);
            }
            continue;
        }
    }
    // if we have found a smaller t we now we intersected if we did not we no we missed
    if (ray.t < t) {
        return true;
    } else {
        return false;
    }
}

// for the cost evaluation of the candidate position for splittin the node
float BoundingVolumeHierarchy::costSAH(BVHNode& node, int axis, float pos) const
{

    glm::vec3 aabbMin = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
    glm::vec3 aabbMax = { std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min() };
    AxisAlignedBox aabbleft { aabbMin, aabbMax };
    AxisAlignedBox aabbright { aabbMin, aabbMax };
    int leftCount = 0;
    int rightCount = 0;

    // keep growing the left adn the right box untill all the corersponding triangle centroids
    //  are in there respective boxes
    for (int i = 0; i < node.triCount; i++) {
        int triId = triIdx[node.leftChild + i];
        glm::uvec3 triangle = alltriangles[triId];
        int meshid = meshpointer[triId];
        Mesh mesh = allmeshes[meshid];

        glm::vec3 centroid = (1.0f / 3.0f) * mesh.vertices[triangle.x].position + (1.0f / 3.0f) * mesh.vertices[triangle.y].position
            + (1.0f / 3.0f) * mesh.vertices[triangle.z].position;
        if (centroid[axis] < pos) {
            leftCount++;
            growAABB(aabbleft, mesh.vertices[triangle.x].position);
            growAABB(aabbleft, mesh.vertices[triangle.y].position);
            growAABB(aabbleft, mesh.vertices[triangle.z].position);
        } else {
            rightCount++;
            growAABB(aabbright, mesh.vertices[triangle.x].position);
            growAABB(aabbright, mesh.vertices[triangle.y].position);
            growAABB(aabbright, mesh.vertices[triangle.z].position);
        }
    }
    // calculate the cost of having the AABBS be this
    float cost = leftCount * areaAABB(aabbleft) + rightCount * areaAABB(aabbright);
    return cost > 0 ? cost : std::numeric_limits<float>::max();
}

void BoundingVolumeHierarchy::growAABB(AxisAlignedBox& aabb, glm::vec3 p) const
{
    // if the vertex is smaller thent he current lower
    // replace the lower with this vertex
    aabb.lower[0] = fminf(aabb.lower[0], p[0]);
    aabb.lower[1] = fminf(aabb.lower[1], p[1]);
    aabb.lower[2] = fminf(aabb.lower[2], p[2]);

    // same idea for the upper of the current aabb
    aabb.upper[0] = fmaxf(aabb.upper[0], p[0]);
    aabb.upper[1] = fmaxf(aabb.upper[1], p[1]);
    aabb.upper[2] = fmaxf(aabb.upper[2], p[2]);
}

float BoundingVolumeHierarchy::areaAABB(AxisAlignedBox& aabb) const
{
    glm::vec3 extent = aabb.upper - aabb.lower;
    return extent[0] * extent[1] + extent[1] * extent[2] + extent[2] * extent[0];
}