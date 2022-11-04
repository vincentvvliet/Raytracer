#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;
struct BVHNode {
    glm::vec3 aabbMin;
    glm::vec3 aabbMax;
    int leftChild;
    int firsttri, triCount;
    bool isLeaf() const { return triCount>0;}
};

struct aabSAHHelper {
    glm::vec3 aabbMin = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
    glm::vec3 aabbMax = { std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min() };
    void grow(glm::vec3 p) {
        aabbMin[0] = fminf(aabbMin[0], p[0]);
        aabbMin[1] = fminf(aabbMin[1], p[1]);
        aabbMin[2] = fminf(aabbMin[2], p[2]);

        aabbMax[0] = fmaxf(aabbMax[0], p[0]);
        aabbMax[1] = fmaxf(aabbMax[1], p[1]);
        aabbMax[2] = fmaxf(aabbMax[2], p[2]);
    }

    float area() {
        glm::vec3 extent = aabbMax - aabbMin;
        return extent[0] * extent[1] + extent[1] * extent[2] + extent[2] * extent[0];
    }
};



class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    

    
    void LevelNodes(int NodeId, std::vector<BVHNode>& resultarray, int currentlevel, int level);
    int tree_height(int NodeId);
    

    void UpdateNodeBounds(int NodeId);
    void subdivide(int NodeId, int axis, bool SAHBinning);
    float EvaluateSAH(BVHNode node, int axis1, float pos);
    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;
    


    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    bool intersectBVH(Ray& ray, HitInfo& hitInfo, int NodeId, Features features) const;



private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<BVHNode> nodes;
    std::vector<int> triIdx;
    std::vector<glm::uvec3> alltriangles;
    int nodesUsed = 1;

    std::vector<Mesh> allmeshes;
    std::vector<int> meshpointer;

    int level;
};