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



class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    

    
    void LevelNodes(int NodeId, std::vector<BVHNode>& resultarray, int currentlevel, int level);
    int tree_height(int NodeId);
    

    void UpdateNodeBounds(int NodeId);
    void subdivide(int NodeId, int axis);
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