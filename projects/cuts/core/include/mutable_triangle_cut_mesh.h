#pragma once

#include "geometry.h"

#include <array>
#include <unordered_set>

// Forward declarations
struct MutableCutHalfedge;

struct MutableCutVertex {

    MutableCutVertex(Vector3 p_);

    Vector3 position;
    MutableCutHalfedge* halfedge = nullptr;
    bool isBoundary = false;
    
    // Not always populated
    size_t index;
    double desiredVertexDensity; // used for adaptive remeshing
    bool touched = false;
    
    Vector3 normal();
    bool isOnCut();
    int degree();
};

struct MutableCutFace {

    MutableCutHalfedge* halfedge = nullptr; // NOTE: Not guaranteed to be real FORNOW, unlike HalfedgeMesh
                                         // (can be done, but not implemented)

    Vector3 normal();
    Vector3 areaNormal();
    Vector3 barycenter();
    Vector3 circumcenter();
    double desiredDensity();
};

struct MutableCutHalfedge {

    MutableCutHalfedge* twin = nullptr;
    MutableCutHalfedge* next = nullptr;
    MutableCutVertex* vertex = nullptr;
    MutableCutFace* face = nullptr;
    bool isCut = false;
    
    bool isOnBoundary();
    bool isReal();

    Vector3 vector();
    double oppositeAngle();
};


class MutableCutMesh {

    public:

        MutableCutMesh(Geometry<Euclidean>* geometry);
        ~MutableCutMesh();

        // Parameters
        double normalThreshold = PI/2.0; // how much a modification is allowed to change a face normal
        double skinnyThresh = 0.01; // relative threshold for deciding if a face is skinny

        // Modifications
        bool tryFlipEdge(MutableCutHalfedge* he1);
        bool tryCollapseEdge(MutableCutHalfedge* he1);
        // bool trySplitEdge(MutableCutHalfedge* he1);
        bool tryMoveVertex(MutableCutVertex* v, Vector3 newPos);

        Vector3 splineMidpoint(Vector3 p1, Vector3 n1, Vector3 p2, Vector3 n2);

        Geometry<Euclidean>* toHalfedgeMesh();
        void indexVertices();
        void computeSizingField(double absoluteEpsilon);
        void setVerticesUntouched();

        void checkDegree();

        // Members
        std::unordered_set<MutableCutVertex*> vertices;
        std::unordered_set<MutableCutFace*> faces;
        std::unordered_set<MutableCutHalfedge*> halfedges;
    


};

