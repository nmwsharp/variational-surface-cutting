#pragma once

#include "fast_triangle_soup.h"


class SurfacePatch {

public:

    SurfacePatch(const std::vector<Tri>& triangles_, 
                 const std::vector<Vector3>& vertices_, 
                 const std::vector<TriBool>& boundaryEdges_, 
                 const std::vector<size_t>& parentVertex_,
                 const std::vector<size_t>& parentFace_);

    ~SurfacePatch();


    FastTriangleSoup soupMesh;
    size_t nVert, nTri;
    std::vector<size_t> parentVertex; // maps each vertex to one on cut mesh
    std::vector<size_t> parentFace;   // maps each face to one on cut mesh
    int iRegion = -1;
    size_t iComp;

    // == Boundary
    size_t nBoundaryVertices;
    std::vector<size_t> boundaryVertexIndices; // list of boundary vertices
    std::vector<size_t> bInd; // indexer for boundary vertices (maps meshInd --> boundaryInd)

    // == Scale factor
    std::vector<double> localScaleFactor;
    double globalScaleFactor;
    
    // Smooth boundary values around sliver elements
    bool haveBoundaryValueSmoother = false;
    std::vector<double> bIdFrac;
    GC::SparseMatrix<double> bValSmoother;
    void buildBoundaryValueSmoother();
    std::vector<double> smoothOverSliverBoundaryElements(const std::vector<double>& vals);


    // == Halfedge Version
    Geometry<Euclidean>* geometry = nullptr;
    void generateGeometry();


    // == Solve systems
    void extendFromBoundary(std::vector<double>& vals);


    // == Compute data
    double area;
    void computeArea();

    Vector3 meanNormal;
    std::vector<double> normalDeviation; // per-face
    void computeNormalDeviation();

    std::vector<double> distortion; // defined everywhere
    std::vector<double> yamabeDuDn; // defined on boundary
    void solveYamabeProblem();

    // == Compute gradients
    std::vector<double> boundaryGradient; // outward-positive convention, defined on boundary
    void clearBoundaryGradient();
    void addDirichletDistortionGradient(double weight, bool useLocalScaling=false);
    void addHenckyDistortionGradient(double weight);
    void addAreaGradient(double weight, double targetArea = -1.0);
    void addNormalDeviationGradient(double weight, bool useLocalScaling=false);

    // == Compute energy
    double computeBoundaryLengthEnergy();
    double computeLocalScaledBoundaryLengthEnergy();
    double computeVisibilityEnergy(const std::vector<double>& cutMeshVisibility);
    double computeDirichletDistortionEnergy();
    double computeLocalScaledDirichletDistortionEnergy();
    double computeHenckyDistortionEnergy();
    double computeAreaEnergy(double targetArea = -1.0);
    double computeNormalDeviationEnergy();

private:

    size_t INVALID_IND; 

    // Mainly for debugging
    void findBoundaryChain();
    std::vector<size_t> boundaryChain; // vertex indices in a connected loop around the boundary
                                       // (any boundary if multiple)
    std::vector<double> boundaryChainDist;  // distance along the boundary

};