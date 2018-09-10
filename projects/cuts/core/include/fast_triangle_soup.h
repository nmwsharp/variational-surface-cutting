#pragma once

#include <cstdlib>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <array>

#include <vector3.h>
#include <utilities.h>
#include <geometry.h>

#include "halfedge_mesh.h"
#include "fast_cholesky.h"

typedef std::array<size_t, 3> Tri;
typedef std::array<bool, 3> TriBool;

class FastTriangleSoup {

public:
    FastTriangleSoup();
    FastTriangleSoup(Geometry<Euclidean>* geometry);
    FastTriangleSoup(const std::vector<Tri>& triangles_, const std::vector<Vector3>& vertices_, const std::vector<TriBool>& isBoundaryEdge_);

    ~FastTriangleSoup();

    // Convert to halfedges mesh (allocates new objects)
    void toHalfedgeMesh(HalfedgeMesh*& mesh, Geometry<Euclidean>*& geometry);

    // Mesh data
    std::vector<Tri> triangles;
    std::vector<Vector3> vertices;
    std::vector<TriBool> isBoundaryEdge;
    std::vector<bool> isBoundaryVertex;
    size_t nBoundaryVertex;
    std::vector<double> faceCurvature;
    
    bool geometryCached = false;
    void cacheGeometry();

    void buildInteriorLaplacian();
    void solveLaplaceDirichlet(std::vector<double>& bVals); // non-finte bVal means interior
    void extendFromBoundary(std::vector<double>& vals); // vals is defined on all vertices, input on boundary, output is overwritten on interior
    std::vector<double> solveInteriorPoisson(const std::vector<double>& pointwiseDensity, bool lumpedRHS=false); // return vector holds pointwise u on the interior and pointwise du/dn on the boundary
    std::vector<double> solveInteriorPoissonFace(const std::vector<double>& pointwiseFaceDensity); // return vector holds pointwise u on the interior and pointwise du/dn on the boundary
    std::vector<double> solveYamabeProblem(); // return vector holds pointwise u on the interior and pointwise du/dn on the boundary

    // Remeshing
    void improveQuality(int nIterations);
    void collapseBoundaryEdges();
    std::vector<double> applyCollapseMap(const std::vector<double>& vals); // transfer values from collapsed mesh to uncollapsed mesh
    std::vector<double> invertCollapseMap(const std::vector<double>& vals); // transfer values from uncollapsed mesh to collapsed mesh
    std::vector<size_t> collapseVertInd; // map from index on original vertex list to index in new vertex list

    // Useful cached geometry things
    size_t nVert;
    size_t nTri;
    std::vector<std::array<double,3>> cotan;
    std::vector<double> integratedGaussCurvature;
    std::vector<double> dualArea;
    std::vector<double> triangleArea;
    std::vector<Vector3> faceNormal;
    std::vector<double> boundaryLength;

    std::vector<Vector3> computeVertexNormals();

private:

    const size_t INVALID_IND = std::numeric_limits<size_t>::max();

    // Cache the interior Laplacian
    bool haveInteriorLaplacian = false;
    FastCholesky* interiorLaplacian = nullptr;
    size_t nInterior = 0;
    vector<size_t> interiorInd;

};



