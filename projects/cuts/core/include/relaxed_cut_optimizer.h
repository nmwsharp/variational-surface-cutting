#pragma once

#include <vector>
#include <cmath>
#include <array>

#include "utilities.h"
#include "geometry.h"
#include "fast_cholesky.h"
#include "fast_triangle_soup.h"


class RelaxedCutOptimizer{
    
public:

    RelaxedCutOptimizer(Geometry<Euclidean> *geometry);
    ~RelaxedCutOptimizer();


    // Parameters
    double alphaReg = 1.0;
    double stepSize = 1.0;
    size_t iIter = 0;

    // High level control
    void doStep();

    // State
    VertexData<double> u, y;

    // Initialization
    void initialize();

    // Updates
    VertexData<double> buildGradY();
    void projectY();
    void updateU();
    void implicitRegularizerStep();


    // Visualization
    VertexData<double> getLaplacianY();

private:

    HalfedgeMesh *mesh;
    Geometry<Euclidean> *geometry;

    void cacheGeometry();
    size_t nVert;
    size_t nEdge;
    size_t nFace;
    double surfaceArea = -1;
    double lengthScale = -1;
    size_t INVALID_IND;
    VertexData<size_t> vInd;
    FaceData<size_t> fInd;
    EdgeData<size_t> eInd;
    std::vector<HalfedgePtr> halfedgeList;  // need to shuffle
    EdgeData<double> edgeLengths;
    FaceData<double> faceAreas;
    DualFaceData<double> dualArea;

    double Ksum;
    VertexData<double> integratedK;

    GC::SparseMatrix<double> laplacian;

    // Cached implicit step
    double cachedCoef = -1;
    GC::SparseMatrix<double> cachedImplicitStep;
    
};