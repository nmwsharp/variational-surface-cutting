#include "relaxed_cut_optimizer.h"

#include <csignal>
#include <queue>
#include <tuple>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "timing.h"

#define PLOT_ENERGY (true)

using namespace GC;
using std::vector;

RelaxedCutOptimizer::RelaxedCutOptimizer(Geometry<Euclidean> *geometry_) : 
    mesh(geometry_->getMesh()), geometry(geometry_)
{
    cacheGeometry();

    initialize();
}

RelaxedCutOptimizer::~RelaxedCutOptimizer() {
}

void RelaxedCutOptimizer::cacheGeometry() {

    // Indices
    INVALID_IND = std::numeric_limits<size_t>::max();
    vInd = mesh->getVertexIndices();
    fInd = mesh->getFaceIndices();
    eInd = mesh->getEdgeIndices();
    nVert = mesh->nVertices();
    nEdge = mesh->nEdges();
    nFace = mesh->nFaces();

    geometry->getEdgeLengths(edgeLengths);
    geometry->getFaceAreas(faceAreas);
    geometry->getDualFaceAreas(dualArea);
    surfaceArea = geometry->totalArea();
    lengthScale = std::sqrt(surfaceArea);

    halfedgeList.reserve(mesh->nHalfedges());
    for(HalfedgePtr he : mesh->halfedges()) {
        halfedgeList.push_back(he);
    }


    // Get curvatures
    geometry->getVertexAngleDefects(integratedK); 
    Ksum = 0.0;
    for(VertexPtr v : mesh->vertices()) {
        Ksum += integratedK[v] * integratedK[v]; // TODO this almost certainly doesn't make sense discretization-wise
    }

    // Build the laplacian
    laplacian = cotanMatrix<double>(geometry, vInd);
    laplacian.shift(1e-4); 
}


void RelaxedCutOptimizer::initialize() {

    iIter = 0;

    // Randomly sample y
    y = VertexData<double>(mesh);
    for(VertexPtr v : mesh->vertices()) {
        y[v] = randomReal(-1.0, 1.0);
    }


    projectY(); 
    u = VertexData<double>(mesh, 0.0);
    updateU();
}


void RelaxedCutOptimizer::doStep() {

    // Build the gradient
    VertexData<double> gradY = buildGradY();

    // Take a gradient step (excluding the regularizer)
    y -= stepSize * gradY;

    // Take an implicit step of the regularizer
    implicitRegularizerStep();

    // Project y
    projectY();

    // Update u
    updateU();

}


void RelaxedCutOptimizer::projectY() {

    // Project to unit ball
    double ySum = 0;
    for(VertexPtr v : mesh->vertices()) {
        ySum = dualArea[v.dual()] * y[v] * y[v];
    }
    double fac = std::sqrt(ySum);
    for(VertexPtr v : mesh->vertices()) {
        y[v] /= fac;
    }
    
    // Project to unit curvature ball 
    // double ySum = 0;
    // VertexData<double> ky = y * integratedK;
    // for(VertexPtr v : mesh->vertices()) {
    //     ySum = ky[v] * ky[v];
    // }
    // double fac = std::sqrt(ySum);
    // for(VertexPtr v : mesh->vertices()) {
    //     y[v] /= fac;
    // }
}

void RelaxedCutOptimizer::updateU() {

    DenseVector<double> rhs = (y * integratedK).toVector();
    DenseVector<double> uVec;

    solvePositiveDefinite(laplacian, uVec, rhs);
    u.fromVector(uVec);
}

VertexData<double> RelaxedCutOptimizer::buildGradY() {

    // Accumulate the gradient
    VertexData<double> gradY(mesh, 0.0);

    // Add u terms
    gradY += u * integratedK;

    // Add y regularizer terms
    // VertexData<double> lapY(mesh, laplacian * y.toVector());
    // gradY += 2.0 * alphaReg * lapY;

    return gradY;
}

void RelaxedCutOptimizer::implicitRegularizerStep() {

    double currCoef = (2.0 * alphaReg * stepSize);
    if(currCoef != cachedCoef) {
        cachedCoef = currCoef;
        cachedImplicitStep = SparseMatrix<double>::identity(mesh->nVertices()) + currCoef * laplacian;
    }
    
    DenseVector<double> rhs = y.toVector();
    DenseVector<double> yVec;

    solvePositiveDefinite(cachedImplicitStep, yVec, rhs);
    y.fromVector(yVec);

}

VertexData<double> RelaxedCutOptimizer::getLaplacianY() {
    return VertexData<double>(mesh, laplacian * y.toVector());
}


#undef PLOT_ENERGY