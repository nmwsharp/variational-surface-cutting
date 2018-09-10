#pragma once

#include "halfedge_mesh.h"
#include "geometry.h"
#include "utilities.h"

#include <Eigen/Core>

#include <array>

struct DevelopableApproximation {


    DevelopableApproximation(Geometry<Euclidean>* geometry_);

    void solve(const VertexData<char>& vertexFixed); // result is in devGeom member below

    Geometry<Euclidean>* devGeom; // WARNING: User needs to free this, since it is usually considered the output of the program

private:

    // Parameters
    double alphaPreserveAngles = 0.01;

    // Methods 
    double computeEnergy();
    VertexData<Vector3> computeGradient();
    double computeInteriorCurvature(); // away from fixed vertices

    void buildFilterMatrix();
    void filterGradient(VertexData<Vector3>& gradient);
    GC::SparseMatrix<double> filterMatrix;

    HalfedgeMesh* mesh;
    Geometry<Euclidean>* initGeometry;
    size_t nVert;
    VertexData<size_t> vInd;

    VertexData<char> vFixed;
    HalfedgeData<double> initAngles;
};