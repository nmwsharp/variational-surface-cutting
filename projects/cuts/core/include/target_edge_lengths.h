#pragma once

#include "halfedge_mesh.h"
#include "geometry.h"
#include "utilities.h"

#include <Eigen/Core>

#include <array>

struct PositionSolver {


    PositionSolver(Geometry<Euclidean>* geometry_);

    VertexData<Vector3> solve(const EdgeData<double>& targetEdgeLengths, const VertexData<char>& vertexFixed);
    void buildMatrices(const EdgeData<double>& targetEdgeLengths);

    // The interface for the LBFGS solver
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);

private:

    HalfedgeMesh* mesh;
    Geometry<Euclidean>* geometry;
    double bendingWeight = 1e-13;
    double globalScale = 1e6;
    // double bendingWeight = .0000;
    // double bendingWeight = 1.0;

    size_t nVert, nEdge;
    VertexData<size_t> vInd;
    EdgeData<size_t> eInd;
    std::vector<double> edgeAreas;
    std::vector<std::array<size_t, 2>> edges;

    Eigen::VectorXd targetLengths;
    std::vector<bool> vCut;


    // Operators
    GC::SparseMatrix<double> hodge0, hodge1, hodge1Inv, hodge0Inv, hodge0Bar, hodge0BarInv;
    GC::SparseMatrix<double> d0, d0T, d1, d1T;
    GC::SparseMatrix<double> eyeV;

    GC::SparseMatrix<double> zeroFormLaplacian;
    GC::SparseMatrix<double> zeroFormWeakLaplacian;
    GC::SparseMatrix<double> oneFormLaplacian;
    
    GC::SparseMatrix<double> B;
    GC::SparseMatrix<double> BBT;

    Eigen::SparseMatrix<double> Beigen;
    Eigen::SparseMatrix<double> BBTeigen;
};