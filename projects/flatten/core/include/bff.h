// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "parameterization.h"

class BoundaryFirstFlattening: public Parameterization {
public:
    // Constructor
    BoundaryFirstFlattening(Geometry<Euclidean> *geometry_);
    
    // Flatten
    bool flatten() override;
    
protected:
    // Setup - performs Cholesky decomposition of cotan matrix
    void setup(bool assignVertexIndices = false);
    
    // Computes wedge angle sum
    double computeWedgeAngleSum(HalfedgePtr h) const;
    
    // Computes boundary scale factors
    void computeBoundaryScaleFactors(GC::DenseMatrix<double>& boundaryScaleFactors);
    
    // Modifies lengths to ensure gamma closes
    void closeLengths(GC::DenseMatrix<double>& lengths,
                      const GC::DenseMatrix<Complex>& tangents) const;
    
    // Computes curve approximation
    void computeCurveApproximation(GC::DenseMatrix<double>& gammaReal,
                                   GC::DenseMatrix<double>& gammaImag,
                                   GC::DenseMatrix<double>& lengths,
                                   const GC::DenseMatrix<double>& curvatures) const;
    
    // Solves dirichlet laplace problem
    void solveDirichletLaplace(GC::DenseMatrix<double>& x, GC::DenseMatrix<double>& rhs);
    
    // Computes harmonic extension
    void computeHarmonicExtension(const GC::DenseMatrix<double>& boundaryData);
    
    // Solves neumann laplace problem
    void solveNeumannLaplace(GC::DenseMatrix<double>& x, GC::DenseMatrix<double>& rhs);
    
    // Computes harmonic conjugate
    void computeHarmonicConjugate();
    
    // Computes flattening with target lengths
    void computeTargetLengthFlattening();
    
    // Extends real and imaginary components of boundary curve
    void extendBoundaryCurve(const GC::DenseMatrix<double>& gammeReal,
                             const GC::DenseMatrix<double>& gammaImag);
    
    // Computes flattening with target angles
    void computeTargetAngleFlattening();
    
    // Computes flattening with prescribed cone angles
    void computeConeFlattening();
    
    // Sets uvs
    void setUVs();
    
    // Member variables
    bool performSetup;
    int wc, iwc, bwc;
    GC::DenseMatrix<double> u, v;
    CornerData<size_t> indices;
    GC::SparseMatrix<double> A, Aii, Aib, Abb;
#ifndef HAVE_SUITESPARSE
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solverLn;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solverLd;
#endif
};
