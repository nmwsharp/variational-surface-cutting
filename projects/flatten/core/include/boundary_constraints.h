// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "parameterization.h"

class BoundaryConstraints {
public:
    // Constructor
    BoundaryConstraints(Parameterization *param_, bool setConstraintsOnCuts_ = true);
    
    // Sets target curvatures for square boundary
    void setSquareConstraints();

    // Sets target curvatures for L shaped boundary
    void setLShapedConstraints();

    // Sets target curvatures for star boundary
    void setStarConstraints();

    // Sets target curvatures for uniform boundary
    void setUniformConstraints();
    
    // Sets target lengths for sinusoid boundary
    void setSinusoidConstraints();
    
    // Sets target lengths for gaussian boundary
    void setGaussianConstraints();
    
    // Sets target lengths for isometric boundary
    void setIsometricConstraints();

protected:
    // Member variables
    Parameterization *param;
    Geometry<Euclidean> *geometry;
    HalfedgeMesh *mesh;
    const CornerData<size_t>& wIndices;
    bool setConstraintsOnCuts;
    int boundarySize;
};
