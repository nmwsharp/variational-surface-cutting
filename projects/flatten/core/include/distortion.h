// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "parameterization.h"

class Distortion {
public:
    // Constructor
    Distortion(Parameterization *param_);
    
    // Compute quasi-conformal distortion
    double computeQcDistortion();
    
    // Compute area distortion
    double computeAreaDistortion();
    
    // Compute angle distortion
    double computeAngleDistortion();
    
    // Compute angle sum distortion
    double computeAngleSumDistortion();
    
    // Compute cross length ratio distortion
    double computeClDistortion();
    
protected:
    // Member variables
    Parameterization *param;
    Geometry<Euclidean> *geometry;
    HalfedgeMesh *mesh;
    const CornerData<Vector2>& uvs;
};
