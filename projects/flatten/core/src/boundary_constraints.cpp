// Author: Rohan Sawhney
// Fall 2016

#include "boundary_constraints.h"
#include "flatten_utils.h"

BoundaryConstraints::BoundaryConstraints(Parameterization *param_, bool setConstraintsOnCuts_):
param(param_),
geometry(param->geometry),
mesh(param->mesh),
wIndices(param->wIndices),
setConstraintsOnCuts(setConstraintsOnCuts_),
boundarySize(0)
{
    if (setConstraintsOnCuts) for (HalfedgePtr h: mesh->cutBoundary()) boundarySize++;
    else boundarySize = mesh->nImaginaryHalfedges();
}

void BoundaryConstraints::setSquareConstraints()
{
    int i = 0;
    int nBy4 = boundarySize/4;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        i++;
        if (i % nBy4 == 0 && param->targetCurvatures.size() < 4) {
            param->targetCurvatures[wIndices[wedge(h)]] = M_PI_2;
        }
    }
}

void BoundaryConstraints::setLShapedConstraints()
{
    int i = 0, j = 0;
    int nBy6 = boundarySize/6;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        i++;
        if (i % nBy6 == 0 && param->targetCurvatures.size() < 6) {
            j++;
            if (j == 5) param->targetCurvatures[wIndices[wedge(h)]] = -M_PI_2;
            else param->targetCurvatures[wIndices[wedge(h)]] = M_PI_2;
        }
    }
}

void BoundaryConstraints::setStarConstraints()
{
    int i = 0, j = 0;
    int nBy10 = boundarySize/10;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        i++;
        if (i % nBy10 == 0 && param->targetCurvatures.size() < 10) {
            j++;
            if (j % 2 == 0) param->targetCurvatures[wIndices[wedge(h)]] = 2.0*M_PI/3.0;
            else param->targetCurvatures[wIndices[wedge(h)]] = -48.0*M_PI/180.0;
        }
    }
}

void BoundaryConstraints::setUniformConstraints()
{
    double curvature = 2.0*M_PI/boundarySize;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        param->targetCurvatures[wIndices[wedge(h)]] = curvature;
    }
}

void BoundaryConstraints::setSinusoidConstraints()
{
    double S = 0.0;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        S += geometry->length(h.edge());
    }
    
    double A = S/boundarySize;
    double B = 1.5;
    double f = 6.0;
    double s = 0.0;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        double theta = 2.0*M_PI*f*s/S;
        param->targetLengths[wIndices[wedge(h)]] = A*(sin(theta) + B);
        s += geometry->length(h.edge());
    }
}

void BoundaryConstraints::setGaussianConstraints()
{
    double A = 1.0;
    double B = 2.5;
    double S = 0.0;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        S += geometry->length(h.edge());
    }
    
    double s = 0.0;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        double theta = 2.0*M_PI*s/S;
        param->targetLengths[wIndices[wedge(h)]] = max(A*exp(-sqr(M_PI - theta)/B),
                                                       numeric_limits<double>::min());
        s += geometry->length(h.edge());
    }
}

void BoundaryConstraints::setIsometricConstraints()
{
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!setConstraintsOnCuts && h.isReal()) continue;
        param->targetLengths[wIndices[wedge(h)]] = geometry->length(h.edge());
    }
}
