// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "geometry.h"
#include <stack>
#include <unordered_map>

using namespace placeholders;

// Helper function to return a boundary halfedge's previous halfedge
inline HalfedgePtr previous(HalfedgePtr h)
{
    bool noCut = true;
    HalfedgePtr he = h.twin().next();
    do {
        if (he.edge().isCut()) {
            noCut = false;
            break;
        }
        he = he.twin().next();
    } while (he.isReal());
    
    return noCut ? he.prev() : he.twin();
}

// Helper function to return the wedge on a boundary halfedge's tail 
inline CornerPtr wedge(HalfedgePtr h)
{
    return h.twin().prev().corner();
}

// Helper function to determine if vertex is on the cut boundary
inline bool onCutBoundary(VertexPtr v)
{
    for (HalfedgePtr h: v.outgoingHalfedges()) {
        EdgePtr e = h.edge();
        if (e.isBoundary() || e.isCut()) return true;
    }
    
    return false;
}

// Implementation of Clausen integral
inline double Cl2(double x)
{
    if (x == 0.0) return 0.0;
    x = remainder(x, 2*M_PI);
    if (x == 0.0) return 0.0;
    
    if (fabs(x) <= 2.0944) {
        double xx = x * x;
        return ((((((((((((2.3257441143020875e-22 * xx
                           + 1.0887357368300848e-20) * xx
                           + 5.178258806090624e-19) * xx
                           + 2.5105444608999545e-17) * xx
                           + 1.2462059912950672e-15) * xx
                           + 6.372636443183181e-14) * xx
                           + 3.387301370953521e-12) * xx
                           + 1.8978869988971e-10) * xx
                           + 1.1482216343327455e-8) * xx
                           + 7.873519778281683e-7) * xx
                           + 0.00006944444444444444) * xx
                           + 0.013888888888888888) * xx
                           - log(fabs(x)) + 1.0) * x;
    }
    
    x += ((x > 0.0) ? - M_PI : M_PI);
    double xx = x * x;
    return ((((((((((((3.901950904063069e-15 * xx
                       + 4.566487567193635e-14) * xx
                       + 5.429792727596476e-13) * xx
                       + 6.5812165661369675e-12) * xx
                       + 8.167010963952222e-11) * xx
                       + 1.0440290284867003e-9) * xx
                       + 1.3870999114054669e-8) * xx
                       + 1.941538399871733e-7) * xx
                       + 2.927965167548501e-6) * xx
                       + 0.0000496031746031746) * xx
                       + 0.0010416666666666667) * xx
                       + 0.041666666666666664) * xx
                       + log(0.5)) * x;
}

namespace UnconstraintedMinimization {
    struct Handle {
        // Typedefs
        typedef function<void(double&, const GC::DenseMatrix<double>&)> ComputeEnergy;
        typedef function<void(GC::DenseMatrix<double>&, const GC::DenseMatrix<double>&)> ComputeGradient;
        typedef function<void(GC::SparseMatrix<double>&, const GC::DenseMatrix<double>&)> ComputeHessian;
        
        // Constructor
        Handle() {}
        
        // Member variables
        ComputeEnergy computeEnergy;
        ComputeGradient computeGradient;
        ComputeHessian computeHessian;
    };
    
    // Implements newton's method
    bool solve(GC::DenseMatrix<double>& x, Handle& handle);
}

namespace Layout {
    class AngleLengthLayout {
    public:
        // Contructor
        AngleLengthLayout(HalfedgeMesh *mesh_, const CornerData<size_t>& wIndices_,
                          const CornerData<double>& angles_, const EdgeData<double>& lengths_);
        
        // Perform layout
        void performLayout(vector<Vector2>& x);
        
    protected:
        // Performs face layout
        void performFaceLayout(HalfedgePtr h, const Vector2& dir, vector<Vector2>& x,
                               unordered_map<int, bool>& visited, stack<EdgePtr>& stack);
        
        // Member variables
        HalfedgeMesh *mesh;
        FaceData<size_t> fIndices;
        const CornerData<size_t>& wIndices;
        const CornerData<double>& angles;
        const EdgeData<double>& lengths;
    };
}
