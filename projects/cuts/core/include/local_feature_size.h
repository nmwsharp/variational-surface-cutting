#pragma once

#include <vector>

#include "geometry.h"
#include "dense_matrix.h"
#include "sparse_matrix.h"

/*

    Compute the "local feature size", a measure of the length-scale of the geometry at each point.

*/

// Use Laplacian smoothing to smooth out the feature size.
// This is simple, but has a downside that regions with a small feature size get "dampened" by nearby
// regions with larger feature size.
VertexData<double> computeLocalFeatureSize_smooth(Geometry<Euclidean>* geometry, GC::SparseMatrix<double>& zeroFormLaplacian, double smoothing);

// Use an FMM-like expanding search to smooth out feature size.
// This is a bit less straightforward, but has the nice property that the smoothed value is always
// no less than the pointwise value.
VertexData<double> computeLocalFeatureSize_eikonal(Geometry<Euclidean>* geometry);