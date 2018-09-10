#pragma once

#include "geometry.h"


std::vector<Vector3> getIcosphereVertices(Vector3 center, double radius, size_t atLeastNPoints);

VertexData<double> computeVisibilityFraction(Geometry<Euclidean>* geometry, size_t nSpherePointsMin = 500, bool headOn = false);