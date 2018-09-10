#pragma once

#include "geometry.h"
#include "eulerian_shape_optimizer.h"

// Initialization strategies
VertexData<LabelVec> randomMSDF(Geometry<Euclidean>* geometry);
VertexData<LabelVec> singleVertex(Geometry<Euclidean>* geometry, int startingVert = -1);
VertexData<LabelVec> xAxisSplitMSDF(Geometry<Euclidean>* geometry);
VertexData<LabelVec> normalClusterMSDF(Geometry<Euclidean>* geometry, double angleThreshold=PI/4);
VertexData<LabelVec> distanceClusterMSDF(Geometry<Euclidean>* geometry, int nCluster=10);
VertexData<LabelVec> meanCurvatureLevelSetMSDF(Geometry<Euclidean>* geometry, double scale=1.0);
VertexData<LabelVec> dChartsMSDF(Geometry<Euclidean>* geometry);
VertexData<LabelVec> cubeMSDF(Geometry<Euclidean>* geometry);
VertexData<LabelVec> cubeCircleMSDF(Geometry<Euclidean>* geometry);
VertexData<LabelVec> sphereAtMSDF(Geometry<Euclidean>* geometry, Vector3 center, double rad);
VertexData<LabelVec> volleyballMSDF(Geometry<Euclidean>* geometry);

VertexData<LabelVec> transferMSDF(Geometry<Euclidean>* sourceGeometry, VertexData<LabelVec>& sourceMSDF, Geometry<Euclidean>* targetGeometry);

// Helpers
VertexData<LabelVec> msdfFromLabels(Geometry<Euclidean>* geometry, VertexData<int> labels);
VertexData<int> colorWithLabels(HalfedgeMesh* mesh, VertexData<int> initLabels);
