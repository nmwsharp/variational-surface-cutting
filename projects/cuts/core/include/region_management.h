#pragma once

#include "geometry.h"
#include "detect_symmetry.h"

#include <vector>

// Utility functions for processing regions

double cutRegionToDisk(Geometry<Euclidean>* regionGeom, SymmetryResult* symmetry); // uses setCut() to mark cuts, returns length of cut made

void canonicalizeParameterization(Geometry<Euclidean>* geom, CornerData<Vector2>& param);

void packParameterization(Geometry<Euclidean>* geom, CornerData<Vector2>& param, double marginFactor=0.01); 

FaceData<int> computeFaceComponents(Geometry<Euclidean>* geom, const CornerData<Vector2>& param);

Geometry<Euclidean>* makePatchExtrinsicDevelopable(Geometry<Euclidean>* geom); 
Geometry<Euclidean>* makeShapeExtrinsicDevelopable(Geometry<Euclidean>* geom, EdgeData<double> targetL, bool stretchOnly = false); 

VertexData<double> yamabeScaleFactors(Geometry<Euclidean>* geom);
VertexData<double> stretchOnlyAdjustmentFactors(Geometry<Euclidean>* geom);

std::vector<std::array<Vector3, 2>> getGeodesicCuts(Geometry<Euclidean>* geom);
std::vector<std::vector<HalfedgePtr>> getDirectedCutEdges(HalfedgeMesh* mesh);
Geometry<Euclidean>* remeshAlongCutGeodesics(Geometry<Euclidean>* geom); 

void computeCylinderParameterization(Geometry<Euclidean>* geom, CornerData<Vector2>& param);

void writeBoundarySVG(std::string filename, Geometry<Euclidean>* geom, CornerData<Vector2>& param);
