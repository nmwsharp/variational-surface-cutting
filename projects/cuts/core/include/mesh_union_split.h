#pragma once

#include <geometry.h>

#include <vector>

Geometry<Euclidean>* unionMeshes(std::vector<Geometry<Euclidean>*> geometries);
Geometry<Euclidean>* unionMeshesDuplicateEdges(std::vector<Geometry<Euclidean>*> geometries);
Geometry<Euclidean>* unionMeshesDuplicateEdgesAndCuts(std::vector<Geometry<Euclidean>*> geometries);
Geometry<Euclidean>* newGeomDuplicateEdgesAndCuts(Geometry<Euclidean>* inputGeometry);
std::vector<Geometry<Euclidean>*> splitAlongCuts(Geometry<Euclidean>* inputGeometry);

Geometry<Euclidean>* improveMeshPreserveCut(Geometry<Euclidean>* geom); 
