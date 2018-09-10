// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "geometry.h"
#include "polygon_soup_mesh.h"

class Subdivision {
public:
	// Constructor
	Subdivision(Geometry<Euclidean> *geometry_);

	// Implements loop subdivision
	Geometry<Euclidean>* subdivide();

protected:
    // Updates current vertices
    Vector3 updateVertexPosition(VertexPtr v);
    
    // Creates a new vertex
    Vector3 createNewVertex(HalfedgePtr he);
    
    // Creates new faces
    void createNewFaces(const vector<unsigned int>& indices);
    
    // Member variables
    Geometry<Euclidean> *geometry;
    HalfedgeMesh *mesh;
    PolygonSoupMesh polygonSoup;
};
