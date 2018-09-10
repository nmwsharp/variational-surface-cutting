// Author: Rohan Sawhney
// Fall 2016

#pragma once

#include "flatten_utils.h"

class Parameterization {
public:
	// Constructor
	Parameterization(Geometry<Euclidean> *geometry_);

    // Destructor
    virtual ~Parameterization() {}
    
	// Flatten
	virtual bool flatten() = 0;
    
    // Assign wedge indices
    void assignWedgeIndices(CornerData<size_t>& indices, bool separateBoundaries = false);
    
    // Member variables
    Geometry<Euclidean> *geometry;
    HalfedgeMesh *mesh;
    CornerData<Vector2> uvs;
    vector<size_t> pinnedWedges;
    
    VertexData<size_t> vIndices;
    FaceData<size_t> fIndices;
    CornerData<size_t> wIndices; // wedge indices
    
    unordered_map<unsigned int, double> targetCurvatures;
    unordered_map<unsigned int, double> targetLengths;
    unordered_map<unsigned int, double> coneAngles;
    
    FaceData<Vector3> qcDistortion;
    FaceData<Vector3> areaDistortion;
    CornerData<Vector3> angleDistortion;
    EdgeData<Vector3> angleSumDistortion;
    EdgeData<Vector3> clDistortion;

protected:
    // Pin wedges
    void pinWedges();
    
    // Checks if wedge is pinned
    bool isPinnedWedge(const unsigned int& wIndex, unsigned int& shift, Vector2& pinnedPosition) const;
    
    // Normalize
    void normalize();

    // Member variables
    vector<Vector2> pinnedPositions;
    int wedgeCount;
    int interiorWedgeCount;
    int boundaryWedgeCount;
};
