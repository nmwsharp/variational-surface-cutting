// Author: Rohan Sawhney
// Fall 2016

#include "parameterization.h"

Parameterization::Parameterization(Geometry<Euclidean> *geometry_):
geometry(geometry_),
mesh(geometry->getMesh()),
uvs(mesh),
vIndices(mesh->getVertexIndices()),
fIndices(mesh->getFaceIndices()),
wIndices(mesh),
qcDistortion(mesh),
areaDistortion(mesh),
angleDistortion(mesh),
angleSumDistortion(mesh),
clDistortion(mesh),
wedgeCount(0),
interiorWedgeCount(0),
boundaryWedgeCount(0)
{
    if (mesh->nBoundaryLoops() != 0 || mesh->hasCut()) assignWedgeIndices(wIndices);
}

void Parameterization::assignWedgeIndices(CornerData<size_t>& indices,
                                          bool separateBoundaries)
{

    interiorWedgeCount = 0;
    wedgeCount = 0;
    for (VertexPtr v: mesh->vertices()) {
        // Grab an outgoing interior halfedge
        HalfedgePtr he = v.halfedge();
        while (!he.isReal()) he = he.twin().next();
        
        // Count interior cut edges and set he to a cut edge if it exists
        int cuts = 0;
        for (HalfedgePtr h: v.outgoingHalfedges()) {
            EdgePtr e = h.edge();
            if (!e.isBoundary() && e.isCut()) {
                he = h;
                cuts++;
            }
        }
        
        // Loop over outgoing halfedges, starting from he and set corner indices, incrementing
        // wedgeCount every time everytime a cut is encountered
        int cut = 0;
        HalfedgePtr h = he;
        do {
            if (h.isReal()) {
                if (separateBoundaries) {
                    if (!v.isBoundary() && cuts == 0) indices[h.next().corner()] = interiorWedgeCount;
                        
                } else {
                    indices[h.next().corner()] = wedgeCount;
                }
            }
            
            // Visit halfedges in CCW order, h.twin().next() visits halfedges in CW order
            h = h.prev().twin();
            if (h.isReal() && h.edge().isCut()) {
                // The check ensures that wedgeCount isn't incremented twice, once in here on
                // the "last" cut in the one ring, and once again in the outer loop
                cut++;
                if (cut < cuts) wedgeCount++;
                
            } else if (v.isBoundary() && !h.isReal() && cuts > 0) {
                wedgeCount++;
            }
            
        } while (h != he);
        
        if (!v.isBoundary() && cuts == 0) interiorWedgeCount++; // Interior vertex
        wedgeCount++;
    }

    // Assign boundary indices
    boundaryWedgeCount = 0;
    // int nB = mesh->nBoundaryLoops() == 0 && mesh->hasCut() ? 1 : (int)mesh->nBoundaryLoops();
    int nB = 1;
    for (int b = 0; b < nB; b++) {
        for (HalfedgePtr h: mesh->cutBoundary(b)) {
            if (separateBoundaries) {
                HalfedgePtr he = h.twin().next();
                do {
                    indices[he.next().corner()] = interiorWedgeCount + boundaryWedgeCount;
                    
                    if (he.edge().isCut()) break;
                    he = he.twin().next();
                } while (he.isReal());
            }
            
            boundaryWedgeCount++;
        }
    }
}

void Parameterization::pinWedges()
{
    pinnedWedges.resize(2);
    pinnedPositions = {Vector2{0, 0}, Vector2{1, 0}};
    
    // Pin diameter vertices on longest boundary loop
    double max = 0.0;
    for (HalfedgePtr h1: mesh->cutBoundary()) {
        const Vector3& p1(geometry->position(h1.vertex()));
        for (HalfedgePtr h2: mesh->cutBoundary()) {
            const Vector3& p2(geometry->position(h2.vertex()));
            double l = norm2(p2-p1);
            
            if (l > max) {
                max = l;
                pinnedWedges[0] = wIndices[wedge(h1)];
                pinnedWedges[1] = wIndices[wedge(h2)];
            }
        }
    }
}

bool Parameterization::isPinnedWedge(const unsigned int& wIndex, unsigned int& shift,
                                     Vector2& pinnedPosition) const
{
    shift = 0;
    for (size_t i = 0; i < pinnedWedges.size(); i++) {
        if (pinnedWedges[i] == wIndex) {
            pinnedPosition = pinnedPositions[i];
            return true;
        }
        
        if (pinnedWedges[i] < wIndex) shift += 1;
    }
    
    return false;
}

double uvArea(FacePtr& f, const CornerData<Vector2>& uvs)
{
    Vector3 AN{0, 0, 0};
    for (CornerPtr c: f.adjacentCorners()) {
        AN += cross(uvs[c], uvs[c.next()]);
    }
    
    return norm(AN)/2.0;
}

Vector2 uvBarycenter(FacePtr& f, const CornerData<Vector2>& uvs)
{
    Vector2 barycenter{0, 0};
    for (CornerPtr c: f.adjacentCorners()) {
        barycenter += uvs[c];
    }
    
    return barycenter/3.0;
}

void Parameterization::normalize()
{
    // Compute center
    double totalArea = 0.0;
    Vector2 center{0, 0};
    for (FacePtr f: mesh->faces()) {
        double area = uvArea(f, uvs);
        
        center += area*uvBarycenter(f, uvs);
        totalArea += area;
    }
    center /= totalArea;
    
    // Shift 
    double r = 0.0;
    for (CornerPtr c: mesh->corners()) {
        uvs[c] -= center;
        r = max(r, norm2(uvs[c]));
    }

    // Scale
    r = sqrt(r);
    for (CornerPtr c: mesh->corners()) {
        uvs[c] /= r;
    }
    
    // Anchor uvs
    HalfedgePtr h = mesh->halfedge(0);
    Vector2 e1 = uvs[h.next().corner()] - uvs[h.corner()]; e1.normalize();
    Vector2 e2{-e1.y, e1.x};
    
    for (CornerPtr c: mesh->corners()) {
        Vector2& uv = uvs[c];
        uv = Vector2{dot(e1, uv), dot(e2, uv)};
    }
}
