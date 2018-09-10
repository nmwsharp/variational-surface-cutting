// Author: Rohan Sawhney
// Fall 2016

#include "subdivision.h"

Subdivision::Subdivision(Geometry<Euclidean> *geometry_):
geometry(geometry_),
mesh(geometry->getMesh())
{
    
}

Vector3 Subdivision::updateVertexPosition(VertexPtr v)
{
    Vector3 coordinate;
    if (v.isBoundary()) {
        coordinate = 3.0*geometry->position(v)/4.0;
        for (HalfedgePtr h: v.incomingHalfedges()) {
            if (!h.isReal()) {
                coordinate += (geometry->position(h.vertex()) +
                               geometry->position(h.next().next().vertex()))/8.0;
                break;
            }
        }
        
    } else {
        int n = 0;
        Vector3 sum{0, 0, 0};
        for (HalfedgePtr h: v.incomingHalfedges()) {
            sum += geometry->position(h.vertex());
            n++;
        }
        
        double beta = n == 3 ? 3.0/16.0 : 3.0/(8.0*n);
        coordinate = (1 - n*beta)*geometry->position(v) + beta*sum;
    }
    
    return coordinate;
}

Vector3 Subdivision::createNewVertex(HalfedgePtr he)
{
    Vector3 vertex;
    if (he.twin().isReal()) {
        vertex = (geometry->position(he.vertex())*3.0 +
                  geometry->position(he.twin().vertex())*3.0 +
                  geometry->position(he.next().next().vertex()) +
                  geometry->position(he.twin().next().next().vertex()))/8.0;
        
    } else {
        vertex = (geometry->position(he.vertex()) +
                  geometry->position(he.twin().vertex()))/2.0;
    }
    
    return vertex;
}

void Subdivision::createNewFaces(const vector<unsigned int>& indices)
{
    vector<size_t> face;
    face.push_back(indices[0]);
    face.push_back(indices[1]);
    face.push_back(indices[5]);
    polygonSoup.polygons.push_back(face);
    
    face.clear();
    face.push_back(indices[5]);
    face.push_back(indices[1]);
    face.push_back(indices[3]);
    polygonSoup.polygons.push_back(face);
    
    face.clear();
    face.push_back(indices[3]);
    face.push_back(indices[1]);
    face.push_back(indices[2]);
    polygonSoup.polygons.push_back(face);
    
    face.clear();
    face.push_back(indices[5]);
    face.push_back(indices[3]);
    face.push_back(indices[4]);
    polygonSoup.polygons.push_back(face);
}

Geometry<Euclidean>* Subdivision::subdivide()
{
    EdgeData<unsigned int> visited(mesh);
    EdgeData<unsigned int> edgeVertexMap(mesh);
    
    // Insert updated position of existing vertices
    for (VertexPtr v: mesh->vertices()) {
        polygonSoup.vertexCoordinates.push_back(updateVertexPosition(v));
    }
    
    // Create new vertices and faces
    VertexData<size_t> vIndices = mesh->getVertexIndices();
    polygonSoup.polygons.reserve(4*mesh->nFaces());
    for (FacePtr f: mesh->faces()) {
        
        vector<unsigned int> indices;
        for (HalfedgePtr h: f.adjacentHalfedges()) {
            indices.push_back(vIndices[h.vertex()]);
            
            EdgePtr e = h.edge();
            if (!visited[e]) {
                edgeVertexMap[e] = (unsigned int)polygonSoup.vertexCoordinates.size();
                polygonSoup.vertexCoordinates.push_back(createNewVertex(h));
                
                visited[e] = true;
            }
            indices.push_back(edgeVertexMap[e]);
        }
        createNewFaces(indices);
    }
   
    // Create new mesh and geometry
    Geometry<Euclidean> *subdividedGeometry = nullptr;
    new HalfedgeMesh(polygonSoup, subdividedGeometry);
    return subdividedGeometry;
}
