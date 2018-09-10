#include "mesh_union_split.h"

#include "geometry.h"
#include "polygon_soup_mesh.h"
#include "mutable_triangle_cut_mesh.h"
#include "combining_hash_functions.h"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <tuple>

Geometry<Euclidean>* unionMeshes(std::vector<Geometry<Euclidean>*> geometries) {

    std::unordered_map<Vector3, size_t> vertInd;
    std::vector<Vector3> vertPos;
    std::vector<std::vector<size_t>> faces; 
    size_t nVert = 0;

    // Build a map of cut edges so we can restore them later
    typedef std::pair<Vector3, Vector3> PointPair;
    std::unordered_set<PointPair> cutEdges;
    for(auto geom : geometries) {
        HalfedgeMesh* mesh = geom->getMesh();
        for(EdgePtr e : mesh->edges()) {

            // Add the edge to the lookup set of cuts
            if(e.isCut() || e.isBoundary()) {
                HalfedgePtr he = e.halfedge();
                PointPair p = std::make_pair(geom->position(he.vertex()), geom->position(he.twin().vertex()));
                cutEdges.insert(p);
            }

            // If it's a cut, we need to add it with both orders to ensure we can find it later
            if(e.isCut()) {
                HalfedgePtr he = e.halfedge().twin();
                PointPair p = std::make_pair(geom->position(he.vertex()), geom->position(he.twin().vertex()));
                cutEdges.insert(p);
            }
        }
    }

    // Index vertices
    for(auto geom : geometries) {
        HalfedgeMesh* mesh = geom->getMesh();
        for(VertexPtr v : mesh->vertices()) {
            Vector3 pos = geom->position(v);
            if(vertInd.find(pos) == vertInd.end()) {
                vertInd[pos] = nVert++;
                vertPos.push_back(pos);
            }
        }
    }

    // Build faces
    for(auto geom : geometries) {
        HalfedgeMesh* mesh = geom->getMesh();

        for(FacePtr f : mesh->faces()) {

            std::vector<size_t> face;

            for(VertexPtr v : f.adjacentVertices()) {
                face.push_back(vertInd[geom->position(v)]);
            }

            faces.push_back(face);
        }
    }

    PolygonSoupMesh pSoup(faces, vertPos);

    Geometry<Euclidean>* newGeom;
    HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, newGeom);


    // Mark cut edges
    for(EdgePtr e : newMesh->edges()) {
        if(e.isBoundary()) continue;
       
        // Is this a cut edge?
        PointPair p = std::make_pair(newGeom->position(e.halfedge().vertex()), newGeom->position(e.halfedge().twin().vertex()));
        if(cutEdges.find(p) != cutEdges.end()) {
            e.markCut(true);
        } 
    }

    return newGeom;
} 

Geometry<Euclidean>* unionMeshesDuplicateEdges(std::vector<Geometry<Euclidean>*> geometries) {

    std::vector<Vector3> vertPos;
    std::vector<std::vector<size_t>> faces; 
    size_t nVert = 0;

    for(auto geom : geometries) {

        // Index vertices
        size_t vertOffset = vertPos.size();
        HalfedgeMesh* mesh = geom->getMesh();
        for(VertexPtr v : mesh->vertices()) {
            Vector3 pos = geom->position(v);
            vertPos.push_back(pos);
        }

        // Build faces
        VertexData<size_t> vInd = mesh->getVertexIndices();
        for(FacePtr f : mesh->faces()) {

            std::vector<size_t> face;

            for(VertexPtr v : f.adjacentVertices()) {
                face.push_back(vInd[v] + vertOffset);
            }

            faces.push_back(face);
        }
    }

    PolygonSoupMesh pSoup(faces, vertPos);

    Geometry<Euclidean>* newGeom;
    HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, newGeom);

    return newGeom;
}

Geometry<Euclidean>* unionMeshesDuplicateEdgesAndCuts(std::vector<Geometry<Euclidean>*> geometries) {

    std::vector<Vector3> vertPos;
    std::vector<std::vector<size_t>> faces; 
    size_t nVert = 0;
    size_t INVALID_IND = std::numeric_limits<size_t>::max();

    for(auto geom : geometries) {


        // Index vertices
        size_t vertOffset = vertPos.size();
        HalfedgeMesh* mesh = geom->getMesh();
        HalfedgeData<size_t> distinctVertInd(mesh);
        for(VertexPtr v : mesh->vertices()) {

            // Check if this vertex is adjacent to any cuts
            bool hasCut = false;
            for(EdgePtr e : v.adjacentEdges()) {
                if(e.isCut()) {
                    hasCut = true;
                }
            }

            // Hard case where vertex has adjacent cuts
            if(hasCut) {
               
                // Find any cut edge to start with
                HalfedgePtr heStart;
                for(HalfedgePtr he : v.outgoingHalfedges()) {
                    if(he.edge().isCut()) {
                        heStart = he;
                        break;
                    }
                }

                // Orbit vertex, creating new entries on all sides of cut
                HalfedgePtr currHe = heStart;
                size_t currInd = INVALID_IND;
                do {

                    // If we're about to cross a cut or a boundary, need a new vertex entry
                    if(currHe.edge().isCut() || !currHe.twin().isReal()) {
                        currInd = INVALID_IND; 
                    }

                    // Cross
                    currHe = currHe.twin().next();

                    // Create a new vertex, if needed
                    if(currInd == INVALID_IND) {
                        Vector3 pos = geom->position(v);
                        vertPos.push_back(pos);
                        currInd = vertPos.size()-1; 
                    }

                    // Label this face-vertex
                    distinctVertInd[currHe] = currInd;

                } while (currHe != heStart);

            }
            // Easy case where vertex has no adjacent cuts
            else {
    
                Vector3 pos = geom->position(v);
                vertPos.push_back(pos);

                for(HalfedgePtr he : v.outgoingHalfedges()) {
                    distinctVertInd[he] = vertPos.size()-1;
                }
            }
        }

        // Build faces
        VertexData<size_t> vInd = mesh->getVertexIndices();
        for(FacePtr f : mesh->faces()) {

            std::vector<size_t> face;

            for(HalfedgePtr he : f.adjacentHalfedges()) {
                face.push_back(distinctVertInd[he]);
            }

            faces.push_back(face);
        }
    }

    PolygonSoupMesh pSoup(faces, vertPos);

    Geometry<Euclidean>* newGeom;
    HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, newGeom);

    return newGeom;
}


Geometry<Euclidean>* newGeomDuplicateEdgesAndCuts(Geometry<Euclidean>* inputGeometry) {

    std::vector<Vector3> vertPos;
    std::vector<std::vector<size_t>> faces; 
    size_t nVert = 0;
    size_t INVALID_IND = std::numeric_limits<size_t>::max();

    // Index vertices
    size_t vertOffset = vertPos.size();
    HalfedgeMesh* mesh = inputGeometry->getMesh();
    HalfedgeData<size_t> distinctVertInd(mesh);
    for(VertexPtr v : mesh->vertices()) {

        // Check if this vertex is adjacent to any cuts
        bool hasCut = false;
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) {
                hasCut = true;
            }
        }

        // Hard case where vertex has adjacent cuts
        if(hasCut) {
            
            // Find any cut edge to start with
            HalfedgePtr heStart;
            for(HalfedgePtr he : v.outgoingHalfedges()) {
                if(he.edge().isCut()) {
                    heStart = he;
                    break;
                }
            }

            // Orbit vertex, creating new entries on all sides of cut
            HalfedgePtr currHe = heStart;
            size_t currInd = INVALID_IND;
            do {

                // If we're about to cross a cut or a boundary, need a new vertex entry
                if(currHe.edge().isCut() || !currHe.twin().isReal()) {
                    currInd = INVALID_IND; 
                }

                // Cross
                currHe = currHe.twin().next();

                // Create a new vertex, if needed
                if(currInd == INVALID_IND) {
                    Vector3 pos = inputGeometry->position(v);
                    vertPos.push_back(pos);
                    currInd = vertPos.size()-1; 
                }

                // Label this face-vertex
                distinctVertInd[currHe] = currInd;

            } while (currHe != heStart);

        }
        // Easy case where vertex has no adjacent cuts
        else {

            Vector3 pos = inputGeometry->position(v);
            vertPos.push_back(pos);

            for(HalfedgePtr he : v.outgoingHalfedges()) {
                distinctVertInd[he] = vertPos.size()-1;
            }
        }
    }

    // Build faces
    for(FacePtr f : mesh->faces()) {

        std::vector<size_t> face;

        for(HalfedgePtr he : f.adjacentHalfedges()) {
            face.push_back(distinctVertInd[he]);
        }

        faces.push_back(face);
    }

    PolygonSoupMesh pSoup(faces, vertPos);

    Geometry<Euclidean>* newGeom;
    HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, newGeom);

    return newGeom;
}

std::vector<Geometry<Euclidean>*> splitAlongCuts(Geometry<Euclidean>* inputGeometry) {

    using std::vector;

    cout << endl << endl << "Splitting along cuts" << endl;

    HalfedgeMesh* mesh = inputGeometry->getMesh();
    size_t INVALID_IND = std::numeric_limits<size_t>::max();

    // Assign faces to components via greedy expansion
    FaceData<size_t> faceComponentIndex(mesh, INVALID_IND);
    size_t nComponents = 0;
    for(FacePtr f : mesh->faces()) {

        // Skip already assigned faces
        if(faceComponentIndex[f] != INVALID_IND) continue;

        // Start a new component
        std::vector<FacePtr> facesToCheck = {f};
        size_t iComponent = nComponents;
        nComponents++;
        while(!facesToCheck.empty()) {

            // Grab next from queue
            FacePtr currF = facesToCheck.back();
            facesToCheck.pop_back();

            // Skip if already assigned
            if(faceComponentIndex[currF] != INVALID_IND) {
                continue;
            }

            // Assign this face
            faceComponentIndex[currF] = iComponent;

            // Add neighbors to search
            for(HalfedgePtr he : currF.adjacentHalfedges()) {

                // Can't cross cuts or boundaries
                if(he.edge().isBoundary() || he.edge().isCut()) continue;

                // Add to search
                FacePtr neighF = he.twin().face();
                if(faceComponentIndex[neighF] == INVALID_IND) {
                    facesToCheck.push_back(neighF);
                }
            }
        }
    }

    cout << "Found " << nComponents << " components while splitting along cuts" << endl;

    // Create lists of vertices and faces for each component
    vector<vector<Vector3>> componentVertPos(nComponents);
    vector<vector<vector<size_t>>> componentFaces(nComponents); 

    // Index vertices 
    HalfedgeData<size_t> distinctVertInd(mesh);
    for(VertexPtr v : mesh->vertices()) {

        // Check if this vertex is adjacent to any cuts
        bool hasCut = false;
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) {
                hasCut = true;
            }
        }

        // Hard case where vertex has adjacent cuts
        if(hasCut) {
            
            // Find any cut edge to start with
            HalfedgePtr heStart;
            for(HalfedgePtr he : v.outgoingHalfedges()) {
                if(he.edge().isCut()) {
                    heStart = he;
                    break;
                }
            }

            // Orbit vertex, creating new entries on all sides of cut
            HalfedgePtr currHe = heStart;
            size_t currInd = INVALID_IND;
            do {

                // If we're about to cross a cut or a boundary, need a new vertex entry
                if(currHe.edge().isCut() || !currHe.twin().isReal()) {
                    currInd = INVALID_IND; 
                }

                // Cross
                currHe = currHe.twin().next();

                // Create a new vertex, if needed
                if(currHe.isReal() && currInd == INVALID_IND) {
                    size_t iComp = faceComponentIndex[currHe.face()];
                    Vector3 pos = inputGeometry->position(v);
                    componentVertPos[iComp].push_back(pos);
                    currInd = componentVertPos[iComp].size()-1; 
                }

                // Label this face-vertex
                distinctVertInd[currHe] = currInd;

            } while (currHe != heStart);

        }
        // Easy case where vertex has no adjacent cuts
        else {

            size_t iComp = faceComponentIndex[v.halfedge().face()];

            Vector3 pos = inputGeometry->position(v);
            componentVertPos[iComp].push_back(pos);

            for(HalfedgePtr he : v.outgoingHalfedges()) {
                distinctVertInd[he] = componentVertPos[iComp].size()-1;
            }
        }
    }

    // Build faces
    for(FacePtr f : mesh->faces()) {

        std::vector<size_t> face;

        for(HalfedgePtr he : f.adjacentHalfedges()) {
            face.push_back(distinctVertInd[he]);
        }

        componentFaces[faceComponentIndex[f]].push_back(face);
    }


    // Build the output meshes
    std::vector<Geometry<Euclidean>*> outGeoms;
    for(size_t iComp = 0; iComp < nComponents; iComp++) {
        PolygonSoupMesh pSoup(componentFaces[iComp], componentVertPos[iComp]);
    
        Geometry<Euclidean>* newGeom;
        HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, newGeom);
        outGeoms.push_back(newGeom);
    }
    
    return outGeoms;
}


Geometry<Euclidean>* improveMeshPreserveCut(Geometry<Euclidean>* geom) {

    // Create mutable mesh
    MutableCutMesh mutMesh(geom); 

    // Remesh
    double angleThresh = 0.0872665; // 5 degrees
    mutMesh.skinnyThresh = -1;
    mutMesh.normalThreshold = PI/2.0;
    int nCollapse = 0;
    int nFlip = 0;
    for(int iPass = 0; iPass < 5; iPass++) {


        // === Collapse
        // Edges to possibly collapse        
        std::vector<MutableCutHalfedge*> candidates(mutMesh.halfedges.begin(), mutMesh.halfedges.end());
        for(MutableCutHalfedge* he : candidates) {
    
            // Check if the halfedge still exists
            if(mutMesh.halfedges.find(he) == mutMesh.halfedges.end()) continue;

            // Skip boundary halfedges
            if(he->isOnBoundary()) continue;
            
            MutableCutVertex* vA = he->vertex;
            MutableCutVertex* vB = he->next->vertex;
            MutableCutVertex* vC = he->next->next->vertex;

            // Only collapse around cut
            if(vA->isOnCut() || vB->isOnCut() || vC->isOnCut()) {

                // Collapse skinny triangles
                double oppAngle = he->oppositeAngle();
                if(oppAngle < angleThresh) {
                    bool success = mutMesh.tryCollapseEdge(he);
                    if(success) {
                        nCollapse++;
                    }
                }

            }

        }

        // === Flip
        for(MutableCutHalfedge* he : mutMesh.halfedges) {
    
            // Can't flip boundary edges
            if(he->isOnBoundary()) continue;
            MutableCutVertex* vA = he->vertex;
            MutableCutVertex* vB = he->next->vertex;
            MutableCutVertex* vC = he->next->next->vertex;
            MutableCutVertex* vD = he->twin->next->next->vertex;
   
            
            // Only flip around cut
            if(vA->isOnCut() || vB->isOnCut() || vC->isOnCut() || vD->isOnCut()) {
    
                // Check Delaunay flipping criterion
                double thetaC = angle(vA->position - vC->position, vB->position - vC->position);
                double thetaD = angle(vB->position - vD->position, vA->position - vD->position);
        
                if(thetaC + thetaD > PI) {
                    bool flipped = mutMesh.tryFlipEdge(he);
                    if(flipped) {
                        nFlip++;
                    }
                }
    
            }
    
        }
                

    }

    cout << "Collapsed " << nCollapse << " edges in mesh." << endl;
    cout << "Flipped " << nFlip<< " edges in mesh." << endl;

    return mutMesh.toHalfedgeMesh();
}
