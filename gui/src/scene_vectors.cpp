#include<scene_vectors.h>

#include <iostream>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <limits>
 
#include <shaders.h>
#include <shaders/shiny_shaders.h>
#include <direction_fields.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

SceneVectors::SceneVectors(Viewer &parent_)
    :   SceneObject(parent_),
        vectorProgram(&SHINY_COLOR_VERT_SHADER, &SHINY_COLOR_FRAG_SHADER, DrawMode::IndexedTriangles)
{
}


void SceneVectors::draw() {

    if(!enabled) {
        return;
    }

    // Set uniforms
    glm::mat4 viewMat = parent.camera.getViewMatrix();
    vectorProgram.setUniform("u_viewMatrix", glm::value_ptr(viewMat));

    glm::mat4 projMat = parent.camera.getPerspectiveMatrix();
    vectorProgram.setUniform("u_projMatrix", glm::value_ptr(projMat));

    Vector3 eyePos = parent.camera.getCameraWorldPosition();
    vectorProgram.setUniform("u_eye", eyePos);

    Vector3 lightPos = parent.camera.getLightWorldPosition();
    vectorProgram.setUniform("u_light", lightPos);

    // Draw
    vectorProgram.draw();

}


void SceneVectors::setVectors(const std::vector<Vector3> &basePoints, const std::vector<Vector3> &vectors, double radius) {

    // TODO ways to make this faster if its slow:
    // - Use a geometry shader
    // - The basis and theta computations are unecessarily repeated, cache them
       
    double totalVectorLength = 0;
    for(const Vector3 &v : vectors) {
        totalVectorLength += norm(v);
    }
    double meanVectorLength = totalVectorLength / vectors.size();

    // Constants and parameters
    const unsigned int nSidesPerArrow = 8;
    const unsigned int nVertsPerArrow = 2*nSidesPerArrow + 1;
    const unsigned int nFacesPerArrow = 3*nSidesPerArrow;
    const double shaftFraction = 0.8;
    if(radius == -1.0) {
        radius = 0.1 * meanVectorLength;
    }
    const Vector3 arbitraryDirection = unit(Vector3 {0.12987391283, -.7089237918, .589712378934}); // This is simple and almost definitely works

    std::vector<Vector3> positions(nVertsPerArrow*vectors.size());
    std::vector<Vector3> normals(nVertsPerArrow*vectors.size());
    std::vector<Vector3> colors(nVertsPerArrow*vectors.size());
    std::vector<uint3> indices(nFacesPerArrow*vectors.size());
     

    // Walk the faces building arrays of positions/normals/colors
    unsigned int iPoint = 0;
    unsigned int iIndex = 0;
    for(unsigned int iArrow = 0; iArrow < vectors.size(); iArrow++) {

        // Useful shared values
        Vector3 v = vectors[iArrow];
        Vector3 tail = basePoints[iArrow];
        Vector3 tip = tail + v; 

        Vector3 basisY = unit(cross(v, arbitraryDirection));
        Vector3 basisX = unit(cross(basisY, v));

        double deltaTheta = 2*PI / nSidesPerArrow;

        // Compute the positions of the points we need for each arrow
        for(unsigned int iTheta = 0; iTheta < nSidesPerArrow; iTheta++) {
            double theta = iTheta * deltaTheta;
       
            Vector3 outwardVec = (cos(theta) * basisX + sin(theta) * basisY);
            Vector3 bottomPoint = tail + radius * outwardVec;
            Vector3 topPoint = bottomPoint + shaftFraction * v;

            positions[iPoint + 2*iTheta] = bottomPoint;
            positions[iPoint + 2*iTheta + 1] = topPoint;
            
            normals[iPoint + 2*iTheta] = outwardVec;
            normals[iPoint + 2*iTheta + 1] = outwardVec;
            
            colors[iPoint + 2*iTheta] = arrowColor;
            colors[iPoint + 2*iTheta + 1] = arrowColor;
        }

        // One extra point at the tip
        positions[iPoint + 2*nSidesPerArrow] = tip;
        normals[iPoint + 2*nSidesPerArrow] = unit(v);
        colors[iPoint + 2*nSidesPerArrow] = arrowColor;


        // Compute the indices needed for each arrow
        for(unsigned int iSide = 0; iSide < nSidesPerArrow; iSide++) {

            unsigned int nSide = (iSide+1) % nSidesPerArrow;

            // Lower face triangle
            indices[iIndex + 3*iSide + 0] = uint3{{iPoint + 2*iSide, iPoint + 2*nSide, iPoint + 2*iSide + 1}};

            // Upper face triangle
            indices[iIndex + 3*iSide + 1] = uint3{{iPoint + 2*nSide + 1, iPoint + 2*iSide + 1, iPoint + 2*nSide}};
            
            // Tip triangle
            indices[iIndex + 3*iSide + 2] = uint3{{iPoint + 2*iSide + 1, iPoint + 2*nSide + 1, iPoint + 2*nSidesPerArrow}};

        }


        iPoint += nVertsPerArrow;
        iIndex += nFacesPerArrow;
    }
       

    // Store the values in GL buffers
    vectorProgram.setAttribute("a_position", positions);
    vectorProgram.setAttribute("a_normal", normals);
    vectorProgram.setAttribute("a_color", colors);
    vectorProgram.setIndex(indices);

}


double SceneVectors::computeNiceMeshLengthScale(Geometry<Euclidean>* geometry) {
   HalfedgeMesh* mesh = geometry->getMesh();

    double totalEdgeLength = 0;
    for(EdgePtr e : mesh->edges()) {
        totalEdgeLength += norm(geometry->vector(e.halfedge()));
    }
    double meanEdgeLength = totalEdgeLength / mesh->nEdges();

    return 0.5 * meanEdgeLength;
}

void SceneVectors::setVectors(Geometry<Euclidean>* geometry, const VertexData<Vector3> &vertexVectors, bool autoscale) {
   HalfedgeMesh* mesh = geometry->getMesh();

    // Pick a nice scale for the data
    double scaleFactor = 1.0;
    if(autoscale) {
        double totalVectorLength = 0;
        for(VertexPtr v : mesh->vertices()) { 
            totalVectorLength += norm(vertexVectors[v]);
        }
        double meanVectorLength = totalVectorLength / mesh->nVertices();
        scaleFactor = computeNiceMeshLengthScale(geometry) / meanVectorLength;
    }

    VertexData<size_t> vertInd = mesh->getVertexIndices();
    std::vector<Vector3> bases(mesh->nVertices());
    std::vector<Vector3> scaledVectors(mesh->nVertices());
    for(VertexPtr v : mesh->vertices()) {
        bases[vertInd[v]] = geometry->position(v);
        scaledVectors[vertInd[v]] = vertexVectors[v] * scaleFactor;
    }


    setVectors(bases, scaledVectors);
}

void SceneVectors::setVectors(Geometry<Euclidean>* geometry, const FaceData<Vector3> &faceVectors, bool autoscale) {
   HalfedgeMesh* mesh = geometry->getMesh();


    // Pick a nice scale for the data
    double scaleFactor = 1.0;
    if(autoscale) {
        double totalVectorLength = 0;
        for(FacePtr f : mesh->faces()) {
            totalVectorLength += norm(faceVectors[f]);
        }
        double meanVectorLength = totalVectorLength / mesh->nFaces();
        scaleFactor = computeNiceMeshLengthScale(geometry) / meanVectorLength;
    }



    FaceData<size_t> faceInd = mesh->getFaceIndices();
    std::vector<Vector3> bases(mesh->nFaces());
    std::vector<Vector3> scaledVectors(mesh->nFaces());
    for(FacePtr f : mesh->faces()) {

        // Use the barycenter as the base position
        int n = 0;
        Vector3 barySum {0,0,0};
        for(VertexPtr v : f.adjacentVertices()) {
            barySum += geometry->position(v);
            n++;
        }
        Vector3 barycenter = barySum / n;

        bases[faceInd[f]] = barycenter; 
        scaledVectors[faceInd[f]] = faceVectors[f] * scaleFactor;
    }


    setVectors(bases, scaledVectors);
}

void SceneVectors::setSymmetricVectors(Geometry<Euclidean>* geometry, const VertexData<Complex> &directionFieldComplex, int nSym, bool autoscale) {
    HalfedgeMesh* mesh = geometry->getMesh();

    VertexData<Vector3> directionField = convertComplexDirectionsToR3Vectors(geometry, directionFieldComplex, nSym);

    setSymmetricVectors(geometry, directionField, nSym, autoscale);
}

void SceneVectors::setSymmetricVectors(Geometry<Euclidean>* geometry, const VertexData<Vector3> &vertexVectors, int nSym, bool autoscale) {
   HalfedgeMesh* mesh = geometry->getMesh();

    // Pick a nice scale for the data
    double scaleFactor = 1.0;
    if(autoscale) {
        double totalVectorLength = 0;
        for(VertexPtr v : mesh->vertices()) { 
            totalVectorLength += norm(vertexVectors[v]);
        }
        double meanVectorLength = totalVectorLength / mesh->nVertices();
        scaleFactor = computeNiceMeshLengthScale(geometry) / meanVectorLength;
    }

    VertexData<size_t> vertInd = mesh->getVertexIndices();
    std::vector<Vector3> bases(mesh->nVertices() * nSym);
    std::vector<Vector3> scaledVectors(mesh->nVertices() * nSym);
    for(VertexPtr v : mesh->vertices()) {

        // Project the vector to the tangent plane
        Vector3 tangentVector = geometry->projectToTangentSpace(v, vertexVectors[v]);

        for(int iSym = 0; iSym < nSym; iSym++) {

            size_t ind = nSym * vertInd[v] + iSym;

            // Rotate the vector in tangent space
            double angle = 2.0 * PI * iSym / (double) nSym;
            Vector3 rotVec = tangentVector.rotate_around(geometry->normal(v), angle);

            bases[ind] = geometry->position(v);
            scaledVectors[ind] = rotVec * scaleFactor;
        }
    }


    setVectors(bases, scaledVectors);
}

// Visualize vector data which is perpindicular to every halfedge
void SceneVectors::setEdgePerpVectors(Geometry<Euclidean>* geometry, const EdgeData<double> &edgeVectors, bool autoscale) {
   HalfedgeMesh* mesh = geometry->getMesh();


    // Pick a nice scale for the data
    double scaleFactor = 1.0;
    double radius = -1;
    if(autoscale) {
        double totalVectorLength = 0;
        double maxVectorLength = 0;
        for(EdgePtr e : mesh->edges()) {
            totalVectorLength += std::abs(edgeVectors[e]);
            maxVectorLength = std::max(maxVectorLength, std::abs(edgeVectors[e]));
        }
        double meanVectorLength = totalVectorLength / mesh->nEdges();
        double targetLen = 0.5 * (meanVectorLength + maxVectorLength);
        scaleFactor = 2.0 * computeNiceMeshLengthScale(geometry) / maxVectorLength;
        // scaleFactor = 0.3 * computeNiceMeshLengthScale(geometry) / meanVectorLength;
    }
    radius = 0.1 * computeNiceMeshLengthScale(geometry);



    std::vector<Vector3> bases(mesh->nEdges());
    std::vector<Vector3> scaledVectors(mesh->nEdges());
    size_t iVec = 0;
    for(EdgePtr e : mesh->edges()) {

        // Normal plane is the average of normal planes at adjacent faces
        Vector3 normal{0.0, 0.0, 0.0};
        if(e.halfedge().face().isReal()) normal += geometry->normal(e.halfedge().face());
        if(e.halfedge().twin().face().isReal()) normal += geometry->normal(e.halfedge().twin().face());
        normal = unit(normal);

        // Center point is the halfway point along the edge
        Vector3 centerPoint = 0.5 * (geometry->position(e.halfedge().vertex()) +
                                     geometry->position(e.halfedge().twin().vertex()));

        // Direction is the edge vector, rotated PI/2
        Vector3 edgeV = unit(geometry->vector(e.halfedge()));
        Vector3 vectorDir = edgeV.rotate_around(normal, PI/2.0);
        double vectorLength = scaleFactor * edgeVectors[e];

        bases[iVec] = centerPoint; 
        scaledVectors[iVec] = vectorDir * vectorLength;
        iVec++;
    }
    
    setVectors(bases, scaledVectors, radius);
}

void SceneVectors::setEdgeVectors(Geometry<Euclidean>* geometry, const EdgeData<double> &form, bool autoscale) {
   
    HalfedgeMesh* mesh = geometry->getMesh();

    // Pick a nice scale for the data
    double scaleFactor = 1.0;
    double radius = -1;
    if(autoscale) {
        double totalVectorLength = 0;
        double maxVectorLength = 0;
        for(EdgePtr e : mesh->edges()) {
            double l = std::abs(form[e] * norm(geometry->vector(e.halfedge())));
            totalVectorLength += l;
            maxVectorLength = std::max(maxVectorLength, l);
        }
        double meanVectorLength = totalVectorLength / mesh->nEdges();
        double targetLen = 0.5 * (meanVectorLength + maxVectorLength);
        scaleFactor = 0.3 * computeNiceMeshLengthScale(geometry) / meanVectorLength;
    }
    radius = 0.1 * computeNiceMeshLengthScale(geometry);


    std::vector<Vector3> bases(mesh->nEdges());
    std::vector<Vector3> scaledVectors(mesh->nEdges());
    size_t iVec = 0;
    for(EdgePtr e : mesh->edges()) {

        // Center point is the halfway point along the edge
        Vector3 centerPoint = 0.5 * (geometry->position(e.halfedge().vertex()) +
                                     geometry->position(e.halfedge().twin().vertex()));

        Vector3 vectorRoot = geometry->vector(e.halfedge());
        double vectorLength = scaleFactor * form[e];

        bases[iVec] = centerPoint; 
        scaledVectors[iVec] = vectorRoot * vectorLength;
        iVec++;
    }
    
    setVectors(bases, scaledVectors, radius);

}

void SceneVectors::setOneForm(Geometry<Euclidean>* geometry, const EdgeData<double> &form, bool autoscale, bool zeroSparse) {
    HalfedgeMesh* mesh = geometry->getMesh();
    using namespace Eigen;

    // == Approximate per-face vectors
    FaceData<Vector3> faceVectors(mesh);

    // Useful geometry stuff
    FaceData<Vector3> faceNormal(mesh);
    geometry->getFaceNormals(faceNormal);
    HalfedgeData<Vector3> halfedgeVector(mesh);
    geometry->getHalfedgeVectors(halfedgeVector);
    HalfedgeData<double>orientationSign(mesh);
    for(HalfedgePtr he : mesh->halfedges()) {
        if(he == he.edge().halfedge()) {
            orientationSign[he] = 1.0;
        } else {
            orientationSign[he] = -1.0;
        }
    }

    for(FacePtr f : mesh->faces()) {

        // Build a local basis
        Vector3 basisX = unit(halfedgeVector[f.halfedge()]);
        Vector3 basisY = basisX.rotate_around(faceNormal[f], PI / 2.0);

        // Compute a local map from 1-forms to face vector
        Vector3 result;
        if(!zeroSparse) {
            Matrix<double, 3, 2> vectorMat;
            unsigned int i = 0;
            for(HalfedgePtr he : f.adjacentHalfedges()) {
                Vector3 v = halfedgeVector[he];
                vectorMat(i, 0) = dot(v, basisX);
                vectorMat(i, 1) = dot(v, basisY);
                i++;
            } 
            auto solver = vectorMat.jacobiSvd(ComputeFullU | ComputeFullV);
    
            Matrix<double, 3, 1> oneFormVector;
    
            i = 0;
            HalfedgePtr firstHe = f.halfedge();
            HalfedgePtr currHe = firstHe;
            do {
                oneFormVector(i) = form[currHe.edge()] * orientationSign[currHe];
                i++;
                currHe = currHe.next();
            } while (currHe != firstHe);
    
            Vector2d faceVec = solver.solve(oneFormVector);
            result = faceVec[0] * basisX + faceVec[1] * basisY;
        } 
        // Handle sparse data with multiple cases 
        else {

            // Count nonzeros
            int nnz = 0;
            bool nonzero[3];
            unsigned int j = 0;
            for(HalfedgePtr he : f.adjacentHalfedges()) {
                if(form[he.edge()] != 0) {
                    nnz++;
                }
                nonzero[j] = form[he.edge()] != 0;
                j++;
            }

            // Branch to an approriate solve based on how many edges have nonzeros 
            switch (nnz) {
                case 0: {
                    result = Vector3{0,0,0};
                } break;
                
                case 1: {

                    Vector3 sampleV;
                    double formVal;
                    for(HalfedgePtr he : f.adjacentHalfedges()) {
                        if(nonzero[j]) {
                            Vector3 v = halfedgeVector[he];
                            sampleV = unit(v);
                            formVal = form[he.edge()] * orientationSign[he];
                        }
                        j++;
                    }

                    result = sampleV * formVal;

                } break;
                
                case 2: {
                    
                    Matrix<double, 2, 2> vectorMat;
                    unsigned int i = 0;
                    j = 0;
                    for(HalfedgePtr he : f.adjacentHalfedges()) {
                        if(nonzero[j]) {
                            Vector3 v = halfedgeVector[he];
                            vectorMat(i, 0) = dot(v, basisX);
                            vectorMat(i, 1) = dot(v, basisY);
                            i++;
                        }
                        j++;
                    } 
                    auto solver = vectorMat.jacobiSvd(ComputeFullU | ComputeFullV);
                    Matrix<double, 2, 1> oneFormVector;
                    i = 0;
                    j = 0;
                    HalfedgePtr firstHe = f.halfedge();
                    HalfedgePtr currHe = firstHe;
                    do {
                        if(nonzero[j]) {
                            oneFormVector(i) = form[currHe.edge()] * orientationSign[currHe];
                            i++;
                        }
                        j++;
                        currHe = currHe.next();
                    } while (currHe != firstHe);
            
                    Vector2d faceVec = solver.solve(oneFormVector);
                    result = faceVec[0] * basisX + faceVec[1] * basisY;

                } break;
                
                case 3: {

                    Matrix<double, 3, 2> vectorMat;
                    unsigned int i = 0;
                    for(HalfedgePtr he : f.adjacentHalfedges()) {
                        Vector3 v = halfedgeVector[he];
                        vectorMat(i, 0) = dot(v, basisX);
                        vectorMat(i, 1) = dot(v, basisY);
                        i++;
                    } 
                    auto solver = vectorMat.jacobiSvd(ComputeFullU | ComputeFullV);
            
                    Matrix<double, 3, 1> oneFormVector;
            
                    i = 0;
                    HalfedgePtr firstHe = f.halfedge();
                    HalfedgePtr currHe = firstHe;
                    do {
                        oneFormVector(i) = form[currHe.edge()] * orientationSign[currHe];
                        i++;
                        currHe = currHe.next();
                    } while (currHe != firstHe);
            
                    Vector2d faceVec = solver.solve(oneFormVector);
                    result = faceVec[0] * basisX + faceVec[1] * basisY;

                } break;
            }


        }


        faceVectors[f] = result;
    }

    setVectors(geometry, faceVectors, autoscale);
}