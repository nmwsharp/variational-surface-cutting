#include "fast_triangle_soup.h"

#include "polygon_soup_mesh.h"
#include "dense_matrix.h"
#include "fast_cholesky.h"

#include <string>
#include <sstream>
#include <set>
#include <map>

using std::vector;

FastTriangleSoup::FastTriangleSoup()
{
    
}


FastTriangleSoup::FastTriangleSoup(const std::vector<Tri>& triangles_, const std::vector<Vector3>& vertices_, const std::vector<TriBool>& isBoundaryEdge_)
    : triangles(triangles_), vertices(vertices_), isBoundaryEdge(isBoundaryEdge_)
{
    nVert = vertices.size();
    nTri = triangles.size();

    // Find boundary vertices
    isBoundaryVertex = std::vector<bool>(vertices.size(), false);
    for(size_t iTri = 0; iTri < triangles.size(); iTri++) {
        const Tri& tri = triangles[iTri];
        const TriBool& isB = isBoundaryEdge[iTri];
        for(size_t i = 0; i < 3; i++) {
            if(isB[i]) {
                isBoundaryVertex[tri[i]] = true;
                isBoundaryVertex[tri[(i+1)%3]] = true;
            }
        }
    }

}

FastTriangleSoup::FastTriangleSoup(Geometry<Euclidean>* geometry) {

    HalfedgeMesh* mesh = geometry->getMesh();

    for(VertexPtr v : mesh->vertices()) {
        vertices.push_back(geometry->position(v));
    }

    VertexData<size_t> vInd = mesh->getVertexIndices();
    for(FacePtr f : mesh->faces()) {

        HalfedgePtr he0 = f.halfedge();
        HalfedgePtr he1 = he0.next();
        HalfedgePtr he2 = he1.next();

        triangles.push_back({{vInd[he0.vertex()], vInd[he1.vertex()], vInd[he2.vertex()]}});
        isBoundaryEdge.push_back({{he0.edge().isBoundary(), he1.edge().isBoundary(), he2.edge().isBoundary()}});
    }


    nVert = vertices.size();
    nTri = triangles.size();

    // Find boundary vertices
    isBoundaryVertex = std::vector<bool>(vertices.size(), false);
    for(size_t iTri = 0; iTri < triangles.size(); iTri++) {
        const Tri& tri = triangles[iTri];
        const TriBool& isB = isBoundaryEdge[iTri];
        for(size_t i = 0; i < 3; i++) {
            if(isB[i]) {
                isBoundaryVertex[tri[i]] = true;
                isBoundaryVertex[tri[(i+1)%3]] = true;
            }
        }
    }
    
}

FastTriangleSoup::~FastTriangleSoup() {
    delete interiorLaplacian;
}


void FastTriangleSoup::cacheGeometry() {

    // Allocate/initialize space
    cotan.resize(nTri);
    // pointwiseGaussCurvature.resize(nVert);
    integratedGaussCurvature.resize(nVert);
    std::fill(integratedGaussCurvature.begin(), integratedGaussCurvature.end(), 0.0);
    dualArea.resize(nVert);
    std::fill(dualArea.begin(), dualArea.end(), 0.0);
    boundaryLength.resize(nVert);
    std::fill(boundaryLength.begin(), boundaryLength.end(), 0.0);
    triangleArea.resize(nTri);
    std::fill(triangleArea.begin(), triangleArea.end(), 0.0);
    faceNormal.resize(nTri);
    std::fill(faceNormal.begin(), faceNormal.end(), Vector3{0.0,0.0,0.0});

    // Compute face-wise quantities
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri tri = triangles[iTri];
        TriBool triB = isBoundaryEdge[iTri];

        // Area
        Vector3 v1d = vertices[tri[1]] - vertices[tri[0]];
        Vector3 v2d = vertices[tri[2]] - vertices[tri[0]];
        double area = 0.5 * norm(cross(v1d,v2d));
        faceNormal[iTri] = unit(cross(v1d,v2d));
        triangleArea[iTri] = area;
        for(int j = 0; j < 3; j++) {
            dualArea[tri[j]] += area / 3.0;
        }

        // Boundary length
        for(int j = 0; j < 3; j++) {
            size_t vA = tri[j];
            size_t vB = tri[(j+1)%3];
            if(triB[j]) {
                double bLen = norm(vertices[vA] - vertices[vB]);
                boundaryLength[vA] += bLen / 2.0;
                boundaryLength[vB] += bLen / 2.0;
            }
        }

        // Angles and cotans
        for(int j = 0; j < 3; j++) {
            Vector3 v1 = vertices[tri[j]] - vertices[tri[(j+2)%3]];
            Vector3 v2 = vertices[tri[(j+1)%3]] - vertices[tri[(j+2)%3]];

            double theta = angle(v1, v2);
            integratedGaussCurvature[tri[(j+2)%3]] += theta; // accumulate here temporarily
            cotan[iTri][j] = dot(v1, v2) / norm(cross(v1, v2));
        }
    }

    // Finish computing curvature
    for(size_t iVert = 0; iVert < nVert; iVert++) {
        if(isBoundaryVertex[iVert]) {
            integratedGaussCurvature[iVert] = 0.0;
        } else {
            integratedGaussCurvature[iVert] = (2.0*PI - integratedGaussCurvature[iVert]);
        }
    }

    geometryCached = true;
}


void FastTriangleSoup::solveLaplaceDirichlet(vector<double>& inputVals) {

    if(!geometryCached) {
        cacheGeometry();
    }


    // Index interior vertices
    vector<size_t> ind(nVert, INVALID_IND);
    size_t nInterior = 0;
    for(size_t i = 0; i < nVert; i++) {
        if(!std::isfinite(inputVals[i])) {
            ind[i] = nInterior++;
        } 
    }

    std::vector<int> vCount(nVert, 0);
    for(const Tri& t : triangles) {
        for(int j = 0; j < 3; j++) {
            vCount[t[j]]++;
        }
    }

    // Estimate the degree of each interior vertex
    std::vector<size_t> degreeEst(nInterior, 2);
    for(const Tri& t : triangles) {
        for(int j = 0; j < 3; j++) {
            size_t vA = t[j];
            size_t vB = t[(j+1)%3];

            if(ind[vA] != INVALID_IND) {
                degreeEst[ind[vA]]++;
            }
            if(ind[vB] != INVALID_IND) {
                degreeEst[ind[vB]]++;
            }
        }
    }

    // Allocate the Cholesky matrix and right hand side
    FastCholesky L(nInterior);
    L.reserveColumns(degreeEst);
    GC::DenseVector<double> rhs(nInterior);

    // Build the Poisson problem
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];


        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];
        
            if(ind[indRoot] == INVALID_IND) continue;
        
            double cotanA = cotan[iTri][iR];
            double cotanB = cotan[iTri][iB];
        
            // Integral against neighbor A 
            if(ind[indA] != INVALID_IND) {
                L.addValue(ind[indRoot], ind[indA], -0.5 * cotanA);
            } else {
                rhs(ind[indRoot]) += 0.5 * cotanA * inputVals[indA];
            }
        
            // Integral against neighbor B
            if(ind[indB] != INVALID_IND) {
                L.addValue(ind[indRoot], ind[indB], -0.5 * cotanB);
            } else {
                rhs(ind[indRoot]) += 0.5 * cotanB * inputVals[indB];
            }
        
            // Integral against central vertex
            L.addValue(ind[indRoot], ind[indRoot], 0.5 * (cotanA + cotanB));
        }
    }

    // Solve the Poisson problem
    L.shiftDiagonal(1e-4);
    L.factor();
    GC::DenseVector<double> x = L.solve(rhs);

    // Store the result 
    for(size_t i = 0; i < nVert; i++) {
        if(ind[i] != INVALID_IND) {
            inputVals[i] = x(ind[i]);
        }
    }

}

vector<double> FastTriangleSoup::solveYamabeProblem() {

    std::vector<double> negCurvature = faceCurvature;
    for(double& x : negCurvature) x *= -1.0;

    return solveInteriorPoissonFace(negCurvature);
}



void FastTriangleSoup::buildInteriorLaplacian() {

    if(!geometryCached) {
        cacheGeometry();
    }

    // Index interior vertices
    interiorInd = vector<size_t>(nVert, INVALID_IND);
    nInterior = 0;
    for(size_t i = 0; i < nVert; i++) {
        if(!isBoundaryVertex[i]) {
            interiorInd[i] = nInterior++;
        }
    }

    std::vector<int> vCount(nVert, 0);
    for(const Tri& t : triangles) {
        for(int j = 0; j < 3; j++) {
            vCount[t[j]]++;
        }
    }

    // Upper bound on degree of each vertex is # neighboring faces + 1.
    // Add one to get number of entries in row.
    std::vector<size_t> degreeEst(nInterior, 2);
    for(const Tri& t : triangles) {
        for(int j = 0; j < 3; j++) {
            size_t vA = t[j];
            size_t vB = t[(j+1)%3];

            if(interiorInd[vA] != INVALID_IND) {
                degreeEst[interiorInd[vA]]++;
            }
            if(interiorInd[vB] != INVALID_IND) {
                degreeEst[interiorInd[vB]]++;
            }
        }
    }

    // Allocate the Cholesky matrix and right hand side
    interiorLaplacian = new FastCholesky(nInterior);
    interiorLaplacian->reserveColumns(degreeEst);

    // Build the Poisson problem
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];


        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];
        
            if(isBoundaryVertex[indRoot]) continue;
        
            double cotanA = cotan[iTri][iR];
            double cotanB = cotan[iTri][iB];
        
            // Integral against neighbor A 
            if(!isBoundaryVertex[indA]) {
                interiorLaplacian->addValue(interiorInd[indRoot], interiorInd[indA], -0.5 * cotanA);
            }
        
            // Integral against neighbor B
            if(!isBoundaryVertex[indB]) {
                interiorLaplacian->addValue(interiorInd[indRoot], interiorInd[indB], -0.5 * cotanB);
            }
        
            // Integral against central vertex
            interiorLaplacian->addValue(interiorInd[indRoot], interiorInd[indRoot], 0.5 * (cotanA + cotanB));
        
        }
    }

    // Solve the Poisson problem
    interiorLaplacian->shiftDiagonal(1e-4);
    interiorLaplacian->factor();

    haveInteriorLaplacian = true;
}


vector<double> FastTriangleSoup::solveInteriorPoisson(const std::vector<double>& pointwiseDensity, bool lumpedRHS) {

    // Ensure we have the Laplacian
    if(!haveInteriorLaplacian) {
        buildInteriorLaplacian();
    }


    // Build the RHS
    GC::DenseVector<double> rhs(nInterior);

    // Build the Poisson problem
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];


        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];
        
            if(isBoundaryVertex[indRoot]) continue;
        
        
            // RHS forcing term
            double fRoot = pointwiseDensity[indRoot];
            double fA = pointwiseDensity[indA];
            double fB = pointwiseDensity[indB];

            if(lumpedRHS) {
              rhs(interiorInd[indRoot]) += triangleArea[iTri] * fRoot / 3.0;
            } else {
              rhs(interiorInd[indRoot]) += triangleArea[iTri] / 12.0 * (2*fRoot + fA + fB);
            }

        }
    }

    // Solve
    GC::DenseVector<double> x = interiorLaplacian->solve(rhs);

    // Store the result in interior vertices
    vector<double> result(nVert, 0.0);
    for(size_t i = 0; i < nVert; i++) {
        if(!isBoundaryVertex[i]) {
            result[i] = x(interiorInd[i]);
        }
    }


    // Compute the boundary normal a boundary vertices
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];

        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];

            if(!isBoundaryVertex[indRoot]) continue;
            
            double cotanA = cotan[iTri][iR];
            double cotanB = cotan[iTri][iB];

            if(!isBoundaryVertex[indA]) {
                result[indRoot] += 0.5 * cotanA * result[indA];
            } 
            if(!isBoundaryVertex[indB]) {
                result[indRoot] += 0.5 * cotanB * result[indB];
            } 

            double fRoot = pointwiseDensity[indRoot];
            double fA = pointwiseDensity[indA];
            double fB = pointwiseDensity[indB];
            if(lumpedRHS) {
              result[indRoot] += triangleArea[iTri] * fRoot / 3.0;
            } else {
              result[indRoot] += triangleArea[iTri] / 12.0 * (2*fRoot + fA + fB);
            }
        }
    }
    

    // We want the pointwise boundary normal
    for(size_t i = 0; i < nVert; i++) {
        if(isBoundaryVertex[i]) {
            result[i] /= boundaryLength[i];
        }
    }

    return result;
}

vector<double> FastTriangleSoup::solveInteriorPoissonFace(const std::vector<double>& pointwiseFaceDensity) {

    // Ensure we have the Laplacian
    if(!haveInteriorLaplacian) {
        buildInteriorLaplacian();
    }


    // Build the RHS
    GC::DenseVector<double> rhs(nInterior);

    // Build the Poisson problem
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];


        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];
        
            if(isBoundaryVertex[indRoot]) continue;
        
        
            // RHS forcing term
            double fFace = pointwiseFaceDensity[iTri];
            rhs(interiorInd[indRoot]) += (triangleArea[iTri] / 3.0) * fFace;
        }
    }

    // Solve
    GC::DenseVector<double> x = interiorLaplacian->solve(rhs);

    // Store the result in interior vertices
    vector<double> result(nVert, 0.0);
    for(size_t i = 0; i < nVert; i++) {
        if(!isBoundaryVertex[i]) {
            result[i] = x(interiorInd[i]);
        }
    }


    // Compute the boundary normal a boundary vertices
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];

        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];

            if(!isBoundaryVertex[indRoot]) continue;
            
            double cotanA = cotan[iTri][iR];
            double cotanB = cotan[iTri][iB];

            if(!isBoundaryVertex[indA]) {
                result[indRoot] += 0.5 * cotanA * result[indA];
            } 
            if(!isBoundaryVertex[indB]) {
                result[indRoot] += 0.5 * cotanB * result[indB];
            } 

            // RHS forcing term
            double fFace = pointwiseFaceDensity[iTri];
            result[indRoot] += (triangleArea[iTri] / 3.0) * fFace;
        }
    }
    

    // We want the pointwise boundary normal
    for(size_t i = 0; i < nVert; i++) {
        if(isBoundaryVertex[i]) {
            result[i] /= boundaryLength[i];
        }
    }

    return result;
}

void FastTriangleSoup::extendFromBoundary(std::vector<double>& vals) {

    // Ensure we have the Laplacian
    if(!haveInteriorLaplacian) {
        buildInteriorLaplacian();
    }

    for(double x : vals) {
        if(!std::isfinite(x)) {
            cout << "Non finite intput to extend" << endl;
        }
    }    

    // Build the RHS
    GC::DenseVector<double> rhs(nInterior);

    // Build the Poisson problem
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri& t = triangles[iTri];

        for(size_t iRoot = 0; iRoot < 3; iRoot++) {

            size_t iR = iRoot;
            size_t iA = (iRoot+1)%3;
            size_t iB = (iRoot+2)%3;
            
            size_t indRoot = t[iR];
            size_t indA = t[iA];
            size_t indB = t[iB];
        
            if(interiorInd[indRoot] == INVALID_IND) continue;
        
            double cotanA = cotan[iTri][iR];
            double cotanB = cotan[iTri][iB];
        
            // Integral against neighbor A 
            if(interiorInd[indA] == INVALID_IND) {
                rhs(interiorInd[indRoot]) += 0.5 * cotanA * vals[indA];
            }
        
            // Integral against neighbor B
            if(interiorInd[indB] == INVALID_IND) {
                rhs(interiorInd[indRoot]) += 0.5 * cotanB * vals[indB];
            }
        }
    }

    // Solve
    GC::DenseVector<double> x = interiorLaplacian->solve(rhs);

    // DEBUG
    for(size_t i = 0; i < x.nRows(); i++) {
        if(!std::isfinite(x(i))) {
            cout << "Bad output from extend Laplacian :(" << endl;
        }
    }

    // Store the result in interior vertices
    for(size_t i = 0; i < nVert; i++) {
        if(!isBoundaryVertex[i]) {
            vals[i] = x(interiorInd[i]);
        }
    }

}


void FastTriangleSoup::improveQuality(int nIterations) {

    const double timestep = 0.1;

    // Hacky method by Nick that seems to work quite well: Move each vertex in the mean direction of its halfedge
    // vectors, towards the average of its vertex star positions. This is nice because it doesn't require any 
    // extra connectivity information which is expensive to get.


    // Identify sharp vertices, which we will not move during the optimization
    // To avoid needed extra connectivity, instead of using dihedral angle, we measure
    // the maximum deviation of any face normal from the vertex average.
    double sharpThreshold = PI / 6.0;
    vector<Vector3> vertexNormals(nVert, Vector3{0.0, 0.0, 0.0});
    for(Tri& tri : triangles) {
        Vector3 vecA = vertices[tri[1]] - vertices[tri[0]];
        Vector3 vecB = vertices[tri[2]] - vertices[tri[0]];
        Vector3 areaNormal = 0.5 * cross(vecA, vecB);
        for(int i = 0; i < 3; i++) {
            vertexNormals[tri[i]] += areaNormal;
        }
    }
    for(size_t i = 0; i < nVert; i++) {
        vertexNormals[i] = unit(vertexNormals[i]);
    }
    vector<double> maxNormalDeviation(nVert, 0.0);
    for(Tri& tri : triangles) {
        Vector3 vecA = vertices[tri[1]] - vertices[tri[0]];
        Vector3 vecB = vertices[tri[2]] - vertices[tri[0]];
        Vector3 normal = unit(cross(vecA, vecB));
        for(int i = 0; i < 3; i++) {
            double theta = angle(vertexNormals[tri[i]], normal);
            maxNormalDeviation[tri[i]] = std::max(maxNormalDeviation[tri[i]], theta);
        }
    }
    vector<char> canMoveVertex(nVert, true);
    for(size_t i = 0; i < nVert; i++) {
        if(isBoundaryVertex[i] || maxNormalDeviation[i] > sharpThreshold) {
            canMoveVertex[i] = false;
        }
    }
    

    
    // Initialize energy gradient to zero
    vector<Vector3> energyGradient(nVert, Vector3{0.0, 0.0, 0.0});

    for(int iIter = 0; iIter < nIterations; iIter++) {

        // Clear out vertex normals and energy gradient 
        std::fill(energyGradient.begin(), energyGradient.end(), Vector3{0.0, 0.0, 0.0});
        std::fill(vertexNormals.begin(), vertexNormals.end(), Vector3{0.0, 0.0, 0.0});
        
        // Compute gradient and update vertex normals
        // Wrong for boundary verts, but we don't use those values
        for(Tri& tri : triangles) {
        
            // Start gradient
            for(int i = 0; i < 3; i++) {
                Vector3 vec = vertices[tri[(i+1)%3]] - vertices[tri[i]];
                energyGradient[tri[i]] += -vec;
            }

            // Boundary term gradient
            // (tries to shrink skinny triangles near boundary)
            /*
            for(int i = 0; i < 3; i++) {

                int iB1 = i;
                int iB2 = (i+1)%3;
                int iTip = (i+2)%3;

                if(isBoundaryVertex[tri[iB1]] && isBoundaryVertex[tri[iB2]]) {
                    Vector3 edgeMidpoint = 0.5 * (vertices[tri[iB1]] + vertices[tri[iB2]]);
                    Vector3 vec = vertices[tri[iTip]] - edgeMidpoint;
                    energyGradient[tri[iTip]] += 2.0 * vec;
                }
            }
            */
    
            // Vertex normals 
            Vector3 vecA = vertices[tri[1]] - vertices[tri[0]];
            Vector3 vecB = vertices[tri[2]] - vertices[tri[0]];
            Vector3 areaNormal = 0.5 * cross(vecA, vecB);
            for(int i = 0; i < 3; i++) {
                vertexNormals[tri[i]] += areaNormal;
            }
        }
    
        // Normalize vertex normals 
        for(size_t i = 0; i < nVert; i++) {
            vertexNormals[i] = unit(vertexNormals[i]);
        }
        
        
        // Apply gradient flow in the tangential direction
        for(size_t i = 0; i < nVert; i++) {
            Vector3 g = timestep * energyGradient[i];
        
        
            if (canMoveVertex[i] && g.isDefined()) {
        
                // TODO would be better to still allow sliding
                // long sharp edge, rather than freezing entirely
    
                // Project out the normal component
                Vector3 N = vertexNormals[i];
                if (N.isFinite() && N.isDefined()) {

                    g -= dot(g, N) * N;
    
                    // Clamp if the move would be too large TODO?
                    /* 
                    double maxDot = 0;
                    for(HalfedgePtr he : v.outgoingHalfedges()) {
                        double d = std::abs(dot(g, unit(geometry->vector(he)))) / geometry->length(he.edge());
                        maxDot = std::max(maxDot, d);
                    }
                    if(maxDot > 0.25) {
                        g *= (0.25 / maxDot);
                    }
                    */
    
    
                    // Integrate using forward Euler
                    vertices[i] -= g;
                }
            }
        } 
    }

    geometryCached = false;
}

// Helper methods for managing a disjoint set array
size_t find(std::vector<size_t>& arr, size_t i) {
    if(arr[i] == i) {
        return i;
    }
    size_t res = find(arr, arr[i]);
    arr[i] = res; // path compression
    return res;
}
// Unions element 'i' to be a part of group 'j'
void unionTo(std::vector<size_t>& arr, size_t i, size_t j) {
    arr[i] = j;
}

void FastTriangleSoup::collapseBoundaryEdges() {

    throw std::runtime_error("Not currently implemented");

    /*
    const double collapseThresh = 5.0;

    // Map of original vertices to collapsed vertices
    // (this is basically an inline union-find datastructure)
    std::vector<size_t> collapseMap(nVert);
    for(size_t i = 0; i < nVert; i++) {
        collapseMap[i] = i;
    }

    // New list of faces
    std::vector<Tri> newTriangles;
    int nCollapse = 0;
    for(Tri t : triangles) {

        bool collapsed = false;

        for(int i = 0; i < 3; i++) {

            size_t indB1 = find(collapseMap, t[i]);
            size_t indB2 = find(collapseMap, t[(i+1)%3]);
            size_t indTip = find(collapseMap, t[(i+2)%3]);

            if(isBoundaryVertex[indB1] && isBoundaryVertex[indB2]) {
            // if(isBoundaryVertex[t[(i+1)%3]] && isBoundaryVertex[t[(i+2)%3]]) {

                // Compute the skinniness ratio for the edge
                Vector3 edge = (vertices[indB1] - vertices[indB2]);
                Vector3 edgeMidpoint = 0.5 * (vertices[indB1] + vertices[indB2]);
                Vector3 vec = vertices[indTip] - edgeMidpoint;
                double ratio = norm(vec) / norm(edge);

                if(ratio < collapseThresh) continue;

                // Collapse the edge
                nCollapse++;
                collapsed = true;

                size_t indBKeep = indB1;
                size_t indBRemove = indB2;

                // Move the remaining vertex to the center
                vertices[indBKeep] = edgeMidpoint;
                pointwiseCurvature[indBKeep] = 0.5 * (pointwiseCurvature[indBKeep] + pointwiseCurvature[indBRemove]);

                //  References to the removed vertex now point to the kept vertex
                unionTo(collapseMap, indBRemove, indBKeep);

                break;
            }
        }


        if(!collapsed) {
            newTriangles.push_back(t);
        }
    }

    // Copy over used vertices
    vector<Vector3> newVertices;
    vector<char> newIsBoundary;
    vector<double> newCurvature;
    vector<size_t> newVertInd(nVert, INVALID_IND);
    size_t nNewVert = 0;
    for(size_t i = 0; i < nVert; i++) {
        if(find(collapseMap, i) == i) {
            newVertices.push_back(vertices[i]);
            newIsBoundary.push_back(isBoundaryVertex[i]); 
            newCurvature.push_back(pointwiseCurvature[i]); 
            newVertInd[i] = nNewVert++;
        }
    }
    
     
    collapseVertInd.resize(nVert);
    for(size_t i = 0; i < nVert; i++) {
        collapseVertInd[i] = newVertInd[find(collapseMap, i)];
    }
    
    // Fix up triangle indices
    for(Tri &t : newTriangles) {
        for(int j = 0; j < 3; j++) {
            t[j] = newVertInd[find(collapseMap, t[j])];
        }
    }

    // Store new values
    triangles = newTriangles;
    vertices = newVertices;
    // isBoundaryVertex = newIsBoundary;
    pointwiseCurvature = newCurvature;
    nTri = triangles.size();
    nVert = vertices.size();

    geometryCached = false;

    cout << "Collapsed " << nCollapse << " skinny boundary faces." << endl;
    */
}
    
std::vector<double> FastTriangleSoup::applyCollapseMap(const std::vector<double>& vals) {

    vector<double> result(collapseVertInd.size(), 33);
    for(size_t i = 0; i < result.size(); i++) {
        result[i] = vals[collapseVertInd[i]];
    }

    return result;
}


std::vector<double> FastTriangleSoup::invertCollapseMap(const std::vector<double>& vals) {

    // Collapsed value is overwritten and thus depends on loop order, would be better to take average
    vector<double> result(vertices.size());
    for(size_t i = 0; i < result.size(); i++) {
        result[collapseVertInd[i]] = vals[i];
    }

    return result;
}

void FastTriangleSoup::toHalfedgeMesh(HalfedgeMesh*& mesh, Geometry<Euclidean>*& geometry) { 

    // Convert to polygon soup
    std::vector<std::vector<size_t>> pInds;
    for(Tri t : triangles) {
        pInds.push_back({t[0], t[1], t[2]});
    }
    PolygonSoupMesh pSoup(pInds, vertices);

    mesh = new HalfedgeMesh(pSoup, geometry);
}


std::vector<Vector3> FastTriangleSoup::computeVertexNormals() {

    std::vector<Vector3> vertNormals(vertices.size());
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        Tri tri = triangles[iTri];

        // Area
        Vector3 v1d = vertices[tri[1]] - vertices[tri[0]];
        Vector3 v2d = vertices[tri[2]] - vertices[tri[0]];
        Vector3 areaNormal = cross(v1d,v2d);

        for(size_t iVert : tri) {
            vertNormals[iVert] += areaNormal;
        }
    }

    // Normalize
    for(Vector3& v : vertNormals) {
        v = unit(v);
    }

    return vertNormals;
}
