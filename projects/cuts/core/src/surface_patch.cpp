#include "surface_patch.h"

#include "sparse_matrix.h"

SurfacePatch::SurfacePatch( const std::vector<Tri>& triangles_, 
                            const std::vector<Vector3>& vertices_,
                            const std::vector<TriBool>& boundaryEdges_,
                            const std::vector<size_t>& parentVertex_,
                            const std::vector<size_t>& parentFace_)

    : soupMesh(triangles_, vertices_, boundaryEdges_), 
      parentVertex(parentVertex_), parentFace(parentFace_)
{

    INVALID_IND = std::numeric_limits<size_t>::max();

    nVert = soupMesh.vertices.size();
    nTri = soupMesh.triangles.size();

    // Index boundary vertices in the soup mesh
    nBoundaryVertices = 0;
    bInd = std::vector<size_t>(nVert, INVALID_IND);
    for(size_t i = 0; i < nVert; i++) {
        if(soupMesh.isBoundaryVertex[i]) {
            boundaryVertexIndices.push_back(i);
            bInd[i] = nBoundaryVertices;
            nBoundaryVertices++;
        }
    }
}



SurfacePatch::~SurfacePatch() {
    if(geometry != nullptr) {
        delete geometry->getMesh();
        delete geometry;
    } 
}


void SurfacePatch::extendFromBoundary(std::vector<double>& vals) {
    soupMesh.extendFromBoundary(vals);
}

void SurfacePatch::generateGeometry() {
    HalfedgeMesh* mesh_;
    soupMesh.toHalfedgeMesh(mesh_, geometry);
}

void SurfacePatch::buildBoundaryValueSmoother() {

    // == First, compute an aspect ratio per boundary vertex, which is the worst of the aspect
    //    ratios at the two adjacent boundary elements
    std::vector<double> boundaryVertexAspect(nBoundaryVertices, 0.0);
    std::vector<double> dualLength(nBoundaryVertices, 0.0);
    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];
        const TriBool& triB = soupMesh.isBoundaryEdge[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            if(!triB[i]) continue; // only process boundary edges

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];
            size_t ind2 = tri[(i+2)%3];

            Vector3 p0 = soupMesh.vertices[ind0];
            Vector3 p1 = soupMesh.vertices[ind1];
            Vector3 p2 = soupMesh.vertices[ind2];

            double width = norm(p0 - p1);
            double height = norm(p2 - 0.5*(p0+p1));
            double aspect = height / width; // aspect ratio against this base
            
            size_t bInd0 = bInd[ind0];
            size_t bInd1 = bInd[ind1];

            // vert 0
            boundaryVertexAspect[bInd0] = std::max(boundaryVertexAspect[bInd0], aspect);
            dualLength[bInd0] += width; // Accumulate dual length while we're at it

            // vert 1
            boundaryVertexAspect[bInd1] = std::max(boundaryVertexAspect[bInd1], aspect);
            dualLength[bInd1] += width;

        }
    }

    // == Determine how much each boundary vertex wants to retain its value, according to a logistic function
    // 'aspectThresh' is essentially the threshold of aspect ratio at which we stop trusting a triangle to give to
    //    good boundary normals. At that point we weight it by 50%, and the weight drops rapidly beyond.
    // 'logistcSlope' should probably stay at 1, unless you're really picky.
    const double aspectThresh = 3.0;
    const double logisticSlope = 1.5;
    auto idFrac = [=](double aspect) {return 1.0 - 1.0 / (1.0 + std::exp(-logisticSlope * (aspect - aspectThresh)));};
    bIdFrac = std::vector<double>(nBoundaryVertices);;
    for(size_t i = 0; i < nBoundaryVertices; i++) {
        bIdFrac[i] = idFrac(boundaryVertexAspect[i]);
    }

    // == Build linear Laplacian smoother
    // convex combination: (a I + (1-a) L)x = a y
    bValSmoother = GC::SparseMatrix<double>(nBoundaryVertices, nBoundaryVertices);

    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];
        const TriBool& triB = soupMesh.isBoundaryEdge[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            if(!triB[i]) continue; // only process boundary edges

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];
            size_t ind2 = tri[(i+2)%3];

            Vector3 p0 = soupMesh.vertices[ind0];
            Vector3 p1 = soupMesh.vertices[ind1];
            Vector3 p2 = soupMesh.vertices[ind2];

            double width = norm(p0 - p1);
            double height = norm(p2 - 0.5*(p0+p1));

            size_t bInd0 = bInd[ind0];
            size_t bInd1 = bInd[ind1];

            double w0 = (1.0 - bIdFrac[bInd0]) * clamp(dualLength[bInd0] / width, 0.0, 1000.0); // inverse distance weight
            bValSmoother(bInd0, bInd0) += w0;
            bValSmoother(bInd0, bInd1) += -w0;
            
            // vert 1 
            double w1 = (1.0 - bIdFrac[bInd1]) * clamp(dualLength[bInd1] / width, 0.0, 1000.0);
            bValSmoother(bInd1, bInd1) += w1;
            bValSmoother(bInd1, bInd0) += -w1;
        }
    }

    // Add remaining identity component
    for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
        bValSmoother(iB, iB) += bIdFrac[iB];
    }

    bValSmoother.shift(1e-4);
}

std::vector<double> SurfacePatch::smoothOverSliverBoundaryElements(const std::vector<double>& vals) {

    // Make sure we have the smoother
    if(!haveBoundaryValueSmoother) {
        buildBoundaryValueSmoother();
    }

    // Build the right hand side
    GC::DenseVector<double> rhs(nBoundaryVertices);
    for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
        rhs(iB) = bIdFrac[iB] * vals[iB];
    }

    // Solve
    GC::DenseVector<double> result;
    // GC::solvePositiveDefinite(bValSmoother, result, rhs);
    #pragma omp critical
    {
        GC::solve(bValSmoother, result, rhs);
    }

    // Copy
    std::vector<double> output(nBoundaryVertices);
    for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
        output[iB] = result(iB);
    }


    // DEBUG things, uncomment to make plots of quantities along boundary
    /*
    findBoundaryChain();
    cout << "[PLOTBOUNDARYSMOOTHING" << iComp << "]#TITLE:Boundary du/dn Smoothing" << endl;
    cout << "[PLOTBOUNDARYSMOOTHING" << iComp << "]#NAMES:Original,Smoothed" << endl;
    cout << "[PLOTBOUNDARYSMOOTHING" << iComp << "]#XLABEL:Boundary" << endl;
    cout << "[PLOTBOUNDARYSMOOTHING" << iComp << "]#YLABEL:Value" << endl;
    cout << "[PLOTBOUNDARYSMOOTHING" << iComp << "]#HASXVAL:True" << endl;
    for(size_t iC = 0; iC < boundaryChain.size(); iC++) {
        size_t iVert = boundaryChain[iC];
        double t = boundaryChainDist[iC];
        double initVal = vals[bInd[iVert]];
        double smoothedVal = output[bInd[iVert]];
        double frac = bIdFrac[bInd[iVert]];

        cout << "[PLOTBOUNDARYSMOOTHING" << iComp << "]"
             << t << "," 
             << initVal << "," 
             << smoothedVal << ","
             << frac 
             << endl;
    }
    */

    return output;
}

void SurfacePatch::computeArea() {

    area = 0;

    if(!soupMesh.geometryCached) {
        soupMesh.cacheGeometry();
    }

    area = std::accumulate(soupMesh.triangleArea.begin(), soupMesh.triangleArea.end(), 0.0);
}

void SurfacePatch::computeNormalDeviation() {

    // FIXME Right now we're using extrinsic averaging across the board for normals, but that's certainly not the best
    // metric for S^2
    
    // Compute mean normal
    meanNormal = Vector3{0.0, 0.0, 0.0};
    for(const Tri& tri : soupMesh.triangles) {
        Vector3 areaNormal = 0.5 * cross(soupMesh.vertices[tri[1]] - soupMesh.vertices[tri[0]], soupMesh.vertices[tri[2]] - soupMesh.vertices[tri[0]]);
        meanNormal += areaNormal;
    }
    meanNormal = unit(meanNormal);

    // Compute deviation
    normalDeviation.resize(nTri); 
    for(size_t iTri = 0; iTri < nTri; iTri++) {
        const Tri& tri = soupMesh.triangles[iTri];
        Vector3 areaNormal = 0.5 * cross(soupMesh.vertices[tri[1]] - soupMesh.vertices[tri[0]], soupMesh.vertices[tri[2]] - soupMesh.vertices[tri[0]]);
        Vector3 normal = unit(areaNormal);

        double err = norm(normal - meanNormal); 
        normalDeviation[iTri] = err;
    }
}

void SurfacePatch::clearBoundaryGradient() {
    boundaryGradient = std::vector<double>(nBoundaryVertices, 0.0); 
}    

// Penalizes negative area if targetArea == -1,
// deviation from target otherwise
void SurfacePatch::addAreaGradient(double weight, double targetArea) {

    double weightParam = weight * std::pow(globalScaleFactor, -4);

    if(targetArea == -1.0) {
        for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
            boundaryGradient[iB] += -weightParam * 2.0 * area;
        }
    } else {
        for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
            boundaryGradient[iB] += weightParam * 2.0 * (area - targetArea);
        }
    }
}

void SurfacePatch::addNormalDeviationGradient(double weight, bool useLocalScaling) {

    std::vector<Vector3> vertexNormals = soupMesh.computeVertexNormals(); 
    double weightParam = weight * std::pow(globalScaleFactor, -2);

    for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
        size_t ind = boundaryVertexIndices[iB];

        Vector3 vNormal = vertexNormals[ind];
        double err = norm(meanNormal - vNormal);

        boundaryGradient[iB] += weightParam * 2.0 * err;
    }
}
    
void SurfacePatch::addDirichletDistortionGradient(double weight, bool useLocalScaling) {

    // Smooth du/dn
    std::vector<double> smoothedYamabeDuDn = smoothOverSliverBoundaryElements(yamabeDuDn);

    if(useLocalScaling) {

        // Solve the adjoint problem
        std::vector<double> pRHS(nVert, 0.0);

        for (size_t iTri = 0; iTri < nTri; iTri++) {
          for (size_t iV = 0; iV < 3; iV++) {

            // = Add Laplace contributions from edge

            // Indices
            size_t iTail = soupMesh.triangles[iTri][iV];
            size_t iTip = soupMesh.triangles[iTri][(iV + 1) % 3];

            // Weights
            double scaleWeight = 0.5 * (std::pow(localScaleFactor[iTail], 3) +
                                        std::pow(localScaleFactor[iTip], 3));
            double cotanWeight = 0.5 * soupMesh.cotan[iTri][iV];
            double edgeWeight = scaleWeight * cotanWeight;

            // Values
            double vTail = distortion[iTail];
            double vTip = distortion[iTip];
            double dualAreaTail = soupMesh.dualArea[iTail];
            double dualAreaTip = soupMesh.dualArea[iTip];

            pRHS[iTip] += -2.0 * edgeWeight * (vTip - vTail) / dualAreaTip;
            pRHS[iTail] += -2.0 * edgeWeight * (vTail - vTip) / dualAreaTail;
          }
        }
       
        // Boundary contribution
        for(size_t iV = 0; iV < nVert; iV++) {
          if(soupMesh.isBoundaryVertex[iV]) {
            double scaleWeight = std::pow(localScaleFactor[iV], 3);
            double dualArea = soupMesh.dualArea[iV];
            pRHS[iV] += -2.0 * scaleWeight * yamabeDuDn[bInd[iV]] * soupMesh.boundaryLength[iV] / dualArea;
          }
        }

        std::vector<double> adjoint = soupMesh.solveInteriorPoisson(pRHS, true);

        //for(size_t iV = 0; iV < nVert; iV++) {
          //if(soupMesh.isBoundaryVertex[iV]) {
            //cout << "boundaryV  -2*dudn = " << -2*yamabeDuDn[bInd[iV]] << "  dpdn = " << adjoint[iV] << endl;
          //} else {
            //cout << "interiorV  -2*u = " << -2*distortion[iV] << "  adjoint = " << adjoint[iV] << endl;
          //}
        //}

        // Extract boundary values
        std::vector<double> dpdn(nBoundaryVertices);
        for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
            size_t ind = boundaryVertexIndices[iB];
            dpdn[iB] = adjoint[ind];
        }

        // Smooth boundary values
        std::vector<double> smoothedDpDn = smoothOverSliverBoundaryElements(dpdn);

        // Compute gradient
        for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
            // theta = -l * du/dn*du/dn - du/dn * dp/dn
            double duduTerm = std::pow(localScaleFactor[boundaryVertexIndices[iB]], 3) * smoothedYamabeDuDn[iB] * smoothedYamabeDuDn[iB];
            double dpduTerm = smoothedDpDn[iB] * smoothedYamabeDuDn[iB];

            //cout << "dudn = " << smoothedYamabeDuDn[iB] << "  dpdn = " << smoothedDpDn[iB] << endl;

            boundaryGradient[iB] += weight * std::pow(globalScaleFactor, -3) * (-duduTerm - dpduTerm);
        }


    } else { 

        for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
            // theta = (du/dn)^2
            boundaryGradient[iB] += weight * smoothedYamabeDuDn[iB] * smoothedYamabeDuDn[iB];
        }

    }
}


void SurfacePatch::addHenckyDistortionGradient(double weight) {

    // Smooth du/dn
    std::vector<double> smoothedYamabeDuDn = smoothOverSliverBoundaryElements(yamabeDuDn);

    // Solve the adjoint problem
    std::vector<double> jPrime(nVert);
    for(size_t iVert = 0; iVert < nVert; iVert++) {
        jPrime[iVert] = weight * std::pow(globalScaleFactor, -2) * 2.0 * distortion[iVert];
    }
    std::vector<double> adjoint = soupMesh.solveInteriorPoisson(jPrime);

    // Extract boundary values
    std::vector<double> dpdn(nBoundaryVertices);
    for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
        size_t ind = boundaryVertexIndices[iB];
        dpdn[iB] = adjoint[ind];
    }

    // Smooth boundary values
    std::vector<double> smoothedDpDn = smoothOverSliverBoundaryElements(dpdn);

    // Compute gradient
    for(size_t iB = 0; iB < nBoundaryVertices; iB++) {
        // theta = l * (j + du/dn * dp/dn) = l * (0 + du/dn * dp/dn)
        boundaryGradient[iB] += smoothedYamabeDuDn[iB] * smoothedDpDn[iB];
    }

}


void SurfacePatch::solveYamabeProblem() {

    // Do the hard work
    distortion = soupMesh.solveYamabeProblem();

    // Move boundary values
    yamabeDuDn.resize(nBoundaryVertices);
    for(size_t bI : boundaryVertexIndices) {
        yamabeDuDn[bInd[bI]] = distortion[bI];
        distortion[bI] = 0.0; // boundary condition
    }
    
}
    
double SurfacePatch::computeBoundaryLengthEnergy() {

    double E = 0;

    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];
        const TriBool& triB = soupMesh.isBoundaryEdge[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            if(!triB[i]) continue; // only process boundary edges

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];

            double L = norm(soupMesh.vertices[ind0] - soupMesh.vertices[ind1]);

            E += L;
        }
    }

    return E / globalScaleFactor;
}

double SurfacePatch::computeLocalScaledBoundaryLengthEnergy() {

    double E = 0;

    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];
        const TriBool& triB = soupMesh.isBoundaryEdge[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            if(!triB[i]) continue; // only process boundary edges

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];

            double meanScale = 0.5 * (localScaleFactor[ind0] + localScaleFactor[ind1]);

            double L = norm(soupMesh.vertices[ind0] - soupMesh.vertices[ind1]) / meanScale;

            E += L;
        }
    }

    return E;
}

double SurfacePatch::computeVisibilityEnergy(const std::vector<double>& cutMeshVisibility) {

    double E = 0;

    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];
        const TriBool& triB = soupMesh.isBoundaryEdge[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            if(!triB[i]) continue; // only process boundary edges

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];

            double meanViz = 0.5 * (cutMeshVisibility[parentVertex[ind0]] + cutMeshVisibility[parentVertex[ind1]]);

            double L = norm(soupMesh.vertices[ind0] - soupMesh.vertices[ind1]) * meanViz;

            E += L * meanViz;
        }
    }

    return E / globalScaleFactor;
}

double SurfacePatch::computeDirichletDistortionEnergy() {

    double E = 0;

    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];

            double cotan = soupMesh.cotan[iTri][i]; 

            double gradDistortion = (distortion[ind0] - distortion[ind1])*(distortion[ind0] - distortion[ind1]) * cotan * 0.5;
            E += gradDistortion;
        }
    }

    return E;    
}

double SurfacePatch::computeLocalScaledDirichletDistortionEnergy() {

    double E = 0;

    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];

            double cotan = soupMesh.cotan[iTri][i];
            double meanScale = 0.5 * (std::pow(localScaleFactor[ind0],3) + std::pow(localScaleFactor[ind1],3));

            double gradDistortion = meanScale * (distortion[ind0] - distortion[ind1])*(distortion[ind0] - distortion[ind1]) * cotan * 0.5;
            E += gradDistortion;
        }
    }


    return E * std::pow(globalScaleFactor, -3);    
}

double SurfacePatch::computeHenckyDistortionEnergy() {

    double E = 0;

    for(size_t iVert = 0; iVert < nVert; iVert++) {

        double A = soupMesh.dualArea[iVert];
        double u = distortion[iVert];

        E += A * u*u;
    }

    return E / (globalScaleFactor * globalScaleFactor);
}

// Penalizes negative size if targetArea == -1, deviation from size otherwise
double SurfacePatch::computeAreaEnergy(double targetArea) {

    if(targetArea == -1) {
        // TODO negative energy makes plots ugly
        return -area * area * std::pow(globalScaleFactor, -4);
    } else {    
        return (area-targetArea) * (area-targetArea) * std::pow(globalScaleFactor, -4);
    }
}

double SurfacePatch::computeNormalDeviationEnergy() {

    double E = 0;

    for(size_t iTri = 0; iTri < nTri; iTri++) {
        const Tri& tri = soupMesh.triangles[iTri];
        Vector3 areaNormal = 0.5 * cross(soupMesh.vertices[tri[1]] - soupMesh.vertices[tri[0]], soupMesh.vertices[tri[2]] - soupMesh.vertices[tri[0]]);
        Vector3 normal = unit(areaNormal);
        double area = norm(areaNormal);

        double err = norm(normal - meanNormal); 
        E += err * err * area;
    }

    return E * std::pow(globalScaleFactor, -2);
}


void SurfacePatch::findBoundaryChain() {

    // Build two-sided lookup
    std::vector<std::vector<double>> neighDist(nVert);
    std::vector<std::vector<size_t>> neighInd(nVert);
    for(size_t iTri = 0; iTri < nTri; iTri++) {

        const Tri& tri = soupMesh.triangles[iTri];
        const TriBool& triB = soupMesh.isBoundaryEdge[iTri];

        // Loop over sides of the triangles
        for(int i = 0; i < 3; i++) {

            if(!triB[i]) continue; // only process boundary edges

            size_t ind0 = tri[i];
            size_t ind1 = tri[(i+1)%3];

            Vector3 p0 = soupMesh.vertices[ind0];
            Vector3 p1 = soupMesh.vertices[ind1];

            double width = norm(p0 - p1);

            neighInd[ind0].push_back(ind1);
            neighDist[ind0].push_back(width);

            neighInd[ind1].push_back(ind0);
            neighDist[ind1].push_back(width);

        }
    }

    
    std::vector<bool> usedVert(nVert, false);
    boundaryChain.clear();

    size_t lastVert = boundaryVertexIndices[0];
    boundaryChain.push_back(lastVert);
    boundaryChainDist.push_back(0);
    usedVert[lastVert] = true;

    // Search until chain complete
    while(true) {

        bool found = false;
        for(size_t iNeigh = 0; iNeigh < 2; iNeigh++) {

            if(usedVert[neighInd[lastVert][iNeigh]]) continue;

            found = true;
            boundaryChainDist.push_back(boundaryChainDist.back() + neighDist[lastVert][iNeigh]);
            lastVert = neighInd[lastVert][iNeigh];
            boundaryChain.push_back(lastVert);
            usedVert[lastVert] = true;
            break;
        }

        if(!found) break;
    }


    cout << "Found boundary chain of length " << boundaryChain.size() << endl;
}
