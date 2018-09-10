#include "eulerian_shape_optimizer.h"

#include <csignal>
#include <queue>
#include <tuple>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "discrete_operators.h"
#include "solvers.h"
#include "timing.h"
#include "fast_marching_method.h"
#include "local_feature_size.h"
#include "polygon_soup_mesh.h"
#include "region_management.h"
#include "combining_hash_functions.h"
#include "compute_visibility.h"
#include "sdf_initialization.h"
#include "mesh_union_split.h"

#include "bff.h"
#include "boundary_constraints.h"

#define PLOT_ENERGY (true)
#define SAFETY_CHECKS (true)

using namespace GC;

EulerianShapeOptimizer::EulerianShapeOptimizer(Geometry<Euclidean> *geometry_) : 
    mesh(geometry_->getMesh()), geometry(geometry_)
{
    cacheGeometry();
    initializeData();
    buildOperators();
    initializeTimestep();
    setState(randomMSDF(geometry));
}

EulerianShapeOptimizer::~EulerianShapeOptimizer() {
    clearCaches(); // frees memory
}

void EulerianShapeOptimizer::cacheGeometry(void) {

    // Indices
    INVALID_IND = std::numeric_limits<size_t>::max();
    vInd = mesh->getVertexIndices();
    fInd = mesh->getFaceIndices();
    eInd = mesh->getEdgeIndices();
    nVert = mesh->nVertices();
    nEdge = mesh->nEdges();
    nFace = mesh->nFaces();

    geometry->getHalfedgeCotans(cotanWeights);
    geometry->getEdgeLengths(edgeLengths);
    geometry->getDualFaceAreas(vertexArea);
    geometry->getFaceAreas(faceArea);
    geometry->getFaceNormals(faceNormal);
    geometry->getHalfedgeVectors(halfedgeVector);
    geometry->getHalfedgeAngles(oppAngles);
    geometry->getEdgeCotanWeights(cotanEdgeWeights);

    // Mean edge length
    meanEdgeLength = 0;
    for(EdgePtr e : mesh->edges()) meanEdgeLength += edgeLengths[e];
    meanEdgeLength /= nEdge;
    EPS_POSITION = meanEdgeLength * 1e-6;

    surfaceArea = geometry->totalArea();
    lengthScale = std::sqrt(surfaceArea);

    // Build orientation sign helper
    orientationSign = HalfedgeData<double>(mesh);
    for(HalfedgePtr he : mesh->halfedges()) {
        if(he == he.edge().halfedge()) {
            orientationSign[he] = 1.0;
        } else {
            orientationSign[he] = -1.0;
        }
    }

    // Build edge area
    edgeArea = EdgeData<double>(mesh, 0.0);
    for(EdgePtr e : mesh->edges()) {
        if(e.halfedge().face().isReal()) {
            edgeArea[e] += faceArea[e.halfedge().face()] / 3.0;
        }
        if(e.halfedge().twin().face().isReal()) {
            edgeArea[e] += faceArea[e.halfedge().twin().face()] / 3.0;
        }
    }

    // Build face curvature
    VertexData<double> angleDefect;
    faceCurvature = FaceData<double>(mesh);
    // if(spaceFillingCurveMode) {
    if(false) {

        // Set the curvature to be that of a sphere with the same area
        double sphereArea = std::sqrt(surfaceArea / (4 * PI));
        double sphereCurvature = 1 / (sphereArea * sphereArea);
        faceCurvature.fill(sphereCurvature);

    } else {

        geometry->getVertexAngleDefects(angleDefect);
        for(FacePtr f : mesh->faces()) {
            double cSum = 0;
            for(HalfedgePtr he : f.adjacentHalfedges()) {
                double oppAngle = geometry->angle(he);
                double oppScale = (2*PI) / (2*PI - angleDefect[he.next().next().vertex()]);
                if(he.next().next().vertex().isBoundary()) {
                    oppScale = 1.0;
                }
                double newAngle = oppAngle * oppScale;
                cSum += newAngle;
            }
    
            double faceC = (cSum - PI) / faceArea[f];
            faceCurvature[f] = faceC;
        }

    }

    cout << "Length scale = " << lengthScale << endl;
}
    

void EulerianShapeOptimizer::initializeData() {

    nConnectedComponents = 0;

    // Plot printing
    if(PLOT_ENERGY) {
        cout << "[ENERGYPLOT]#TITLE:Energy during optimization" << endl;
        cout << "[ENERGYPLOT]#NAMES:"
             << "Length" << ","
             << "Distortion (Dirichlet)" << ","
             << "Distortion (Hencky)" << ","
             << "Area" << ","
             << "Visibility" << ","
             << "Normal Deviation"
             << endl;
        cout << "[ENERGYPLOT]#XLABEL:Iteration" << endl;
        cout << "[ENERGYPLOT]#YLABEL:Energy Contribution" << endl;
    }
}

void EulerianShapeOptimizer::buildOperators() {

    cout << "Building operators..." << endl;

    hodge0 = buildHodge0<double>(geometry);
    hodge1 = buildHodge1<double>(geometry);
    hodge2 = buildHodge2<double>(geometry);
    d0 = buildDerivative0<double>(mesh);
    d0T = d0.transpose();
    d1 = buildDerivative1<double>(mesh);
    d1T = d1.transpose();

    // hodge0Inv = clampedDiagonalInverse(hodge0, EPS_POSITION * EPS_POSITION);
    hodge0Inv = clampedDiagonalInverse(hodge0, -1);
    // hodge1Inv = clampedDiagonalInverse(hodge1, 1e-3);
    hodge1Inv = clampedDiagonalInverse(hodge1, -1);

    eyeV = SparseMatrix<double>::identity(nVert);
    zeroFormWeakLaplacian = d0T * hodge1 * d0;
    zeroFormLaplacian = hodge0Inv * zeroFormWeakLaplacian;


    cout << "  ..done building operators." << endl;
}

void EulerianShapeOptimizer::setFakeFaceCuvature(VertexData<double> scalarFunction) {

    // Set the curvature to be relative to that of a sphere with the same area
    double sphereArea = std::sqrt(surfaceArea / (4 * PI));
    double sphereCurvature = 1 / (sphereArea * sphereArea);

    for(FacePtr f : mesh->faces()) {

        faceCurvature[f] = 1/3.0 * sphereCurvature * (
                                scalarFunction[f.halfedge().vertex()] +
                                scalarFunction[f.halfedge().next().vertex()] + 
                                scalarFunction[f.halfedge().next().next().vertex()] ); 

    }

    cout << "sphere area = " << sphereArea << endl;
    cout << "sphere curvature = " << sphereArea << endl;
}

// Do various custom initializations
/*
void EulerianShapeOptimizer::initializeRegions() {
        
    VertexData<LabelVec> customPhi(mesh, LabelVec::Zero());

    // Initialize the SDF with random distances
    for(VertexPtr v : mesh->vertices()) {
        double x = geometry->position(v).x;
        customPhi[v][0] = x;
        customPhi[v][1] = -x;
    }

    // Passthrough to the main version
    initializeRegions(customPhi);
}
*/

void EulerianShapeOptimizer::setState(VertexData<LabelVec> initialPhi) {

    // Use initial phi vaules
    phi = initialPhi;
    if(hasSymmetry) projectSymmetric();
    regionLabels = VertexData<int>(mesh);
    connectedComponent = VertexData<int>(mesh);
    inferVertexRegions();

    boundaryEdgeCrossing = EdgeData<char>(mesh, false);
    boundaryFaceCrossing = FaceData<char>(mesh, false);
    crossingInnerSegment = HalfedgeData<BoundarySegment*>(mesh, nullptr);
    crossingOuterSegment = HalfedgeData<BoundarySegment*>(mesh, nullptr);
    faceTripleSegments = FaceData<std::array<BoundarySegment*, 3>>(mesh);
    boundarySegments.resize(K_REGIONS);

    clearCaches();
}

void EulerianShapeOptimizer::computeFeatureSize() {
    featureSize = computeLocalFeatureSize_smooth(geometry, zeroFormLaplacian, 0.02);
    // featureSize = computeLocalFeatureSize_eikonal(geometry);
  
    // featureSize = VertexData<double>(mesh, 1.0);
}

void EulerianShapeOptimizer::initializeTimestep() {
    
    // Compute feature size
    computeFeatureSize();

    // localTimestep = VertexData<double>(mesh);
    // for(VertexPtr v : mesh->vertices()) {
    //     localTimestep[v] = std::pow(featureSize[v], 2);
    // } 

}


void EulerianShapeOptimizer::doStep() {

    cerr << endl << endl << "=== Stepping iteration " << iIter << endl;

    double initEnergy = computeEnergy(true);

    buildGradientAtBoundary();
    VertexData<LabelVec> surfaceGradient = extendGradientToSurface();
    //double slope = computeGradientSlope();
    takeGradientStep(surfaceGradient, stepSizeParam);

    double energy = computeEnergy(true);
    double deltaT = stepSizeParam * lengthScale * lengthScale * lengthScale / 1000;
    //double observedSlope = (initEnergy - energy) / deltaT;

    //cout << "Observed slope = "  << observedSlope << " predicted slope = " << slope << endl;



    // for space filling
    if(spaceFillingCurveMode) {
        // weightDirichletDistortion *= 1.00391968; 
        // weightHenckyDistortion *= 1.00391968; 
        weightDirichletDistortion *= 1.005;
        weightHenckyDistortion *= 1.005;
    }


    iIter++;
}


int EulerianShapeOptimizer::doStepLineSearch() {

    cerr << endl << endl << "=== Stepping iteration " << iIter << " (with line search)" << endl;

    double initEnergy = computeEnergy(); // DEBUG test
    VertexData<LabelVec> initState = phi;

    buildGradientAtBoundary();
    VertexData<LabelVec> surfaceGradient = extendGradientToSurface();
    double slope = computeGradientSlope();
    double deltaT = stepSizeParam * lengthScale * lengthScale * lengthScale / 1000;


    // Backtracking line search
    double alphaStep = 1.0;
    double tauShortenFactor = 0.5;
    double cExpectedDecrease = 0.2;
    int nMaxReductions = 10;
    int iSearchIter = 0;
    while(true) {

        // Take a candidate step
        takeGradientStep(surfaceGradient, alphaStep * stepSizeParam);
        double newEnergy = computeEnergy(true);

        double energyDecrease = initEnergy - newEnergy; 
        double expectedDecrease = alphaStep * deltaT * slope;

        cerr << endl << "alpha = " << alphaStep << endl;
        cerr << "decrease = " << energyDecrease << endl;
        cerr << "expected decrease = " << expectedDecrease << endl;

        // Accept the step
        if(energyDecrease > cExpectedDecrease * expectedDecrease ||
           iSearchIter >= nMaxReductions) {
            
            break;
        } 
        // Try a smaller step
        else {

            // Shorten the step size
            alphaStep *= tauShortenFactor;

            // Roll back the iteration
            setState(initState);
        }

        iSearchIter++;
    }

    cerr << "Line search reduced step size " << iSearchIter << " times to " << (alphaStep * stepSizeParam) << endl;

    iIter++;

    return iSearchIter;
}

void EulerianShapeOptimizer::projectSymmetric() {

    if(!hasSymmetry) return;

    for(VertexPtr v : symmetry.canonicalVertices) {
        for(VertexPtr otherV : symmetry.symmetrySet[v]) {
            phi[otherV] = phi[v];
        }
    }
}


void EulerianShapeOptimizer::ensureHaveBoundaryGeometry() {

    if(haveBoundaryGeometry) return;

    inferVertexRegions();
    inferBoundaryGeometry();
    snapSDF();
    recomputeSDFFromBoundary();
    inferVertexRegions();
    inferBoundaryGeometry();

    haveBoundaryGeometry = true;
}
    
void EulerianShapeOptimizer::ensureHaveTriangleSoupCutMeshes() {

    if(haveTriangleSoupCutMeshes) return;

    // Dependencies:
    ensureHaveBoundaryGeometry();
    ensureHaveConnectedComponents();

    extractCutMeshes();
    buildPatchMeshes();
    copyFeatureSizeToPatchMeshes();

    haveTriangleSoupCutMeshes = true;
    
    initializeBoundaryGradient(); // needs to be allocated after cut mesh is extracted,
}

void EulerianShapeOptimizer::ensureHaveYamabeSolution() {

    if(haveYamabeSolution) return;
    
    // Dependencies:
    ensureHaveTriangleSoupCutMeshes();

    solveYamabeProblems(); 

    haveYamabeSolution = false;
}
    
    
void EulerianShapeOptimizer::ensureHaveVisibilityFraction() {

    if(haveVisibilityFraction) return;
    
    // Dependencies:
    ensureHaveTriangleSoupCutMeshes();

    // Raycast to compute visibility. Only needs to be done once since shape does not
    // change.
    if(!visibilityCasted) {
        if(pinPatchMode) {
            visibilityFraction = generateLocationTerm();
            visibilityCasted = true;
        } else {
            // visibilityFraction = computeVisibilityFraction(geometry);
            visibilityFraction = computeVisibilityFraction(geometry, 20, true);
            visibilityCasted = true;
        }
    }
    
    cutMeshVisibility = interpolateToCutMesh(visibilityFraction);

    haveVisibilityFraction = true;
}
    
void EulerianShapeOptimizer::ensureHaveConnectedComponents() {

    if(haveConnectedComponents) return;

    // Dependencies
    ensureHaveBoundaryGeometry();

    findConnectedComponents();

    haveConnectedComponents = true;
}

void EulerianShapeOptimizer::ensureHavePatchArea() {

    if(havePatchArea) return;
    
    // Dependencies:
    ensureHaveTriangleSoupCutMeshes();
    ensureHaveConnectedComponents();

    computePatchArea();

    havePatchArea = true;
}

void EulerianShapeOptimizer::ensureHavePatchNormals() {

    if(havePatchNormals) return;
    
    // Dependencies:
    ensureHaveTriangleSoupCutMeshes();
    ensureHaveConnectedComponents();

    computePatchNormals();

    havePatchNormals = true;
}

void EulerianShapeOptimizer::ensureHaveCutHalfedgeMesh() {

    if(haveCutHalfedgeMesh) return;

    // Dependencies:
    ensureHaveTriangleSoupCutMeshes();

    constructCutHalfedgeMesh();

    haveCutHalfedgeMesh = true;
}


void EulerianShapeOptimizer::ensureHavePatchGeometries() {

    if(havePatchGeometries) return;

    // Dependencies:
    ensureHaveTriangleSoupCutMeshes();

    constructPatchGeometries();

    havePatchGeometries = true;
}


void EulerianShapeOptimizer::ensureHaveCoarseParameterizedMesh() {

    if(haveCoarseParameterizedMesh) return;

    // Dependencies:
    ensureHaveCutHalfedgeMesh();
    ensureHavePatchGeometries();

    computeCoarseFlattenedGeometry();

    haveCoarseParameterizedMesh = true;
}

void EulerianShapeOptimizer::ensureHaveGeodesicParameterizedMesh() {

    if(haveGeodesicParameterizedMesh) return;

    // Dependencies:
    ensureHaveCutHalfedgeMesh();
    ensureHavePatchGeometries();

    computeGeodesicFlattenedGeometry();

    haveGeodesicParameterizedMesh = true;
}

void EulerianShapeOptimizer::ensureHaveGeodesicParameterizedUVPosMesh() {

    if(haveGeodesicParameterizedUVPosMesh) return;

    // Dependencies:
    ensureHaveCutHalfedgeMesh();
    ensureHavePatchGeometries();

    computeGeodesicFlattenedGeometryUVPos();

    haveGeodesicParameterizedUVPosMesh = true;
}

void EulerianShapeOptimizer::ensureHaveExtrinsicDevelopableMesh() {

    if(haveExtrinsicDevelopableMesh) return;

    // Dependencies:
    ensureHaveGeodesicParameterizedMesh();

    computeExtrinsicDevelopableMesh();

    haveExtrinsicDevelopableMesh = true;

}

void EulerianShapeOptimizer::inferVertexRegions() {

    for(VertexPtr v : mesh->vertices()) {
        int myRegion;
        phi[v].minCoeff(&myRegion); // the most-negative distance is the label
                                    // ideally exactly one label should be negative, but
                                    // the representation is not exact
        regionLabels[v] = myRegion; 
    }

}

void EulerianShapeOptimizer::inferBoundaryGeometry() {
    
    START_TIMING(bGeo)


    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
        boundarySegments[iRegion].clear();
    }
    crossingInnerSegment.fill(nullptr);
    crossingOuterSegment.fill(nullptr);
    faceTripleSegments.fill({{nullptr, nullptr, nullptr}});

    // Exhaustively look for boundary segments in every triangle 
    for(FacePtr f : mesh->faces()) {
    
        // Keep track of triplet points so that we can resolve their center points later
        vector<BoundarySegment*> tripleSegments;

        for(HalfedgePtr he : f.adjacentHalfedges()) {
    
            for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {

                VertexPtr v0 = he.vertex();
                VertexPtr v1 = he.next().vertex();
                VertexPtr v2 = he.next().next().vertex();
                
                // Check for a triplet boundary chunk
                if(regionLabels[v1] == iRegion && regionLabels[v0] != iRegion && regionLabels[v2] != iRegion && regionLabels[v0] != regionLabels[v2]) {

                    BoundarySegment b;
                    b.type = BType::TRIPLE;

                    {   // Start crossing
                        b.heStart = he;
                        int otherLabel = regionLabels[v0];
                        b.tStart = 1.0 - computeCrossing(phi[v1][iRegion], phi[v1][otherLabel], phi[v0][otherLabel], phi[v0][iRegion]);
                    }

                    {   // End crossing
                        b.heEnd = he.next();
                        int otherLabel = regionLabels[v2];
                        b.tEnd = computeCrossing(phi[v1][iRegion], phi[v1][otherLabel], phi[v2][otherLabel], phi[v2][iRegion]);
                    }


                    boundarySegments[iRegion].push_back(b);
                    tripleSegments.push_back(&boundarySegments[iRegion].back());
                }
                // Check for a triangular boundary chunk
                else if(regionLabels[v1] == iRegion && regionLabels[v0] != iRegion && regionLabels[v2] != iRegion) {

                    BoundarySegment b;
                    b.type = BType::TRI;

                    {   // Start crossing
                        b.heStart = he;
                        int otherLabel = regionLabels[v0];
                        b.tStart = 1.0 - computeCrossing(phi[v1][iRegion], phi[v1][otherLabel], phi[v0][otherLabel], phi[v0][iRegion]);
                    }

                    {   // End crossing
                        b.heEnd = he.next();
                        int otherLabel = regionLabels[v2];
                        b.tEnd = computeCrossing(phi[v1][iRegion], phi[v1][otherLabel], phi[v2][otherLabel], phi[v2][iRegion]);
                    }

                    boundarySegments[iRegion].push_back(b);
                }
                // Check for a trapezoidal boundary chunk
                else if(regionLabels[v1] == iRegion && regionLabels[v0] != iRegion && regionLabels[v2] == iRegion) {

                    BoundarySegment b;
                    b.type = BType::QUAD;

                    {   // Start crossing
                        b.heStart = he;
                        int otherLabel = regionLabels[v0];
                        b.tStart = 1.0 - computeCrossing(phi[v1][iRegion], phi[v1][otherLabel], phi[v0][otherLabel], phi[v0][iRegion]);
                    }

                    {   // End crossing
                        b.heEnd = he.next().next();
                        int otherLabel = regionLabels[v0];
                        b.tEnd = computeCrossing(phi[v2][iRegion], phi[v2][otherLabel], phi[v0][otherLabel], phi[v0][iRegion]);
                    }

                    boundarySegments[iRegion].push_back(b);
                    
                }
    
            }
        }

        // If this face has triplet morphology, resolve the triple point as the barycenter
        if(tripleSegments.size() == 3) { 

            // Compute the center of the crossing points

            // Use the point where all three distances are equal under linear interpolation
            VertexPtr vArr[] = {f.halfedge().vertex(), f.halfedge().next().vertex(), f.halfedge().next().next().vertex()};
            int labelArr[] = {regionLabels[vArr[0]], regionLabels[vArr[1]], regionLabels[vArr[2]]};
            Vector3 phiVals[] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
            for(int iLab = 0; iLab < 3; iLab++) {
                for(int iVert = 0; iVert < 3; iVert++) {
                    phiVals[iLab][iVert] = phi[vArr[iVert]][labelArr[iLab]];
                } 
            }

            /* HACK: compute as center of crossings instead, to ensure convexity
            Vector3 centerW = triplePointBarycentrics(phiVals[0], phiVals[1], phiVals[2]);
            Vector3 center = centerW[0] * geometry->position(vArr[0]) +
                             centerW[1] * geometry->position(vArr[1]) +
                             centerW[2] * geometry->position(vArr[2]);
            */

            Vector3 center = Vector3{0.0, 0.0, 0.0};
            for(BoundarySegment* b : tripleSegments) {
                center += (startPoint(*b)) / 3.0;
            }

            for(BoundarySegment* b : tripleSegments) {
                b->triplePoint = center;
                // b->tripleBary = centerW;
            }

            faceTripleSegments[f] = {{tripleSegments[0], tripleSegments[1], tripleSegments[2]}};

        } else if(!tripleSegments.empty()) { // only 0 and 3 make sense
            throw std::runtime_error("Something has gone terribly wrong");
        }

    }


    // Mark all boundary crossing edges and faces
    boundaryFaceCrossing.fill(false);
    boundaryEdgeCrossing.fill(false);
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
        for(BoundarySegment &b : boundarySegments[iRegion]) {

            // Mark crossings
            boundaryEdgeCrossing[b.heStart.edge()] = true;
            boundaryEdgeCrossing[b.heEnd.edge()] = true;
            boundaryFaceCrossing[b.heStart.face()] = true;

            BoundarySegment* bPtr = &b;
            crossingOuterSegment[b.heStart] = bPtr;
            crossingInnerSegment[b.heEnd] = bPtr;

        }
    }
    
    
    auto bGeoTime = FINISH_TIMING(bGeo);
    cout << "Inferring boundary geometry took " << pretty_time(bGeoTime) << endl;
}

// Use the fast marching method to compute signed distances from each region boundary for each vertex
// Assumes the current boundary segements are valid
void EulerianShapeOptimizer::recomputeSDFFromBoundary() {


    // Search outward to compute distance for each label, one at a time
    #pragma omp parallel for
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
    
        VertexData<char> onBoundary(mesh, (char)false);
        size_t nOnBoundary = 0;

        // === Establish boundary conditions
        for(BoundarySegment &b : boundarySegments[iRegion]) {
            for(VertexPtr v : b.heStart.face().adjacentVertices()) {
                onBoundary[v] = true;
                nOnBoundary++;
            }
        }

        // If non are on boundary, region doesn't exist
        if(nOnBoundary == 0) {
            cout << "Skipping region " << iRegion << ", not present." << endl;

            for(VertexPtr v : mesh->vertices()) {
                // phi[v][iRegion] = std::numeric_limits<double>::infinity();
                phi[v][iRegion] = 1000 * lengthScale;
            }

            continue;
        }

        // Assemble a sparse list of boundary conditions for FMM algorithm 
        vector<std::pair<VertexPtr, double>> boundaryDistances;
        for(VertexPtr v : mesh->vertices()) {
            if(onBoundary[v]) {
                boundaryDistances.push_back(std::make_pair(v, std::abs(phi[v][iRegion])));
            }
        }
            
        // Run FMM
        START_TIMING(fmm)
        VertexData<double> distResult = FMMDistance(mesh, boundaryDistances, edgeLengths, oppAngles);
        auto fmmTime = FINISH_TIMING(fmm);
        cout << "FMM took " << pretty_time(fmmTime) << endl;

        // Copy over the result (with added sign)
        for(VertexPtr v : mesh->vertices()) {
            if(!onBoundary[v]) {
                if(regionLabels[v] == iRegion) {
                    phi[v][iRegion] = -distResult[v];
                } else {
                    phi[v][iRegion] = distResult[v];
                }
            } 
        }


    }

}

void EulerianShapeOptimizer::snapSDF() {

    // Snap the (multi) signed distance function to a valid configuration.
    // This is identical to the procedure in "Multiple Interacting Liquids"
    // by Losasso et. al.

    for(VertexPtr v : mesh->vertices()) {

        LabelVec& p = phi[v];

        // Find the lowest and second lowest elements
        int smallestInd, secondSmallestInd;
        p.minCoeff(&smallestInd); 
        double smallestVal = p[smallestInd];
        p[smallestInd] = std::numeric_limits<double>::infinity();
        p.minCoeff(&secondSmallestInd); 
        double secondSmallestVal = p[secondSmallestInd];

        // Project
        double mean = 0.5 * (smallestVal + secondSmallestVal);

        // Shift only relevant values
        // smallestVal -= mean;
        // secondSmallestVal -= mean;
        // p[smallestInd] = smallestVal;
        // p[secondSmallestInd] = secondSmallestVal;

        // Shift all values
        p[smallestInd] = smallestVal;
        p -= LabelVec::Ones()*mean;

    }


}

void EulerianShapeOptimizer::findConnectedComponents() {
    
    DisjointSets connectedComponentSets(nVert);

    for(VertexPtr v : mesh->vertices()) {
        for(VertexPtr vn : v.adjacentVertices()) {
            if(regionLabels[v] == regionLabels[vn]) {
                connectedComponentSets.merge(vInd[v], vInd[vn]);
            }
        }
    }

    // Densely re-number connected components
    vector<bool> compIndUsed(nVert, false);
    componentRegion.clear();
    for(VertexPtr v : mesh->vertices()) {
        compIndUsed[connectedComponentSets.find(vInd[v])] = true;
    }
    vector<size_t> compIndNew(nVert);
    size_t nComps = 0;
    for(size_t i = 0; i < nVert; i++) {
        if(compIndUsed[i]) {
            compIndNew[i] = nComps++;
            componentRegion.push_back(regionLabels[mesh->vertex(i)]);
        }
    }

    for(VertexPtr v : mesh->vertices()) {
        connectedComponent[v] = compIndNew[connectedComponentSets.find(vInd[v])];
    }
    nConnectedComponents = nComps;
}

void EulerianShapeOptimizer::computePatchArea() {

    double sumArea = 0.0;
    for(SurfacePatch* p : patches) {
        p->computeArea();
        sumArea += p->area;
    }

    // Safety check
    if(!approxEqualsAbsolute(sumArea, surfaceArea)) {
        throw std::runtime_error("Sum of component areas does not match surface area");
    }
}

void EulerianShapeOptimizer::computePatchNormals() {
    for(SurfacePatch* p : patches) {
        p->computeNormalDeviation();
    }
}

void EulerianShapeOptimizer::extractCutMeshes() {

    using std::vector;

    // === Find and index all vertices in the output mesh
    // These will be
    //   - Vertices in the input mesh
    //   - Edges crossed by a cut
    //   - Faces which contain a triple point

    START_TIMING(extract)

    vector<Vector3> customVertPos;
    cutMeshVertexParent.clear();
    cutMeshRegionBoundary.clear();
    size_t nCustomVert = 0;

    // Used to lookup vertices when building triangles later in the method
    VertexData<size_t> origVertexToCustomVertexInd(mesh, INVALID_IND);
    EdgeData<size_t> origEdgeToCustomVertexInd(mesh, INVALID_IND);
    FaceData<size_t> origFaceToCustomVertexInd(mesh, INVALID_IND);


    // Index vertices
    for(VertexPtr v : mesh->vertices()) {
        customVertPos.push_back(geometry->position(v));
        cutMeshVertexParent.push_back(v);
        cutMeshRegionBoundary.push_back(false);
        origVertexToCustomVertexInd[v] = nCustomVert;
        nCustomVert++;
    } 
    
    // Index crossing edges and faces from boundary segments
    // Note: Since we track segments on both sides of the boundary,
    //       each edge is crossed by several segments. Track the ones we have
    //       already indexed to avoid overcounting.
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
        for(BoundarySegment& b : boundarySegments[iRegion]) {

            // Edge crossings
            if(b.type == BType::TRI || b.type == BType::QUAD) {

                { // Start of segment
                    EdgePtr startEdge = b.heStart.edge();
                    if(origEdgeToCustomVertexInd[startEdge] == INVALID_IND) {
                        customVertPos.push_back(startPoint(b));
                        cutMeshVertexParent.push_back(VertexPtr());
                        cutMeshRegionBoundary.push_back(true);
                        origEdgeToCustomVertexInd[startEdge] = nCustomVert;
                        nCustomVert++;
                    }
                }
                
                { // End of segment
                    EdgePtr endEdge = b.heEnd.edge();
                    if(origEdgeToCustomVertexInd[endEdge] == INVALID_IND) {
                        customVertPos.push_back(endPoint(b));
                        cutMeshVertexParent.push_back(VertexPtr());
                        cutMeshRegionBoundary.push_back(true);
                        origEdgeToCustomVertexInd[endEdge] = nCustomVert;
                        nCustomVert++;
                    }
                }
            }
            // Triple points in faces
            else if(b.type == BType::TRIPLE) {

                { // Start of segment
                    EdgePtr startEdge = b.heStart.edge();
                    if(origEdgeToCustomVertexInd[startEdge] == INVALID_IND) {
                        customVertPos.push_back(startPoint(b));
                        cutMeshVertexParent.push_back(VertexPtr());
                        cutMeshRegionBoundary.push_back(true);
                        origEdgeToCustomVertexInd[startEdge] = nCustomVert;
                        nCustomVert++;
                    }
                }
                
                { // End of segment
                    EdgePtr endEdge = b.heEnd.edge();
                    if(origEdgeToCustomVertexInd[endEdge] == INVALID_IND) {
                        customVertPos.push_back(endPoint(b));
                        cutMeshVertexParent.push_back(VertexPtr());
                        cutMeshRegionBoundary.push_back(true);
                        origEdgeToCustomVertexInd[endEdge] = nCustomVert;
                        nCustomVert++;
                    }
                }

                { // Center point
                    FacePtr face = b.heStart.face();
                    if(origFaceToCustomVertexInd[face] == INVALID_IND) {
                        customVertPos.push_back(b.triplePoint);
                        cutMeshVertexParent.push_back(VertexPtr());
                        cutMeshRegionBoundary.push_back(true);
                        origFaceToCustomVertexInd[face] = nCustomVert;
                        nCustomVert++;
                    }
                }
            }
        }
    } 


    // === Build lists of faces, for the main mesh as well as each of the sub meshes

    // Vectors to fill (all indices point to the single master vertex list)
    vector<Tri> cutMeshTris;
    vector<TriBool> cutMeshIsBoundaryEdge;
    vector<double> customFaceCurvature; 
    cutMeshIsComponentBoundaryEdge.clear();
    cutMeshFacePatchIndex.clear();

    // Faces on the interior
    for(FacePtr f : mesh->faces()) {
        if(!boundaryFaceCrossing[f]) {

            size_t ind0 = origVertexToCustomVertexInd[f.halfedge().vertex()];
            size_t ind1 = origVertexToCustomVertexInd[f.halfedge().next().vertex()];
            size_t ind2 = origVertexToCustomVertexInd[f.halfedge().next().next().vertex()];
            Tri inds = {{ind0, ind1, ind2}};

            TriBool triBoundary = {{
                f.halfedge().edge().isBoundary(),
                f.halfedge().next().edge().isBoundary(),
                f.halfedge().next().next().edge().isBoundary()
            }};

            // Add to the main mesh
            cutMeshTris.push_back(inds);
            customFaceCurvature.push_back(faceCurvature[f]);
            cutMeshFacePatchIndex.push_back(connectedComponent[f.halfedge().vertex()]);
            cutMeshIsComponentBoundaryEdge.push_back(triBoundary);
            cutMeshIsBoundaryEdge.push_back(triBoundary);
        }
    }

    // Faces on the boundary
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
        for(BoundarySegment& b : boundarySegments[iRegion]) {

            switch (b.type) {

                // Each triangular boundary segment contributes 1 face
                case BType::TRI: {
                    
                    // Initial indices
                    size_t indRoot = origVertexToCustomVertexInd[b.heEnd.vertex()];
                    size_t indLeft = origEdgeToCustomVertexInd[b.heStart.edge()];
                    size_t indRight = origEdgeToCustomVertexInd[b.heEnd.edge()];

                    // Add a new vertex in the center
                    /*
                    Vector3 cVert = (startPoint(b) + endPoint(b) + geometry->position(b.heEnd.vertex())) / 3.0;
                    size_t indCenter = nCustomVert;
                    customVertPos.push_back(cVert);
                    customVertIsBoundary.push_back(false);
                    customVertIsRegionBoundary.push_back(false);
                    nCustomVert++;

                    // Triangles
                    std::array<Tri,3> tris = {{ {{indRoot, indRight, indCenter}}, 
                                                {{indRoot, indCenter, indLeft}}, 
                                                {{indCenter, indRight, indLeft}} }};
                    // Add to the meshes
                    for(Tri t : tris) {
                        cutMeshTris.push_back(t);
                        customRegionMeshTris[iRegion].push_back(t);
                    }
                    */
                   
                    std::array<Tri,1> tris = {{ {{indRoot, indRight, indLeft}} }};
                    std::array<TriBool,1> triBoundaries = {{ {{b.heEnd.edge().isBoundary(), false, b.heStart.edge().isBoundary()}} }};
                    std::array<TriBool,1> triCompBoundaries = {{ {{b.heEnd.edge().isBoundary(), true, b.heStart.edge().isBoundary()}} }};
                    
                    // Add to the meshes
                    for(size_t i = 0; i < tris.size(); i++) {
                        cutMeshTris.push_back(tris[i]);
                        customFaceCurvature.push_back(faceCurvature[b.heStart.face()]);
                        cutMeshFacePatchIndex.push_back(connectedComponent[b.heEnd.vertex()]);
                        cutMeshIsBoundaryEdge.push_back(triBoundaries[i]);
                        cutMeshIsComponentBoundaryEdge.push_back(triCompBoundaries[i]);
                    }

                } break;

                // Each quad boundary segment contributes 4 faces
                case BType::QUAD: {
                    
                    // Initial indices (counter-clockwise A-D)
                    size_t indA = origEdgeToCustomVertexInd[b.heStart.edge()];
                    size_t indB = origVertexToCustomVertexInd[b.heStart.twin().vertex()];
                    size_t indC = origVertexToCustomVertexInd[b.heEnd.vertex()];
                    size_t indD = origEdgeToCustomVertexInd[b.heEnd.edge()];
                    
                    // Add a new vertex in the center
                    Vector3 cVert = (startPoint(b) + endPoint(b) + geometry->position(b.heStart.next().vertex()) + geometry->position(b.heEnd.vertex())) / 4.0;
                    size_t indCenter = nCustomVert;
                    customVertPos.push_back(cVert);
                    cutMeshVertexParent.push_back(VertexPtr());
                    cutMeshRegionBoundary.push_back(false);
                    nCustomVert++;

                    std::array<Tri,4> tris = {{ {{indCenter, indA, indB}}, 
                                                {{indCenter, indB, indC}}, 
                                                {{indCenter, indC, indD}}, 
                                                {{indCenter, indD, indA}} }};
                    
                    std::array<TriBool,4> triCompBoundaries = {{{{false, b.heStart.edge().isBoundary(), false}}, 
                                                                {{false, b.heStart.next().edge().isBoundary(), false}}, 
                                                                {{false, b.heEnd.edge().isBoundary(), false}}, 
                                                                {{false, true, false}} }};
                                                            
                    std::array<TriBool,4> triBoundaries = {{{{false, b.heStart.edge().isBoundary(), false}}, 
                                                            {{false, b.heStart.next().edge().isBoundary(), false}}, 
                                                            {{false, b.heEnd.edge().isBoundary(), false}}, 
                                                            {{false, false, false}} }};
                    // Add to the meshes
                    for(size_t i = 0; i < tris.size(); i++) {
                        cutMeshTris.push_back(tris[i]);
                        customFaceCurvature.push_back(faceCurvature[b.heStart.face()]);
                        cutMeshFacePatchIndex.push_back(connectedComponent[b.heEnd.vertex()]);
                        cutMeshIsBoundaryEdge.push_back(triBoundaries[i]);
                        cutMeshIsComponentBoundaryEdge.push_back(triCompBoundaries[i]);
                    }

                } break;

                // Each triple face boundary segment contributes 4 faces
                case BType::TRIPLE: {

                    // Indices (counter-clockwise A-D)
                    size_t indA = origEdgeToCustomVertexInd[b.heStart.edge()];
                    size_t indB = origVertexToCustomVertexInd[b.heEnd.vertex()];
                    size_t indC = origEdgeToCustomVertexInd[b.heEnd.edge()];
                    size_t indD = origFaceToCustomVertexInd[b.heStart.face()];
                    
                    // Add a new vertex in the center
                    Vector3 cVert = (startPoint(b) + endPoint(b) + geometry->position(b.heStart.next().vertex()) + geometry->position(b.heEnd.vertex())) / 4.0;
                    size_t indCenter = nCustomVert;
                    customVertPos.push_back(cVert);
                    cutMeshVertexParent.push_back(VertexPtr());
                    cutMeshRegionBoundary.push_back(false);
                    nCustomVert++;

                    std::array<Tri,4> tris = {{ {{indCenter, indA, indB}}, 
                                                {{indCenter, indB, indC}}, 
                                                {{indCenter, indC, indD}}, 
                                                {{indCenter, indD, indA}} }};
                    
                    std::array<TriBool,4> triCompBoundaries = {{{{false, b.heStart.edge().isBoundary(), false}}, 
                                                                {{false, b.heEnd.edge().isBoundary(), false}}, 
                                                                {{false, true, false}}, 
                                                                {{false, true, false}} }};
                    
                    std::array<TriBool,4> triBoundaries = {{{{false, b.heStart.edge().isBoundary(), false}}, 
                                                            {{false, b.heEnd.edge().isBoundary(), false}}, 
                                                            {{false, false, false}}, 
                                                            {{false, false, false}} }};
                    // Add to the meshes
                    for(size_t i = 0; i < tris.size(); i++) {
                        cutMeshTris.push_back(tris[i]);
                        customFaceCurvature.push_back(faceCurvature[b.heStart.face()]);
                        cutMeshFacePatchIndex.push_back(connectedComponent[b.heEnd.vertex()]);
                        cutMeshIsBoundaryEdge.push_back(triBoundaries[i]);
                        cutMeshIsComponentBoundaryEdge.push_back(triCompBoundaries[i]);
                    }
                    
                } break;

            }
        }
    }

    // === Constuct the meshes and build associations


    // The whole underlying mesh
    cutMesh = new FastTriangleSoup(cutMeshTris, customVertPos, cutMeshIsBoundaryEdge);
    cutMesh->faceCurvature = customFaceCurvature;
    
    // cutMeshGood = std::vector<double>(cutMesh->nVert, 0.0);
    // for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
    //     if(cutMeshRegionBoundary[i]) {
    //         cutMeshGood[i] = 1.0;
    //     }
    // }


    // Associations for the whole underlying mesh
    origMeshToCutMesh = origVertexToCustomVertexInd;
    
    // Data on mesh
    cutMeshInterpolatedPhi = interpolateToCutMesh(phi);

    cout << "Cut mesh extraction took " << pretty_time(FINISH_TIMING(extract)) << endl;
}

void EulerianShapeOptimizer::buildPatchMeshes() {

    const std::vector<Tri>& triangles = cutMesh->triangles;
    const std::vector<Vector3>& vertices = cutMesh->vertices;
    const std::vector<TriBool>& isBoundaryEdge = cutMeshIsComponentBoundaryEdge;

    // Sort the triangles in to buckets for eacn component
    std::vector<std::vector<Tri>> componentTriangles(nConnectedComponents);
    std::vector<std::vector<TriBool>> componentIsBoundaryEdge(nConnectedComponents);
    std::vector<std::vector<size_t>> cutMeshTriInd(nConnectedComponents); // triangle index in cut mesh
    std::vector<std::vector<double>> componentFaceCurvatures(nConnectedComponents);
    for(size_t i = 0; i < triangles.size(); i++) {
        componentTriangles[cutMeshFacePatchIndex[i]].push_back(triangles[i]);
        componentIsBoundaryEdge[cutMeshFacePatchIndex[i]].push_back(isBoundaryEdge[i]);
        cutMeshTriInd[cutMeshFacePatchIndex[i]].push_back(i);
        componentFaceCurvatures[cutMeshFacePatchIndex[i]].push_back(cutMesh->faceCurvature[i]);
    }

    // Build a patch mesh for each component
    std::vector<size_t> vertPatchIndex(vertices.size(), INVALID_IND);
    for(size_t iComp = 0; iComp < nConnectedComponents; iComp++) {

        size_t nPatchVert = 0;

        // Vectors of data for this patch
        std::vector<Tri> patchTriangles; // has indices translated to this mesh
        std::vector<Vector3> patchVertices;
        std::vector<size_t> cutMeshVertexInd; // vertex index in cut mesh

        // Loop over all triangles that are a part of this patch, assembling necessary values
        for(const Tri& tri : componentTriangles[iComp]) {

            // If any over the vertices are being seen for the first time, index the vertex
            // and copy data to vectors for this patch
            for(size_t iVert : tri) {
                if(vertPatchIndex[iVert] != INVALID_IND) continue; // seen before

                // Index it
                vertPatchIndex[iVert] = nPatchVert;
                patchVertices.push_back(vertices[iVert]);

                // Copy relevant data
                cutMeshVertexInd.push_back(iVert);

                nPatchVert++;
            }

            // Copy the triangle with new indices
            patchTriangles.push_back({{
                vertPatchIndex[tri[0]],
                vertPatchIndex[tri[1]],
                vertPatchIndex[tri[2]]
            }});

        }

        // Build the actual patch
        SurfacePatch* newPatch = new SurfacePatch(patchTriangles, patchVertices, componentIsBoundaryEdge[iComp], cutMeshVertexInd, cutMeshTriInd[iComp]);
        newPatch->iComp = iComp;
        newPatch->iRegion = componentRegion[iComp]; 
        newPatch->soupMesh.faceCurvature = componentFaceCurvatures[iComp];
        newPatch->globalScaleFactor = lengthScale;
        patches.push_back(newPatch);


        // Unmark seen vertices
        for(size_t ind : cutMeshVertexInd) {
            vertPatchIndex[ind] = INVALID_IND;
        }
    }
}

void EulerianShapeOptimizer::copyFeatureSizeToPatchMeshes() {

    GC::DenseVector<double> cutMeshFeatureSize = getCutMeshScaleVector(true);

    // Copy the global and local scales to each patch
    for(SurfacePatch* p : patches) {


        p->localScaleFactor.resize(p->nVert);
        for(size_t iVert = 0; iVert < p->nVert; iVert++) {
            p->localScaleFactor[iVert] = cutMeshFeatureSize(p->parentVertex[iVert]);
        }
    }
}

void EulerianShapeOptimizer::solveYamabeProblems() {

    #pragma omp parallel for
    for(size_t iP = 0; iP < patches.size(); iP++) {
        SurfacePatch* p = patches[iP];
        p->solveYamabeProblem();

#ifdef SAFETY_CHECKS
        for(size_t i = 0; i < p->soupMesh.vertices.size(); i++) {
            double u = p->distortion[i];
            if(!std::isfinite(u)) {
                invalidValuePanic("solveYamabeProblems()");
            }
        }
#endif

    }

}


void EulerianShapeOptimizer::initializeBoundaryGradient() {

    // Clear boundary gradient in patch meshes
    for(SurfacePatch* p : patches) {
        p->clearBoundaryGradient();
    }

}

void EulerianShapeOptimizer::buildGradientAtBoundary() {

    initializeBoundaryGradient();

    addBoundaryGradientTermLengthRegularization();
    addBoundaryGradientTermDirichletDistortion();
    addBoundaryGradientTermHenckyDistortion();
    addBoundaryGradientTermArea();
    addBoundaryGradientTermVisibility();
    addBoundaryGradientTermNormalDeviation();


#ifdef SAFETY_CHECKS
    for(SurfacePatch* p : patches) {
        for(size_t iB = 0; iB < p->nBoundaryVertices; iB++) {
            size_t ind = p->boundaryVertexIndices[iB];
            double bGrad = p->boundaryGradient[iB];
            if(!std::isfinite(bGrad)) {
                invalidValuePanic("buildGradientAtBoundary()");
            }
        }
    }
#endif 

}

void EulerianShapeOptimizer::addBoundaryGradientTermLengthRegularization() {
    
    if(weightLengthRegularization == 0.0) return;

    ensureHaveBoundaryGeometry();

    // Stepped implicitly in takeGradientStep()

    // == Grad l term with local scaling
    if(localScaleLengthRegularization) {
    
        ensureHaveTriangleSoupCutMeshes();

        GC::DenseVector<double> cutMeshScaleInv = getCutMeshScaleVector(true, -1);
    
        // This derivative has two terms, dg / dn, which is computed here,
        // and h * g, which is stepped implicitly in takeGradientStep
    
        // Note: We compute this quantity here, rather than in the patch because
        // here we can evaluate gradients on both sides of the patch
    
        // Make sure we have triangle area
        if(!cutMesh->geometryCached) {
            cutMesh->cacheGeometry();
        }
    
        // The dg/dn term can be computed easily at edges, remap to vertices
        // using area weighting
        double weightParam = weightLengthRegularization;
        std::vector<LabelVec> dldnSum(cutMesh->vertices.size(), LabelVec::Zero());
        std::vector<double> dldnWeight(cutMesh->vertices.size(), 0.0);
        for(size_t iTri = 0; iTri < cutMesh->triangles.size(); iTri++) {
            const Tri& tri = cutMesh->triangles[iTri];
    
            size_t edgeArea = cutMesh->triangleArea[iTri] / 3.0;
    
            for(int triInd = 0; triInd < 3; triInd++) {
    
                size_t vertI = tri[triInd];
                size_t vertJ = tri[(triInd+1)%3];
    
                double edgeLength = norm(cutMesh->vertices[vertI] - cutMesh->vertices[vertJ]);
    
                double lengthScaleGrad = (cutMeshScaleInv(vertJ) - cutMeshScaleInv(vertI)) / edgeLength;
                LabelVec phiGrad = (cutMeshInterpolatedPhi[vertJ] - cutMeshInterpolatedPhi[vertI]) / edgeLength;
                double weight = edgeLength;
    
                // Add gradient for regions adjacent to boundary vertex
                dldnSum[vertI] -= lengthScaleGrad * phiGrad * weight;
                dldnWeight[vertI] += weight;
    
                dldnSum[vertJ] -= lengthScaleGrad * phiGrad * weight;
                dldnWeight[vertJ] += weight;
    
            }
        }
    
        // The actual boundary motion gradient comes from dividing out the area average
        std::vector<LabelVec> dldn(cutMesh->vertices.size());
        for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
            dldn[i] += dldnSum[i] / dldnWeight[i];
        }
    
    
        // Copy the value up to each patch
        for(SurfacePatch* p : patches) {
    
            for(size_t iB = 0; iB < p->nBoundaryVertices; iB++) {
                size_t ind = p->boundaryVertexIndices[iB];
    
                // Get boundary gradient and its index on the cut mesh
                size_t cutInd = p->parentVertex[ind];
                int patchRegion = p->iRegion;

#ifdef SAFETY_CHECKS
                if(!std::isfinite(dldn[cutInd][patchRegion])) {
                    invalidValuePanic("addBoundaryGradientTermLengthRegularization()");
                }
#endif 
    
                // Add the appropariate value to the gradient
                p->boundaryGradient[iB] += weightParam * dldn[cutInd][patchRegion];
            }
        }

        
    }

}

void EulerianShapeOptimizer::addBoundaryGradientTermDirichletDistortion() {

    if(weightDirichletDistortion == 0.0) return;

    ensureHaveYamabeSolution();

    if(pinPatchMode) {

        // #pragma omp parallel for
        for(size_t iP = 0; iP < patches.size(); iP++) {
            SurfacePatch* p = patches[iP];
            if(p->iRegion == 0) {
                p->addDirichletDistortionGradient(weightDirichletDistortion, localScaleDirichletDistortion);
            }
        }

    } else {
        #pragma omp parallel for
        for(size_t iP = 0; iP < patches.size(); iP++) {
            SurfacePatch* p = patches[iP];
            p->addDirichletDistortionGradient(weightDirichletDistortion, localScaleDirichletDistortion);
        }
    }
}


void EulerianShapeOptimizer::addBoundaryGradientTermHenckyDistortion() {

    if(weightHenckyDistortion == 0.0) return;

    ensureHaveYamabeSolution();

    // #pragma omp parallel for
    for(size_t iP = 0; iP < patches.size(); iP++) {
        SurfacePatch* p = patches[iP];
        p->addHenckyDistortionGradient(1000 * weightHenckyDistortion);
    }
}



void EulerianShapeOptimizer::addBoundaryGradientTermArea() {

    if(weightArea == 0.0) return;

    ensureHavePatchArea();


    if(patchZeroTargetArea == -1) {
        for(SurfacePatch* p : patches) {
            p->addAreaGradient(weightArea);
        }
    } else {
        for(SurfacePatch* p : patches) {
            if(p->iRegion == 0) {
                p->addAreaGradient(weightArea, patchZeroTargetArea);
            }
        }
    }

}



void EulerianShapeOptimizer::addBoundaryGradientTermVisibility() {
    
    if(weightVisibility == 0.0) return;
    
    ensureHaveVisibilityFraction();
    ensureHaveTriangleSoupCutMeshes();

    // Remap visibility fraction and phi to cut mesh


    // This derivative has two terms, dg / dn, which is computed here,
    // and h * g, which is stepped implicitly in takeGradientStep

    // Note: We compute this quantity here, rather than in the patch because
    // here we can evaluate gradients on both sides of the patch

    // Make sure we have triangle area
    if(!cutMesh->geometryCached) {
        cutMesh->cacheGeometry();
    }
   

    // The dg/dn term can be computed easily at edges, remap to vertices
    // using area weighting
    double weightParam = weightVisibility / lengthScale;
    std::vector<LabelVec> dgdnSum(cutMesh->vertices.size(), LabelVec::Zero());
    std::vector<double> dgdnWeight(cutMesh->vertices.size(), 0.0);
    for(size_t iTri = 0; iTri < cutMesh->triangles.size(); iTri++) {
        const Tri& tri = cutMesh->triangles[iTri];

        double edgeArea = cutMesh->triangleArea[iTri] / 3.0;

        for(int triInd = 0; triInd < 3; triInd++) {

            size_t vertI = tri[triInd];
            size_t vertJ = tri[(triInd+1)%3];

            double edgeLength = norm(cutMesh->vertices[vertI] - cutMesh->vertices[vertJ]);

            double visibilityGrad = (cutMeshVisibility[vertJ] - cutMeshVisibility[vertI]) / edgeLength;
            LabelVec phiGrad = (cutMeshInterpolatedPhi[vertJ] - cutMeshInterpolatedPhi[vertI]) / edgeLength;
            double weight = edgeArea;

            // Add gradient for regions adjacent to boundary vertex
            dgdnSum[vertI] -= visibilityGrad * phiGrad * weight;
            dgdnWeight[vertI] += weight;

            dgdnSum[vertJ] -= visibilityGrad * phiGrad * weight;
            dgdnWeight[vertJ] += weight;

        }
    }

    // The actual boundary motion gradient comes from dividing out the area average
    std::vector<LabelVec> dgdn(cutMesh->vertices.size());
    for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
        dgdn[i] = dgdnSum[i] / dgdnWeight[i];
    }


    // Copy the value up to each patch
    for(SurfacePatch* p : patches) {

        for(size_t iB = 0; iB < p->nBoundaryVertices; iB++) {
            size_t ind = p->boundaryVertexIndices[iB];

            // Get boundary gradient and its index on the cut mesh
            size_t cutInd = p->parentVertex[ind];
            int patchRegion = p->iRegion;

#ifdef SAFETY_CHECKS
            if(!std::isfinite(dgdn[cutInd][patchRegion])) {
                invalidValuePanic("addBoundaryGradientTermVisibility()");
            }
#endif 

            // Add the appropariate value to the gradient
            p->boundaryGradient[iB] += -weightParam * dgdn[cutInd][patchRegion];
        }
    }

}

void EulerianShapeOptimizer::addBoundaryGradientTermNormalDeviation() {
    
    if(weightNormalDeviation == 0.0) return;
    
    ensureHavePatchNormals();
    ensureHaveTriangleSoupCutMeshes();

    for(SurfacePatch* p : patches) {
        p->addNormalDeviationGradient(weightNormalDeviation);
    }
}

VertexData<LabelVec> EulerianShapeOptimizer::extendGradientToSurface() {

    VertexData<LabelVec> gradient(mesh, LabelVec::Zero());
    
    if(haveTriangleSoupCutMeshes) {
        extendedGradient = std::vector<LabelVec>(cutMesh->vertices.size());
    }

    // Early-out if the boundary gradient is all zero (meaning there is nothing to extend)
    bool boundaryGradientExists = false; 
    for(SurfacePatch* p : patches) {
        for(double x : p->boundaryGradient) {
            if(x != 0.0) {
                boundaryGradientExists = true;
            }
            if(!std::isfinite(x)) {
                cout << "Non finite gradient :(" << endl;
            }
        }
    }
    if(!boundaryGradientExists) {
        cout << "No boundary gradient to extend." << endl;
        return gradient;        
    }


    cout << "Extending boundary motion to SDF derivative..." << endl;

    // === Convert the boundary velocity to an SDF update by applying advection equation

    // Initialize with zero on boundary, NaN elsewhere
    std::vector<LabelVec> motionVals(cutMesh->vertices.size(), LabelVec::Zero()*std::numeric_limits<double>::quiet_NaN());
    for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
        if(cutMeshRegionBoundary[i]) {
            motionVals[i] = LabelVec::Zero(); 
        }
    }


    // Accumulate values from each patch
    for(SurfacePatch* p : patches) {
        for(size_t iB = 0; iB < p->nBoundaryVertices; iB++) {
            size_t ind = p->boundaryVertexIndices[iB];

            // Get boundary gradient and its index on the cut mesh
            size_t cutInd = p->parentVertex[ind];
            double bGrad = p->boundaryGradient[iB];
            int patchRegion = p->iRegion;

            // Convert pointwise boundary motion to change in SDF (derived from advection equation)
            for(int jRegion = 0; jRegion < K_REGIONS; jRegion++) {
                if(jRegion == patchRegion) {
                    motionVals[cutInd][jRegion] += -bGrad;
                } else {
                    motionVals[cutInd][jRegion] += bGrad;
                }
            }
            
        }
    }

    // Handle pin constraint
    // if(pinPatchMode) {
    //     for(SurfacePatch* p : patches) {
    //         if(p->iRegion == 0) {
    //             for(size_t iV = 0; iV < p->soupMesh.vertices.size(); iV++) {

    //                 // Find the pinned vertex
    //                 if(p->parentVertex[iV] == pinVert) {

    //                     double currDist = phi[mesh->vertex(pinVert)][0];
    //                     double targetDist = targetBoundaryDist;

    //                     double weightParam = weightArea * std::pow(lengthScale, -3);

    //                     // Add a gradient term that encourages the vertex to be the target distance from the boundary
    //                     // motionVals[iV][0] = -2*(currDist - targetDist);
    //                     // for(int iRegion = 1; iRegion < K_REGIONS; iRegion++) {

    //                     // }

    //                 }
    //             }
    //         }
    //     }
    // }

    // Accumulate the extended result

    // Extend within each patch
    for(SurfacePatch* p : patches) {

        // Extend each coordinate independently
        for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {

            // NOTE A hack to handle boundaries here; better to do this
            // by building a separate Laplacian with natural BCs at shape
            // boundary

            // Build a vector of values for this region
            std::vector<double> extendedVals(p->soupMesh.vertices.size(), std::numeric_limits<double>::quiet_NaN());
            // std::vector<double> extendedVals(p->soupMesh.vertices.size(), 0);
            for(size_t bInd : p->boundaryVertexIndices) {
                size_t cutInd = p->parentVertex[bInd];
                extendedVals[bInd] = motionVals[cutInd][iRegion];


                // NOTE hacky boundary handling
                if(cutMeshVertexParent[cutInd] != VertexPtr() &&
                   cutMeshVertexParent[cutInd].isBoundary()) {
                        extendedVals[bInd] = std::numeric_limits<double>::quiet_NaN();
                   }
            }
    
            // Harmonic extension
            p->soupMesh.solveLaplaceDirichlet(extendedVals);
            // p->extendFromBoundary(extendedVals);
    
            // Store values for this region
            for(size_t i = 0; i < p->soupMesh.vertices.size(); i++) {
                extendedGradient[p->parentVertex[i]][iRegion] = extendedVals[i];

#ifdef SAFETY_CHECKS
                if(!std::isfinite(extendedVals[i])) {
                    invalidValuePanic("extendGradientToSurface()");
                }
#endif
            }
        }

    }

    // Extend across whole surface
    /*
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {

        std::vector<double> vals(cutMesh->nVert, std::numeric_limits<double>::quiet_NaN());

        // Copy known values at region boundaries in to vals for extension
        for(SurfacePatch* p : patches) {

            // Only process patches from this region
            if(p->iRegion != iRegion) {
                continue; 
            }     

            // Copy boundary value for all boundary verts
            for(size_t vInd : p->boundaryVertexIndices) {
                size_t cutMeshInd = p->parentVertex[vInd];
                vals[cutMeshInd] = motionVals[cutMeshInd][iRegion];
            }
        }

        // Solve the harmonic extension Poisson problem
        // (requires factor :( )
        cutMesh->solveLaplaceDirichlet(vals);

        // Copy result back
        for(size_t iV = 0; iV < cutMesh->nVert; iV++) {
            extendedGradient[iV][iRegion] = vals[iV];
        }
    }
    */


    // Copy result to a vertexdata object on the original mesh
    for(VertexPtr v : mesh->vertices()) {
        gradient[v] = extendedGradient[origMeshToCutMesh[v]];

        // NaN safety... (arbitrarily bad meshes can do terrible things to you)
        for(int i = 0; i < K_REGIONS; i++) {
            if(!std::isfinite(gradient[v][i])) {
                gradient[v][i] = 0.0;
            }
        }
    } 

    return gradient;
}


// This is necessary because of Eulerian represntation and implicit terms -- we can't accurately compute it directly from the
// gradient we use to timestep.
double EulerianShapeOptimizer::computeGradientSlope() {

    // Accumulate the gradient here (use the SDF representation for accumulation) 
    std::vector<LabelVec> cutMeshGradient(cutMesh->nVert, LabelVec::Zero()*std::numeric_limits<double>::quiet_NaN());
    for(size_t i = 0; i < cutMesh->nVert; i++) {
        if(cutMeshRegionBoundary[i]) {
            cutMeshGradient[i] = LabelVec::Zero(); 
        }
    }

    // Accumulate the boundary gradient values from each patch
    for(SurfacePatch* p : patches) {
        for(size_t iB = 0; iB < p->nBoundaryVertices; iB++) {
            size_t ind = p->boundaryVertexIndices[iB];

            // Get boundary gradient and its index on the cut mesh
            size_t cutInd = p->parentVertex[ind];
            double bGrad = p->boundaryGradient[iB];
            int patchRegion = p->iRegion;

            // Convert normal velocity to SDF gradient
            cutMeshGradient[cutInd][p->iRegion] += -bGrad;
        }
    }


    // Explicit gradient contribution
    // This term is normally stepped implicitly, so we need to manually add it to our gradient

    // Compute Laplacian on original mesh and interpolate to cut
    updateExplicitOperator();
    VertexData<LabelVec> LPhi(mesh);
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {

        DenseVector<double> regionPhi(nVert);
        for(VertexPtr v : mesh->vertices()) {
            regionPhi(vInd[v]) = phi[v][iRegion];
        }

        GC::DenseVector<double> LPhiVec = explicitMat * regionPhi;

        // Store
        for(VertexPtr v : mesh->vertices()) {
            LPhi[v][iRegion] = LPhiVec(vInd[v]) / 2.0; // the value in explicitMat includes contributions from both side
        }
    }

    std::vector<LabelVec> lPhiAllCut = interpolateToCutMesh(LPhi);

    // Copy only those laplacianPhiCut values for the regions that
    // actually border that cut
    std::vector<LabelVec> lPhiCut(cutMesh->nVert, LabelVec::Zero());
    for(SurfacePatch* p : patches) {

        for(size_t iTri = 0; iTri < p->nTri; iTri++) {
    
            Tri& t = p->soupMesh.triangles[iTri];
            TriBool& tEdge = p->soupMesh.isBoundaryEdge[iTri];
    
            for(size_t iEdge = 0; iEdge < 3; iEdge++) {

                if(!tEdge[iEdge]) continue; // only process boundary edges

                size_t iA = iEdge;
                size_t iB = (iEdge+1)%3;
                size_t indA = p->parentVertex[t[iA]];
                size_t indB = p->parentVertex[t[iB]];

                lPhiCut[indA][p->iRegion] = lPhiAllCut[indA][p->iRegion];
                lPhiCut[indB][p->iRegion] = lPhiAllCut[indB][p->iRegion];
            }
        }
    }

    // Add the implicit contribution to the accumulated gradient
    for(size_t i = 0; i < cutMesh->nVert; i++) {
        if(cutMeshRegionBoundary[i]) {
            cutMeshGradient[i] += lPhiCut[i];
        }
    }


    // Sum values on oppopsite sides of every dual segment
    std::vector<LabelVec> cutMeshGradientSum(cutMesh->nVert, LabelVec::Zero()*std::numeric_limits<double>::quiet_NaN());
    for(size_t i = 0; i < cutMesh->nVert; i++) {
        if(cutMeshRegionBoundary[i]) {

            // Initialize each region with the value on that side
            LabelVec val = cutMeshGradient[i];

            // Count neighboring regions
            int nNeigh = 0;
            for(int i = 0; i < K_REGIONS; i++) {
                if(val[i] != 0) {
                    nNeigh++;
                }
            }

            // Sum across dual boundary
            // (note that the nNeigh count makes us do the right thing at the triple point,
            //  where the dual boundary is half of one region and half of another)
            double sum = 0;
            for(int iReg = 0; iReg < K_REGIONS; iReg++) {
                if(val[iReg] != 0) {
                    for(int jReg = 0; jReg < K_REGIONS; jReg++) {
                        if(iReg == jReg) continue;
                        val[iReg] -= cutMeshGradient[i][jReg] / (nNeigh - 1); 
                    }
                }
            }

            cutMeshGradientSum[i] = val;
        }
    }

            
    // Motion comes from sum on both sides of patch
    double slope = 0.0;
    for(SurfacePatch* p : patches) {

        if(!p->soupMesh.geometryCached) {
            p->soupMesh.cacheGeometry();
        }

        for(size_t iB = 0; iB < p->nBoundaryVertices; iB++) {
            size_t ind = p->boundaryVertexIndices[iB];

            // Get boundary gradient and its index on the cut mesh
            size_t cutInd = p->parentVertex[ind];

            if(!cutMeshRegionBoundary[cutInd]) continue; // needed for meshes with boundary

            double grad = cutMeshGradientSum[cutInd][p->iRegion];
            double dualLen = p->soupMesh.boundaryLength[ind];

            slope += grad*grad*dualLen / 2.0;  // divide by two because we are computing on both sides
                                               // of every boundary segment
        }
    }


    return slope;
}

void EulerianShapeOptimizer::updateExplicitOperator() {

    // NOTE: When there is no visibility term, we can pull out the diagonal component and solve an SPD
    // problem in implicit step. However, I don't see any way to apply this transformation when there is a
    // visibility term, so it seems it would create cases in the code.

    if(weightLengthRegularization != cachedLengthCoef || weightVisibility != cachedVisibilityCoef) {
        cachedLengthCoef = weightLengthRegularization;
        cachedVisibilityCoef = weightVisibility;

        // Length term
        SparseMatrix<double> lengthMatExplicit;
        SparseMatrix<double> bilapMatExplicit;
        if(localScaleLengthRegularization) {
            // lengthMatExplicit = weightLengthRegularization / std::pow(lengthScale, 4) * getScaleVector(true, 3).asDiagonal() * zeroFormLaplacian;
            lengthMatExplicit = weightLengthRegularization * getScaleVector(true, -1).asDiagonal() * zeroFormLaplacian;
        } else {
            lengthMatExplicit = weightLengthRegularization / lengthScale * zeroFormLaplacian;
        }
        bilapMatExplicit = weightBilapRegularization * lengthScale * zeroFormLaplacian * zeroFormLaplacian;

        // Visibility term
        SparseMatrix<double> visibilityMatExplicit;
        if(weightVisibility == 0.0) {
            visibilityMatExplicit = SparseMatrix<double>(nVert, nVert);
        } else {
            visibilityMatExplicit = weightVisibility / lengthScale * visibilityFraction.toVector().asDiagonal() * zeroFormLaplacian;
        }

        explicitMat = 2.0 * (lengthMatExplicit + bilapMatExplicit + visibilityMatExplicit);
        
        // Invalidate implicitMat cache
        implicitMatCache.clear();
    }
      
}

void EulerianShapeOptimizer::takeGradientStep(VertexData<LabelVec> gradient, double deltaTParam) {

    // Backward Euler
    double deltaT = deltaTParam * lengthScale * lengthScale * lengthScale / 1000;

    
    // If desired, shrink the stepsize to ensure the steps aren't too large
    // Note: ignores implicit regularizer
    // Note: "too large" is determined by mean edge length 
    if(spaceFillingCurveMode) {

        double maxStepRate = 0;
        for(VertexPtr v : mesh->vertices()) {
            for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
                double stepRate = std::abs(gradient[v][iRegion]);
                maxStepRate = std::max(maxStepRate, stepRate);  
            }
        }

        // Decrease by factors of two, so the factorization caching still works
        if(deltaT * maxStepRate > stepSizeFactor * meanEdgeLength) {
            cout << "Step would be too large at " << (deltaT * maxStepRate) << ".";
            cout << " Clamping deltaT from " << deltaT << " to ";
            while(deltaT * maxStepRate > stepSizeFactor * meanEdgeLength) {
                deltaT /= 2;
            }
            cout << deltaT << "." << endl;
        }
    }


    // Recompute implicit matrix only if needed, so the factorization gets cached
    updateExplicitOperator();

    // Make sure we have the implicitMat
    // (we cache it to re-use factorizations)
    if(implicitMatCache.find(deltaT) == implicitMatCache.end()) {

        // Operator for the implicit step
        SparseMatrix<double> newImplicitMat = SparseMatrix<double>::identity(nVert) + deltaT*explicitMat;
        implicitMatCache[deltaT] = newImplicitMat;
    }
    SparseMatrix<double>& implicitMat = implicitMatCache[deltaT];

    // Weird hack: GC::SparseMatrix solvers aren't threadsafe on factorization,
    // but they are for solves. Force a single-threaded factorization here
    // so the solves below are safe.
    DenseVector<double> temp1(nVert);
    DenseVector<double> temp2(nVert);
    solveSquare(implicitMat, temp1, temp2);
    


    // Take an implicit step of the regularizer (for each region)
    #pragma omp parallel for
    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {

        // Implicit stepping
        // Take an explicit step of the optimizer to compute the RHS for just this region
        DenseVector<double> regionRHS(nVert);
        for(VertexPtr v : mesh->vertices()) {
            regionRHS(vInd[v]) = phi[v][iRegion] - deltaT * gradient[v][iRegion];
        }

        // Solve
        DenseVector<double> regionPhi(nVert);
        // solve(implicitMat, regionPhi, regionRHS);
        solveSquare(implicitMat, regionPhi, regionRHS);

        // Store
        for(VertexPtr v : mesh->vertices()) {
            phi[v][iRegion] = regionPhi(vInd[v]);
        }

        // Explicit stepping
        /*
        DenseVector<double> regionPhi(nVert);
        for(VertexPtr v : mesh->vertices()) {
            regionPhi(vInd[v]) = phi[v][iRegion];
        }
        DenseVector<double> regStep = explicitMat * regionPhi;

        for(VertexPtr v : mesh->vertices()) {
            phi[v][iRegion] = regionPhi(vInd[v]) - deltaT * (regStep(vInd[v]) + gradient[v][iRegion]);
        }
        */
    }
   
    // Handle symmetry
    if(hasSymmetry) {
        projectSymmetric();
    }

    // Handle pin constraints
    if(pinPatchMode) {
        VertexPtr pinV = mesh->vertex(pinVert); 
       
        // Sloppy projection
        int minRegion;
        phi[pinV].minCoeff(&minRegion);
        if(minRegion != 0) {
            phi[pinV][0] = 2*phi[pinV][minRegion];
        }
    }

    clearCaches();
}

GC::DenseVector<double> EulerianShapeOptimizer::getScaleVector(bool local, int pow) {

    GC::DenseVector<double> scaleVec;

    if(local) {
        scaleVec = featureSize.toVector();
        if(pow != 1) {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = std::pow(scaleVec(i), pow);
            }
        }
    } else {
        scaleVec = GC::DenseVector<double>(nVert); 
        if(pow == 1) {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = lengthScale;
            }
        } else {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = std::pow(lengthScale, pow);
            }
        }
    }

    return scaleVec;
}

GC::DenseVector<double> EulerianShapeOptimizer::getCutMeshScaleVector(bool local, int pow) {

    GC::DenseVector<double> scaleVec;

    if(local) {

        // Extend featuresize to the cut mesh
        std::vector<double> cutMeshFeatureSize = interpolateToCutMesh(featureSize);

        scaleVec = GC::DenseVector<double>(cutMesh->vertices.size()); 
        if(pow == 1) {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = cutMeshFeatureSize[i];
            }
        } else {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = std::pow(cutMeshFeatureSize[i], pow);
            }
        }

    } else {

        scaleVec = GC::DenseVector<double>(cutMesh->vertices.size()); 

        if(pow == 1) {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = lengthScale;
            }
        } else {
            for(size_t i = 0; i < scaleVec.nRows(); i++) {
                scaleVec(i) = std::pow(lengthScale, pow);
            }
        }

    }

    return scaleVec;
}

void EulerianShapeOptimizer::clearCaches() {

    // Everything is now stale
    haveBoundaryGeometry = false;
    haveTriangleSoupCutMeshes = false;
    haveYamabeSolution = false;
    haveVisibilityFraction = false;
    haveConnectedComponents = false;
    havePatchArea = false;
    havePatchNormals = false;
    haveCutHalfedgeMesh = false;
    havePatchGeometries = false;
    haveCoarseParameterizedMesh = false;
    haveGeodesicParameterizedMesh = false;
    haveGeodesicParameterizedUVPosMesh = false;
    haveExtrinsicDevelopableMesh = false;

    // Really go the extra mile and free memory
    safeDelete(cutMesh);
    deleteGeometryAndMesh(cutMeshGeometry);
    deleteGeometryAndMesh(geodesicCutMeshGeometry);
    deleteGeometryAndMesh(geodesicCutMeshUVPosGeometry);
    deleteGeometryAndMesh(extrinsicDevelopableMeshGeometry);

    for(size_t i = 0; i < patches.size(); i++) {
        safeDelete(patches[i]);
    }
    patches.clear();
    nConnectedComponents = 0;
}

void EulerianShapeOptimizer::constructCutHalfedgeMesh() {
    HalfedgeMesh* temp;
    cutMesh->toHalfedgeMesh(temp, cutMeshGeometry);
}

void EulerianShapeOptimizer::constructPatchGeometries() {
    for(SurfacePatch* patch : patches) {
        patch->generateGeometry();
    }
}


// This version does not do geodesic smoothing on additional cuts, but conveniently returns its results on the cutMeshGeometry
void EulerianShapeOptimizer::computeCoarseFlattenedGeometry() {

    typedef std::pair<Vector3, Vector3> PointPair;

    HalfedgeMesh* cutHalfedgeMesh = cutMeshGeometry->getMesh();


    // Build a halfedge hashmap to transfer data
    std::unordered_map<PointPair, HalfedgePtr> halfedgeLookup;
    for(HalfedgePtr he : cutHalfedgeMesh->halfedges()) {
        PointPair p = std::make_pair(cutMeshGeometry->position(he.vertex()), cutMeshGeometry->position(he.twin().vertex()));
        halfedgeLookup[p] = he;
    }

    // Flatten each patch and transfer the result back to the whole mesh
    coarseGeometryCoords = CornerData<Vector2>(cutHalfedgeMesh);
    for(SurfacePatch* p : patches) {

        Geometry<Euclidean>* geom = p->geometry;

        // Further cut the region meshes to ensure that they are disks
        double addedLength = cutRegionToDisk(geom, &symmetry);
        CornerData<Vector2> param = parameterize(geom);
        // CornerData<Vector2> param(geom->getMesh(), Vector2{0.0, 0.0});

        // Transfer to the original mesh
        for(HalfedgePtr he : geom->getMesh()->halfedges()) {
            PointPair p = std::make_pair(geom->position(he.vertex()), geom->position(he.twin().vertex()));
            HalfedgePtr origHe = halfedgeLookup[p];

            if(he.edge().isCut() || (he.edge().isBoundary() && !origHe.edge().isBoundary())) {
                origHe.edge().markCut(true);
            }

            coarseGeometryCoords[origHe.corner()] = param[he.corner()];
        }
    }
    
    // Translate/flip the parameterizations to pack injectively
    // packParameterization(cutMeshGeometry, coarseGeometryCoords); 
}


// Helper function to union meshes

// This version DOES do geodesic smoothing on additional cuts, but inconveniently returns its results on a mesh which is new and only appears here
void EulerianShapeOptimizer::computeGeodesicFlattenedGeometry() {

    geodesicCutMeshGeometry = generateParameterizedNiceGeometry(geodesicGeometryCoords);
}

// Very similar to above version, except returns a mesh that has UVs in the x,y position coordates
void EulerianShapeOptimizer::computeGeodesicFlattenedGeometryUVPos() {

    
    geodesicCutMeshUVPosGeometry = generateParameterizedNiceGeometry(geodesicGeometryCoordsUVPos, true);

    // At this point each vertex should have a single set of texture coordinates, store them as positions
    for(VertexPtr v : geodesicCutMeshUVPosGeometry->getMesh()->vertices()) {
        Vector2 coord = geodesicGeometryCoordsUVPos[v.halfedge().next().corner()];
        (*geodesicCutMeshUVPosGeometry)[v] = Vector3{coord.x, coord.y, 0};
    }

}

Geometry<Euclidean>* EulerianShapeOptimizer::generateParameterizedNiceGeometry(CornerData<Vector2>& paramCoords, bool splitPatchSeams) {

    // Steps:
    //   1) Cut patches to disk
    //   2) Union patches to single mesh
    //   3) Remesh (maybe)
    //   4) Split to components
    //   5) Parameterize each component
    //   6) Transfer parameterization back to mesh

    // == 1 == Smooth out the added cuts geodesics and remesh accordingly
    std::vector<Geometry<Euclidean>*> geodesicCutGeoms;
    for(SurfacePatch* p : patches) {

        Geometry<Euclidean>* geom = p->geometry;

        // Further cut the region meshes to ensure that they are disks
        double addedLength = cutRegionToDisk(geom, &symmetry);
        // double addedLength = 0;

        // Straighten the cuts and get the new corresponding mesh
        Geometry<Euclidean>* newGeom = remeshAlongCutGeodesics(geom);

        geodesicCutGeoms.push_back(newGeom);
    }

    // == 2 == Union patches to single mesh
    Geometry<Euclidean>* initUnionGeom;
    if(splitPatchSeams) {
        initUnionGeom = unionMeshesDuplicateEdgesAndCuts(geodesicCutGeoms);
    } else {
        initUnionGeom = unionMeshes(geodesicCutGeoms);
    }

    // == 3 == Remesh to eliminate super bad triangles
    // NOTE: not doing this for now
    // Geometry<Euclidean>* remeshedUnionGeom = improveMeshPreserveCut(initUnionGeom);
    Geometry<Euclidean>* remeshedUnionGeom = initUnionGeom;

    // == 4 == Split to components
    std::vector<Geometry<Euclidean>*> splitGeoms = splitAlongCuts(remeshedUnionGeom);

    // Build map to transfer paramereizations
    typedef std::pair<Vector3, Vector3> PointPair;
    std::unordered_map<PointPair, HalfedgePtr> halfedgeLookup;
    for(HalfedgePtr he : remeshedUnionGeom->getMesh()->halfedges()) {
        PointPair p = std::make_pair(remeshedUnionGeom->position(he.vertex()), remeshedUnionGeom->position(he.twin().vertex()));
        halfedgeLookup[p] = he;
    }

    // Parameterize and transfer result
    paramCoords = CornerData<Vector2>(remeshedUnionGeom->getMesh());
    for(Geometry<Euclidean>* splitGeom : splitGeoms) {

        // == 5 == Parameterize each component
        CornerData<Vector2> param = parameterize(splitGeom);
    
        
        // == 6 == Transfer parameterization back
        for(HalfedgePtr he : splitGeom->getMesh()->halfedges()) {
            PointPair p = std::make_pair(splitGeom->position(he.vertex()), splitGeom->position(he.twin().vertex()));
            
            HalfedgePtr origHe = halfedgeLookup[p];
            if(he.edge().isCut() || (he.edge().isBoundary() && !origHe.edge().isBoundary())) {
                origHe.edge().markCut(true);
            }

            paramCoords[origHe.corner()] = param[he.corner()];
        }
    }
    
    packParameterization(remeshedUnionGeom, paramCoords); 

    // Leak slightly less memory
    if(initUnionGeom != remeshedUnionGeom) {
        // Conditionally delete because some configurations of the above method
        // actually perform remeshing and allocate a new mesh, whereas for others
        // the "remeshed" geom is just a reference to the same mesh.
        deleteGeometryAndMesh(initUnionGeom);
    }
    for(auto x : geodesicCutGeoms) {
        deleteGeometryAndMesh(x);
    }
    for(auto x : splitGeoms) {
        deleteGeometryAndMesh(x);
    }
    
    return remeshedUnionGeom;
}


void EulerianShapeOptimizer::computeExtrinsicDevelopableMesh() { 

    // Smooth out the added cuts geodesics and remesh accordingly
    // shares work with computeGeodesicFlattenedGeometry()

    // == 1 == Smooth out the added cuts geodesics and remesh accordingly
    std::vector<Geometry<Euclidean>*> geodesicCutGeoms;
    for(SurfacePatch* p : patches) {

        Geometry<Euclidean>* geom = p->geometry;

        // Further cut the region meshes to ensure that they are disks
        double addedLength = cutRegionToDisk(geom, &symmetry);

        // Straighten the cuts and get the new corresponding mesh
        Geometry<Euclidean>* newGeom = remeshAlongCutGeodesics(geom);

        geodesicCutGeoms.push_back(newGeom);
    }

    // == 2 == Union patches to single mesh
    Geometry<Euclidean>* initUnionGeom = unionMeshes(geodesicCutGeoms);

    // == 3 == Remesh to eliminate super bad triangles
    Geometry<Euclidean>* remeshedUnionGeom = improveMeshPreserveCut(initUnionGeom);
    // Geometry<Euclidean>* remeshedUnionGeom = initUnionGeom;

    // == 4 == Split to components
    std::vector<Geometry<Euclidean>*> splitGeoms = splitAlongCuts(remeshedUnionGeom);

    // Build map to transfer lengths
    typedef std::pair<Vector3, Vector3> PointPair;
    std::unordered_map<PointPair, HalfedgePtr> halfedgeLookup;
    for(HalfedgePtr he : remeshedUnionGeom->getMesh()->halfedges()) {
        PointPair p = std::make_pair(remeshedUnionGeom->position(he.vertex()), remeshedUnionGeom->position(he.twin().vertex()));
        halfedgeLookup[p] = he;
    }

    // Parameterize and transfer result
    CornerData<Vector2> paramCoords(remeshedUnionGeom->getMesh());
    for(Geometry<Euclidean>* splitGeom : splitGeoms) {

        // == 5 == Parameterize each component
        CornerData<Vector2> param = parameterize(splitGeom);
    
        
        // == 6 == Transfer parameterization back
        for(HalfedgePtr he : splitGeom->getMesh()->halfedges()) {
            PointPair p = std::make_pair(splitGeom->position(he.vertex()), splitGeom->position(he.twin().vertex()));
            
            HalfedgePtr origHe = halfedgeLookup[p];
            if(he.edge().isCut() || (he.edge().isBoundary() && !origHe.edge().isBoundary())) {
                origHe.edge().markCut(true);
            }

            paramCoords[origHe.corner()] = param[he.corner()];
        }
    }
    
    packParameterization(remeshedUnionGeom, paramCoords); 
    extrinsicDevelopableMeshGeometry = remeshedUnionGeom;

    EdgeData<double> flatLengths(extrinsicDevelopableMeshGeometry->getMesh());
    for(EdgePtr e : extrinsicDevelopableMeshGeometry->getMesh()->edges()) {
        
        HalfedgePtr he = e.halfedge();

        Vector2 c1 = paramCoords[he.next().corner()];
        Vector2 c2 = paramCoords[he.next().next().corner()];
        double l = norm(c1 - c2);

        // Compute on both sides of the edge to account for the fact that at boundary edges,
        // the parameterization scheme probably did not yield an exactly isometric output, 
        // and thus the values differ on either side
        if(!e.isBoundary()) {

            HalfedgePtr het = he.twin();
            Vector2 c1 = paramCoords[het.next().corner()];
            Vector2 c2 = paramCoords[het.next().next().corner()];
            double l2 = norm(c1 - c2);
            l = 0.5 * (l + l2);
        }

        flatLengths[e] = l;
    }
    

    // Make extrinsic developable 
    makeShapeExtrinsicDevelopable(extrinsicDevelopableMeshGeometry, flatLengths);

    // Leak slightly less memory
    if(initUnionGeom != remeshedUnionGeom) {
        // Conditionally delete because some configurations of the above method
        // actually perform remeshing and allocate a new mesh, whereas for others
        // the "remeshed" geom is just a reference to the same mesh.
        deleteGeometryAndMesh(initUnionGeom);
    }
    for(auto x : geodesicCutGeoms) {
        deleteGeometryAndMesh(x);
    }
    for(auto x : splitGeoms) {
        deleteGeometryAndMesh(x);
    }
    
    /*
    std::vector<Geometry<Euclidean>*> subGeoms;

    // Build a halfedge hashmap to transfer data
    typedef std::pair<Vector3, Vector3> PointPair;
    std::unordered_map<PointPair, double> targetLen;

    for(SurfacePatch* p : patches) {

        Geometry<Euclidean>* geom = p->geometry;

        // Further cut the region meshes to ensure that they are disks
        double addedLength = cutRegionToDisk(geom, &symmetry);

        // Straighten the cuts and get the new corresponding mesh
        Geometry<Euclidean>* newGeom = remeshAlongCutGeodesics(geom);

        // Find the target length in the flattened domain
        BoundaryFirstFlattening bff(newGeom);
        BoundaryConstraints constraints(&bff);
        constraints.setIsometricConstraints();
        bff.flatten();
        CornerData<Vector2> param = bff.uvs;
        canonicalizeParameterization(newGeom, param);
        
        for(EdgePtr e : newGeom->getMesh()->edges()) {
    
            HalfedgePtr he = e.halfedge();
    
            Vector2 c1 = param[he.next().corner()];
            Vector2 c2 = param[he.next().next().corner()];
            double l = norm(c1 - c2);

            // Store target length to transfer 
            PointPair p1 = std::make_pair(newGeom->position(he.vertex()), newGeom->position(he.twin().vertex()));
            targetLen[p1] = l;
            PointPair p2 = std::make_pair(newGeom->position(he.twin().vertex()), newGeom->position(he.vertex()));
            targetLen[p2] = l;
        }

        // Make the patch developable
        // Geometry<Euclidean>* devGeom = makePatchExtrinsicDevelopable(newGeom);
        // subGeoms.push_back(devGeom);
        // deleteGeometryAndMesh(newGeom);

        subGeoms.push_back(newGeom);

    }


    // Build a union mesh of all of the geometries
    // extrinsicDevelopableMeshGeometry = unionMeshesDuplicateEdges(subGeoms);
    extrinsicDevelopableMeshGeometry = unionMeshes(subGeoms);

    // Restore lenghts
    EdgeData<double> flatLengths(extrinsicDevelopableMeshGeometry->getMesh());
    for(EdgePtr e : extrinsicDevelopableMeshGeometry->getMesh()->edges()) {
        HalfedgePtr he = e.halfedge();
        PointPair p1 = std::make_pair(extrinsicDevelopableMeshGeometry->position(he.vertex()), extrinsicDevelopableMeshGeometry->position(he.twin().vertex()));
        flatLengths[e] = targetLen[p1];
    }

    // Make extrinsic developable 
    makeShapeExtrinsicDevelopable(extrinsicDevelopableMeshGeometry, flatLengths);

    // Free old meshes
    for(auto g : subGeoms) {
        deleteGeometryAndMesh(g);
    }
    */
}

// Precondition: All values (u, regions, etc) must be accurate
double EulerianShapeOptimizer::computeEnergy(bool print) {

    double E_length = weightLengthRegularization * evaluateEnergyTermLengthRegularization();
    double E_dirichletDistortion = weightDirichletDistortion * evaluateEnergyTermDirichletDistortion();
    double E_henckyDistortion = 1000 * weightHenckyDistortion * evaluateEnergyTermHenckyDistortion();
    double E_area = weightArea * evaluateEnergyTermArea();
    double E_visibility = weightVisibility * evaluateEnergyTermVisibility();
    double E_normalDeviation = weightNormalDeviation * evaluateEnergyTermNormalDeviation();

    double energy = E_length + E_dirichletDistortion + E_henckyDistortion + 
                    E_area + E_visibility + E_normalDeviation;

    // Write actual energy values
    // if(PLOT_ENERGY && print) {
    //     cout << "[ENERGYPLOT]"
    //          << E_length << ","
    //          << E_dirichletDistortion << ","
    //          << E_henckyDistortion << ","
    //          << E_area << ","
    //          << E_visibility << ","
    //          << E_normalDeviation 
    //          << endl;
    // }

    // Write raw values
    if(PLOT_ENERGY && print) {
        cout << "[ENERGYPLOT]"
             << E_length/weightLengthRegularization << ","
             << E_dirichletDistortion/weightDirichletDistortion << ","
             << E_henckyDistortion/weightHenckyDistortion << ","
             << E_area/weightArea << ","
             << E_visibility/weightVisibility << ","
             << E_normalDeviation/weightNormalDeviation 
             << endl;
    }
    

    return energy;
}
    
double EulerianShapeOptimizer::evaluateEnergyTermLengthRegularization() {

    if(weightLengthRegularization == 0.0) return 0.0;
    
    ensureHaveBoundaryGeometry();
    ensureHaveTriangleSoupCutMeshes();
    
    double E = 0;
    for(SurfacePatch* p : patches) {
        if(localScaleLengthRegularization) {
            E += p->computeLocalScaledBoundaryLengthEnergy();
        } else {
            E += p->computeBoundaryLengthEnergy();
        }
    }

    return E;
}

double EulerianShapeOptimizer::evaluateEnergyTermDirichletDistortion() {
    
    if(weightDirichletDistortion == 0.0) return 0.0;
    
    ensureHaveYamabeSolution();

    double E = 0;
    for(SurfacePatch* p : patches) {
        if(localScaleDirichletDistortion) {
            E += p->computeLocalScaledDirichletDistortionEnergy();
        } else {
            E += p->computeDirichletDistortionEnergy();
        }
    }

    return E;
}

double EulerianShapeOptimizer::evaluateEnergyTermHenckyDistortion() {
    
    if(weightHenckyDistortion == 0.0) return 0.0;
    
    ensureHaveYamabeSolution();

    double E = 0;
    for(SurfacePatch* p : patches) {
        E += p->computeHenckyDistortionEnergy();
    }

    return E;
}

double EulerianShapeOptimizer::evaluateEnergyTermArea() {
    
    if(weightArea == 0.0) return 0.0;

    ensureHavePatchArea();

    double E = 0;
    if(patchZeroTargetArea == -1) {
        for(SurfacePatch* p : patches) {
            E += p->computeAreaEnergy();
        }
    } else {
        for(SurfacePatch* p : patches) {
            if(p->iRegion == 0) {
                E += p->computeAreaEnergy(patchZeroTargetArea);
            }
        }
    }
    
    return E;
}

double EulerianShapeOptimizer::evaluateEnergyTermVisibility() {
    
    if(weightVisibility == 0.0) return 0.0;
    
    ensureHaveVisibilityFraction();
    ensureHaveTriangleSoupCutMeshes();

    double E = 0; 
    for(SurfacePatch* p : patches) {
        E += p->computeVisibilityEnergy(cutMeshVisibility);
    }

    return E;
}

double EulerianShapeOptimizer::evaluateEnergyTermNormalDeviation() {
    
    if(weightNormalDeviation == 0.0) return 0.0;
    
    ensureHavePatchNormals();
    ensureHaveTriangleSoupCutMeshes();

    double E = 0;
    for(SurfacePatch* p : patches) {
        E += p->computeNormalDeviationEnergy();
    }

    return E;
}


CornerData<Vector2> EulerianShapeOptimizer::parameterize(Geometry<Euclidean>* geom) {

    cout << "Flattening..." << endl;

    BoundaryFirstFlattening bff(geom);
    BoundaryConstraints constraints(&bff);
    constraints.setIsometricConstraints();
    bff.flatten();
    CornerData<Vector2> param = bff.uvs;

    /*
    LeastSquaresConformalMaps lscm(geom);
    lscm.flatten();
    CornerData<Vector2> param = lscm.uvs;
    */
   
    /*
    ConformalEquivalence cetm(geom);
    cetm.flatten();
    CornerData<Vector2> param = cetm.uvs;
    */

    canonicalizeParameterization(geom, param); 

    return param;
}

void EulerianShapeOptimizer::saveToFile(std::string filename) {

    cout << "Saving state to file " << filename << endl;

    std::ofstream outFile(filename);

    // Header indicating number of vertices and number of regions
    outFile << nVert << " " << K_REGIONS << endl;

    // Regions for every vertex
    for(VertexPtr v : mesh->vertices()) {
        for(size_t i = 0; i < K_REGIONS; i++) {

            if(std::isfinite(phi[v][i])) {
                outFile << phi[v][i];
            } else {
                // If we print 'inf', it's hard to read back in
                outFile << 1e5 * lengthScale;
            }

            if(i == K_REGIONS-1) {
                outFile << endl;
            } else {
                outFile << " ";
            }
        }
    }

    outFile.close();    
}

void EulerianShapeOptimizer::loadFromFile(std::string filename) {

    cout << "Reading state from file " << filename << endl;

    std::ifstream inFile(filename);

    // Check that the file exists
    if(!inFile) {
        throw std::runtime_error("Cannot open file " + filename);
    }

    // Parse the header
    size_t fileVert, fileRegion;
    inFile >> fileVert;
    inFile >> fileRegion;

    // Check if the header matches up
    if(fileVert != nVert) {
        throw std::runtime_error("# vertices in file (" + std::to_string(fileVert) + ") does not match number in current mesh (" + std::to_string(nVert) + ")");
    } 
    if(fileRegion > K_REGIONS) {
        throw std::runtime_error("# regions in file (" + std::to_string(fileRegion) + ") is greater than number supported (" + std::to_string(K_REGIONS) + ")");
    } 

    // Read the data
    VertexData<LabelVec> filePhi(mesh, LabelVec::Zero());
    for(VertexPtr v : mesh->vertices()) {
        for(size_t i = 0; i < fileRegion; i++) {
            inFile >> filePhi[v][i];
        }
    }
    inFile.close(); 

   // Prepare various other state fo4 the optimizer
   setState(filePhi);
}

void EulerianShapeOptimizer::saveDistortionHistogram(std::string filename) {

    ensureHaveGeodesicParameterizedMesh();

    cout << "Writing distortion histogram to file " << filename << endl;

    std::ofstream outFile(filename);
        
    outFile << "#TITLE:Length distortion" << endl;
    outFile << "#XLABEL:Distortion factor $e^u$" << endl;


    double maxDistortion = -std::numeric_limits<double>::infinity();
    double minDistortion = std::numeric_limits<double>::infinity();
    double maxOldDistortion = -std::numeric_limits<double>::infinity();
    double minOldDistortion = std::numeric_limits<double>::infinity();
    double l2StretchSum = 0;
    double l2StretchWeightSum = 0;
    for(FacePtr f : geodesicCutMeshGeometry->getMesh()->faces()) {

        // UV coordinates
        Vector2 p1 = geodesicGeometryCoords[f.halfedge().corner()];
        Vector2 p2 = geodesicGeometryCoords[f.halfedge().next().corner()];
        Vector2 p3 = geodesicGeometryCoords[f.halfedge().next().next().corner()];

        // Spatial coordinates
        Vector3 q1 = geodesicCutMeshGeometry->position(f.halfedge().next().next().vertex());
        Vector3 q2 = geodesicCutMeshGeometry->position(f.halfedge().vertex());
        Vector3 q3 = geodesicCutMeshGeometry->position(f.halfedge().next().vertex());

        // Area in original space
        double origArea = geodesicCutMeshGeometry->area(f);        

        // Distortion as the square root of area
        // Area in param space
        double paramArea = 0.5 * std::abs(cross(p2 - p1, p3 - p1).z);
        double oldLengthDistortion = std::sqrt(paramArea / origArea);
        maxOldDistortion = std::max(maxOldDistortion, oldLengthDistortion);
        minOldDistortion = std::min(minOldDistortion, oldLengthDistortion);


        // L2 "Texture Stretch Metric" from "Texture Mapping Progressive Meshes"
        double A = ((p2.x-p1.x)*(p3.y-p1.y) - (p3.x-p1.x)*(p2.y-p1.y)) / 2.0;
        Vector3 Ss = (q1*(p2.y-p3.y) + q2*(p3.y-p1.y) + q3*(p1.y-p2.y)) / (2 * A);
        Vector3 St = (q1*(p3.x-p2.x) + q2*(p1.x-p3.x) + q3*(p2.x-p1.x)) / (2 * A);
        double a = dot(Ss,Ss);
        double b = dot(Ss,St);
        double c = dot(St,St);
        double Gamma = std::sqrt(((a+c) + std::sqrt((a-c)*(a-c) + 4*b*b) ) / 2.0);
        double gamma = std::sqrt(((a+c) - std::sqrt((a-c)*(a-c) + 4*b*b) ) / 2.0);
        double L2stretch = std::sqrt((Gamma*Gamma + gamma*gamma)/2.0);
        l2StretchSum += L2stretch * L2stretch * origArea;
        l2StretchWeightSum += origArea;

        double lengthDistortion = L2stretch;

        // cout << "Old dist = " << oldLengthDistortion << " new dist = " << lengthDistortion << endl;

        // Write to file
        outFile << lengthDistortion << "," << origArea << endl;
        maxDistortion = std::max(maxDistortion, lengthDistortion);
        minDistortion = std::min(minDistortion, lengthDistortion);
    }

    /* Estimate distortion from the scale factors
    for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
        outFile << std::exp(distortion[i]) << "," << cutMesh->dualArea[i] << endl;
    }
    */

    double l2Stretch = std::sqrt(l2StretchSum/l2StretchWeightSum);
    cout << "L2 stretch total norm = " << l2Stretch << endl;
    outFile << "# L2 stretch total norm: " << l2Stretch << endl;

    cout << "Max AREA-ROOT distortion = " << maxOldDistortion << "  Min AREA-ROOT distortion = " << minOldDistortion << endl;
    cout << "Max distortion = " << maxDistortion << "  Min distortion = " << minDistortion << endl;
    cout << "Finished writing distortion histogram to file " << filename << endl;
    outFile.close();    
}

void EulerianShapeOptimizer::saveCutLines(std::string filename) {

    cout << "Writing cut lines to file " << filename << endl;
    
    // Note: This gets the lines from a parameterization rather than
    // just using the getBoundaryLines() method in order to ensure that
    // it picks up any extra geodesic cuts that might have been added
    // to make regions disk-like.
    
    // Get the lines 
    CornerData<Vector2> param;
    Geometry<Euclidean>* flatGeom = getGeodesicFlattenedGeometry(param);
    // Geometry<Euclidean>* flatGeom = getCoarseFlattenedGeometry(param);
    HalfedgeMesh* flatMesh = flatGeom->getMesh();

    // Create a list of points and lines
    std::vector<Vector3> positions;
    size_t nPos = 0;
    VertexData<size_t> inds(flatMesh);
    for(VertexPtr v : flatMesh->vertices()) {

        bool isOnCut = false;
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) {
                isOnCut = true;
            }
        }

        if(isOnCut) {
            inds[v] = nPos++;
            positions.push_back(flatGeom->position(v));
        }
    }

    std::vector<std::array<size_t, 2>> lines;
    for(EdgePtr e : flatMesh->edges()) {
        if(e.isCut()) {
            lines.push_back({{inds[e.halfedge().vertex()], inds[e.halfedge().twin().vertex()]}});
        }
    }


    // Write the file
    std::ofstream outFile(filename);
        
    outFile << "# points: " << nPos << endl;
    for(Vector3 p : positions) {
        outFile << p << endl;
    }
    
    outFile << "# lines: " << lines.size() << endl;
    for(auto l : lines) {
        outFile << l[0] << "," << l[1] << endl;
    }

    outFile.close();    

}


VertexData<double> EulerianShapeOptimizer::getPhiK(int iRegion) {

    ensureHaveBoundaryGeometry();

    VertexData<double> p(mesh);
    double thresh = 8 * meanEdgeLength;
    for(VertexPtr v : mesh->vertices()) {
        p[v] = clamp(phi[v][iRegion], -thresh, thresh);
    }

    return p;
}

VertexData<Eigen::VectorXd> EulerianShapeOptimizer::getRegionValues() {
    ensureHaveBoundaryGeometry();
    return fixed2Dynamic(phi);
}



FaceData<double> EulerianShapeOptimizer::getPatchArea() {
    
    ensureHavePatchArea();
    ensureHaveCutHalfedgeMesh();

    HalfedgeMesh* cutMesh = cutMeshGeometry->getMesh();
    FaceData<double> result(cutMesh);
    FaceData<size_t> cutMeshFInd = cutMesh->getFaceIndices();
    for(SurfacePatch* p : patches) {
        for(size_t iTri = 0; iTri < p->soupMesh.triangles.size(); iTri++) {
            result[cutMesh->face(p->parentFace[iTri])] = p->area;
        }
    }

    return result;
}

FaceData<double> EulerianShapeOptimizer::getPatchNormalDeviation() {
    
    ensureHavePatchNormals();
    ensureHaveCutHalfedgeMesh();
    
    HalfedgeMesh* cutMesh = cutMeshGeometry->getMesh();
    FaceData<double> result(cutMesh);
    FaceData<size_t> cutMeshFInd = cutMesh->getFaceIndices();
    for(SurfacePatch* p : patches) {
        for(size_t iTri = 0; iTri < p->soupMesh.triangles.size(); iTri++) {
            result[cutMesh->face(p->parentFace[iTri])] = p->normalDeviation[iTri];
        }
    }

    return result;
}
    
VertexData<double> EulerianShapeOptimizer::getVisibility() {
    ensureHaveVisibilityFraction();
    return visibilityFraction;
}

VertexData<double> EulerianShapeOptimizer::getDistortion() {

    ensureHaveYamabeSolution();

    // Get the mesh
    HalfedgeMesh* cutMesh = getCustomGeometry()->getMesh();

    // Transfer to a VertexData object
    double maxDist = -std::numeric_limits<double>::infinity();
    double minDist = std::numeric_limits<double>::infinity();
    VertexData<double> uVals(cutMesh, 0.0);
    for(SurfacePatch* p : patches) {
        for(size_t i = 0; i < p->soupMesh.vertices.size(); i++) {
            uVals[cutMesh->vertex(p->parentVertex[i])] = p->distortion[i];

            maxDist = std::max(maxDist, p->distortion[i]);
            minDist = std::min(minDist, p->distortion[i]);
        }
    }

    cout << "u limits: " << minDist << "  -  " << maxDist << endl; 

    return uVals;
}

VertexData<double> EulerianShapeOptimizer::getScaledDistortion() {

    ensureHaveYamabeSolution();

    // Get the mesh
    HalfedgeMesh* cutMesh = getCustomGeometry()->getMesh();

    // Transfer to a VertexData object
    double maxDist = -std::numeric_limits<double>::infinity();
    double minDist = std::numeric_limits<double>::infinity();
    VertexData<double> uVals(cutMesh, 0.0);
    for(SurfacePatch* p : patches) {
        for(size_t i = 0; i < p->soupMesh.vertices.size(); i++) {
            double scaledU = std::pow(p->localScaleFactor[i], 3) * p->distortion[i];
            uVals[cutMesh->vertex(p->parentVertex[i])] = scaledU;

            maxDist = std::max(maxDist, scaledU);
            minDist = std::min(minDist, scaledU);
        }
    }

    cout << "scaled u limits: " << minDist << "  -  " << maxDist << endl; 

    return uVals;
}
    
VertexData<double> EulerianShapeOptimizer::getRegionShapeOptimizerMotion(int iRegion) {

    iRegion = clamp(iRegion, 0, K_REGIONS-1);

    ensureHaveCutHalfedgeMesh();
    HalfedgeMesh* cutHalfedgeMesh = getCustomGeometry()->getMesh();

    buildGradientAtBoundary();
    VertexData<LabelVec> surfaceGradient = extendGradientToSurface();

    VertexData<double> val(cutHalfedgeMesh);
    VertexData<size_t> vInd = cutHalfedgeMesh->getVertexIndices();
    for(VertexPtr v : cutHalfedgeMesh->vertices()) {
        val[v] = extendedGradient[vInd[v]][iRegion];
    }

    return val;
}

VertexData<double> EulerianShapeOptimizer::getShapeOptimizerMotion() {

    buildGradientAtBoundary();
    VertexData<LabelVec> surfaceGradient = extendGradientToSurface();

    VertexData<double> val(mesh);
    for(VertexPtr v : mesh->vertices()) {
        val[v] = std::abs(surfaceGradient[v][regionLabels[v]]);
    }

    return val;
}

VertexData<double> EulerianShapeOptimizer::getBoundaryDistance() {

    ensureHaveBoundaryGeometry();

    VertexData<double> result(mesh);
    for(VertexPtr v : mesh->vertices()) {
        result[v] = -phi[v].minCoeff();
    }
    return result;
}

Geometry<Euclidean>* EulerianShapeOptimizer::getCustomGeometry() {

    ensureHaveCutHalfedgeMesh();
    return cutMeshGeometry;
}


VertexData<double> EulerianShapeOptimizer::getDebug() {

    /*
    Geometry<Euclidean>* cutGeom = getCustomGeometry();
    VertexData<double> result(cutGeom->getMesh(), 0);
    for(size_t i = 0; i < cutMesh->nVert; i++) {
        result[cutGeom->getMesh()->vertex(i)] = debugVal[i];
    }
    
    return result;
    */

    return debugValVert;
}

Geometry<Euclidean>* EulerianShapeOptimizer::getCoarseFlattenedGeometry(CornerData<Vector2>& paramCoords) {

    ensureHaveCoarseParameterizedMesh();

    paramCoords = coarseGeometryCoords;
    return cutMeshGeometry;
}

Geometry<Euclidean>* EulerianShapeOptimizer::getGeodesicFlattenedGeometry(CornerData<Vector2>& paramCoords) {

    ensureHaveGeodesicParameterizedMesh();

    paramCoords = geodesicGeometryCoords;
    return geodesicCutMeshGeometry;
}

Geometry<Euclidean>* EulerianShapeOptimizer::getGeodesicFlattenedUVPosGeometry(CornerData<Vector2>& paramCoords) {

    ensureHaveGeodesicParameterizedUVPosMesh();

    paramCoords = geodesicGeometryCoordsUVPos;
    return geodesicCutMeshUVPosGeometry;
}

Geometry<Euclidean>* EulerianShapeOptimizer::getExtrinsicDevelopableGeometry() {

    ensureHaveExtrinsicDevelopableMesh();

    return extrinsicDevelopableMeshGeometry;
}

std::vector<std::array<Vector3, 2>> EulerianShapeOptimizer::getExtraCutLines() {

    ensureHavePatchGeometries();

    std::vector<std::array<Vector3, 2>> lines;

    // Flatten each sub-mesh and transfer the result back to the whole mesh
    for(SurfacePatch* p : patches) {

        Geometry<Euclidean>* geom = p->geometry;

        cutRegionToDisk(geom, &symmetry);
        std::vector<std::array<Vector3, 2>> thisLines = getGeodesicCuts(geom);
        lines.insert(lines.end(), thisLines.begin(), thisLines.end());
    }

    return lines;
}
    


Geometry<Euclidean>* EulerianShapeOptimizer::getPatchGeometry(int iPatch) {

    ensureHavePatchGeometries();
    iPatch = clamp(iPatch, 0, (int)patches.size()-1);
    
    return patches[iPatch]->geometry;
}



std::vector<std::array<Vector3, 2>> EulerianShapeOptimizer::getBoundaryLines() {

    ensureHaveBoundaryGeometry();

    std::vector<std::array<Vector3, 2>> lines;

    for(int iRegion = 0; iRegion < K_REGIONS; iRegion++) {
        for(BoundarySegment &b : boundarySegments[iRegion]) {

            Vector3 vStart = startPoint(b);
            Vector3 vEnd = endPoint(b);

            if(b.type == BType::TRIPLE) {
                lines.push_back({{vStart, b.triplePoint}});
                lines.push_back({{b.triplePoint, vEnd}});
            } else {
                lines.push_back({{vStart, vEnd}});
            }
        }
    }

    return lines;
}

VertexData<Eigen::VectorXd> EulerianShapeOptimizer::fixed2Dynamic(const VertexData<LabelVec>& input) {
    VertexData<Eigen::VectorXd> result(mesh);
    for(VertexPtr v : mesh->vertices()) {
        result[v] = input[v];
    }
    return result; 
}

SparseMatrix<double> EulerianShapeOptimizer::clampedDiagonalInverse(SparseMatrix<double>& A, double thresh) {

    size_t N = A.nColumns(); 
    SparseMatrix<double> Ainv(N, N);
    int count = 0;

    for(size_t i = 0; i < N; i++) {

        double Aval = A(i,i);
        if(std::abs(Aval) < thresh) {
            if(Aval >= 0) {
                Ainv(i, i) = 1.0 / thresh;
            } else {
                Ainv(i, i) = -1.0 / thresh;
            }
            count++;
        } else {
            Ainv(i, i) = 1.0 / Aval;
        }
    }

    if(count > 0) {
        cout << "WARNING: Diagonal inverse clamped " << (count * 100.0 / N) << "% (" << count << " total) of entries." << endl;
    }

    return Ainv;
}

VertexData<double> EulerianShapeOptimizer::generateLocationTerm() {

    // Run FMM to generate distance from point
    vector<std::pair<VertexPtr, double>> boundaryDistances = {std::make_pair(mesh->vertex(pinVert), 0)};
    VertexData<double> distResult = FMMDistance(mesh, boundaryDistances, edgeLengths, oppAngles);
   
    // Generate a falloff function that is large near the point and small away.
    // Make it negative to encourage the patch to stay near here, rather than discourage
    for(VertexPtr v : mesh->vertices()) {
        distResult[v] = -1.0 / (1.0 + 5.0 * targetBoundaryDist * distResult[v]);
    }

    return distResult;
}

double EulerianShapeOptimizer::computeCrossing(double sourceHigh, double sourceLow, double targetHigh, double targetLow) {
    // double t1 = clamp(sourceHigh / (sourceHigh - targetLow) , 0.0, 1.0);
    // double t2 = clamp(sourceLow / (sourceLow - targetHigh) , 0.0, 1.0);
    // return 0.5 * (t1 + t2);
    double t = (sourceHigh - sourceLow) / (sourceHigh - sourceLow + targetHigh - targetLow);
    if(t <= 0 || t >= 1) {
        cout << "BAD CROSSING VALUES t = " << t << endl;
    }
    t = clamp(t, 1e-3, 1.0-1e-3);

    return t;
}

// Given 3 smoothed indicator values on each of 3 vertices of a  triangle, returns the triple point where they are all equal
// (copied from beautiful mathematica, obviously)
// ONEDAY: do actual neat math instead of this
Vector3 EulerianShapeOptimizer::triplePointBarycentrics(Vector3 A, Vector3 B, Vector3 C) {
    double w0 = -((-B[2] *C[1]+A[2]*(-B[1]+C[1])+A[1]*(B[2]-C[2])+B[1]*C[2])/((B[0]-B[2]-C[0]+C[2])*(A[1]-A[2]-C[1]+C[2])+(-A[0]+A[2]+C[0]-C[2])*(B[1]-B[2]-C[1]+C[2])));
    double w1 = (A[2]*(B[0]-C[0])+B[2]*C[0]-B[0]*C[2]+A[0]*(-B[2]+C[2]))/(-B[1]*C[0]+B[2]*C[0]+B[0]*C[1]-B[2]*C[1]+A[2]*(B[0]-B[1]-C[0]+C[1])+A[1]*(-B[0]+B[2]+C[0]-C[2])-B[0]*C[2]+B[1]*C[2]+A[0]*(B[1]-B[2]-C[1]+C[2]));
    double w2 = (A[1]*(B[0]-C[0])+B[1]*C[0]-B[0]*C[1]+A[0]*(-B[1]+C[1]))/(B[1]*C[0]-B[2]*C[0]+A[2]*(-B[0]+B[1]+C[0]-C[1])-B[0]*C[1]+B[2]*C[1]+A[0]*(-B[1]+B[2]+C[1]-C[2])+B[0]*C[2]-B[1]*C[2]+A[1]*(B[0]-B[2]-C[0]+C[2]));

    // Clamping and recularization needed because input might not satisfy assumption, causing the point to lie outside the triangle
    double EPS = 1e-3;
    w0 = clamp(w0, EPS, 1.0 - EPS);
    w1 = clamp(w1, EPS, 1.0 - EPS);
    w2 = clamp(w2, EPS, 1.0 - EPS);
    return Vector3{w0, w1, w2} / (w0 + w1 + w2);
}    

double EulerianShapeOptimizer::shortestDistanceFromLine(Vector3 target, Vector3 p1, Vector3 p2) {

    double lineLen2 = norm2(p2-p1);

    // Check degenerate cases
    if(lineLen2 < EPS_POSITION*EPS_POSITION) {
        return norm(target - 0.5*(p1+p2));
    }

    double t = clamp(dot(target - p1, p2 - p1) / lineLen2, 0.0, 1.0);
    Vector3 proj = p1 + t * (p2 - p1);
    return norm(target - proj);
}

Vector3 EulerianShapeOptimizer::startPoint(const BoundarySegment &b) {
    return geometry->position(b.heStart.vertex()) + halfedgeVector[b.heStart]*b.tStart;
}

Vector3 EulerianShapeOptimizer::endPoint(const BoundarySegment &b) {
    return geometry->position(b.heEnd.vertex()) + halfedgeVector[b.heEnd]*b.tEnd;
}
    
double EulerianShapeOptimizer::boundaryLength(const BoundarySegment &b) {
    if(b.type == BType::TRI || b.type == BType::QUAD) {
        return norm(startPoint(b) - endPoint(b));
    } else if(b.type == BType::TRIPLE) {
        return norm(startPoint(b) - b.triplePoint) + norm(b.triplePoint - endPoint(b));
    }
    return -1;
}

    
void EulerianShapeOptimizer::deleteGeometryAndMesh(Geometry<Euclidean>* &geom) {
    if(geom != nullptr) {
        delete geom->getMesh();
        delete geom;
        geom = nullptr;
    }
}
   
// Useful for debugging NaNs
void EulerianShapeOptimizer::invalidValuePanic(std::string sourceName) {

    cout << "ERROR: Invalid value in " << sourceName << endl;
    cout << "Writing state to panic.msdf" << endl;

    saveToFile("panic.msdf");    

    throw std::runtime_error("ERROR: Invalid value in " + sourceName);
}
    
std::vector<double> EulerianShapeOptimizer::interpolateToCutMesh(VertexData<double> vals) {

    std::vector<double> resultSum(cutMesh->vertices.size(), 0.0);
    std::vector<double> resultWeight(cutMesh->vertices.size(), 0.0);
    std::vector<bool> resultKnown(cutMesh->vertices.size(), false);
    size_t nUnknown = cutMesh->nVert;

    // Copy the values that we already know over to the cut mesh
    // NaN marks unknown values
    std::vector<double> interpolatedVals(cutMesh->vertices.size(), std::numeric_limits<double>::quiet_NaN());
    for(VertexPtr v : mesh->vertices()) {
        interpolatedVals[origMeshToCutMesh[v]] = vals[v];
        resultKnown[origMeshToCutMesh[v]] = true;
        nUnknown--;
    }

    // On the new cut vertices, set the value to be the inverse-distance-weighted average
    // of the values at the interior neighbors
    // Note: This weights boundary neighbors by 1/2, but that seems fine
    // Iterate to BFS until all are known
    while(nUnknown > 0) {

        for(size_t iTri = 0; iTri < cutMesh->triangles.size(); iTri++) {
            const Tri& tri = cutMesh->triangles[iTri];
    
            for(int triInd = 0; triInd < 3; triInd++) {
                
                size_t vertI = tri[triInd];
                size_t vertJ = tri[(triInd+1)%3];
    
                double edgeLen = norm(cutMesh->vertices[vertI] - cutMesh->vertices[vertJ]);
                double w = 1.0 / edgeLen;
    
                // Should transfer from j --> i
                if(!resultKnown[vertI] && resultKnown[vertJ]) {
                    resultWeight[vertI] += w;
                    double val = interpolatedVals[vertJ];
                    resultSum[vertI] += val * w;
                }
    
                // Should transfer from i --> j
                if(!resultKnown[vertJ] && resultKnown[vertI]) {
                    resultWeight[vertJ] += w;
                    double val = interpolatedVals[vertI];
                    resultSum[vertJ] += val * w;
                }
            }
        }
    
        bool anyLearned = false;
        for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
            if(!resultKnown[i] && resultWeight[i] > 0) {
                interpolatedVals[i] = resultSum[i] / resultWeight[i];
                resultKnown[i] = true;
                anyLearned = true;
                nUnknown--;
            }
        }

        if(!anyLearned) {
            throw std::runtime_error("Could not interpolate. Some values must be isolated.");
        }

    }

    return interpolatedVals;
}

std::vector<LabelVec> EulerianShapeOptimizer::interpolateToCutMesh(VertexData<LabelVec> vals) {

    std::vector<LabelVec> resultSum(cutMesh->vertices.size(), LabelVec::Zero());
    std::vector<double> resultWeight(cutMesh->vertices.size(), 0.0);
    std::vector<bool> resultKnown(cutMesh->vertices.size(), false);
    size_t nUnknown = cutMesh->nVert;

    // Copy the values that we already know over to the cut mesh
    // NaN marks unknown values
    std::vector<LabelVec> interpolatedVals(cutMesh->vertices.size(), LabelVec::Ones()*std::numeric_limits<double>::quiet_NaN());
    for(VertexPtr v : mesh->vertices()) {
        interpolatedVals[origMeshToCutMesh[v]] = vals[v];
        resultKnown[origMeshToCutMesh[v]] = true;
        nUnknown--;
    }

    // On the new cut vertices, set the value to be the inverse-distance-weighted average
    // of the values at the interior neighbors
    // Note: This weights boundary neighbors by 1/2, but that seems fine
    // Iterate to BFS until all are known
    while(nUnknown > 0) {

        for(size_t iTri = 0; iTri < cutMesh->triangles.size(); iTri++) {
            const Tri& tri = cutMesh->triangles[iTri];
    
            for(int triInd = 0; triInd < 3; triInd++) {
                
                size_t vertI = tri[triInd];
                size_t vertJ = tri[(triInd+1)%3];
    
                double edgeLen = norm(cutMesh->vertices[vertI] - cutMesh->vertices[vertJ]);
                double w = 1.0 / edgeLen;
    
                // Should transfer from j --> i
                if(!resultKnown[vertI] && resultKnown[vertJ]) {
                    resultWeight[vertI] += w;
                    LabelVec val = interpolatedVals[vertJ];
                    resultSum[vertI] += val * w;
                }
    
                // Should transfer from i --> j
                if(!resultKnown[vertJ] && resultKnown[vertI]) {
                    resultWeight[vertJ] += w;
                    LabelVec val = interpolatedVals[vertI];
                    resultSum[vertJ] += val * w;
                }
            }
        }
    
        // Compute the remaining values as the average
        bool anyLearned = false;
        for(size_t i = 0; i < cutMesh->vertices.size(); i++) {
            if(!resultKnown[i] && resultWeight[i] > 0) {
                interpolatedVals[i] = resultSum[i] / resultWeight[i];
                resultKnown[i] = true;
                anyLearned = true;
                nUnknown--;
            }
        }

        if(!anyLearned) {
            throw std::runtime_error("Could not interpolate. Some values must be isolated.");
        }
    }

    return interpolatedVals;
}

#undef PLOT_ENERGY
#undef SAFETY_CHECKS
