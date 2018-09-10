#include "sdf_initialization.h"

#include "d_charts.h"
#include "fast_marching_method.h"

#include <unordered_set>

VertexData<LabelVec> randomMSDF(Geometry<Euclidean>* geometry) {

    HalfedgeMesh* mesh = geometry->getMesh();

    // Randomly choose a region for each vertex
    VertexData<int> randRegionLabels(mesh);
    for(VertexPtr v : mesh->vertices()) {
        randRegionLabels[v] = randomInt(0, K_REGIONS-1);
    }

    return msdfFromLabels(geometry, randRegionLabels);
}

VertexData<LabelVec> singleVertex(Geometry<Euclidean>* geometry, int startingVert) {
    
    HalfedgeMesh* mesh = geometry->getMesh();

    // If no vertex was specified, use the first one
    size_t startingVertInd = 0;
    if(startingVert != -1) {
        startingVertInd = startingVert;
    }

   
    VertexData<int> regionLabels(mesh, 1);
    regionLabels[mesh->vertex(startingVertInd)] = 0;

    return msdfFromLabels(geometry, regionLabels);
}

VertexData<LabelVec> sphereAtMSDF(Geometry<Euclidean>* geometry, Vector3 center, double rad) {
    
    HalfedgeMesh* mesh = geometry->getMesh();


    VertexData<int> regionLabels(mesh, 0);

    for(VertexPtr v : mesh->vertices()) { 

        double dist = norm(geometry->position(v) - center);
        cout << "dist = " << dist << " rad = " << rad << endl;
        if(dist < rad) {
            regionLabels[v] = 1;
        }
    }

    return msdfFromLabels(geometry, regionLabels);
}

/*
VertexData<LabelVec> cubeMSDF(Geometry<Euclidean>* geometry) {

    if(K_REGIONS < 6) {
        throw std::runtime_error("Must have at least 6 regions for cube initialization.");
    }
    
    HalfedgeMesh* mesh = geometry->getMesh();

    auto whichHalf = [&](double val) { 
        if(val < 0) { 
            return 0;
        }
        return 1;
    };


    // "exact" version
    // Assign labels based on the coordinate with maximal absolute value,
    // which ends up looking like the faces on a cube.
    VertexData<LabelVec> cubeLabels(mesh, LabelVec::Zero());
    for(VertexPtr v : mesh->vertices()) {

        double maxCoord = -std::numeric_limits<double>::infinity();
        for(int i = 0; i < 3; i++) {
            if(std::abs(geometry->position(v)[i]) > maxCoord) {
                maxCoord = std::abs(geometry->position(v)[i]);


                int minInd = 2*i;
                if(geometry->position(v)[i] < 0) {
                    minInd++;
                }
               
                // Distance from one set of boundaries is determined by
                // gap between face-determining coordinate and second-largest
                // coordinate
                double faceCoord = std::abs(geometry->position(v)[i]);
                double faceDist = std::numeric_limits<double>::infinity();
                for(int j = 0; j < 3; j++) {
                    if(j != i) {
                        faceDist = std::min(faceDist, std::abs(faceCoord - std::abs(geometry->position(v)[j])));
                    }
                }

                double minDist = faceDist;
                cubeLabels[v] = LabelVec::Ones() * minDist;
                cubeLabels[v][minInd] = -minDist;
            }
        }

    }

    return cubeLabels;
    // return msdfFromLabels(geometry, cubeRegionLabels);
}
*/

VertexData<LabelVec> cubeCircleMSDF(Geometry<Euclidean>* geometry) {

    if(K_REGIONS < 9) {
        throw std::runtime_error("Must have at least 6 regions for cube initialization.");
    }
    
    HalfedgeMesh* mesh = geometry->getMesh();

    auto whichHalf = [&](double val) { 
        if(val < 0) { 
            return 0;
        }
        return 1;
    };


    // "exact" version
    // Assign labels based on the coordinate with maximal absolute value,
    // which ends up looking like the faces on a cube.
    VertexData<LabelVec> cubeLabels(mesh, LabelVec::Zero());
    for(VertexPtr v : mesh->vertices()) {

        double maxCoord = -std::numeric_limits<double>::infinity();
        for(int i = 0; i < 3; i++) {
            if(std::abs(geometry->position(v)[i]) > maxCoord) {
                maxCoord = std::abs(geometry->position(v)[i]);

                Vector2 otherCoords{geometry->position(v)[(i+1)%3],geometry->position(v)[(i+2)%3]};

                int minInd;

                double rad = norm(otherCoords);
                // double cRad = 0.6;
                // bool inCircle = rad < cRad;
                double cRad = 0.0;
                bool inCircle = false;

                if(inCircle) {
                    minInd = 3*i;
                } else {
                    if(otherCoords.x < 0) {
                        minInd = 3*i + 1;
                    } else {
                        minInd = 3*i + 2;
                    }
                }

               
                // Distance from one set of boundaries is determined by
                // gap between face-determining coordinate and second-largest
                // coordinate
                double faceCoord = std::abs(geometry->position(v)[i]);
                double faceDist = std::numeric_limits<double>::infinity();
                for(int j = 0; j < 3; j++) {
                    if(j != i) {
                        faceDist = std::min(faceDist, std::abs(faceCoord - std::abs(geometry->position(v)[j])));
                    }
                }

                double cDist = std::abs(cRad - rad);

                double splitDist = 1000;
                if(!inCircle) {
                    splitDist = std::abs(otherCoords.x);
                }

                double minDist = std::min(std::min(faceDist, cDist), splitDist);
                cubeLabels[v] = LabelVec::Ones() * minDist;
                cubeLabels[v][minInd] = -minDist;
            }
        }

    }

    return cubeLabels;
    // return msdfFromLabels(geometry, cubeRegionLabels);
}


VertexData<LabelVec> volleyballMSDF(Geometry<Euclidean>* geometry) {
  
    if(K_REGIONS < 9) {
        throw std::runtime_error("Must have at least 9 regions for cube initialization.");
    }
    
    HalfedgeMesh* mesh = geometry->getMesh();

    // Use to split in to strips like a volleyball
    double cutoff = 1.0/4.5;
    auto whichThird = [&](double val) { 
        if(val < -cutoff) {
            return 0;
        } else if(val < cutoff) {
            return 1;
        } else {
            return 2;
        }
    };

    /*
    // "discrete" version
    // Assign labels based on the coordinate with maximal absolute value,
    // which ends up looking like the faces on a cube.
    VertexData<int> cubeRegionLabels(mesh);
    for(VertexPtr v : mesh->vertices()) {

        double maxCoord = -std::numeric_limits<double>::infinity();
        for(int i = 0; i < 3; i++) {
            if(std::abs(geometry->position(v)[i]) > maxCoord) {
                maxCoord = std::abs(geometry->position(v)[i]);
                cubeRegionLabels[v] = 2*i*3 + whichThird(geometry->position(v)[(i+1)%3]);
                // cubeRegionLabels[v] = i;
            }
        }

        // FIXME
        // cubeRegionLabels[v] = whichThird(geometry->position(v).y);
    }

    return msdfFromLabels(geometry, colorWithLabels(mesh, cubeRegionLabels));
    // return msdfFromLabels(geometry, cubeRegionLabels);
    */

    // "exact" version
    // Assign labels based on the coordinate with maximal absolute value,
    // which ends up looking like the faces on a cube.
    VertexData<LabelVec> cubeLabels(mesh, LabelVec::Zero());
    for(VertexPtr v : mesh->vertices()) {

        double maxCoord = -std::numeric_limits<double>::infinity();
        for(int i = 0; i < 3; i++) {
            if(std::abs(geometry->position(v)[i]) > maxCoord) {
                maxCoord = std::abs(geometry->position(v)[i]);


                int minInd = i*3 + whichThird(geometry->position(v)[(i+1)%3]);
               
                // Distance from one set of boundaries is determined by
                // gap between face-determining coordinate and second-largest
                // coordinate
                double faceCoord = std::abs(geometry->position(v)[i]);
                double faceDist = std::numeric_limits<double>::infinity();
                for(int j = 0; j < 3; j++) {
                    if(j != i) {
                        faceDist = std::min(faceDist, std::abs(faceCoord - std::abs(geometry->position(v)[j])));
                    }
                }

                // Distance to other boundaries is determined by distance to ladder cutoff line
                double ladderCoord = std::abs(geometry->position(v)[(i+1)%3]);
                double ladderDist = std::abs(ladderCoord- cutoff);

                double minDist = std::min(faceDist, ladderDist);
                cubeLabels[v] = LabelVec::Ones() * minDist;
                cubeLabels[v][minInd] = -minDist;
            }
        }

    }

    return cubeLabels;
    // return msdfFromLabels(geometry, cubeRegionLabels);
}


VertexData<LabelVec> xAxisSplitMSDF(Geometry<Euclidean>* geometry) {

    HalfedgeMesh* mesh = geometry->getMesh();
    VertexData<LabelVec> sdf(mesh);
    double EPS = geometry->lengthScale() * 1e-3;

    // Compute shape center (surface area weighted)
    Vector3 center{0.0, 0.0, 0.0};
    double cWeight = 0;
    for(VertexPtr v : mesh->vertices()) {

        Vector3 p = geometry->position(v);
        double area = geometry->area(v.dual());

        center += p * area;
        cWeight += area;
    }
    center /= cWeight;

    for(VertexPtr v : mesh->vertices()) {

        Vector3 p = geometry->position(v) - center;
        double val = p.x + EPS;
        // double val = p.x + EPS - 0.7;

        LabelVec s = LabelVec::Zero();
        s[0] = val;
        s[1] = -val;

        sdf[v] = s;
    }

    return sdf;
}


VertexData<LabelVec> normalClusterMSDF(Geometry<Euclidean>* geometry, double angleThreshold) {

    HalfedgeMesh* mesh = geometry->getMesh();

    VertexData<Vector3> normal;
    geometry->getVertexNormals(normal);

    // === Sample clusters until all vertices have been assigned
    VertexData<int> clusterLabel(mesh, -1);
    size_t iUnusedVert = 0;
    int iCluster = 0;
    while(true) {

        // Find any vertex which has not been labeled
        // TODO Use most distant?
        VertexPtr seedVert;
        for(; iUnusedVert < mesh->nVertices(); iUnusedVert++) {
            if(clusterLabel[mesh->vertex(iUnusedVert)] == -1) {
                seedVert = mesh->vertex(iUnusedVert);
                break;
            }
        }
        if(seedVert == VertexPtr()) {
            break; // we have assigned all to a cluster!
        }

        // Grow the cluster as much as possible
        Vector3 clusterNormal = normal[seedVert];
        std::vector<VertexPtr> vertsToCheck = {seedVert};
        while(!vertsToCheck.empty()) {

            VertexPtr currVert = vertsToCheck.back();
            vertsToCheck.pop_back();

            // If we've already assigned this vertex since adding it to the queue, skip
            if(clusterLabel[currVert] != -1) continue;

            // Test if this vertex is sufficiently close to the cluster
            double relAngle = angle(clusterNormal, normal[currVert]);
            if(relAngle > angleThreshold) continue;

            // Accept the vertex
            clusterLabel[currVert] = iCluster; 

            // Consider neighbors
            for(VertexPtr vNeigh : currVert.adjacentVertices()) {
                if(clusterLabel[vNeigh] == -1) {
                    vertsToCheck.push_back(vNeigh);
                }
            }
        }


        iCluster++;
    }
    int nClusters = iCluster;

    // Color to remap to proper number of labels
    VertexData<int> regionLabels = colorWithLabels(mesh, clusterLabel);

    return msdfFromLabels(geometry, regionLabels);
}

VertexData<LabelVec> distanceClusterMSDF(Geometry<Euclidean>* geometry, int nCluster) {

    HalfedgeMesh* mesh = geometry->getMesh();


    // Generate cluster center
    std::vector<VertexPtr> clusterCenters;
    {
        std::vector<std::pair<VertexPtr, double>> clusterCentersBCond;
        VertexData<double> distances(mesh, std::numeric_limits<double>::infinity());
        for(int iC = 0; iC < nCluster; iC++) {
            
            // Find most distanct point
            double maxDist = 0;
            VertexPtr maxDistV;
            for(VertexPtr v : mesh->vertices()) {
                if(distances[v] > maxDist) {
                    maxDist = distances[v];
                    maxDistV = v;
                }
            }
    
            // Set it as the new cluster center
            clusterCenters.push_back(maxDistV);
            clusterCentersBCond.push_back(std::make_pair(maxDistV, 0.0));
    
            // Compute new distances
            distances = GC::FMMDistance(geometry, clusterCentersBCond);
        }
    }

    // Assign each vertex to closest cluster
    VertexData<double> minDistances(mesh, std::numeric_limits<double>::infinity());
    VertexData<int> clusterLabel(mesh);
    for(int iC = 0; iC < nCluster; iC++) {

        // Find distance to this cluster
        VertexData<double> distances = GC::FMMDistance(geometry, {std::make_pair(clusterCenters[iC], 0.0)});

        // Note minimum
        for(VertexPtr v : mesh->vertices()) {
            if(distances[v] < minDistances[v]) {
                minDistances[v] = distances[v];
                clusterLabel[v] = iC;
            }
        }
        

    }

    // Color to remap to proper number of labels
    VertexData<int> regionLabels = colorWithLabels(mesh, clusterLabel);

    return msdfFromLabels(geometry, regionLabels);
}

VertexData<int> colorWithLabels(HalfedgeMesh* mesh, VertexData<int> initLabels) {

    // === Attempt to color clusters with region numbers
    int nLabels = 0;
    for(VertexPtr v : mesh->vertices()) {
        nLabels = std::max(nLabels, initLabels[v]);
    }
    nLabels++;

    // Detect adjacent clusters
    std::vector<std::unordered_set<int>> clusterNeighbors(nLabels);
    for(VertexPtr v : mesh->vertices()) {

        int myCluster = initLabels[v];

        for(VertexPtr vNeigh : v.adjacentVertices()) {

            int neighCluster = initLabels[vNeigh];
            if(myCluster == neighCluster) continue;

            // Add to neighbors
            clusterNeighbors[myCluster].insert(neighCluster);
        }
    }

    // Attempt a very lazy greedy coloring

    // Sort vertices by decreasing neighbor degree
    std::vector<int> clusterOrder(nLabels);
    for(int i = 0; i < nLabels; i++) clusterOrder[i] = i;
    std::sort(clusterOrder.begin(), clusterOrder.end(),
        [&]( const int &a, const int &b ) { return clusterNeighbors[a].size() > clusterNeighbors[b].size(); } );

    // Color greedily
    bool coloringSucceeded = true;
    std::vector<int> clustersToRegions(nLabels, -1);
    int regionToTry = 0;  // cycle this to encourage diversity even when coloring is easy
    for(int iC : clusterOrder) {
            
        // Try candidate reigon labels to find the one with minimal conflicts
        int currBestRegion = -1;
        int nBestConflicts = 100000;
        // Cycle through possible region labels, starting with 'regionToTry'
        for(int iLabelCand = regionToTry; iLabelCand != (regionToTry+K_REGIONS-1)%K_REGIONS; iLabelCand = (iLabelCand+1)%K_REGIONS) {

            // Count conflicts
            int nConflict = 0;
            for(int neighCluster : clusterNeighbors[iC]) {
                if(clustersToRegions[neighCluster] == iLabelCand) {
                    nConflict++;
                }
            }

            // Check if this is the new best
            if(nConflict < nBestConflicts) {
                nBestConflicts = nConflict;
                currBestRegion = iLabelCand;
            }
        }

        // If we managed to use the target region, increment it to next (to encourage diversity)
        if(currBestRegion == regionToTry) {
            regionToTry = (regionToTry+1) % K_REGIONS;
        }

        // If we did not find a zero-conflict labelling, then this will unintentionally merge regions
        if(nBestConflicts != 0) {
            coloringSucceeded = false;
        }

        // Save the labelling
        clustersToRegions[iC] = currBestRegion;
    }

    // Remap clusters labels to region labels
    VertexData<int> regionLabels(mesh);
    for(VertexPtr v : mesh->vertices()) {
        regionLabels[v] = clustersToRegions[initLabels[v]];
    }

    // Warn if coloring failed
    if(!coloringSucceeded) {
        cout << "Warning: Greedy coloring failed and unintentionally merged clusters." << endl;
    }


    return regionLabels;
}


VertexData<LabelVec> dChartsMSDF(Geometry<Euclidean>* geometry) {

    HalfedgeMesh* mesh = geometry->getMesh();

    // Run simplified version of d-charts
    DCharts dCharts(geometry);
    dCharts.generateCharts();
    FaceData<int> faceLabels = dCharts.chartLabels;
    int nCharts = dCharts.numCharts;

    // Label vertices from faces
    VertexData<int> vertexLabels(mesh, -1);
    for(VertexPtr v : mesh->vertices()) {

        std::vector<int> labelCounts(nCharts, 0);
        for(FacePtr f : v.adjacentFaces()) {
            labelCounts[faceLabels[f]]++;
        }
        int maxCount = 0;
        for(int iC = 0; iC < nCharts; iC++) {
            if(labelCounts[iC] > maxCount) {
                maxCount = labelCounts[iC];
                vertexLabels[v] = iC;
            }
        }
    }

    // Color the vertex label
    VertexData<int> coloredLabels = colorWithLabels(mesh, vertexLabels);

    return msdfFromLabels(geometry, coloredLabels);
}

VertexData<LabelVec> msdfFromLabels(Geometry<Euclidean>* geometry, VertexData<int> labels) {

    HalfedgeMesh* mesh = geometry->getMesh();

    VertexData<LabelVec> sdf(mesh);

    for(VertexPtr v : mesh->vertices()) {
        
        // Initialize with large values
        int myLabel = labels[v];
        LabelVec f = LabelVec::Ones() * std::numeric_limits<double>::infinity();
        f[myLabel] *= -1;

        // Set distances in the MSDF assuming that crossings occur at the midpoint
        // between vertices of different labels
        for(HalfedgePtr he : v.incomingHalfedges()) {

            double len = geometry->length(he.edge());
            int otherLabel = labels[he.vertex()];

            if(otherLabel == myLabel) continue;

            f[myLabel] = std::max(f[myLabel], -len/2.0);
            f[otherLabel] = std::min(f[otherLabel], len/2.0);

        }

        // If any value didn't get set, set it arbitrarily to 1. Such values correspond to
        // interior vertices, and their values will be overwritten as soon as we do an FMM
        // projection anyway.
        for(int i = 0; i < K_REGIONS; i++) {
            if(!std::isfinite(f[i])) {
                if(i == myLabel) {
                    f[i] = -1.0;
                } else {
                    f[i] = 1.0;
                }
            }
        }

        sdf[v] = f;
    }

    return sdf;
}

VertexData<LabelVec> meanCurvatureLevelSetMSDF(Geometry<Euclidean>* geometry, double scale) {

    HalfedgeMesh* mesh = geometry->getMesh();

    // Compute mean curvature at each vertex
    EdgeData<double> cotanWeights;
    geometry->getEdgeCotanWeights(cotanWeights);
    DualFaceData<double> dualAreas;
    geometry->getDualFaceAreas(dualAreas);
    VertexData<double> meanCurvature(mesh);
    for(VertexPtr v : mesh->vertices()) {

        Vector3 H{0.0, 0.0, 0.0};

        for(int i = 0; i < 3; i++) {
            for(HalfedgePtr he : v.incomingHalfedges()) {
                H[i] += cotanWeights[he.edge()] * (geometry->position(he.vertex())[i] - geometry->position(v)[i]);
            }
        } 

        double hMag = norm(H);
        meanCurvature[v] = hMag / dualAreas[v.dual()];
    }


    // Find level sets
    VertexData<int> baseLabels(mesh);
    for(VertexPtr v : mesh->vertices()) {
        int label = (int) std::round(meanCurvature[v]*scale);

        int octant = 0;
        Vector3 p = geometry->position(v);
        if(p.x > 0) octant += 4;
        if(p.y > 0) octant += 2;
        if(p.z > 0) octant += 1;

        baseLabels[v] = 8*label + octant;
    }

    // Re-label for connectect components
    VertexData<int> labels(mesh, -1);
    int currLab = 0;
    for(VertexPtr v : mesh->vertices()) {
        if(labels[v] != -1) continue;

        // Greedily label neighbors
        std::vector<VertexPtr> toProc{v};
        while(!toProc.empty()) {
            VertexPtr currV = toProc.back();
            toProc.pop_back();
            if(labels[currV] != -1) continue;
            labels[currV] = currLab;

            for(VertexPtr neighV : currV.adjacentVertices()) {
                if(baseLabels[neighV] == baseLabels[currV]) {
                    toProc.push_back(neighV);
                }
            }
        }


        currLab++;
    }


    return msdfFromLabels(geometry, colorWithLabels(mesh, labels));
}



VertexData<LabelVec> transferMSDF(Geometry<Euclidean>* sourceGeometry, VertexData<LabelVec>& sourceMSDF, Geometry<Euclidean>* targetGeometry)  {
  // TODO implement
  return sourceMSDF;
}
