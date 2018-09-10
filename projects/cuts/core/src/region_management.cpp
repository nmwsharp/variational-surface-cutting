#include "region_management.h"

#include "disjoint_sets.h"
#include "surface_line.h"
#include "polygon_soup_mesh.h"
#include "target_edge_lengths.h"
#include "developable_approximation.h"
#include "fast_cholesky.h"
#include "bff.h"
#include "boundary_constraints.h"
#include "fast_marching_method.h"
#include "mesh_union_split.h"
#include "combining_hash_functions.h"
#include "fast_triangle_soup.h"

// Rectangle bin packing
#include "Rect.h"
#include "ShelfBinPack.h"
#include "GuillotineBinPack.h"

#include <tuple>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

double cutRegionToDisk(Geometry<Euclidean>* regionGeom, SymmetryResult* symmetry) {

    HalfedgeMesh* mesh = regionGeom->getMesh();

    // Unmark all cut edges, if some happened to already be marked
    for(EdgePtr e : mesh->edges()) {
        e.markCut(false);
    }


    // If the region is disk-like, there is nothing to be done
    // TODO verify that this check is correct
    int eulerChar = mesh->nVertices() - mesh->nEdges() + mesh->nFaces();
    if(eulerChar == 1) {
        return 0.0;
    }


    // == Find the minimum spanning tree of edges that connects all cuts
    size_t nVert = mesh->nVertices();
    VertexData<size_t> vInd = mesh->getVertexIndices();
    MarkedDisjointSets vertexSets(nVert);
    
    // Datastructures for an expanding search outward
    VertexData<HalfedgePtr> parent(mesh, HalfedgePtr());
    typedef std::tuple<double, VertexPtr, HalfedgePtr> VEntry;
    std::priority_queue<VEntry, vector<VEntry>, std::greater<VEntry>> PQ;

    // Merge each loop to initialize
    for(FacePtr bl : mesh->boundaryLoops()) {
        HalfedgePtr prevHe = HalfedgePtr();
        for(HalfedgePtr he : bl.adjacentHalfedges()) {
            if(prevHe != HalfedgePtr()) {
                vertexSets.merge(vInd[he.vertex()], vInd[prevHe.vertex()]);
            }
            prevHe = he;

            // Add initial neighbors to the search
            VertexPtr v = he.vertex();
            for(HalfedgePtr vHe : v.outgoingHalfedges()) {
                if(vHe.edge().isBoundary()) continue;

                VertexPtr otherV = vHe.twin().vertex();
                PQ.push(std::make_tuple(regionGeom->length(vHe.edge()), otherV, vHe));
            }
        }
        
        // Mark the set to indicate it is a boundary set
        vertexSets.mark(vInd[prevHe.vertex()]);
    }


    // Search to find the shortest connection
    int nBoundariesToUnion = mesh->nBoundaryLoops()-1;
    // while(!PQ.empty() && nBoundariesToUnion > 0) {
    while(!PQ.empty()) {
        
        VEntry e = PQ.top();
        PQ.pop();

        double currDist = std::get<0>(e);
        VertexPtr currV = std::get<1>(e);
        HalfedgePtr prevHe = std::get<2>(e);

        // If we've already unioned these two vertices, skip it
        VertexPtr oldV = prevHe.vertex();
        if(vertexSets.find(vInd[currV]) == vertexSets.find(vInd[oldV])) {
            continue;
        }
        // If unioning across this edge would join to marked regions, then we are connecting
        // holes. Note this and mark the chain of edges that led us here as cut edges.
        if(vertexSets.isMarked(vInd[currV]) && vertexSets.isMarked(vInd[oldV])) {
            nBoundariesToUnion--; 

            // One side
            HalfedgePtr currHe = parent[currV];
            while(currHe != HalfedgePtr()) {
                currHe.edge().markCut(true);
                currHe = parent[currHe.vertex()];
            }

            // Bridge
            prevHe.edge().markCut(true);

            // Other side 
            currHe = parent[oldV];
            while(currHe != HalfedgePtr()) {
                currHe.edge().markCut(true);
                currHe = parent[currHe.vertex()];
            }

        }

        // Accept this edge and note the parent (this may be overwriting a previous setting, but we want that)
        parent[currV] = prevHe;


        // Accept this vertex, unioning the sets and adding neighbors for consideration
        vertexSets.merge(vInd[currV], vInd[oldV]);
        for(HalfedgePtr outHe : currV.outgoingHalfedges()) {
            PQ.push(std::make_tuple(regionGeom->length(outHe.edge()) + currDist, outHe.twin().vertex(), outHe));
        }

    }
    
    // FIXME Need to cut along homology basis to handle genus


    // If a symmetry specification is given, mark all symmetric edges
    // Note that his may mean we cut the surface in to more than one disk
    if(symmetry != nullptr && symmetry->symmetryFound) {

        // Redo symmetry detection on this patch (rather than trying to align to the old symmetry result)
        SymmetryResult patchSym = detectSymmetryAuto(regionGeom);
        if(patchSym.symmetryFound) {

            // Index symmetry pairs for easy symmetry tests
            size_t iSym = 0;
            VertexData<size_t> symInd(mesh);
            for(VertexPtr v : patchSym.canonicalVertices) {
                symInd[v] = iSym;
                for(VertexPtr vs : patchSym.symmetrySet[v]) {
                    symInd[vs] = iSym;
                }
                iSym++;
            }
    
            // Build symmetric edge sets
            std::vector<std::vector<EdgePtr>> edgeSets; // warning: because of the way we build this from halfedges, each set actually appears twice... oh well
            for(VertexPtr v : patchSym.canonicalVertices) {
                for(HalfedgePtr he : v.incomingHalfedges()) {
    
                    std::vector<EdgePtr> edgeSet = {he.edge()};
    
                    // Check all possible symmetric edges to build the set
                    for(VertexPtr vs : patchSym.symmetrySet[v]) {
                        for(HalfedgePtr hes : vs.incomingHalfedges()) {
                            if(symInd[he.vertex()] == symInd[hes.vertex()]) {
                                edgeSet.push_back(hes.edge());
                            }
                        }
                    }

                    edgeSets.push_back(edgeSet);
                }
            }
    
            // Copy cuts across edge sets
            for(std::vector<EdgePtr> edgeSet : edgeSets) {
                bool anyCut = false;
                for(EdgePtr e : edgeSet) {
                    if(e.isCut()) {
                        anyCut = true;
                    }
                }
    
                if(anyCut) {
                    for(EdgePtr e : edgeSet) {
                        e.markCut(true);
                    }
                }
            }

        } else {
            cout << "Whole mesh has symmetry, but could not find symmetry for patch. Skipping." << endl;
        }

    }


    // Add up the length of cut edges
    double cutEdgeLength = 0;
    size_t nCutEdges = 0;
    for(EdgePtr e : mesh->edges()) {
        if(e.isCut()) {
            cutEdgeLength += regionGeom->length(e);
            nCutEdges++;
        }
    }
    
    cout << "Added " << nCutEdges << " cut edges with length " << cutEdgeLength << " to region cut to disk." << endl;
    return cutEdgeLength;
}

// Assumes cuts have already been set, and that cuts never terminate inside of the surface
std::vector<std::vector<HalfedgePtr>> getDirectedCutEdges(HalfedgeMesh* mesh) {

    // === Find a directed graph of cuts

    // Any vertex that has > 2 cut edges touching it, or has any cut edges and is on the boundary
    VertexData<char> isCutVertex(mesh, (char)false);
    int nCutVert = 0;
    for(VertexPtr v : mesh->vertices()) {

        int cutCount = 0;
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) {
                cutCount++;
            }
        }

        if(cutCount > 2 || (cutCount > 0 && v.isBoundary())) {
            isCutVertex[v] = true;
            nCutVert++;
        }

    }


    // Find directed paths of edges connecting the cut vertices
    std::vector<std::vector<HalfedgePtr>> directedCutEdges;
    EdgeData<char> edgeUsed(mesh, (char)false);
    for(VertexPtr v : mesh->vertices()) {
        if(!isCutVertex[v]) continue;

        // Look for the start of any edges
        for(HalfedgePtr startHe : v.outgoingHalfedges()) {
            if(startHe.edge().isCut() && !edgeUsed[startHe.edge()]) {

                // Trace out a halfedge path
                HalfedgePtr currHe = startHe;
                std::vector<HalfedgePtr> thisCutEdge;
                while(true) {

                    thisCutEdge.push_back(currHe);
                    edgeUsed[currHe.edge()] = true;

                    VertexPtr nextVertex = currHe.twin().vertex();
                    if(isCutVertex[nextVertex]) break;

                    
                    for(HalfedgePtr nextHe : nextVertex.outgoingHalfedges()) {
                        if(nextHe.edge().isCut() && !edgeUsed[nextHe.edge()]) {
                            currHe = nextHe;
                            break;
                        }
                    }
                }

                directedCutEdges.push_back(thisCutEdge);
            }
        }
    }

    if(directedCutEdges.size() > 0) {
        cout << "Found " << directedCutEdges.size() << " directed cut edges." << endl;
    }
    return directedCutEdges;

}


std::vector<std::array<Vector3, 2>> getGeodesicCuts(Geometry<Euclidean>* geom) {

    HalfedgeMesh* mesh = geom->getMesh();

    std::vector<std::vector<HalfedgePtr>> directedCutEdges = getDirectedCutEdges(mesh);
    
    // Get the relaxed lines for each of the cuts
    std::vector<std::array<Vector3, 2>> allLines;

    SurfaceLine surfaceLine(geom);
    for(std::vector<HalfedgePtr>& cutEdge : directedCutEdges) {

        surfaceLine.clearLine();
        surfaceLine.constructFromHalfedgePath(cutEdge);
        surfaceLine.geodesicRelax();

        std::vector<std::array<Vector3, 2>> thisLines = surfaceLine.getLineSegments();

        allLines.insert(allLines.end(), thisLines.begin(), thisLines.end());
    }

    return allLines;
}

namespace {
    struct array_hash {
    inline std::size_t operator()(const std::array<size_t,2>& arr) const {
        return std::hash<size_t>{}(std::get<0>(arr)) ^ (std::hash<size_t>{}(std::get<1>(arr)) << 2);
    }
    };
}


Geometry<Euclidean>* remeshAlongCutGeodesics(Geometry<Euclidean>* geom) {

    HalfedgeMesh* mesh = geom->getMesh();

    // === Relax all of the edges and get the segments where they cross each face.
    // Note that this code assumes that no two segments will ever cross one another inside of a face

    std::vector<std::vector<HalfedgePtr>> directedCutEdges = getDirectedCutEdges(mesh);

    size_t INVALID_IND = std::numeric_limits<size_t>::max();
    struct Crossing {

        Crossing(HalfedgePtr he_, double t_, Crossing* otherEnd_) :
            he(he_), t(t_), otherEnd(otherEnd_)
        {
            index = std::numeric_limits<size_t>::max();
            firstSeenInd = std::numeric_limits<size_t>::max();
        }

        HalfedgePtr he;
        double t;
        Crossing* otherEnd;
        size_t index;
        size_t firstSeenInd;
    };

    HalfedgeData<std::vector<Crossing*>> halfedgeCrossings(mesh);
    SurfaceLine surfaceLine(geom);
    for(std::vector<HalfedgePtr>& cutEdge : directedCutEdges) {

        surfaceLine.clearLine();
        surfaceLine.constructFromHalfedgePath(cutEdge);
        surfaceLine.geodesicRelax();

        std::vector<std::tuple<HalfedgePtr, double, HalfedgePtr, double>> crossings = surfaceLine.getSegmentsInFaces();

        // Special case for degenerate edge-aligned geodesic -- no cuts are needed
        if(crossings.size() == 1) {
            continue;
        }

        for(auto c : crossings) {

            HalfedgePtr h1 = std::get<0>(c);
            double t1 = std::get<1>(c);
            HalfedgePtr h2 = std::get<2>(c);
            double t2 = std::get<3>(c);
            Crossing* c1 = new Crossing(h1, t1, nullptr);
            Crossing* c2 = new Crossing(h2, t2, nullptr);
            c1->otherEnd = c2;
            c2->otherEnd = c1;
            halfedgeCrossings[h1].push_back(c1);
            halfedgeCrossings[h2].push_back(c2);

            // Some sanity checks
            if(h1 == h2) {
                throw std::runtime_error("Halfedge crossings are on same halfedge");
            }
            if(h1.face() != h2.face()) {
                throw std::runtime_error("Halfedge crossings are not in shared face");
            }
        }
    }

    // === Construct a new mesh by triangulating the polygons induced by the relaxed lines

    // Sort the crossings along each halfedge
    for(HalfedgePtr he : mesh->halfedges()) {
        std::sort(halfedgeCrossings[he].begin(), halfedgeCrossings[he].end(), [](const Crossing* lhs, const Crossing* rhs)
        {
            return lhs->t < rhs->t;
        });
    }

    // Build triangles out of the crossings
    vector<Vector3> positions;
    vector<vector<size_t>> faces;

    // Create vertices

    // Original vertices
    size_t nVert = mesh->nVertices();
    VertexData<size_t> vInd = mesh->getVertexIndices();
    for(VertexPtr v: mesh->vertices()) {
        positions.push_back(geom->position(v));
    }

    // New vertices induced by line
    for(EdgePtr e : mesh->edges()) {

        HalfedgePtr he = e.halfedge();
        size_t nCross = halfedgeCrossings[he].size();
        for(size_t i = 0; i < nCross; i++) {

            // This crossing terminates at the vertex
            if(halfedgeCrossings[he][i]->t == 0.0) {
                halfedgeCrossings[he][i]->index = vInd[he.vertex()];
            } 
            // Typical crossing (we know there will be a complement on the other side of the edge)
            else {
                size_t thisInd = nVert++;
                halfedgeCrossings[he][i]->index = thisInd;

                // Find and mark the twin
                // (can't directly index in to twin because terminating halfedges to not exist symmetrically)
                for(size_t j = 0; j < halfedgeCrossings[he.twin()].size(); j++) {
                    if(approxEqualsAbsolute(halfedgeCrossings[he.twin()][j]->t, 1.0 - halfedgeCrossings[he][i]->t, 1e-7)) { // should be exact?
                        halfedgeCrossings[he.twin()][j]->index = thisInd;
                    }
                }
                positions.push_back(geom->position(he.vertex()) + geom->vector(he)*halfedgeCrossings[he][i]->t);
            }
        }

        // Check for terminating crossings on the opposite halfedge
        he = e.halfedge().twin();
        nCross = halfedgeCrossings[he].size();
        for(size_t i = 0; i < nCross; i++) {
            if(halfedgeCrossings[he][i]->t == 0.0) {
                halfedgeCrossings[he][i]->index = vInd[he.vertex()];
            } 
        }
    }

    // Safety check
    for(HalfedgePtr he : mesh->halfedges()) {
        std::vector<Crossing*> crossings = halfedgeCrossings[he];
        std::vector<Crossing*> twinCrossings = halfedgeCrossings[he.twin()];
        for(Crossing* c : crossings) {
            if(c->index == INVALID_IND) {
                throw std::runtime_error("Unnumbered crossing escaped");
            }
        }
    }

    // Create faces
    std::unordered_set<std::array<size_t,2>, array_hash> cutEdges; // track which edges correspond to cuts
    for(FacePtr f : mesh->faces()) {

        // Walk around the face, pushing points as we see them and popping once we
        // see an edge the second time, until we reach the vertex where we saw it the first time
        std::vector<size_t> pointStack; 

        for(HalfedgePtr he : f.adjacentHalfedges()) {

            // Push a point for this corner
            pointStack.push_back(vInd[he.vertex()]);

            for(Crossing* c : halfedgeCrossings[he]) {
                
                // Push the vertex for the crossing
                c->otherEnd->firstSeenInd = c->index;
                // Don't need to push a second point for a cut edge sharing a vertex
                if(c->t != 0.0) {
                    pointStack.push_back(c->index);
                }
                

                // Pop off a polygon if this is the second time we've seen the edge
                if(c->firstSeenInd != INVALID_IND) {

                    std::vector<size_t> newPoly;

                    // Pop points to create a polygon
                    while(pointStack.back() != c->firstSeenInd) {
                        size_t p = pointStack.back();
                        newPoly.push_back(p);
                        pointStack.pop_back();
                    }
                    size_t p = pointStack.back();
                    newPoly.push_back(p);
                    pointStack.pop_back();
                    std::reverse(std::begin(newPoly), std::end(newPoly)); // popping winds it the wrong way
                    faces.push_back(newPoly);
                }

                // Push second point
                pointStack.push_back(c->index);
            }

        }

        // Pop the final polygon
        faces.push_back(pointStack);

    }

    // Build the mesh
    PolygonSoupMesh pSoup(faces, positions);
    pSoup.triangulate();

    Geometry<Euclidean>* newGeom;
    HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, newGeom);

    // Mark the vertices along the cut path
    VertexData<char> isCutVertex(newMesh, (char)false);
    for(size_t i = mesh->nVertices(); i < newMesh->nVertices(); i++) {
        isCutVertex[newMesh->vertex(i)] = true; // all new vertices are along cuts
    }
    // The first and last vertex of each directed cut chain are cut vertices
    for(std::vector<HalfedgePtr> cutChain : directedCutEdges) {
        isCutVertex[newMesh->vertex(vInd[cutChain[0].vertex()])] = true;
        isCutVertex[newMesh->vertex(vInd[cutChain.back().twin().vertex()])] = true;
    }

    // Mark cut edges
    for(EdgePtr e : newMesh->edges()) {
        if(isCutVertex[e.halfedge().vertex()] && isCutVertex[e.halfedge().twin().vertex()]) {
            e.markCut(true);
        }
    }


    // Free all of the crossing objects we created
    for(HalfedgePtr he : mesh->halfedges()) {
        for(Crossing* c : halfedgeCrossings[he]) {
            delete c;
        }
    }

    return newGeom;
}



void canonicalizeParameterization(Geometry<Euclidean>* geom, CornerData<Vector2>& param) {

    HalfedgeMesh* mesh = geom->getMesh();

    // === Rescale the flattening such that it has unit scale along the boundary
    double realBoundaryLength = 0;
    double paramBoundaryLength = 0;
    for(HalfedgePtr he : mesh->halfedges()) {
        if(!he.twin().isReal()) {
            realBoundaryLength += norm(geom->vector(he));

            Vector2 p1 = param[he.next().corner()];
            Vector2 p2 = param[he.next().next().corner()];
            paramBoundaryLength += norm(p1 - p2);
        }
    }

    double ratio = realBoundaryLength / paramBoundaryLength;
    for(HalfedgePtr he : mesh->halfedges()) {
        param[he.corner()] *= ratio;
    }

    // === Center at the origin
    Vector2 meanCoord{0.0, 0.0};
    double totalArea = 0.0;
    for(HalfedgePtr he : mesh->halfedges()) {
        double area = geom->area(he.face()) / 3.0;
        totalArea += area;
        meanCoord += param[he.corner()]*area;
    }
    meanCoord /= totalArea;
    for(HalfedgePtr he : mesh->halfedges()) {
        param[he.corner()] -= meanCoord;
    }

    // === Rotate to align principal axis with x-axis
    double Ixx = 0.0;
    double Iyy = 0.0;
    double Ixy = 0.0;
    for(HalfedgePtr he : mesh->halfedges()) {
        double area = geom->area(he.face()) / 3.0;
        Vector2 p = param[he.corner()];
        Ixx += area * p.x * p.x;
        Iyy += area * p.y * p.y;
        Ixy += -area * p.x * p.y;
    }

    // Compute major axis
    double theta = -std::atan(2 * Ixy / (Ixx - Iyy)) / 2.0;
    double matXX = std::cos(theta);
    double matXY = -std::sin(theta);
    double matYX = std::sin(theta);
    double matYY = std::cos(theta);
    for(HalfedgePtr he : mesh->halfedges()) {
        Vector2 p = param[he.corner()];
        Vector2 newP = {matXX * p.x + matXY * p.y, matYX * p.x + matYY * p.y};
        param[he.corner()] = newP;
    }
    
}


// Helper function
namespace {

    bool paramContinuousAcrossHalfedge(const CornerData<Vector2>& param, HalfedgePtr he) {
        if(he.edge().isBoundary()) return false;
        CornerPtr thisWedgeA = he.next().corner();
        CornerPtr otherWedgeA = he.twin().next().next().corner();
        CornerPtr thisWedgeB = he.next().next().corner();
        CornerPtr otherWedgeB = he.twin().next().corner();

        // Check if the parameterization is continuous across the shared edge
        if(param[thisWedgeA] == param[otherWedgeA] && param[thisWedgeB] == param[otherWedgeB]) {
            return true;
        }
        return false;
    }
    
};

// Given a discontinuous texture such that the mesh has multiple connected components in texture space,
// translate the connected components such that none overlap.
void packParameterization(Geometry<Euclidean>* geom, CornerData<Vector2>& param, double marginFactor) { 
    HalfedgeMesh* mesh = geom->getMesh();

    // === Start by finding connected componnets
    FaceData<size_t> connComp(mesh);
    size_t nConnComp = 0;
    FaceData<size_t> fInd = mesh->getFaceIndices();

    {
        // Union to find connected components
        DisjointSets faceConnComp(mesh->nFaces());
        for(FacePtr f : mesh->faces()) {
            for(HalfedgePtr he : f.adjacentHalfedges()) {

                if(he.edge().isBoundary()) continue;

                CornerPtr thisWedgeA = he.next().corner();
                CornerPtr otherWedgeA = he.twin().next().next().corner();
                CornerPtr thisWedgeB = he.next().next().corner();
                CornerPtr otherWedgeB = he.twin().next().corner();
    
                // Check if the parameterization is continuous across the shared edge
                if(param[thisWedgeA] == param[otherWedgeA] && param[thisWedgeB] == param[otherWedgeB]) {
                    faceConnComp.merge(fInd[f], fInd[he.twin().face()]);
                }
            }
        }
    
        // Desnsely re-number the connected components
        std::vector<long long int> compInd(mesh->nFaces(), -1);
        for(FacePtr f : mesh->faces()) {
            size_t ind = faceConnComp.find(fInd[f]);
            if(compInd[ind] < 0) {
                compInd[ind] = nConnComp++;
            }
            
            connComp[f] = compInd[ind];
        }
    }


    // Find bounding boxes for each connected component in parameter space
    std::vector<Vector2> minBounds(nConnComp, Vector2{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()});
    std::vector<Vector2> maxBounds(nConnComp, Vector2{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()});
    
    for(FacePtr f : mesh->faces()) {
        size_t comp = connComp[f];

        for(CornerPtr c : f.adjacentCorners()) {
            Vector2 p = param[c];
            minBounds[comp] = componentwiseMin(minBounds[comp], p); 
            maxBounds[comp] = componentwiseMax(maxBounds[comp], p); 
        }
    }

    // Initial pass to compute area
    double initTotalArea = 0.0; // before margin expansion
    for(size_t i = 0; i < nConnComp; i++) {

        Vector2 initMin = minBounds[i];
        Vector2 initMax = maxBounds[i];
        initTotalArea += (initMax.x - initMin.x) * (initMax.y - initMin.y);
    }
    double bufferWidth = std::sqrt(initTotalArea) * marginFactor;

    // Expand the bounds to include a buffer margin
    // Compute some other useful values while we're at it
    Vector2 globalMin = Vector2{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    std::vector<Vector2> originalPatchCenters(nConnComp); // we will need these later
    double totalArea = 0.0; // after margin expansion
    for(size_t i = 0; i < nConnComp; i++) {

        Vector2 initMin = minBounds[i];
        Vector2 initMax = maxBounds[i];

        Vector2 c = 0.5 * (initMin + initMax);
        originalPatchCenters[i] = c;
        Vector2 extent = initMax - initMin;
        Vector2 newExtent = extent + Vector2{1.0, 1.0}*bufferWidth;

        Vector2 newMin = c - newExtent/2;
        Vector2 newMax = c + newExtent/2;


        minBounds[i] = newMin;
        maxBounds[i] = newMax;
        totalArea += (newMax.x - newMin.x) * (newMax.y - newMin.y);
        globalMin = componentwiseMin(globalMin, newMin);
    }

    // Quantitize the problem, as we conver the data to rectangles
    // since apparently all bin-packing libraries use integer coordinates
    double approxTotalSideLen = std::sqrt(totalArea);
    double unitsPerInt = approxTotalSideLen / 10000.0;
    globalMin = globalMin - Vector2{10*unitsPerInt, 10*unitsPerInt}; // numerical shift
    std::vector<rbp::Rect> rects;
    int sumDim = 0;
    for(size_t i = 0; i < nConnComp; i++) {

        Vector2 minF = minBounds[i];
        Vector2 maxF = maxBounds[i];

        // Shift by the global min, so all coordinates are positive and the lower left is near 0
        minF -= globalMin;
        maxF -= globalMin;

        // Round outwards, so the bounding box is always conservatively large
        int minX = static_cast<int>(std::floor(minF.x / unitsPerInt));
        int maxX = static_cast<int>(std::ceil(maxF.x / unitsPerInt));
        int minY = static_cast<int>(std::floor(minF.y / unitsPerInt));
        int maxY = static_cast<int>(std::ceil(maxF.y / unitsPerInt));

        int width = maxX - minX;
        int height = maxY - minY;
        rects.push_back(rbp::Rect{minX, minY, width, height});
        sumDim += std::max(width, height);
    }

    // Store the result here
    std::vector<Vector2> newCenters(nConnComp);    
    std::vector<bool> isFlipped(nConnComp);

    // Split in to two different algorithms, depending on how many pieces there area.
    if(nConnComp > 200) {

        // Run a simple and efficient shelf-based bin packing.
        
        // Double the size of the box until it succeeds
        int currSize = static_cast<int>(std::ceil(std::sqrt(totalArea)/unitsPerInt));
        bool packingFailed;
        do {
            packingFailed = false;
    
            newCenters.clear();
            isFlipped.clear();
    
            // Try with this sized box
            rbp::ShelfBinPack packer(currSize, currSize, false);
            for(size_t i = 0; i < nConnComp; i++) {
        
                rbp::Rect& rect = rects[i];

                // Insert
                rbp::Rect packRect = packer.Insert(rect.width, rect.height, rbp::ShelfBinPack::ShelfChoiceHeuristic::ShelfFirstFit);
    
                // Check for failure
                if(packRect.width == 0) {
                    packingFailed = true;
                    break;
                }
    
                // Check if it was flipped
                if(packRect.width == rect.width) {
                    isFlipped[i] = false;
                } else {
                    isFlipped[i] = true;
                }
        
                // Note the new center
                Vector2 newC{packRect.x + packRect.width/2.0, packRect.y + packRect.height/2.0};
                newC *= unitsPerInt; // undo quantitization scaling
                newCenters[i] = newC;
            }
        
            // Increase size for the next iteration (if needed)
            currSize = static_cast<int>(std::ceil(currSize * 1.5));
        } while (packingFailed);

    } else {

        // Run an expensive but superior Guillotine scheme
        
        // Increase the size of the box until it succeeds
        int currSize = static_cast<int>(std::ceil(std::sqrt(totalArea)/unitsPerInt));
        bool packingFailed;
        do {
            packingFailed = false;
    
            newCenters.clear();
            isFlipped.clear();
    
            // Try with this sized box
            rbp::GuillotineBinPack packer(currSize, currSize);
            for(size_t i = 0; i < nConnComp; i++) {
        
                rbp::Rect& rect = rects[i];
        
                // Insert
                rbp::Rect packRect = packer.Insert(rect.width, rect.height, true,
                                                   rbp::GuillotineBinPack::FreeRectChoiceHeuristic::RectBestAreaFit,
                                                   rbp::GuillotineBinPack::GuillotineSplitHeuristic::SplitMinimizeArea);
                
                // Check for failure
                if(packRect.width == 0) {
                    packingFailed = true;
                    break;
                }
    
                // Check if it was flipped
                if(packRect.width == rect.width) {
                    isFlipped[i] = false;
                } else {
                    isFlipped[i] = true;
                }
        
                // Note the new center
                Vector2 newC{packRect.x + packRect.width/2.0, packRect.y + packRect.height/2.0};
                newC *= unitsPerInt; // undo quantitization scaling
                newCenters[i] = newC;
            }
        
            // Increase size for the next iteration (if needed)
            currSize = static_cast<int>(std::ceil(currSize * 1.2));
        } while (packingFailed);

    }



    // Transform parameterization according to the results of the packing
    for(FacePtr f : mesh->faces()) {
        size_t iComp = connComp[f];

        bool flip = isFlipped[iComp];
        Vector2 origC = originalPatchCenters[iComp];
        Vector2 targetC = newCenters[iComp];

        for(CornerPtr c : f.adjacentCorners()) {
            Vector2 p = param[c];

            // Shift to the origin
            p -= origC;

            // Rotate, if needed
            if(flip) {
                p = Vector2{-p.y, p.x};
            }

            // Shift to new center
            p += targetC;

            param[c] = p;
        }
    }
}



Geometry<Euclidean>* makePatchExtrinsicDevelopable(Geometry<Euclidean>* geom) {

    HalfedgeMesh* mesh = geom->getMesh();
        
    // Vertices at cuts and boundaries are fixed
    VertexData<char> fixed(mesh, (char)false);
    for(VertexPtr v : mesh->vertices()) {
        // fixed[v] = (char)(v.isBoundary());
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) {
                fixed[v] = (char)true;
            }
        }
    }
    

    VertexData<Vector3> newPositions;
    for(int iter = 0; iter < 1; iter++) {
    
        VertexData<double> scaleFactors = yamabeScaleFactors(geom);
    
        // Read compute target edge lengths from scales
        /*
        EdgeData<double> targetL(mesh);
        for(EdgePtr e : mesh->edges()) {
    
            double initL = geom->length(e);
            double meanScale = 0.5 * (scaleFactors[e.halfedge().vertex()] + scaleFactors[e.halfedge().twin().vertex()]);
    
            targetL[e] = initL * meanScale;
            if(!std::isfinite(targetL[e])) {
                throw std::runtime_error("non finite target length");
            }
        }
        */
    
        // Flatten patch to find target edge lengths
        BoundaryFirstFlattening bff(geom);
        BoundaryConstraints constraints(&bff);
        constraints.setIsometricConstraints();
        bff.flatten();
        CornerData<Vector2> param = bff.uvs;
        canonicalizeParameterization(geom, param);
        
        EdgeData<double> targetL(mesh);
        double globalScale = 0.0;
        int globalScaleCount = 0;
        for(EdgePtr e : mesh->edges()) {
    
            HalfedgePtr he = e.halfedge();
    
            Vector2 c1 = param[he.next().corner()];
            Vector2 c2 = param[he.next().next().corner()];
            double l = norm(c1 - c2);
    
            targetL[e] = l;
        }
        
    
        // Iterative solve for positions satisfying these lenghts
        PositionSolver p(geom);
        newPositions = p.solve(targetL, fixed);

        for(VertexPtr v : mesh->vertices()) {
            (*geom)[v] = newPositions[v];
        }
    
    }
    
    // Return a new geometry/mesh
    HalfedgeMesh* newMesh = mesh->copy();
    Geometry<Euclidean>* devG = new Geometry<Euclidean>(*newMesh);
    for(size_t i = 0; i < mesh->nVertices(); i++) {
        (*devG)[newMesh->vertex(i)] = newPositions[mesh->vertex(i)];
    }


    // Check just how developable it is
    double maxDefect = -1;
    double meanDefectMag = 0;
    VertexData<double> defects;
    devG->getVertexAngleDefects(defects);
    for(size_t i = 0; i < mesh->nVertices(); i++) {
        VertexPtr oldV = mesh->vertex(i);
        VertexPtr newV = newMesh->vertex(i);
        if(!fixed[oldV]) {
            double defect = std::abs(defects[newV]);
            maxDefect = std::max(maxDefect, defect);
            meanDefectMag += defect;
        }
    }
    meanDefectMag /= mesh->nVertices();
    cout << "Extrinsic developable approximation: Max defect = " << maxDefect << "  mean defect magnitude = " << meanDefectMag << endl;
    

    return devG;
}


Geometry<Euclidean>* makeShapeExtrinsicDevelopable(Geometry<Euclidean>* geom, EdgeData<double> targetL, bool stretchOnly) {
    
        HalfedgeMesh* mesh = geom->getMesh();
            
        // Vertices at cuts and boundaries are cut vertices
        VertexData<char> cutVert(mesh, (char)false);
        for(VertexPtr v : mesh->vertices()) {
            // fixed[v] = (char)(v.isBoundary());
            for(EdgePtr e : v.adjacentEdges()) {
                if(e.isCut()) {
                    cutVert[v] = (char)true;
                }
            }
        }

        if(stretchOnly) {
            VertexData<double> stretchOnlyRatios = stretchOnlyAdjustmentFactors(geom);

            for(EdgePtr e : mesh->edges()) {
                double meanScale = 0.5 * (stretchOnlyRatios[e.halfedge().vertex()] + stretchOnlyRatios[e.halfedge().twin().vertex()]);
                targetL[e] *= meanScale;
            }
        }
        
    
        VertexData<Vector3> newPositions;
        for(int iter = 0; iter < 1; iter++) {
        
            // Iterative solve for positions satisfying these lenghts
            cout << "Solving for target lengths" << endl;
            PositionSolver p(geom);
            newPositions = p.solve(targetL, cutVert);
    
            for(VertexPtr v : mesh->vertices()) {
                (*geom)[v] = newPositions[v];
            }
        
        }
        
        // Return a new geometry/mesh
        HalfedgeMesh* newMesh = mesh->copy();
        Geometry<Euclidean>* devG = new Geometry<Euclidean>(*newMesh);
        for(size_t i = 0; i < mesh->nVertices(); i++) {
            (*devG)[newMesh->vertex(i)] = newPositions[mesh->vertex(i)];
        }
    
    
        // Check just how developable it is
        double maxDefect = -1;
        double meanDefectMag = 0;
        VertexData<double> defects;
        devG->getVertexAngleDefects(defects);
        for(size_t i = 0; i < mesh->nVertices(); i++) {
            VertexPtr oldV = mesh->vertex(i);
            VertexPtr newV = newMesh->vertex(i);
            if(!cutVert[oldV]) {
                double defect = std::abs(defects[newV]);
                maxDefect = std::max(maxDefect, defect);
                meanDefectMag += defect;
            }
        }
        meanDefectMag /= mesh->nVertices();
        cout << "Extrinsic developable approximation: Max defect = " << maxDefect << "  mean defect magnitude = " << meanDefectMag << endl;
        
    
        return devG;
    }


VertexData<double> yamabeScaleFactors(Geometry<Euclidean>* geom) {

    HalfedgeMesh* mesh = geom->getMesh();
    const size_t INVALID_IND = std::numeric_limits<size_t>::max();

    // Boundary or cut vertices
    VertexData<char> isB(mesh, false);
    for(VertexPtr v : mesh->vertices()) {
        if(v.isBoundary()) isB[v] = true;
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) isB[v] = true;
        }
    }

    // Index interior 
    VertexData<size_t> interiorInd(mesh, INVALID_IND);
    size_t nInterior = 0;
    for(VertexPtr v : mesh->vertices()) {
        if(!isB[v]) {
            interiorInd[v] = nInterior++;
        }
    }

    // Upper bound on degree of each vertex is # neighboring faces + 1.
    // Add one to get number of entries in row.
    std::vector<size_t> degreeEst(nInterior);
    for(VertexPtr v : mesh->vertices()) {
        if(!isB[v]) {
            degreeEst[interiorInd[v]] = v.degree() + 1;
        }
    }

    // Allocate the Cholesky matrix and right hand side
    FastCholesky interiorLaplacian(nInterior);
    interiorLaplacian.reserveColumns(degreeEst);
    

    // Useful values
    VertexData<double> k; geom->getVertexAngleDefects(k);
    DualFaceData<double> dualArea; geom->getDualFaceAreas(dualArea);
    for(VertexPtr v : mesh->vertices()) {
        if(isB[v]) {
            k[v] = 0.0;
        } else {
            // k[v] /= dualArea[v.dual()];
        }
    }
    EdgeData<double> cotanWeight; geom->getEdgeCotanWeights(cotanWeight);
    FaceData<double> area; geom->getFaceAreas(area);

    // Build the Poisson problem (Dirichlet-0 BC)
    for(EdgePtr e : mesh->edges()) {

        size_t indA = interiorInd[e.halfedge().vertex()];
        size_t indB = interiorInd[e.halfedge().twin().vertex()];
        double w = cotanWeight[e];
        if(!std::isfinite(w)) {
            throw std::runtime_error("Bad cotan weight");
        }

        if(indA != INVALID_IND) {
            interiorLaplacian.addValue(indA, indA, w);
        }
        if(indB != INVALID_IND) {
            interiorLaplacian.addValue(indB, indB, w);
        }
        if((indA != INVALID_IND) && (indB != INVALID_IND)) {
            interiorLaplacian.addValue(indA, indB, -w);
            interiorLaplacian.addValue(indB, indA, -w);
        }
    }
    interiorLaplacian.shiftDiagonal(1e-4);
    interiorLaplacian.factor();

    // Integrate to get the RHS
    GC::DenseVector<double> rhs(nInterior);
    /*
    for(FacePtr f : mesh->faces()) {

        double coef = area[f] / 12;

        for(VertexPtr va : f.adjacentVertices()) {
            size_t indA = interiorInd[va];
            if(indA == INVALID_IND) continue;
            for(VertexPtr vb : f.adjacentVertices()) {
                if(vb == va) {
                    rhs(indA) += -2.0 * k[va] * coef; 
                } else {
                    rhs(indA) += -k[vb] * coef; 
                }
            }
        }
    }
    */
    

    for(VertexPtr v : mesh->vertices()) {
        size_t i = interiorInd[v];
        if(i != INVALID_IND) {
            rhs(i) = -k[v];
        }
    }

    // Solve
    GC::DenseVector<double> u = interiorLaplacian.solve(rhs);

    VertexData<double> scaleFactors(mesh);
    for(VertexPtr v : mesh->vertices()) {
        if(isB[v]) {
            scaleFactors[v] = 1.0;
        } else {
            if(!std::isfinite(u(interiorInd[v]))) {
                throw std::runtime_error("Non finite scale factor");
            }
            scaleFactors[v] = std::exp(u(interiorInd[v]));
        }
    }

    return scaleFactors;
}

// Compute an amount to scale target edge lengths by such that they always need to stretch to map to the target geometry
VertexData<double> stretchOnlyAdjustmentFactors(Geometry<Euclidean>* geom) {

    // Get the initial Yamabe scale factors (TODO could reuse Laplacian from below)
    HalfedgeMesh* mesh = geom->getMesh();
    VertexData<double> initialYamabeFactors = yamabeScaleFactors(geom);
    for(VertexPtr v : mesh->vertices()) {
        initialYamabeFactors[v] = std::log(initialYamabeFactors[v]); // convert back to log form
    }


    // Boundary or cut vertices
    VertexData<char> isB(mesh, false);
    VertexData<double> bVal(mesh, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> fastTriBvals;
    int adjCount = 0;
    for(VertexPtr v : mesh->vertices()) {
        
        // Typical isometric boundary vertices
        if(v.isBoundary()) {
            isB[v] = true;
            bVal[v] = 0;
        } 
        for(EdgePtr e : v.adjacentEdges()) {
            if(e.isCut()) { 
                isB[v] = true;
                bVal[v] = 0;
            }
        }
       
        // Adjustment vertices, for which we are trying to find the adjustment scale factor
        if(initialYamabeFactors[v] > 0)  {
            isB[v] = true;
            bVal[v] = -initialYamabeFactors[v];
            adjCount++;
        }

        fastTriBvals.push_back(bVal[v]);
    }
    cout << "Adjusting " << adjCount << " scale factors to be positive" << endl;

    // Solve
    FastTriangleSoup fastTri(geom);
    fastTri.solveLaplaceDirichlet(fastTriBvals);

    // Copy over the adjustment factors
    VertexData<size_t> vInd = mesh->getVertexIndices();
    VertexData<double> ratios(mesh);
    for(VertexPtr v : mesh->vertices()) {

        double initialScaleFactor = std::exp(initialYamabeFactors[v]);
        double stretchScaleFactor = std::exp(initialYamabeFactors[v] + fastTriBvals[vInd[v]]);
        double ratio = stretchScaleFactor / initialScaleFactor;
        ratios[v] = ratio;

    }

    return ratios;
}
    

void computeCylinderParameterization(Geometry<Euclidean>* geom, CornerData<Vector2>& param) {

    cout << "Computing disk parameterization" << endl;
    HalfedgeMesh* mesh = geom->getMesh();
    param = CornerData<Vector2>(mesh);

    // === Split in to patches 
    std::vector<Geometry<Euclidean>*> splitGeoms = splitAlongCuts(geom);
    cout << "Split in to " << splitGeoms.size() << " components" << endl;


    // Build map to transfer paramerizations
    typedef std::pair<Vector3, Vector3> PointPair;
    std::unordered_map<PointPair, HalfedgePtr> halfedgeLookup;
    for(HalfedgePtr he : mesh->halfedges()) {
        PointPair p = std::make_pair(geom->position(he.vertex()), geom->position(he.twin().vertex()));
        halfedgeLookup[p] = he;
    }

    // Special case to ensure that if there are two regions, the parameterizations align
    bool alignedParam = (splitGeoms.size() == 2);
    Vector3 startingVertPos;
    double loopLen;
    size_t iGeom = 0;

    // Parameterize and transfer result
    for(Geometry<Euclidean>* splitGeom : splitGeoms) {

        // == 5 == Parameterize each component
        HalfedgeMesh* splitMesh = splitGeom->getMesh();
        CornerData<Vector2> splitParam(splitMesh);
    
        // === Compute the U coordinate as distance from the boundary, using FMM
        std::vector<std::pair<VertexPtr, double>> initialDistances;
        for(VertexPtr v : splitMesh->vertices()) {
            if(v.isBoundary()) {
                initialDistances.push_back(std::make_pair(v, 0));
            }
        }
        VertexData<double> bDist = GC::FMMDistance(splitGeom, initialDistances);
        for(VertexPtr v : splitMesh->vertices()) {
            for(CornerPtr c : v.adjacentCorners()) {
                splitParam[c].x = bDist[v];
            }
        }

        // === Compute the V coordinate by walking the boundary, then using harmonic interpolation

        // Initialize by walking boundary
        double cumDist = 0;
        VertexData<double> vertTheta(splitMesh, -1);

        // Special case to line up coordinates for two segments
        if(alignedParam) {

            // Start the first one anywhere
            if(iGeom == 0) { 
                startingVertPos = splitGeom->position((*(splitMesh->boundaryLoop(0).adjacentHalfedges().begin())).vertex()); // this monstrosity gets the position of the frst vertex in the boundary loop
                for(HalfedgePtr he : splitMesh->boundaryLoop(0).adjacentHalfedges()) {
                    vertTheta[he.vertex()] = cumDist;
                    cumDist += splitGeom->length(he.edge());
                }
                loopLen = cumDist;
            } 
            // Start the second loop where the first one started, and go backwards 
            else {

                // Find the halfedge to start on
                HalfedgePtr startingHe;
                for(HalfedgePtr he : splitMesh->boundaryLoop(0).adjacentHalfedges()) {
                    if(splitGeom->position(he.vertex()) == startingVertPos) {
                        startingHe = he;
                        break;
                    }
                }

                // Walk around the boundary, parameterizing backwards
                HalfedgePtr currHe = startingHe;
                cumDist = loopLen;
                do {
                    vertTheta[currHe.vertex()] = loopLen;
                    loopLen -= splitGeom->length(currHe.edge());
                    currHe = currHe.next();
                } while(currHe != startingHe);
                vertTheta[startingHe.vertex()] = 0.0; // we want coordinates on [0,-), not (0,-]

            }
        } 
        // General case 
        else {
            for(HalfedgePtr he : splitMesh->boundaryLoop(0).adjacentHalfedges()) {
                vertTheta[he.vertex()] = cumDist;
                cumDist += splitGeom->length(he.edge());
            }
        }



        VertexData<double> vertThetaX(splitMesh, 0);
        VertexData<double> vertThetaY(splitMesh, 0);
        for(HalfedgePtr he : splitMesh->boundaryLoop(0).adjacentHalfedges()) {
            VertexPtr v = he.vertex();
            vertTheta[v] /= cumDist;
            // double angle = vertTheta[v] * 2 * PI;
            // vertThetaX[v] = std::cos(angle);
            // vertThetaY[v] = std::sin(angle);
            // vertThetaX[v] = angle;
            // vertThetaY[v] = std::sin(angle);
        }

        // == Extend

        /*
        // Extend with gradient perp to distance function gradient
        size_t nInteriorVertices = 0;
        VertexData<size_t> interiorInd(splitMesh);
        for(VertexPtr v : splitMesh->vertices()) {
            if(!v.isBoundary()) {
                interiorInd[v] = nInteriorVertices++;
            }
        }
        GC::SparseMatrix<Complex> M(splitMesh->nFaces(), nInteriorVertices);
        size_t iFace = 0;
        // GC::DenseMatrix<double> rhsX(splitMesh->nFaces());
        // GC::DenseMatrix<double> rhsY(splitMesh->nFaces());
        GC::DenseMatrix<Complex> rhs(splitMesh->nFaces());
        for(FacePtr f : splitMesh->faces()) {

            // Gather values
            std::array<VertexPtr, 3> vert = {{
                f.halfedge().vertex(),
                f.halfedge().next().vertex(),
                f.halfedge().next().next().vertex()
            }};
                
            Vector3 basisX = unit(splitGeom->position(vert[1]) - splitGeom->position(vert[0]));
            Vector3 basisY = basisX.rotate_around(splitGeom->normal(f), PI/2.0);
            double x1 = dot(basisX, splitGeom->position(vert[1]) - splitGeom->position(vert[0]));
            double x2 = dot(basisX, splitGeom->position(vert[2]) - splitGeom->position(vert[0]));
            double y2 = dot(basisY, splitGeom->position(vert[2]) - splitGeom->position(vert[0]));
            double w0 = bDist[vert[0]];
            double w1 = bDist[vert[1]];
            double w2 = bDist[vert[2]];

            // Build composite matrices
            Eigen::Matrix<double, 2, 2> A; 
            Eigen::Matrix<double, 2, 3> B;
            Eigen::Vector3d w;

            A << 1.0/x1, 0, -x2/(x1*y2), 1/y2;
            B << -1, 1, 0, -1, 0, 1;
            w << w0, w1, w2;

            Eigen::Vector3d c = (A*B*w).transpose()*A*B;
           
            for(int iV = 0; iV < 3; iV++) {

                // Move boundary values to righthand side
                if(vert[iV].isBoundary()) {
                    // rhsX(iFace) -= c[iV] * vertThetaX[vert[iV]];
                    // rhsY(iFace) -= c[iV] * vertThetaY[vert[iV]];
                    Complex cVal = vertThetaX[vert[iV]] + IM_I * vertThetaY[vert[iV]];
                    rhs(iFace) -= c[iV] * cVal; 
                } else {
                    size_t i = interiorInd[vert[iV]];
                    M(iFace, i) += c[iV];
                }

            } 

            iFace++;
        }

        // Solve
        GC::DenseVector<Complex> extendedC;
        GC::solveLeastSquares(M, extendedC, rhs);

        // Store result
        for(VertexPtr v : splitMesh->vertices()) {
            if(!v.isBoundary()) {
                vertThetaX[v] = extendedC(interiorInd[v]).real();
                vertThetaY[v] = extendedC(interiorInd[v]).imag();
            }
        }
       
        // Convert back to coords
        VertexData<double> VCoord(splitMesh);
        for(VertexPtr v : splitMesh->vertices()) {
            double xVal = vertThetaX[v];
            double yVal = vertThetaY[v];
            // cout << "xval = " << xVal << " yVal = " << yVal << endl;
            double angle = PI + std::atan2(yVal, xVal);
            double coord = angle*cumDist / (2*PI);
            VCoord[v] = coord;
        }
        */

        // Brute-force lookup the closest point on the boundary in R3
        VertexData<double> VCoord(splitMesh);
        for(VertexPtr v : splitMesh->vertices()) {
            if(!v.isBoundary()) {
                Vector3 vPos = splitGeom->position(v);

                // Walk around boundary to find closest
                double closestDistance2 = std::numeric_limits<double>::infinity();
                double interpVal;
                for(HalfedgePtr he : splitMesh->boundaryLoop(0).adjacentHalfedges()) {

                    VertexPtr v1 = he.vertex();
                    VertexPtr v2 = he.twin().vertex();
                    Vector3 v1Pos = splitGeom->position(v1);
                    Vector3 v2Pos = splitGeom->position(v2);

                    // Check distance
                    Vector3 d = v2Pos - v1Pos;
                    double t = dot(vPos - v1Pos, d) / norm2(d);
                    t = clamp(t, 0.0, 1.0);
                    Vector3 p = v1Pos + t * d;
                    double dist2 = norm2(p - vPos);

                    if(dist2 < closestDistance2) {
                        closestDistance2 = dist2;

                        // Set the new value
                        double v1Val = vertTheta[v1];
                        double v2Val = vertTheta[v2];

                        // Handle special case where this is the last edge and wraps around
                        bool backwardWinding = (alignedParam && iGeom == 1);
                        if(backwardWinding) {
                            if(v1Val == 0) {
                                v1Val += 1;
                            }
                        } else {
                            if(v2Val == 0) {
                                v2Val += 1;
                            }
                        }

                        interpVal = (1.0 - t) * v1Val + t * v2Val;
                    }
                }

                vertTheta[v] = interpVal;
            }

            VCoord[v] = vertTheta[v];
        }

        
        for(FacePtr f : splitMesh->faces()) {

            // Heuristic to avoid wrapping around inside a triangle: if any coords are more than halfway around, all are
            bool anyGreater = false;
            for(VertexPtr v : f.adjacentVertices()) {
                if(VCoord[v] > .75) {
                    anyGreater = true;
                }
            }

            for(CornerPtr c : f.adjacentCorners()) {
                double coord = VCoord[c.vertex()];
                if(coord < .25 && anyGreater) {
                    coord += 1;
                }
                // Coordinates in world scale
                splitParam[c].y = coord * cumDist; 

                // Put both coords in the [0-1] (ish) range
                // splitParam[c].y = coord;
                // splitParam[c].x /= cumDist;
            }
        }

        
        // == 6 == Transfer parameterization back
        for(HalfedgePtr he : splitGeom->getMesh()->halfedges()) {
            PointPair p = std::make_pair(splitGeom->position(he.vertex()), splitGeom->position(he.twin().vertex()));
            
            HalfedgePtr origHe = halfedgeLookup[p];
            if(he.edge().isCut() || (he.edge().isBoundary() && !origHe.edge().isBoundary())) {
                origHe.edge().markCut(true);
            }

            param[origHe.corner()] = splitParam[he.corner()];
        }

        iGeom++;
    }

}

// Greedy labelling of connected components
FaceData<int> computeFaceComponents(Geometry<Euclidean>* geom, const CornerData<Vector2>& param) {

    HalfedgeMesh* mesh = geom->getMesh();

    FaceData<int> faceComp(mesh, -1);
    int iComp = 0;
    for(FacePtr f : mesh->faces()) {
        if(faceComp[f] >= 0) continue;

        // Grow new component
        std::vector<FacePtr> toProc = {f};
        while(!toProc.empty()) {
            FacePtr currF = toProc.back();
            toProc.pop_back();
            if(faceComp[currF] == -1) {
                faceComp[currF] = iComp;

                // Add adjacent faces not across cut
                for(HalfedgePtr he : currF.adjacentHalfedges()) {
                    if(he.edge().isBoundary()) continue;

                    if(paramContinuousAcrossHalfedge(param, he)) {
                        FacePtr neighF = he.twin().face();
                        toProc.push_back(neighF);
                    }
                }
            }
        }

        iComp++;
    }

    return faceComp;
}

void writeBoundarySVG(std::string filename, Geometry<Euclidean>* geom, CornerData<Vector2>& param) {

    cout << "Writing boundary to SVG at " << filename << endl;

    HalfedgeMesh* mesh = geom->getMesh();

    // Open file for writing
    ofstream out(filename);
    out.precision(numeric_limits<double>::max_digits10);


    // Mark all cuts
    HalfedgeData<char> isCut(mesh);
    for(HalfedgePtr he : mesh->halfedges()) {
        isCut[he] = he.edge().isBoundary() || !paramContinuousAcrossHalfedge(param, he);
    }

    // Write header
    out << "<svg>"  << endl;
    out << "<style type=\"text/css\">" << endl;
    out << ".st0{fill:none;stroke:#000000;stroke-width:0.5;stroke-miterlimit:10;}" << endl;
    out << "</style>" << endl;

    // Walk boundaries writing point coordinates
    HalfedgeData<char> halfedgeWritten(mesh, false);
    for(HalfedgePtr he : mesh->halfedges()) {

        if(halfedgeWritten[he]) continue;

        if(isCut[he]) {

            // Start a new polygon
            out << "<polygon class=\"st0\" points=\" " << endl;

            // Walk boundary
            HalfedgePtr firstHe = he;
            HalfedgePtr currHe = firstHe;
            do {

                // Write out point at tail
                Vector2 pos = param[currHe.next().corner()];
                out << pos.x << "," << pos.y << " " << endl;
                halfedgeWritten[currHe] = true;

                // Advance
                currHe = currHe.next();
                while (!isCut[currHe]) {
                    currHe = currHe.twin().next();
                }

            } while(currHe != firstHe);

            // End polygon
            out << "\"/>" << endl;
        }
    }


    // Write footer
    out << "</svg>"  << endl;
    out.close();
}