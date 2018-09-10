#include "surface_line.h"

#include <utility>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "timing.h"


using std::vector;



// ==== Surface Graph ====

SurfaceLine::SurfaceLine(Geometry<Euclidean> *geometry_, std::vector<HalfedgePtr> initPath) 
    : mesh(geometry_->getMesh()), geometry(geometry_)
{
    cacheGeometry();
    constructFromHalfedgePath(initPath);
}

SurfaceLine::SurfaceLine(Geometry<Euclidean> *geometry_)
    : mesh(geometry_->getMesh()), geometry(geometry_)
{
    cacheGeometry();
}

SurfaceLine::~SurfaceLine() {
    clearLine();
}

void SurfaceLine::cacheGeometry() {

    geometry->getVertexPositions(positions);
    geometry->getAngularCoordinates(angularCoordinates);

}

void SurfaceLine::clearLine() {
    if(firstP != nullptr && lastP != nullptr) {
        deleteSegmentRange(firstP, lastP, true);
        firstP = nullptr;
        lastP = nullptr;
    }
}

void SurfaceLine::constructFromHalfedgePath(std::vector<HalfedgePtr> halfedges) {

    cout << "Constructing from halfedge path with length " << halfedges.size() << endl;

    double EPS = 0.0001;

    // === First, build the sequence of corners that the segments will pass through
    std::vector<FacePtr> facePtrs;
    for(size_t iHe = 1; iHe < halfedges.size(); iHe++) {

        HalfedgePtr currHe = halfedges[iHe];
        HalfedgePtr prevHe = halfedges[iHe-1];

        if(currHe.vertex() != prevHe.twin().vertex()) {
            throw std::runtime_error("Bad halfedge chain");
        }

        // The vertex we will orbit around
        VertexPtr rootVert = currHe.vertex();

        // Collect all of the outgoing halfedge pointers in the 1-ring in to a clockwise list and a ccw list
        std::vector<HalfedgePtr> cwList;
        std::vector<HalfedgePtr> ccwList;
        std::vector<HalfedgePtr>* currList = &cwList;

        {
            HalfedgePtr he = prevHe.next();
            HalfedgePtr firstHe = he;
            do {
                currList->push_back(he);
                if(he == currHe) {
                    currList = &ccwList; 
                }
                he = he.twin().next();
            } while(he != firstHe);
        }


        // Since the mesh is manifold, one of these two corner sequences must avoid the boundary. Check which.
        bool cwListReal = true;
        for(HalfedgePtr he : cwList) {
            if(!he.isReal()) cwListReal = false;
        }
        bool ccwListReal = true;
        for(HalfedgePtr he : ccwList) {
            if(!he.isReal()) ccwListReal = false;
        }

        if(!cwListReal && !ccwListReal) {
            throw std::runtime_error("Neither vertex orbit is real... that's odd.");
        }


        // Use the real list (or just one if both were)
        if(cwListReal) {
            for(HalfedgePtr he : cwList) {
                facePtrs.push_back(he.face());
            }
        } else {
            for(long long int i = ccwList.size()-1; i >= 0; i--) {
                HalfedgePtr he = ccwList[i];
                facePtrs.push_back(he.face());
            }
        }
    }

    // Trim duplicates
    std::vector<FacePtr> distinctFacePtrs;
    for(FacePtr f : facePtrs) {
        if(distinctFacePtrs.size() == 0 || f != distinctFacePtrs.back()) {
            distinctFacePtrs.push_back(f);
        }
    }

    // === Second, stitch the face pointers to form a coherent sequence of segments

    // First endpoint
    firstP = new SegmentEndpoint(HalfedgePtr(), std::numeric_limits<double>::quiet_NaN());
    firstP->vertex = halfedges[0].vertex();
    firstP->isVertexPt = true;
    
    SegmentEndpoint* openEndpoint = firstP;
    for(long long int iFp = 0; iFp < ((long long int)distinctFacePtrs.size())-1; iFp++) {

        FacePtr currFace = distinctFacePtrs[iFp]; 
        FacePtr nextFace = distinctFacePtrs[iFp+1]; 

        HalfedgePtr crossingHe = halfedgeBetween(currFace, nextFace);

        SegmentEndpoint* closingEndpoint = new SegmentEndpoint(crossingHe, 0.5);
        LineSegment* newSegment = new LineSegment(openEndpoint, closingEndpoint);

        openEndpoint = new SegmentEndpoint(crossingHe.twin(), 0.5);
        closingEndpoint->twin = openEndpoint;
        openEndpoint->twin = closingEndpoint;
    }

    // Last endpoint / segment
    lastP = new SegmentEndpoint(HalfedgePtr(), std::numeric_limits<double>::quiet_NaN());
    lastP->vertex = halfedges.back().twin().vertex();
    lastP->isVertexPt = true;
    new LineSegment(openEndpoint, lastP);

}

/*
void SurfaceLine::repairEdgeAlignedSegments(std::vector<LineSegment*> &segments, std::vector<FacePoint*> &facePoints) {

    // Check each segment in the initial list
    int nEdgeAlignedSegments = 0;
    for(size_t i = 0; i < segments.size(); i++) {
        LineSegment* seg = segments[i];
        
        // If this segment terminates inside a face, or its endpoints lie on distinct halfedges, we have
        // no work to do.
        if(seg->p1->facePoint != nullptr || seg->p2->facePoint != nullptr || seg->p1->he != seg->p2->he) {
            continue;
        }

        // === Repair this segment
        nEdgeAlignedSegments++;

        // Create a new point slightly inside the the face
        Vector3 p1Bary = getPositionBary(seg->p1);
        Vector3 p2Bary = getPositionBary(seg->p2);
        Vector3 newPoint = p1Bary + p2Bary + Vector3::constant(0.001);
        newPoint /= (newPoint.x + newPoint.y + newPoint.z);

        // Split the segment and hook it to this new point
        SegmentEndpoint* right = seg->p2;
        FacePoint* newFacepoint = new FacePoint{seg->p1->he.face(), newPoint};
        SegmentEndpoint* newFaceEndpointLeft = new SegmentEndpoint{seg, HalfedgePtr(), nullptr, -777, newFacepoint};
        seg->p2 = newFaceEndpointLeft;

        SegmentEndpoint* newFaceEndpointRight = new SegmentEndpoint{nullptr, HalfedgePtr(), nullptr, -777, newFacepoint};
        LineSegment* newSegment = new LineSegment(newFaceEndpointRight, right);

        segments.push_back(newSegment);
        facePoints.push_back(newFacepoint);
    }

    cout << "Input contained " << nEdgeAlignedSegments << " edge aligned segments. Repaired by perturbing in to containing face." << endl;
}
*/


// Make sure we can hash pair types (stupid C++)
struct pairHash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto x = std::hash<T1>{}(p.first);
        auto y = std::hash<T2>{}(p.second);
        return x ^ y;  
    }
};



void SurfaceLine::deleteSegmentRange(SegmentEndpoint* start, SegmentEndpoint* end, bool inclusive) {

        // Acumulate lists to delete to avoid breaking data structure while traversing it
        vector<LineSegment*> segmentsToDelete;
        vector<SegmentEndpoint*> endpointsToDelete;

        // Travers the segments, marking points for deletion
        SegmentEndpoint* currEndpoint = start->segment->pEnd;
        while(true) {
            segmentsToDelete.push_back(currEndpoint->segment);

            if(currEndpoint == end) break;

            endpointsToDelete.push_back(currEndpoint);
            endpointsToDelete.push_back(currEndpoint->twin);

            currEndpoint = currEndpoint->twin->segment->pEnd;
        }

        // If requested, mark the bounding endpoints for deletion
        if(inclusive) {
            endpointsToDelete.push_back(start);
            endpointsToDelete.push_back(end);
        }

        // Actually delete
        // (if we deleted inline it would break the traversal
        for(LineSegment* l : segmentsToDelete) {
            delete l;
        }
        for(SegmentEndpoint* e : endpointsToDelete) {
            delete e;
        }

}

Vector3 SurfaceLine::halfedgePointToR3(HalfedgePtr he, double t) const {
    return (1.0 - t) * positions[he.vertex()] + t * positions[he.twin().vertex()];
}



Vector3 SurfaceLine::getPositionR3(SegmentEndpoint* endpoint) const {
    if(endpoint->isVertexPt) {
        return geometry->position(endpoint->vertex);
    } else {
        return halfedgePointToR3(endpoint->he, endpoint->t);
    }
}

double SurfaceLine::closestPointOnLine(Vector3 lineStart, Vector3 lineVec, Vector3 pa, Vector3 pb) {

    // TODO I derived this, someone else should check it

    double EPS = 1e-8;

    // Necessary distance from each point to the line
    Vector3 lineEnd = lineStart + lineVec;
    double ta = -dot(lineStart - pa, lineVec) / norm2(lineVec);
    double La = norm(cross(pa - lineStart, pa - lineEnd)) / norm(lineVec);
    double tb = -dot(lineStart - pb, lineVec) / norm2(lineVec);
    double Lb = norm(cross(pb - lineStart, pb - lineEnd)) / norm(lineVec);

    // If the points are both on the line, anywhere works
    if(La < EPS && Lb < EPS) {
        return (ta + tb) / 2.0;
    }

    // How much of the distance the first point should cover for optimality
    double za = La * (tb-ta) / (La + Lb);     

    // Convert this to a point along the line
    return ta + za; 
}



void SurfaceLine::geodesicRelax() {

    // cout << "Geodesic relaxing edge" << endl;
    
    double EPS = 1e-5;

    size_t pathLen = 0;
    { // Compute a reasionable number of iterations from the initial path length
        SegmentEndpoint* currBeginPoint = firstP;
        while(currBeginPoint != nullptr) {
            pathLen++;
            currBeginPoint = currBeginPoint->segment->pEnd->twin;
        }
    }

    size_t ITERMAX = 50*pathLen;


    // Can't relax the trivial path (and it creates a lot of special cases)
    if(pathLen == 1) {
        return;
    }

    // Iteratively relax
    double lastMaxChange = 1;
    size_t iIter = 0;
    while(lastMaxChange > EPS && iIter < ITERMAX) {
       
        // Watch the max change to decide when we are done 
        lastMaxChange = 0;
   
        // === Phase 1: Walk along the edge, relaxing each edge crossing along its edge
        // towards the geodesic minimum
        {
            // Performance heuristic: do this step several times for each topology check
            int nWalks = 20;
            for(int iWalk = 0; iWalk < nWalks; iWalk++) {
                SegmentEndpoint* firstMidpoint = firstP->segment->pEnd;
                SegmentEndpoint* secondMidpoint = firstMidpoint->twin;
                while(secondMidpoint) {
                    
                    // Accumulate the endpoints
                    SegmentEndpoint* firstOuter = firstMidpoint->segment->pStart;
                    SegmentEndpoint* secondOuter = secondMidpoint->segment->pEnd;

                    // Get the points corresponding to each
                    Vector3 pa = getPositionR3(firstOuter);
                    Vector3 pb = getPositionR3(secondOuter);

                    // Get the halfedge line that they share
                    Vector3 lineStart = positions[firstMidpoint->he.vertex()];
                    Vector3 lineVec = positions[secondMidpoint->he.vertex()] - lineStart;
            
                    // Find the closest point along the edge
                    double tClosest = closestPointOnLine(lineStart, lineVec, pa, pb);
                    tClosest = clamp(tClosest, EPS, 1.0 - EPS);

                    // Record the change
                    double myChange = std::abs(tClosest - firstMidpoint->t);
                    lastMaxChange = std::max(lastMaxChange, myChange);

                    // Update the crossing position
                    firstMidpoint->t = tClosest; 
                    secondMidpoint->t = 1.0 - tClosest;

                    // Advance to next segment
                    firstMidpoint = secondOuter;
                    secondMidpoint = firstMidpoint->twin;
                }
            }
        }

        // === Phase 2: Look for topological improvements of the edge, where we need to add/remove 
        // segments or shift them in to neighboring triangles.
        // There are two things to look for here:
        //   - If the edge ever has multiple segments in the same face, we should take a shortcut
        //     and insert a new halfedge directly between 
        //   - If the edge comes very close to a vertex, flip over that vertex (TODO use a better condition
        //     here, maybe angle-based?)

        // === First, resolve multi-segment-in-face shortcuts
        {
            unordered_map<FacePtr, SegmentEndpoint*> lastFaceEndpoint;

            // Iterate over the segments accumulating maps of the last time we cross each face
            SegmentEndpoint* currBeginPoint = firstP;
            while(currBeginPoint != nullptr) {

                LineSegment* currSegment = currBeginPoint->segment; 
                FacePtr f = currBeginPoint->segment->containingFace();

                // Mark the last endpoints always
                lastFaceEndpoint[f] = currSegment->pEnd;

                // Advance to next segment
                currBeginPoint = currSegment->pEnd->twin;
            }

            // The last point is the exit point for all adjacent faces
            for(FacePtr f : lastP->vertex.adjacentFaces()) {
                lastFaceEndpoint[f] = lastP;
            }

            // Make a second pass through the segments, taking shortcuts within faces when appropriate.
            // NOTE: We must do this with a walk through the edge, rather than iterating through the map,
            //     because taking one shortcut might render another shortcut invalid.

            // The first point is a start point for all adjacent faces
            for(FacePtr f : firstP->vertex.adjacentFaces()) {
                if(f != firstP->segment->containingFace() && lastFaceEndpoint.find(f) != lastFaceEndpoint.end()) {

                    lastMaxChange = 1.0;
                
                    SegmentEndpoint* lastEndpoint = lastFaceEndpoint[f];

                    // Delete the intermediate segments which are no longer needed
                    deleteSegmentRange(firstP, lastEndpoint, false);

                    // Create a new segment between the points
                    LineSegment* newSegment = new LineSegment(firstP, lastEndpoint);
                    
                }
            }


            currBeginPoint = firstP->segment->pEnd->twin;
            while(currBeginPoint != nullptr) {

                FacePtr f = currBeginPoint->segment->containingFace();

                // Test if we can take a shortcut here
                if(lastFaceEndpoint[f]->segment != currBeginPoint->segment) {

                    lastMaxChange = 1.0;
                
                    SegmentEndpoint* firstEndpoint = currBeginPoint;
                    SegmentEndpoint* lastEndpoint = lastFaceEndpoint[f];

                    // Delete the intermediate segments which are no longer needed
                    deleteSegmentRange(firstEndpoint, lastEndpoint, false);

                    // Create a new segment between the points
                    LineSegment* newSegment = new LineSegment(firstEndpoint, lastEndpoint);

                }

                // Advance to next segment
                currBeginPoint = currBeginPoint->segment->pEnd->twin;
            }
        }
        
        
        // === Second, resolve vertex-flip shortcuts
        {
            unordered_map<VertexPtr, SegmentEndpoint*> oneRingExits;

            // Iterate over the segments accumulating maps of the last time we exit a 1-ring
            SegmentEndpoint* currBeginPoint = firstP;
            while(currBeginPoint != nullptr) {

                FacePtr f = currBeginPoint->segment->containingFace();

                // Mark the last endpoints
                SegmentEndpoint* laterEndpoint = currBeginPoint->segment->pEnd;
                if(laterEndpoint->isVertexPt) {
                    // Need a special case for the last endpoint 
                    for(VertexPtr v : f.adjacentVertices()) {
                        if(v != laterEndpoint->vertex) {
                            oneRingExits[v] = laterEndpoint;
                        }
                    }
                } else {
                    oneRingExits[laterEndpoint->he.next().next().vertex()] = laterEndpoint;
                }

                // Advance to next segment
                currBeginPoint = currBeginPoint->segment->pEnd->twin;
            }

            // Make a second pass through the segments, taking shortcuts within 1-rings when appropriate.
            // FIXME it is possible for this to fail on squiggly paths if it removes a 1-ring exit which
            // it later tries to utilize
            currBeginPoint = firstP;
            while(currBeginPoint != nullptr) {

                // The 1-rings we might be entering


                // Special case for the first endpoint 
                if(currBeginPoint->isVertexPt) {

                    FacePtr f = currBeginPoint->segment->containingFace();
                    bool flipped = false;
                    for(VertexPtr v : f.adjacentVertices()) {
                        if(v != currBeginPoint->vertex) {

                            // Get the points at which we enter exit the 1-ring
                            SegmentEndpoint* firstEndpoint = currBeginPoint;
                            SegmentEndpoint* lastEndpoint = oneRingExits[v];
                            if(lastEndpoint == nullptr) {
                                break;
                            }
                            
                            // Try to flip over the vertex if appropriate
                            bool vertexFlipped = maybeFlipOneRing(v, firstEndpoint, lastEndpoint); 
                            
                            if(vertexFlipped) {
            
                                lastMaxChange = 1.0;
            
                                // Resume processing from the newly placed segment
                                currBeginPoint = lastEndpoint->twin;
                                flipped = true;
                                break;
                            }
                            
                        }
                    }

                    if(flipped) {
                        continue;                        
                    }

                } else {

                    FacePtr f = currBeginPoint->he.face();
                    VertexPtr oneRingV = currBeginPoint->he.next().next().vertex();

                    // Get the points at which we enter exit the 1-ring
                    SegmentEndpoint* firstEndpoint = currBeginPoint;
                    SegmentEndpoint* lastEndpoint = oneRingExits[oneRingV];
                    if(lastEndpoint == nullptr) {
                        break;
                    }
                    
                    // Try to flip over the vertex if appropriate
                    bool vertexFlipped = maybeFlipOneRing(oneRingV, firstEndpoint, lastEndpoint); 
                    
                    if(vertexFlipped) {
    
                        lastMaxChange = 1.0;
    
                        // Resume processing from the newly placed segment
                        currBeginPoint = lastEndpoint->twin;
                        continue;
                    }

                }

                // Advance to next segment
                currBeginPoint = currBeginPoint->segment->pEnd->twin;
            }
        }


        iIter++;
    }


    if(iIter == ITERMAX) {
        cout << "WARNING: geodesicRelax() terminated on iteration count" << endl;
    }

}
    

// Flip the edge over a given vertex "v". The entryPoint is where the edge enters the one-ring around 
// the vertex, and exitPoint is where it leaves it. Both of these should be the "inner" of the two entry
// points along the boundary.
// TODO Absolutely need to verify the the angular logic in this is correct valid math. Until I come up with a 
// proof, this should be called a heuristic.
bool SurfaceLine::maybeFlipOneRing(VertexPtr v, SegmentEndpoint* entryPoint, SegmentEndpoint* exitPoint) {

    double DELTA = 0.01; // new crossings get put this close to the center vertex

    // If the endpoints are already in the same face, we have no work to do here
    // (fyi: if a collapse were helpful, it would be caught by the face-shortcut check)
    if(entryPoint->segment->containingFace() == exitPoint->segment->containingFace()) {
        return false; 
    }

    // Don't try to flip across the boundary
    if(v.isBoundary()) {
        return false;
    }

    // First, compute the angles of the endpoints in the 1-ring
    double entryAngle = endpointAngle(entryPoint, v);
    double exitAngle = endpointAngle(exitPoint, v);

    // Decide whether or not to flip, based on the angles
    double angularDistance = regularizeAngle(exitAngle - entryAngle);
            
    // Special case for node
    HalfedgePtr nextHe = entryPoint->segment->pEnd->he;

    // The shortest path is clockwise
    if(angularDistance > PI) {

        if(nextHe.vertex() != v) {
            
            // Wrong direction, flip needed!
           
            FacePtr targetFace = exitPoint->segment->containingFace();

            // Delete all of the old segments
            deleteSegmentRange(entryPoint, exitPoint, false);

            // Walk around the vertex, inserting new segments until we reach the target halfedge
            SegmentEndpoint* currPt = entryPoint;
            
            // Special case for first node
            HalfedgePtr currHe;
            if(entryPoint->isVertexPt) {
                HalfedgePtr sharedHe = sharedHalfedge(entryPoint->vertex, v);
                currHe = sharedHe.next();
            }  else {
                currHe = entryPoint->he.next().next();
            }


            do {

                // Create a new endpoint on this side of the edge
                SegmentEndpoint* newEndpoint = new SegmentEndpoint(currHe, DELTA);

                // Create a new segment joining the edges
                LineSegment* newSegment = new LineSegment(currPt, newEndpoint);

                // Create a new endpoint on the far side of the edge
                currPt = new SegmentEndpoint(currHe.twin(), 1.0 - DELTA);
                newEndpoint->twin = currPt;
                currPt->twin = newEndpoint;


                // Advance to the next crossing halfedge
                currHe = currHe.twin().next();
            } while(currPt->he.face() != targetFace);

            // Create a final segment, linking to the exit point
            new LineSegment(currPt, exitPoint);

            return true;
        } else {
            // Correct direction, no flip needed
            return false;
        }
    } 
    // The shortest path is counter-clockwise 
    else {
        if(nextHe.next().vertex() == v) {
            // Correct direction, no flip needed
            return false;
        } else {
            // Wrong direction, flip needed!

            FacePtr targetFace = exitPoint->segment->containingFace();

            // Delete all of the old segments
            deleteSegmentRange(entryPoint, exitPoint, false);

            // Walk around the vertex, inserting new segments until we reach the target halfedge
            SegmentEndpoint* currPt = entryPoint;

            // Special case for first endpoint 
            HalfedgePtr currHe;
            if(entryPoint->isVertexPt) {
                HalfedgePtr sharedHe = sharedHalfedge(v, entryPoint->vertex);
                currHe = sharedHe.next().next();
            }  else {
                currHe = entryPoint->he.next();
            }

            do {

                // Create a new endpoint on this side of the edge
                SegmentEndpoint* newEndpoint = new SegmentEndpoint(currHe, 1.0 - DELTA);

                // Create a new segment joining the edges
                LineSegment* newSegment = new LineSegment(currPt, newEndpoint);

                // Create a new endpoint on the far side of the edge
                currPt = new SegmentEndpoint(currHe.twin(), DELTA);
                newEndpoint->twin = currPt;
                currPt->twin = newEndpoint;

                // Advance to the next crossing halfedge
                currHe = currHe.twin().next().next();
            } while(currPt->he.face() != targetFace);

            // Create a final segment, linking to the exit point
            new LineSegment(currPt, exitPoint);

            return true;
        }

    }

}

std::vector<std::array<Vector3,2>> SurfaceLine::getLineSegments() {

    std::vector<std::array<Vector3, 2>> line;
    
    // Clamp a bit tighter to mitigate degeneracies
    clampEdgeCrossings(0.01); 

    SegmentEndpoint* currBeginPoint = firstP;
    while(currBeginPoint != nullptr) {

        Vector3 pA = getPositionR3(currBeginPoint);
        Vector3 pB = getPositionR3(currBeginPoint->segment->pEnd);
        line.push_back({{pA, pB}});

        currBeginPoint = currBeginPoint->segment->pEnd->twin;
    }

    return line; 
}

// Get all of the segments in the relaxed line.
// All segments are clamped to a finite distance away from vertices, except the beginning and
// terminating endpoints which are indicated by t=0
std::vector<std::tuple<HalfedgePtr, double, HalfedgePtr, double>> SurfaceLine::getSegmentsInFaces() {

    // Clamp a bit tighter to mitigate degeneracies
    clampEdgeCrossings(0.01); 
    typedef std::tuple<HalfedgePtr, double, HalfedgePtr, double> Entry;
    std::vector<Entry> result;


    // Verify that the start and end are not degnerate faces
    if(firstP->segment->pEnd != lastP) {

        HalfedgePtr heFirst = firstP->segment->pEnd->he;
        if(heFirst.vertex() == firstP->vertex || heFirst.twin().vertex() == firstP->vertex) {
            throw std::runtime_error("Degenerate initial line edge");
        }
        
        HalfedgePtr heLast = lastP->segment->pStart->he;
        if(heFirst.vertex() == lastP->vertex || heFirst.twin().vertex() == lastP->vertex) {
            throw std::runtime_error("Degenerate final line edge");
        }

    }

    // Special case for a degenerate single-segment line
    if(firstP->segment->pEnd == lastP) {

        FacePtr containingFace = firstP->segment->containingFace();

        HalfedgePtr startHe;
        double tStart = 0.0;
        HalfedgePtr endHe;
        double tEnd = 0.0;
        for(HalfedgePtr he : containingFace.adjacentHalfedges()) {
            if(he.vertex() == firstP->vertex) {
                startHe = he;
            }
            if(he.vertex() == lastP->vertex) {
                endHe = he;
            }
        }

        result.push_back(std::make_tuple(startHe, tStart, endHe, tEnd));
        return result;
    }

    SegmentEndpoint* currBeginPoint = firstP;
    while(currBeginPoint != nullptr) {

        SegmentEndpoint* pStart = currBeginPoint;
        SegmentEndpoint* pEnd = currBeginPoint->segment->pEnd;

        HalfedgePtr startHe;
        double tStart;
        if(pStart == firstP) {
            startHe = pStart->segment->pEnd->he.next().next();
            tStart = 0.0;
        } else {
            startHe = pStart->he;
            tStart = pStart->t;
        }

        HalfedgePtr endHe;
        double tEnd;
        if(pEnd == lastP) {
            endHe = pEnd->segment->pStart->he.next().next();
            tEnd = 0.0;
        } else {
            endHe = pEnd->he;
            tEnd = pEnd->t;
        }

        result.push_back(std::make_tuple(startHe, tStart, endHe, tEnd));

        currBeginPoint = currBeginPoint->segment->pEnd->twin;
    }
    
    return result;
}

void SurfaceLine::clampEdgeCrossings(double dist) {
    
    SegmentEndpoint* currBeginPoint = firstP;
    while(currBeginPoint != nullptr) {

        currBeginPoint->t = clamp(currBeginPoint->t, dist, 1.0-dist);
        if(currBeginPoint->twin != nullptr) {
            currBeginPoint->twin->t = clamp(currBeginPoint->twin->t, dist, 1.0-dist);
        }

        currBeginPoint = currBeginPoint->segment->pEnd->twin;
    }
}
            
SurfaceLine::SegmentEndpoint::SegmentEndpoint(HalfedgePtr he_, double t_) :
    segment(nullptr), he(he_), twin(nullptr), t(t_), vertex(VertexPtr()), isVertexPt(false)
{

}

SurfaceLine::LineSegment::LineSegment(SegmentEndpoint* p1_, SegmentEndpoint* p2_) 
    : pStart(p1_), pEnd(p2_)
{
    pStart->segment = this;
    pEnd->segment = this;
}

FacePtr SurfaceLine::LineSegment::containingFace() {
    
    if(!pStart->isVertexPt) {
        return pStart->he.face();
    }
    if(!pEnd->isVertexPt) {
        return pEnd->he.face();
    }

    for(HalfedgePtr he : pStart->vertex.outgoingHalfedges()) {
        if(he.twin().vertex() == pEnd->vertex) {
            if(he.isReal()) {
                return he.face();
            } else {
                return he.twin().face();
            }
        }
    }

    throw std::runtime_error("Invalid single segment line");
}
        
double SurfaceLine::endpointAngle(SegmentEndpoint* point, VertexPtr v) {

    // Handle vertex points separately
    if(point->isVertexPt) {
        for(HalfedgePtr he : v.outgoingHalfedges()) {
            if(point->vertex == he.twin().vertex()) {
                return angularCoordinates[he];
            }
        }
    }


    if(point->he.vertex() == v) {
        return angularCoordinates[point->he];
    }
    if(point->he.twin().vertex() == v) {
        return angularCoordinates[point->he.twin()];
    }
    if(point->he.next().twin().vertex() == v) {
        double c1 = angularCoordinates[point->he.next().twin()];
        double c2 = angularCoordinates[point->he.next().next()];
        return circularMean(c1, point->t, c2, 1.0 - point->t);
    }
    throw std::runtime_error("Invalid endpoint angle call");
}
        
double SurfaceLine::circularMean(double theta1, double w1, double theta2, double w2) {
    return std::arg(w1 * std::exp(IM_I * theta1) + 
                    w2 * std::exp(IM_I * theta2));
}        

HalfedgePtr SurfaceLine::sharedHalfedge(VertexPtr start, VertexPtr end) {

    for(HalfedgePtr he : start.outgoingHalfedges()) {
        if(he.twin().vertex() == end) {
            return he;
        }
    }

    throw std::runtime_error("Bad shared halfedge call");
}

HalfedgePtr SurfaceLine::halfedgeBetween(FacePtr f1, FacePtr f2) {

    for(HalfedgePtr he : f1.adjacentHalfedges()) {
        if(he.twin().face() == f2) {
            return he;
        }
    }
    throw std::runtime_error("Bad halfedge between call");
}
