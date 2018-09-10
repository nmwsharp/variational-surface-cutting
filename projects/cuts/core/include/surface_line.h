#pragma once

#include <vector>
#include <set>
#include <cmath>
#include <limits>
#include <array>
#include <tuple>

#include "utilities.h"
#include "geometry.h"


// Represenents an undirected graph embedded in a surface.
//
// Vertices are represented as a point in a face.
//
// Edges are represented by the two vertices they connect, as well as a squence of edge crossings
// giving a path between the two vertices; this path is implicitly (extrinsically) straight within
// each face.


class SurfaceLine {

    public:

        // === Inner classes ===

        // Forward declarations
        struct LineSegment;

        struct SegmentEndpoint {

            SegmentEndpoint(HalfedgePtr he_, double t_);

            LineSegment* segment;

            HalfedgePtr he;             // The edge crossed
            SegmentEndpoint* twin;      // The crossing on he.twin()
            double t;                   // The point [0,1] along the halfedge

            VertexPtr vertex;         // If this is non-null, then the segment terminates at a mesh vertex
            bool isVertexPt;
        };

        struct LineSegment {

            LineSegment(SegmentEndpoint* p1, SegmentEndpoint* p2);

            SegmentEndpoint* pStart;
            SegmentEndpoint* pEnd;

            FacePtr containingFace();

        };


        // === Surface Graph === 


        // Constructor and destructor
        SurfaceLine(Geometry<Euclidean> *geometry, std::vector<HalfedgePtr> initPath);
        SurfaceLine(Geometry<Euclidean> *geometry);
        ~SurfaceLine();


        // Does the actual work of construction and destruction
        void constructFromHalfedgePath(std::vector<HalfedgePtr> edges);
        void clearLine();

        // === Members
        SegmentEndpoint* firstP = nullptr;
        SegmentEndpoint* lastP = nullptr;


        // === High level methods

        // Relax the line to a geodesic 
        void geodesicRelax();

        std::vector<std::array<Vector3,2>> getLineSegments();
        std::vector<std::tuple<HalfedgePtr, double, HalfedgePtr, double>> getSegmentsInFaces();

        // === Utility methods

        // Geomtric helpers
        Vector3 getPositionR3(SegmentEndpoint* endpoint) const;
        Vector3 halfedgePointToR3(HalfedgePtr he, double t) const;

        // Misc data
        HalfedgeMesh *mesh;
        Geometry<Euclidean> *geometry;

    private:

        // Does the hard work for gedoesic-ing
        bool maybeFlipOneRing(VertexPtr v, SegmentEndpoint* entryPoint, SegmentEndpoint* exitPoint);

        // Rotated star and end

        // Utilities
        void deleteSegmentRange(SegmentEndpoint* start, SegmentEndpoint* end, bool inclusive);

        void clampEdgeCrossings(double dist); // epsilon to clamp t values away from vertices

        // Geometry computation helpers
        double closestPointOnLine(Vector3 lineStart, Vector3 lineVec, Vector3 pa, Vector3 pb);
        double endpointAngle(SegmentEndpoint* point, VertexPtr v);
        double circularMean(double theta1, double w1, double theta2, double w2);
        HalfedgePtr sharedHalfedge(VertexPtr start, VertexPtr end);
        HalfedgePtr halfedgeBetween(FacePtr f1, FacePtr f2);

        // Geometry caching
        void cacheGeometry();
        VertexData<Vector3> positions;
        HalfedgeData<double> angularCoordinates;

};