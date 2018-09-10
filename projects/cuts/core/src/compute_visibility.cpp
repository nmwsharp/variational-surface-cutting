#include "compute_visibility.h"

#include "mesh_ray_tracer.h"
#include "combining_hash_functions.h"

#include <unordered_map>


std::vector<Vector3> getIcosphereVertices(Vector3 center, double radius, size_t atLeastNPoints) {

    // Initial icosahedron
    std::vector<Vector3> pts = {
        Vector3{0,-0.525731,0.850651},
        Vector3{0.850651,0,0.525731},
        Vector3{0.850651,0,-0.525731},
        Vector3{-0.850651,0,-0.525731},
        Vector3{-0.850651,0,0.525731},
        Vector3{-0.525731,0.850651,0},
        Vector3{0.525731,0.850651,0},
        Vector3{0.525731,-0.850651,0},
        Vector3{-0.525731,-0.850651,0},
        Vector3{0,-0.525731,-0.850651},
        Vector3{0,0.525731,-0.850651},
        Vector3{0,0.525731,0.850651}
    };

    typedef std::array<size_t, 3> Tri;
    std::vector<Tri> faces = {
        {{2,3,7}},
        {{2,8,3}},
        {{4,5,6}},
        {{5,4,9}},
        {{7,6,12}},
        {{6,7,11}},
        {{10,11,3}},
        {{11,10,4}},
        {{8,9,10}},
        {{9,8,1}},
        {{12,1,2}},
        {{1,12,5}},
        {{7,3,11}},
        {{2,7,12}},
        {{4,6,11}},
        {{6,5,12}},
        {{3,8,10}},
        {{8,2,1}},
        {{4,10,9}},
        {{5,9,1}}
    };
    // Shift to 0-indexing
    for(Tri& t : faces) {
        for(int i = 0; i < 3; i++) t[i]--;
    }


    // Subdivide until we have enough vertices
    while(pts.size() < atLeastNPoints) {

        // === Subdivide

        // Create new vertices
        std::unordered_map<std::pair<size_t, size_t>, size_t> newVertMap;
        for(Tri t : faces) {
            for(int i = 0; i < 3; i++) {
                size_t vA = t[i];
                size_t vB = t[(i+1)%3];
                if(vA > vB) continue;

                // Create a new vert along the edge
                Vector3 newP = (pts[vA] + pts[vB]) / 2.0;
                pts.push_back(newP);
                std::pair<size_t, size_t> key = std::make_pair(vA, vB);
                newVertMap[key] = pts.size()-1;
            }
        }

        // Normalize points
        for(Vector3& v : pts) {
            v = unit(v);
        }

        // Don't both creating last set of faces
        if(pts.size() >= atLeastNPoints) break;

        // Create new faces
        std::vector<Tri> newFaces;
        for(Tri t : faces) {

            // 1-4 subdivision

            // Gather indices
            size_t v0 = t[0];
            size_t v1 = t[1];
            size_t v2 = t[2];
            size_t vn0 = newVertMap[std::make_pair(std::min(v0, v1), std::max(v0, v1))];
            size_t vn1 = newVertMap[std::make_pair(std::min(v1, v2), std::max(v1, v2))];
            size_t vn2 = newVertMap[std::make_pair(std::min(v2, v0), std::max(v2, v0))];

            // Create 4 faces
            newFaces.push_back({{v0, vn0, vn2}});
            newFaces.push_back({{v1, vn1, vn0}});
            newFaces.push_back({{v2, vn2, vn1}});
            newFaces.push_back({{vn0, vn1, vn2}});
        }

        faces = newFaces;
    }


    // Reposition and scale
    for(Vector3& v : pts) {
        v *= radius;
        v += center;
    }
    

    return pts;
}


VertexData<double> computeVisibilityFraction(Geometry<Euclidean>* geometry, size_t nSpherePointsMin, bool headOn) {


    // Build the bvh
    MeshRayTracer tracer(geometry);

    // Get the points at which we will test
    std::vector<Vector3> pts;
    double scale = geometry->lengthScale();
    Vector3 headOnCameraPos = geometry->center() + Vector3{0, 3*scale, -10 * scale}; 
    if(headOn) { 
        double sphereRad = 0.2 * scale;
        pts = getIcosphereVertices(headOnCameraPos, sphereRad, nSpherePointsMin);
    } else {
        double sphereRad = scale * 10;
        Vector3 sphereC = geometry->center();
        pts = getIcosphereVertices(sphereC, sphereRad, nSpherePointsMin);
    }
    size_t nPts = pts.size();

    // Trace a lot of points
    VertexData<double> viewFrac(geometry->getMesh());
    for(VertexPtr v : geometry->getMesh()->vertices()) {

        // Pull slightly in to an adjacent face to avoid ambiguities
        Vector3 vPos = geometry->position(v);
        FacePtr adjFace = v.halfedge().face();
        Vector3 targetP = 0.001 * geometry->barycenter(adjFace) + 0.999 * vPos;

        size_t hitCount = 0; 
        for(Vector3 sampleP : pts) {

            Vector3 vec = targetP - sampleP;
            RayHitResult result = tracer.trace(sampleP, vec);

            if(result.hit && result.face == adjFace) {
                hitCount++;
            }
        }

        double hitFrac = (double)hitCount / nPts;
        viewFrac[v] = hitFrac;
    }

    // If head on, dot with normal with ray to camera (makes function smooth)
    if(headOn) {

        for(VertexPtr v : geometry->getMesh()->vertices()) {

            Vector3 normal = geometry->normal(v);
            Vector3 rayToCamera = unit(headOnCameraPos - geometry->position(v));
            viewFrac[v] *= dot(normal, rayToCamera);
        }
    }

    return viewFrac;
}