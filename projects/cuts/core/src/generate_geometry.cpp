#include "generate_geometry.h"

#include "fast_triangle_soup.h"
#include "combining_hash_functions.h"

#include <unordered_map>


// Icosahedrom mesh with symmetry along "Disdyakis triacontahedron" face canonicalization
Geometry<Euclidean>* generateSymmetricIcosahedron(size_t atLeastNPoints, SymmetryResult& symResult) {

    std::vector<Vector3> initPts = {
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
    std::vector<Vector3> pts = initPts;

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


    // Build an array of barycentric coordinates in the original faces
    // Build an array of which original face each face comes from
    std::vector<std::array<Vector3, 3>> faceBarys;
    std::vector<size_t> faceParent;
    for(size_t i = 0; i < faces.size(); i++) {
        faceBarys.push_back({{
            Vector3{1.0, 0.0, 0.0},
            Vector3{0.0, 1.0, 0.0},
            Vector3{0.0, 0.0, 1.0}
        }});
        faceParent.push_back(i);
    }



    // Subdivide until we have enough vertices
    atLeastNPoints = std::max(atLeastNPoints, (size_t)500); // need enough points for subdivision to generate meaningful symmetry
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

        // Create new faces
        std::vector<Tri> newFaces;
        std::vector<std::array<Vector3, 3>> newFaceBarys;
        std::vector<size_t> newFaceParent;
        for(size_t i = 0; i < faces.size(); i++) {

            Tri& t = faces[i];
            std::array<Vector3,3>& bary = faceBarys[i];
            size_t parent = faceParent[i];

            // 1-4 subdivision

            // Gather indices
            size_t v0 = t[0];
            size_t v1 = t[1];
            size_t v2 = t[2];
            size_t vn0 = newVertMap[std::make_pair(std::min(v0, v1), std::max(v0, v1))];
            size_t vn1 = newVertMap[std::make_pair(std::min(v1, v2), std::max(v1, v2))];
            size_t vn2 = newVertMap[std::make_pair(std::min(v2, v0), std::max(v2, v0))];

            // Gather barycentric coords
            Vector3 c0 = bary[0];
            Vector3 c1 = bary[1];
            Vector3 c2 = bary[2];
            Vector3 cn0 = (c0 + c1) / 2;
            Vector3 cn1 = (c1 + c2) / 2;
            Vector3 cn2 = (c2 + c0) / 2;

            // == Create 4 faces
            newFaces.push_back({{v0, vn0, vn2}});
            newFaceBarys.push_back({{c0, cn0, cn2}});
            newFaceParent.push_back(parent);

            newFaces.push_back({{v1, vn1, vn0}});
            newFaceBarys.push_back({{c1, cn1, cn0}});
            newFaceParent.push_back(parent);
            
            newFaces.push_back({{v2, vn2, vn1}});
            newFaceBarys.push_back({{c2, cn2, cn1}});
            newFaceParent.push_back(parent);

            newFaces.push_back({{vn0, vn1, vn2}});
            newFaceBarys.push_back({{cn0, cn1, cn2}});
            newFaceParent.push_back(parent);
        }

        faces = newFaces;
        faceBarys = newFaceBarys;
        faceParent = newFaceParent;
    }


    // Build a halfedge mesh
    std::vector<std::array<bool, 3>> noBoundaries(faces.size(), {{false, false, false}});
    HalfedgeMesh* mesh;
    Geometry<Euclidean>* geom;
    FastTriangleSoup soup(faces, pts, noBoundaries);
    soup.toHalfedgeMesh(mesh, geom);

    // Build symmetry information
    SymmetryResult sym;
    sym.symmetryFound = true; 
    sym.symmetrySet = VertexData<std::vector<VertexPtr>>(mesh);

    auto sortRevVector3 = [](Vector3 v)
    {
        std::vector<double> coords = {v.x, v.y, v.z};
        std::sort(coords.begin(), coords.end());
        return Vector3{coords[2], coords[1], coords[0]};
    };

    // Record the canonical vertices and build a lookup map
    std::vector<bool> isCanonical(pts.size(), false);
    std::vector<Vector3> canonicalCoord(pts.size());
    std::unordered_map<Vector3, size_t> canonicalIndex;
    for(size_t iFace = 0; iFace < faces.size(); iFace++) {

        for(size_t i = 0; i < 3; i++) {
            size_t iVert = faces[iFace][i];

            Vector3 coord = faceBarys[iFace][i]; 
            Vector3 coordCanonical = sortRevVector3(coord);

            // Canonical coord is the same in all faces, so we can just track one
            // and avoid duplicate loops
            canonicalCoord[iVert] = coordCanonical;

            // Check if this is a canonical vertex 
            if(faceParent[iFace] == 0)  {
                if(coord == coordCanonical) {
                    isCanonical[iVert] = true; 
                    canonicalIndex[coord] = iVert;
                }
            }
        }
    }


    // Build symmetry sets and list of canonical vertices
    for(size_t iVert = 0; iVert < pts.size(); iVert++) {
        if(isCanonical[iVert]) {

            sym.canonicalVertices.push_back(mesh->vertex(iVert));

        } else {

            Vector3 canonCoord = canonicalCoord[iVert];
            size_t canonInd = canonicalIndex[canonCoord];
            VertexPtr canonVert = mesh->vertex(canonInd);
            VertexPtr thisVert = mesh->vertex(iVert);
            sym.symmetrySet[canonVert].push_back(thisVert);

        }
    }

    symResult = sym;


    // Perturb geometry to be bumpy
    VertexData<Vector3> normals;
    geom->getVertexNormals(normals);

    double H = .8; // height factor
    double b = .93;
    double L = .28; // base sphere cap radius
    double r = 1.1; // ratio of lower ramp spheres to tip spheres
    double maxT = .629; // max value x will take

 
    std::function<double(double)> fPerturb = [&](double x) { 

        if(x < b*L) {
            return H*(std::sqrt(L*L - x*x) - std::sqrt(L*L - b*b*L*L)) + fPerturb(b*L);
        } else if(x < (1+r)*b*L) {
            return H*(r*L - std::sqrt(r*r*L*L - std::pow(r*L - (x + L - 2*b*L), 2) ));
        } else {
            return 0.0;
        }
            
    };

    for(VertexPtr v : mesh->vertices()) {

        Vector3 pos = geom->position(v);

        // Parameterize by distance from one of the 12 initial points
        double t = std::numeric_limits<double>::infinity();
        for(Vector3 p : initPts) {
            t = std::min(t, norm(p - pos));
        }
        maxT = std::max(maxT, t);

        geom->position(v) += fPerturb(t) * normals[v];

    }

    cout << "maxT = " << maxT << endl;


    return geom;
}