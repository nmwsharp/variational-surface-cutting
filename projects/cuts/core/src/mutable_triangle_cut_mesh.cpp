#include "mutable_triangle_cut_mesh.h"

#include "polygon_soup_mesh.h"

#include <stack>


MutableCutMesh::MutableCutMesh(Geometry<Euclidean>* geometry)
{
    HalfedgeMesh* mesh = geometry->getMesh();

    // Create vertex objects
    VertexData<MutableCutVertex*>  vertexPtrs(mesh);
    for(VertexPtr v : mesh->vertices()) {
        MutableCutVertex* mV = new MutableCutVertex(geometry->position(v));
        vertices.insert(mV);
        vertexPtrs[v] = mV;
        mV->isBoundary = v.isBoundary();
    }


    // Create triangle objects
    FaceData<MutableCutFace*>  facePtrs(mesh);
    for(FacePtr f : mesh->faces()) {
        MutableCutFace *mT = new MutableCutFace();
        faces.insert(mT);
        facePtrs[f] = mT;
    }

    // Create halfedge objects
    HalfedgeData<MutableCutHalfedge*> halfedgePtrs(mesh);
    for(HalfedgePtr he : mesh->halfedges()) {
        MutableCutHalfedge* mHe = new MutableCutHalfedge();
        halfedges.insert(mHe);
        halfedgePtrs[he] = mHe;
        mHe->isCut = he.edge().isCut();
    }
    for(HalfedgePtr he : mesh->imaginaryHalfedges()) {
        MutableCutHalfedge* mHe = new MutableCutHalfedge();
        halfedges.insert(mHe);
        halfedgePtrs[he] = mHe;
        mHe->isCut = false;
    }

    // Connect all the things
    for(HalfedgePtr he : mesh->halfedges()) {
        halfedgePtrs[he]->twin = halfedgePtrs[he.twin()];
        halfedgePtrs[he]->next = halfedgePtrs[he.next()];
        halfedgePtrs[he]->face = facePtrs[he.face()];
        halfedgePtrs[he]->vertex = vertexPtrs[he.vertex()];
    }
    for(HalfedgePtr he : mesh->imaginaryHalfedges()) {
        halfedgePtrs[he]->twin = halfedgePtrs[he.twin()];
        halfedgePtrs[he]->next = halfedgePtrs[he.next()];
        halfedgePtrs[he]->vertex = vertexPtrs[he.vertex()];
    }
    for(FacePtr f : mesh->faces()) {
        facePtrs[f]->halfedge = halfedgePtrs[f.halfedge()];
    }
    for(VertexPtr v : mesh->vertices()) {
        vertexPtrs[v]->halfedge = halfedgePtrs[v.halfedge()];
    }

}        

MutableCutMesh::~MutableCutMesh() {
    for(MutableCutVertex* v : vertices) {
        delete v;
    }    

    for(MutableCutFace* f : faces) {
        delete f;
    }
    
    for(MutableCutHalfedge* he : halfedges) {
        delete he;
    }
}
        
void MutableCutMesh::indexVertices() {
    size_t i = 0;
    for(MutableCutVertex* v : vertices) {
        v->index = i++;
    }
}
        
        
bool MutableCutMesh::tryFlipEdge(MutableCutHalfedge* ha1) {

    if(ha1->isOnBoundary()) return false; // don't flip boundary edges
    if(ha1->isCut) return false; // don't flip cut edges


    // Get halfedges of first face
    MutableCutHalfedge* ha2 = ha1->next;
    MutableCutHalfedge* ha3 = ha2->next;
    
    // Get halfedges of second face
    MutableCutHalfedge* hb1 = ha1->twin;
    MutableCutHalfedge* hb2 = hb1->next;
    MutableCutHalfedge* hb3 = hb2->next;
    
    // Get vertices and faces
    MutableCutVertex* va = ha1->vertex;
    MutableCutVertex* vb = hb1->vertex;
    MutableCutVertex* vc = ha3->vertex;
    MutableCutVertex* vd = hb3->vertex;
    MutableCutFace* fa = ha1->face;
    MutableCutFace* fb = hb1->face;

    // Abort if any vertices in the diamond have been touched already
    // if(va->touched || vb->touched || vc->touched || vd->touched) {
    //     return false;
    // }
    
    if(va->degree() <= 3 || vb->degree() <= 3) return false; // can't flip a degree-3 vertex

    // Make sure the face isn't nonconvex (or nearly so)
    double angleThresh = 2.79; // 160 degrees
    double angle1 = angle(vd->position - va->position, vc->position - va->position);
    double angle2 = angle(vc->position - vb->position, vd->position - vb->position);
    if(angle1 > angleThresh || angle2 > angleThresh) {
        return false;
    }


    { // Check how much normals would be changed by the flip
      // TODO seems like a great use for error quadrics
        Vector3 nAOrig = cross(vb->position - va->position, vc->position - va->position);
        Vector3 nBOrig = cross(vd->position - va->position, vb->position - va->position);
        Vector3 nANew = cross(vd->position - vc->position, vb->position - vc->position);
        Vector3 nBNew = cross(va->position - vc->position, vd->position - vc->position);

        double maxNormalChange = std::max({
            angle(nAOrig, nANew),
            angle(nAOrig, nBNew),
            angle(nBOrig, nANew),
            angle(nBOrig, nBNew)
        });

        if(maxNormalChange > normalThreshold || !nANew.isFinite() || !nBNew.isFinite()) {
            return false;
        }

        // Make sure we don't create any excpetionally skinny faces
        double meanArea = 0.5 *(norm(nAOrig) + norm(nBOrig));
        if(norm(nANew) < skinnyThresh * meanArea || norm(nBNew) < skinnyThresh * meanArea) { 
            return false;
        }

    }

    {   // Verify that this will not create a nonmanifold edges
        // (vc and vd should not have any shared neighbors other than va and vb)

        // Build map of neighbors of va
        std::unordered_set<MutableCutVertex*> vcNeighbors;
        MutableCutHalfedge* currHe = vc->halfedge;
        do {
            vcNeighbors.insert(currHe->twin->vertex);
            currHe = currHe->twin->next;
        } while (currHe != vc->halfedge);
        
        // Check neighbors of vd
        currHe = vd->halfedge;
        do {
            MutableCutVertex* testVert = currHe->twin->vertex;
            if(testVert != va && testVert != vb && vcNeighbors.find(testVert) != vcNeighbors.end()) {
                return false;
            }
            currHe = currHe->twin->next;
        } while (currHe != vd->halfedge);

    } 

    // Update vertex pointers
    if( va->halfedge == ha1 ) va->halfedge = hb2;
    if( vb->halfedge == hb1 ) vb->halfedge = ha2;
    
    // Update face pointers
    fa->halfedge = ha1;
    fb->halfedge = hb1;
    
    // Update halfedge pointers
    ha1->next = hb3;
    hb3->next = ha2;
    ha2->next = ha1;
    hb1->next = ha3;
    ha3->next = hb2;
    hb2->next = hb1;
    ha1->vertex = vc;
    hb1->vertex = vd;
    ha3->face = fb;
    hb3->face = fa;

    va->touched = true;
    vb->touched = true;
    vc->touched = true;
    vd->touched = true;

    return true;
}
        
bool MutableCutMesh::tryCollapseEdge(MutableCutHalfedge* ha1) {

    if(ha1->isOnBoundary()) {
        // cout << "Can't collapse, boundary" << endl;
        return false; // don't collapse boundary edges
    }
    MutableCutHalfedge* ha2 = ha1->next;
    MutableCutHalfedge* ha3 = ha2->next;
    
    // Get halfedges of second face
    MutableCutHalfedge* hb1 = ha1->twin;
    MutableCutHalfedge* hb2 = hb1->next;
    MutableCutHalfedge* hb3 = hb2->next;
    
    // Get vertices and faces
    MutableCutVertex* va = ha1->vertex;
    MutableCutVertex* vb = hb1->vertex;
    MutableCutVertex* vc = ha3->vertex;
    MutableCutVertex* vd = hb3->vertex;
    MutableCutFace* fa = ha1->face;
    MutableCutFace* fb = hb1->face;

    // Get neighbor halfedges 
    MutableCutHalfedge* hnac = ha3->twin;
    MutableCutHalfedge* hncb = ha2->twin;
    MutableCutHalfedge* hnbd = hb3->twin;
    MutableCutHalfedge* hnda = hb2->twin;
    
    // Abort if any vertices in the diamond have been touched already
    // if(va->touched || vb->touched || vc->touched || vd->touched) {
    //     return false;
    // }
    
    if(vc->degree() <= 3 || vd->degree() <= 3) return false; // can't collapse around a degree-3 vertex
  
    // Don't change topology of cut network
    bool vaIsCut = va->isOnCut();
    bool vbIsCut = vb->isOnCut();
    if(!ha1->isCut && (vaIsCut && vbIsCut)) {
        // cout << "Can't collapse, cut topology" << endl;
        return false;
    }

   
    // If one of the vertices is a cut vertex, don't move it
    Vector3 newAPos; 
    if(vaIsCut && !vbIsCut) {
        newAPos = va->position;
    } else if (vbIsCut && !vaIsCut) {
        newAPos = vb->position;
    } else {
        newAPos = splineMidpoint(va->position, va->normal(), vb->position, vb->normal()); 
    }
    
    {   // Make sure normals would not change too much, and that we don't create an exceptionall skinny triangle
        // (note that the latter check is mainly to circumvent a particluar combination of spilts and collapses
        // which results in a triangle with exactly zero area)
        // TODO seems like a great use for error quadrics
      
        // Loop over all faces adjacent to va and check how much their normals will change
        MutableCutHalfedge* currHe = va->halfedge;
        do {
            MutableCutFace* f = currHe->face;
            if(f != nullptr) {

                Vector3 initNormal = f->areaNormal();
                Vector3 oldPos = va->position;
                va->position = newAPos;
                Vector3 newNormal = f->areaNormal();
                va->position = oldPos;

                if( !newNormal.isFinite() || 
                    angle(initNormal, newNormal) > normalThreshold ||
                    norm(newNormal) < skinnyThresh * norm(initNormal)
                    ) {
                        // cout << "Can't collapse, normal change" << endl;
                        return false;
                    }

            }
            currHe = currHe->twin->next;
        } while (currHe != va->halfedge);
        
        // Loop over all faces adjacent to vb and check how much their normals will change
        currHe = vb->halfedge;
        do {
            MutableCutFace* f = currHe->face;
            if(f != nullptr) {

                Vector3 initNormal = f->normal();
                Vector3 oldPos = vb->position;
                vb->position = newAPos;
                Vector3 newNormal = f->normal();
                vb->position = oldPos;

                if( !newNormal.isFinite() || 
                    angle(initNormal, newNormal) > normalThreshold ||
                    norm(newNormal) < skinnyThresh * norm(initNormal)
                    ) {
                        // cout << "Can't collapse, normal change" << endl;
                        return false;
                    }

            }
            currHe = currHe->twin->next;
        } while (currHe != vb->halfedge);


    }


    {   // Verify that this will not create a nonmanifold edges
        // (va and vb should not have any shared neighbors other than vc and vd)

        // Build map of neighbors of va
        std::unordered_set<MutableCutVertex*> vaNeighbors;
        MutableCutHalfedge* currHe = va->halfedge;
        do {
            vaNeighbors.insert(currHe->twin->vertex);
            currHe = currHe->twin->next;
        } while (currHe != va->halfedge);
        
        // Check neighbors of vb
        currHe = vb->halfedge;
        do {
            MutableCutVertex* testVert = currHe->twin->vertex;
            if(testVert != vc && testVert != vd && vaNeighbors.find(testVert) != vaNeighbors.end()) {
                // cout << "Can't collapse, manifold" << endl;
                return false;
            }
            currHe = currHe->twin->next;
        } while (currHe != vb->halfedge);

    } 



    // Find new vertex halfedges if needed
    while(va->halfedge == ha1 || va->halfedge == hb2) va->halfedge = va->halfedge->twin->next;
    while(vb->halfedge == ha2 || vb->halfedge == hb1) vb->halfedge = vb->halfedge->twin->next;
    while(vc->halfedge == ha3) vc->halfedge = vc->halfedge->twin->next;
    while(vd->halfedge == hb3) vd->halfedge = vd->halfedge->twin->next;


    // Find new halfedge vertices around vB
    MutableCutHalfedge* currHe = vb->halfedge;
    do {
        currHe->vertex = va;
        currHe = currHe->twin->next;
    } while (currHe != vb->halfedge);
    
    // Fix twin pointers
    hnac->twin = hncb;
    hncb->twin = hnac;
    hnda->twin = hnbd;
    hnbd->twin = hnda;
   
    // Copy cut state
    if(hnac->isCut) hnac->twin->isCut = true;
    if(hncb->isCut) hncb->twin->isCut = true;
    if(hnda->isCut) hnda->twin->isCut = true;
    if(hnbd->isCut) hnbd->twin->isCut = true;

    // Update vA to be a mix of vA and vB
    va->position = newAPos;
    va->desiredVertexDensity = 0.5 * (va->desiredVertexDensity + vb->desiredVertexDensity); 

    // Update degree
    if(vc->isBoundary) va->isBoundary = true;
    
    va->touched = true;
    vc->touched = true;
    vd->touched = true;

    // Out with the old
    halfedges.erase(ha1); delete ha1;
    halfedges.erase(ha2); delete ha2;
    halfedges.erase(ha3); delete ha3;
    halfedges.erase(hb1); delete hb1;
    halfedges.erase(hb2); delete hb2;
    halfedges.erase(hb3); delete hb3;
    vertices.erase(vb); delete vb;
    faces.erase(fa); delete fa;
    faces.erase(fb); delete fb;
    

    return true;
}

// bool MutableCutMesh::trySplitEdge(MutableCutHalfedge* ha1) {
    
//     MutableCutHalfedge* ha2 = ha1->next;
//     MutableCutHalfedge* ha3 = ha2->next;
    
//     // Get halfedges of second face
//     MutableCutHalfedge* hb1 = ha1->twin;
//     MutableCutHalfedge* hb2 = hb1->next;
//     MutableCutHalfedge* hb3 = hb2->next;
    
//     // Get vertices and faces
//     MutableCutVertex* va = ha1->vertex;
//     MutableCutVertex* vb = hb1->vertex;
//     MutableCutVertex* vc = ha3->vertex;
//     MutableCutVertex* vd = hb3->vertex;
//     MutableCutFace* fa = ha1->face;
//     MutableCutFace* fb = hb1->face;

//     // Get neighbor halfedges 
//     MutableCutHalfedge* hnac = ha3->twin;
//     MutableCutHalfedge* hncb = ha2->twin;
//     MutableCutHalfedge* hnbd = hb3->twin;
//     MutableCutHalfedge* hnda = hb2->twin;
    
//     // Abort if any vertices in the diamond have been touched already
//     if(va->touched || vb->touched || vc->touched || vd->touched) {
//         return false;
//     }

//     // Add a new vertex at the edge midpoint
//     Vector3 newPos = splineMidpoint(va->position, va->normal(), vb->position, vb->normal()); 
//     MutableCutVertex* vn = new MutableCutVertex(newPos); vertices.insert(vn);
//     vn->desiredVertexDensity = 0.5 * (va->desiredVertexDensity + vb->desiredVertexDensity);
//     vn->touched = true;

//     // Create the new faces
//     MutableCutFace* fan = new MutableCutFace(); faces.insert(fan);
//     MutableCutFace* fbn = new MutableCutFace(); faces.insert(fbn);

//     // Create the new halfedges
//     MutableCutHalfedge* han1 = new MutableCutHalfedge(); halfedges.insert(han1);
//     MutableCutHalfedge* han2 = new MutableCutHalfedge(); halfedges.insert(han2);
//     MutableCutHalfedge* han3 = new MutableCutHalfedge(); halfedges.insert(han3);
//     MutableCutHalfedge* hbn1 = new MutableCutHalfedge(); halfedges.insert(hbn1);
//     MutableCutHalfedge* hbn2 = new MutableCutHalfedge(); halfedges.insert(hbn2);
//     MutableCutHalfedge* hbn3 = new MutableCutHalfedge(); halfedges.insert(hbn3);

//     // Set halfedge->face
//     han1->face = fan; han2->face = fan; han3->face = fan; 
//     hbn1->face = fbn; hbn2->face = fbn; hbn3->face = fbn; 

//     // Set face->halfedge
//     fan->halfedge = han1;
//     fbn->halfedge = hbn1;

//     // Set halfedge->next
//     han1->next = han2; han2->next = han3; han3->next = han1;
//     hbn1->next = hbn2; hbn2->next = hbn3; hbn3->next = hbn1;

//     // Set halfedge->vertex
//     han1->vertex = va;
//     han2->vertex = vn;
//     han3->vertex = vc;
//     hbn1->vertex = vn;
//     hbn2->vertex = va;
//     hbn3->vertex = vd;
//     ha1->vertex = vn;
//     hb2->vertex = vn;

//     // Set halfedge->twin
//     han2->twin = ha3; ha3->twin = han2;
//     hbn3->twin = hb2; hb2->twin = hbn3;
//     hnac->twin = han3; han3->twin = hnac;
//     han1->twin = hbn1; hbn1->twin = han1;
//     hnda->twin = hbn2; hbn2->twin = hnda;

//     // Set vertex->halfedge
//     vn->halfedge = ha1;
//     va->halfedge = hbn2;
//     vb->halfedge = hb1;
//     vc->halfedge = ha3;
//     vd->halfedge = hbn3;

//     // Update degree
//     vn->degree = 4;
//     vc->degree++;
//     vd->degree++;
    
//     va->touched = true;
//     vb->touched = true;
//     vc->touched = true;
//     vd->touched = true;

//     return true;
// }

bool MutableCutMesh::tryMoveVertex(MutableCutVertex* v, Vector3 newPos) {

    // Loop over all faces adjacent to v and check how much their normals will change
    MutableCutHalfedge* currHe = v->halfedge;
    do {
        MutableCutFace* f = currHe->face;
        if(f != nullptr) {

            Vector3 initNormal = f->areaNormal();
            Vector3 oldPos = v->position;
            v->position = newPos;
            Vector3 newNormal = f->areaNormal();
            v->position = oldPos;

            if( !newNormal.isFinite() || 
                angle(initNormal, newNormal) > normalThreshold ||
                norm(newNormal) < skinnyThresh * norm(initNormal)
                ) return false;

        }
        currHe = currHe->twin->next;
    } while (currHe != v->halfedge);

    v->position = newPos;

    return true;
}


Geometry<Euclidean>* MutableCutMesh::toHalfedgeMesh() {

    indexVertices();
    size_t nVert = vertices.size();

    // Build vertex position vector 
    std::vector<Vector3> vPos(nVert);
    int minDegree = std::numeric_limits<int>::max();
    int maxDegree = 0;
    for(MutableCutVertex* v : vertices) {
        vPos[v->index] = v->position;
        minDegree = std::min(minDegree, v->degree());
        maxDegree = std::max(maxDegree, v->degree());
    }

    cout << "Min vertex degree = " << minDegree << "  max vertex degree = " << maxDegree << endl;

    // Convert to polygon soup
    std::vector<std::vector<size_t>> pInds;
    for(MutableCutFace* t : faces) {
        MutableCutHalfedge* he = t->halfedge;
        MutableCutVertex* vA = he->vertex;
        MutableCutVertex* vB = he->next->vertex;
        MutableCutVertex* vC = he->next->next->vertex;
        pInds.push_back({   vA->index,
                            vB->index,
                            vC->index });
    }
    PolygonSoupMesh pSoup(pInds, vPos);

    Geometry<Euclidean>* geom;
    HalfedgeMesh* newMesh = new HalfedgeMesh(pSoup, geom);

    // Copy cuts
    size_t i = 0;
    for(MutableCutFace* mutF : faces) {

        FacePtr f = newMesh->face(i);
        HalfedgePtr he0 = f.halfedge();
        HalfedgePtr he1 = f.halfedge().next();
        HalfedgePtr he2 = f.halfedge().next().next();

        MutableCutHalfedge* mhe0 = mutF->halfedge;
        MutableCutHalfedge* mhe1 = mhe0->next;
        MutableCutHalfedge* mhe2 = mhe1->next;

        if(mhe0->isCut) he0.edge().markCut(true);
        if(mhe1->isCut) he1.edge().markCut(true);
        if(mhe2->isCut) he2.edge().markCut(true);

        i++;
    }

    return geom;
}

void MutableCutMesh::setVerticesUntouched() {
    for(MutableCutVertex* v : vertices) {
        v->touched = false;
    }
}

MutableCutVertex::MutableCutVertex(Vector3 p_) 
    : position(p_)
{}

bool MutableCutVertex::isOnCut() {
    MutableCutHalfedge* currHe = halfedge;
    MutableCutHalfedge* firstHe = halfedge;
    do {

        if(currHe->isCut) return true;

        currHe = currHe->twin->next;
    } while (currHe != firstHe);

    return false;
}

int MutableCutVertex::degree() {
    int d = 0;
    MutableCutHalfedge* currHe = halfedge;
    MutableCutHalfedge* firstHe = halfedge;
    do {
        d++;
        currHe = currHe->twin->next;
    } while (currHe != firstHe);

    return d;
}

bool MutableCutHalfedge::isOnBoundary() {
    return (!isReal()) || (!twin->isReal());
}

bool MutableCutHalfedge::isReal() {
    return face != nullptr; 
}

Vector3 MutableCutFace::normal() {
    return unit(areaNormal());
}

Vector3 MutableCutFace::areaNormal() {

    Vector3 p0 = halfedge->vertex->position;
    Vector3 p1 = halfedge->next->vertex->position;
    Vector3 p2 = halfedge->next->next->vertex->position;

    Vector3 vA = p1 - p0;
    Vector3 vB = p2 - p0;

    return cross(vA, vB) / 2.0;
}

Vector3 MutableCutFace::barycenter() {

    Vector3 p0 = halfedge->vertex->position;
    Vector3 p1 = halfedge->next->vertex->position;
    Vector3 p2 = halfedge->next->next->vertex->position;

    return (p0 + p1 + p2) / 3.0;
}

Vector3 MutableCutFace::circumcenter() {

    Vector3 a = halfedge->vertex->position;
    Vector3 b = halfedge->next->vertex->position;
    Vector3 c = halfedge->next->next->vertex->position;

    return a + (norm2(c - a)*cross(cross(b - a, c - a), b - a) + 
        norm2(b - a)*cross(c - a, cross(b - a, c - a)))/(2.*norm2(
        cross(b - a, c - a)));

}

double MutableCutFace::desiredDensity() {

    double s0 = halfedge->vertex->desiredVertexDensity;
    double s1 = halfedge->next->vertex->desiredVertexDensity;
    double s2 = halfedge->next->next->vertex->desiredVertexDensity;

    return (s0 + s1 + s2) / 3.0;
}

Vector3 MutableCutVertex::normal() {

    // TODO maybe want angle-weighted for behavior around sharp features

    Vector3 normalAvg{0.0, 0.0, 0.0};
    MutableCutHalfedge* currHe = halfedge;
    do {
        MutableCutFace* f = currHe->face;
        if(f != nullptr) {
            Vector3 areaNormal = f->areaNormal();
            normalAvg += areaNormal;
        }
        currHe = currHe->twin->next;
    } while (currHe != halfedge);

    return unit(normalAvg);
}

Vector3 MutableCutHalfedge::vector() {
    Vector3 p0 = vertex->position;
    Vector3 p1 = next->vertex->position;
    return p1 - p0;
}

double MutableCutHalfedge::oppositeAngle() {

    Vector3 pRoot = next->next->vertex->position;
    Vector3 pA = vertex->position;
    Vector3 pB = next->vertex->position;

    return angle(pA - pRoot, pB - pRoot); 
}

Vector3 MutableCutMesh::splineMidpoint(Vector3 p1, Vector3 n1, Vector3 p2, Vector3 n2) {

    Vector3 basisX = unit(p2 - p1);
    Vector3 sharedNormal = unit(n1 + n2);
    Vector3 basisY = unit(sharedNormal - basisX * dot(sharedNormal, basisX));
    Vector3 midpoint = 0.5 * (p1 + p2);

    double slope1 = dot(basisX, -n1) / dot(basisY, n1);
    double slope2 = dot(basisX, -n2) / dot(basisY, n2);
    
    // If both of the normals point in the same direction in this shapred basis, something crazy is
    // going on with the surface here. The spline would try to adapt, but we're probably better
    // off just returning the midpoint.
    if(slope1 * slope2 > 0.0) {
        return midpoint;
    }

    // Cubic spline equation at midpoint
    double yVal = 0.125 * slope1 - 0.125 * slope2;
    yVal = clamp(yVal, -2.0, 2.0); // if normals differ significantly, spline will put the point very far away
    Vector3 pos = midpoint + yVal * basisY * norm(p2-p1);
    
    
    // Point-slope formulation can yield degenerate values for ridiculous geometries
    if(!pos.isFinite()) {
        return midpoint;
    }

    return pos;
}
