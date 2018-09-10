// Author: Rohan Sawhney
// Fall 2016

#include "distortion.h"

Distortion::Distortion(Parameterization *param_):
param(param_),
geometry(param->geometry),
mesh(param->mesh),
uvs(param->uvs)
{
    
}

Vector3 hsv(double h, double s, double v)
{
    double r = 0, g = 0, b = 0;
    
    if (s == 0) {
        r = v;
        g = v;
        b = v;
        
    } else {
        h = (h == 1 ? 0 : h) * 6;
        
        int i = (int)floor(h);
        
        double f = h - i;
        double p = v*(1 - s);
        double q = v*(1 - (s*f));
        double t = v*(1 - s*(1 - f));
        
        switch (i) {
            case 0:
                r = v;
                g = t;
                b = p;
                break;
                
            case 1:
                r = q;
                g = v;
                b = p;
                break;
                
            case 2:
                r = p;
                g = v;
                b = t;
                break;
                
            case 3:
                r = p;
                g = q;
                b = v;
                break;
                
            case 4:
                r = t;
                g = p;
                b = v;
                break;
                
            case 5:
                r = v;
                g = p;
                b = q;
                break;
                
            default:
                break;
        }
    }
    
    return Vector3{r, g, b};
}

Vector3 color(double distortion)
{
    // Clamp to range [1, 1.5]
    distortion = max(1.0, min(1.5, distortion));
    
    // Compute color
    return hsv((2.0 - 4.0*(distortion-1.0))/3.0, 0.7, 0.65);
}

double computeQcDistortion(vector<Vector3>& p, vector<Vector3>& q)
{
    // Compute edge vectors
    Vector3 u1 = p[1] - p[0];
    Vector3 u2 = p[2] - p[0];
    
    Vector3 v1 = q[1] - q[0];
    Vector3 v2 = q[2] - q[0];
    
    // Compute orthonormal bases
    Vector3 e1 = u1; e1.normalize();
    Vector3 e2 = (u2 - dot(u2, e1)*e1); e2.normalize();
    
    Vector3 f1 = v1; f1.normalize();
    Vector3 f2 = (v2 - dot(v2, f1)*f1); f2.normalize();
    
    // Project onto bases
    p[0] = Vector3{0, 0, 0};
    p[1] = Vector3{dot(u1, e1), dot(u1, e2), 0};
    p[2] = Vector3{dot(u2, e1), dot(u2, e2), 0};
    
    q[0] = Vector3{0, 0, 0};
    q[1] = Vector3{dot(v1, f1), dot(v1, f2), 0};
    q[2] = Vector3{dot(v2, f1), dot(v2, f2), 0};
    
    double A = 2.0*norm(cross(u1, u2));
    
    Vector3 Ss = (q[0]*(p[1].y - p[2].y) + q[1]*(p[2].y - p[0].y) + q[2]*(p[0].y - p[1].y)) / A;
    Vector3 St = (q[0]*(p[2].x - p[1].x) + q[1]*(p[0].x - p[2].x) + q[2]*(p[1].x - p[0].x)) / A;
    double a = dot(Ss, Ss);
    double b = dot(Ss, St);
    double c = dot(St, St);
    double det = sqrt(sqr(a-c) + 4.0*b*b);
    double Gamma = sqrt(0.5*(a + c + det));
    double gamma = sqrt(0.5*(a + c - det));
    
    if (Gamma < gamma) swap(Gamma, gamma);
    
    return Gamma/gamma;
}

double Distortion::computeQcDistortion()
{
    double totArea = 0.0, sum = 0.0;
    for (FacePtr f: mesh->faces()) {
        vector<Vector3> p, q;
        for (CornerPtr c: f.adjacentCorners()) {
            p.push_back(geometry->position(c.vertex()));
            q.push_back(Vector3{uvs[c].x, uvs[c].y, 0});
        }
        
        double distortion = ::computeQcDistortion(p, q);
        double area = geometry->area(f);
        sum += area*distortion;
        totArea += area;
        param->qcDistortion[f] = color(distortion);
    }
    
    return sum/totArea;
}

double computeAreaDistortion(const vector<Vector3>& p, const vector<Vector3>& q)
{
    // Compute edge vectors and areas
    Vector3 u1 = p[1] - p[0];
    Vector3 u2 = p[2] - p[0];
    double Area = norm(cross(u1, u2));
    
    Vector3 v1 = q[1] - q[0];
    Vector3 v2 = q[2] - q[0];
    double area = norm(cross(v1, v2));
            
    return log(sqrt(Area/area));
}

double Distortion::computeAreaDistortion()
{
    FaceData<double> distortion(mesh);
    for (FacePtr f: mesh->faces()) {
        vector<Vector3> p, q;
        for (CornerPtr c: f.adjacentCorners()) {
            p.push_back(geometry->position(c.vertex()));
            q.push_back(Vector3{uvs[c].x, uvs[c].y, 0});
        }
        
        distortion[f] = ::computeAreaDistortion(p, q);
        param->areaDistortion[f] = color(distortion[f]);
    }
    
    // Average log scale factors from faces to vertices
    VertexData<size_t> vIndices = mesh->getVertexIndices();
    GC::DenseMatrix<double> u(mesh->nVertices()); u.zero();
    for (VertexPtr v: mesh->vertices()) {
        for (FacePtr f: v.adjacentFaces()) u(vIndices[v]) += distortion[f];
        u(vIndices[v]) /= v.degree();
    }
    
    // Remove shift
    u.removeMean();
    
    // Return Dirichlet energy
    GC::SparseMatrix<double> L = cotanMatrix<double>(geometry, vIndices);
    GC::DenseMatrix<double> uT = u.transpose();
    return 0.25*(uT*(L*u))(0);
}

double computeAngleDistortion(const vector<Vector3>& p, const vector<Vector3>& q)
{
    // Compute edge vectors and angles
    Vector3 u1 = p[1] - p[0];
    Vector3 u2 = p[2] - p[0];
    double Ang = angle(u1, u2);
    
    Vector3 v1 = q[1] - q[0];
    Vector3 v2 = q[2] - q[0];
    double ang = angle(v1, v2);
    
    if (Ang < ang) swap(Ang, ang);
    
    return Ang / ang;
}

double Distortion::computeAngleDistortion()
{
    double sum = 0.0;
    for (CornerPtr c: mesh->corners()) {
        vector<Vector3> p = {geometry->position(c.vertex()),
            geometry->position(c.next().vertex()),
            geometry->position(c.prev().vertex())};
        
        const Vector2& uv1 = uvs[c];
        const Vector2& uv2 = uvs[c.next()];
        const Vector2& uv3 = uvs[c.prev()];
        vector<Vector3> q = {Vector3{uv1.x, uv1.y, 0},
                             Vector3{uv2.x, uv2.y, 0},
                             Vector3{uv3.x, uv3.y, 0}};
        
        double distortion = ::computeAngleDistortion(p, q);
        sum += distortion;
        param->angleDistortion[c] = color(distortion);
    }
    
    return sum/mesh->nCorners();
}

double computeAngleSumDistortion(const vector<Vector3>& p, const vector<Vector3>& q)
{
    // Compute edge vectors and angles
    Vector3 u1 = p[0] - p[2];
    Vector3 u2 = p[1] - p[2];
    double Ang = angle(u1, u2);
    if (p.size() > 3) {
        u1 = p[0] - p[3];
        u2 = p[1] - p[3];
        Ang += angle(u1, u2);
    }
    
    Vector3 v1 = q[0] - q[2];
    Vector3 v2 = q[1] - q[2];
    double ang = angle(v1, v2);
    if (q.size() > 3) {
        v1 = q[0] - q[3];
        v2 = q[1] - q[3];
        ang += angle(v1, v2);
    }
    
    if (Ang < ang) swap(Ang, ang);
    
    return Ang / ang;
}

double Distortion::computeAngleSumDistortion()
{
    double sum = 0.0;
    for (EdgePtr e: mesh->edges()) {
        if (!e.isBoundary() && !e.isCut()) {
            HalfedgePtr h = e.halfedge();
            CornerPtr c = h.next().corner();
            CornerPtr cOp = h.twin().corner();
            vector<Vector3> p = {geometry->position(c.vertex()),
                                 geometry->position(c.next().vertex()),
                                 geometry->position(c.prev().vertex()),
                                 geometry->position(cOp.vertex())};
            
            const Vector2& uv1 = uvs[c];
            const Vector2& uv2 = uvs[c.next()];
            const Vector2& uv3 = uvs[c.prev()];
            const Vector2& uv4 = uvs[cOp];
            vector<Vector3> q = {Vector3{uv1.x, uv1.y, 0},
                                 Vector3{uv2.x, uv2.y, 0},
                                 Vector3{uv3.x, uv3.y, 0},
                                 Vector3{uv4.x, uv4.y, 0}};
            
            double distortion = ::computeAngleSumDistortion(p, q);
            sum += distortion;
            param->angleSumDistortion[e] = color(distortion);
        }
    }
    
    int nB = mesh->nBoundaryLoops() == 0 && mesh->hasCut() ? 1 : (int)mesh->nBoundaryLoops();
    for (int b = 0; b < nB; b++) {
        for (HalfedgePtr h: mesh->cutBoundary(b)) {
            CornerPtr c = h.twin().next().corner();
            vector<Vector3> p = {geometry->position(c.vertex()),
                                 geometry->position(c.next().vertex()),
                                 geometry->position(c.prev().vertex())};
            
            const Vector2& uv1 = uvs[c];
            const Vector2& uv2 = uvs[c.next()];
            const Vector2& uv3 = uvs[c.prev()];
            vector<Vector3> q = {Vector3{uv1.x, uv1.y, 0},
                                 Vector3{uv2.x, uv2.y, 0},
                                 Vector3{uv3.x, uv3.y, 0}};
            
            double distortion = ::computeAngleSumDistortion(p, q);
            sum += distortion;
            param->angleSumDistortion[h.edge()] = color(distortion);
        }
    }
    
    return sum/mesh->nEdges();
}

double computeClDistortion(const vector<Vector3>& p, const vector<Vector3>& q)
{
    // Compute edge vectors and cross lengths
    double uim = norm(p[1] - p[0]);
    double umj = norm(p[2] - p[1]);
    double ujk = norm(p[3] - p[2]);
    double uki = norm(p[0] - p[3]);
    double C = (uim*ujk) / (umj*uki);
    
    double vim = norm(q[1] - q[0]);
    double vmj = norm(q[2] - q[1]);
    double vjk = norm(q[3] - q[2]);
    double vki = norm(q[0] - q[3]);
    double c = (vim*vjk) / (vmj*vki);
    
    if (C < c) swap(C, c);
    
    return sqrt(C / c);
}

double Distortion::computeClDistortion()
{
    int count = 0;
    double sum = 0.0;
    for (EdgePtr e: mesh->edges()) {
        if (!e.isBoundary() && !e.isCut()) {
            HalfedgePtr h = e.halfedge();
            CornerPtr c = h.prev().corner();
            CornerPtr cOp = h.twin().prev().corner();
            vector<Vector3> p = {geometry->position(c.vertex()),
                                 geometry->position(c.next().vertex()),
                                 geometry->position(cOp.vertex()),
                                 geometry->position(cOp.next().vertex())};
            
            const Vector2& uv1 = uvs[c];
            const Vector2& uv2 = uvs[c.next()];
            const Vector2& uv3 = uvs[cOp];
            const Vector2& uv4 = uvs[cOp.next()];
            vector<Vector3> q = {Vector3{uv1.x, uv1.y, 0},
                                 Vector3{uv2.x, uv2.y, 0},
                                 Vector3{uv3.x, uv3.y, 0},
                                 Vector3{uv4.x, uv4.y, 0}};
            
            double distortion = ::computeClDistortion(p, q);
            sum += distortion;
            param->clDistortion[e] = color(distortion);
            count++;
            
        } else {
            param->clDistortion[e] = Vector3{1, 1, 1};
        }
    }
    
    return sum/count;
}
