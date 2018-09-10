#include "target_edge_lengths.h"


#include "LBFGS.h"
#include "discrete_operators.h"

using namespace GC;

PositionSolver::PositionSolver(Geometry<Euclidean>* geometry_) :
    geometry(geometry_)
{

    mesh = geometry->getMesh();

    nVert = mesh->nVertices();
    nEdge = mesh->nEdges();
    vInd = mesh->getVertexIndices();
    eInd = mesh->getEdgeIndices();

    // Build edge tuple array
    edges.resize(nEdge);
    for(EdgePtr e : mesh->edges()) {
        edges[eInd[e]] = {{vInd[e.halfedge().vertex()], vInd[e.halfedge().twin().vertex()]}};
    }
}

GC::SparseMatrix<double> clampedDiagonalInverse(SparseMatrix<double>& A, double thresh) {
    
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

void PositionSolver::buildMatrices(const EdgeData<double>& targetEdgeLengths) {

    cout << "Building operators" << endl;


    // Build angles and area from edge lengths
    FaceData<double> faceAreas(mesh);
    HalfedgeData<double> oppAngles(mesh);
    for(FacePtr f : mesh->faces()) {
        double a = targetEdgeLengths[f.halfedge().edge()];
        double b = targetEdgeLengths[f.halfedge().next().edge()];
        double c = targetEdgeLengths[f.halfedge().next().next().edge()];
        std::array<double, 3> l = {{a, b, c}};

        // Edge lengths
        double p = (a + b + c) / 2;
        double A = std::sqrt(p * (p-a) * (p-b) * (p-c));
        // cout << a << " " << b << " " << " " << c << endl;
        if(!std::isfinite(A)) {
            cout << "non finite area" << endl;
            // throw std::runtime_error("boo");
            A = 0.0; // not sure why this is happening? is BFF always exactly isometric?
        }
        faceAreas[f] = A;

        // Angles
        std::array<double, 3> oAngles;
        for(int j = 0; j < 3; j++) {

            double myL = l[j];
            double L1 = l[(j+1)%3]; 
            double L2 = l[(j+2)%3]; 

            double cosTheta = (L1*L1 + L2*L2 - myL*myL) / (2*L1*L2);
            double theta = std::acos(cosTheta);

            if(!std::isfinite(theta)) {
                cout << "non finite theta" << endl;
                // throw std::runtime_error("boo2");
                theta = .9 * PI / 2.0; // not sure why this is happening? is BFF always exactly isometric?
            }
            
            oAngles[j] = theta;
        }
        oppAngles[f.halfedge()] = oAngles[0];
        oppAngles[f.halfedge().next()] = oAngles[1];
        oppAngles[f.halfedge().next().next()] = oAngles[2];
    }

    // Build Hodge0
    hodge0 = SparseMatrix<double>(nVert, nVert);
    for(VertexPtr v : mesh->vertices()) {
        size_t i = vInd[v];
        double area = 0;
        for(FacePtr f : v.adjacentFaces()) {
            area += faceAreas[f] / 3;
        }
        hodge0(i,i) = area;
    }

    // Build Hodge1
    hodge1 = SparseMatrix<double>(nEdge, nEdge);
    for(EdgePtr e : mesh->edges()) {
        double cotanW = 0;
        size_t i = eInd[e];
        if(e.halfedge().isReal()) {
            cotanW += 1.0 / std::tan(oppAngles[e.halfedge()]);
        }
        if(e.halfedge().next().isReal()) {
            cotanW += 1.0 / std::tan(oppAngles[e.halfedge().next()]);
        }
        hodge1(i,i) = cotanW;
    }

    // hodge0 = buildHodge0<double>(geometry);
    // hodge1 = buildHodge1<double>(geometry);
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

    hodge0Bar = hodge0;
    hodge0BarInv = hodge0Inv;
    for(size_t i  = 0; i < nVert; i++) {
        if(vCut[i]) {
            hodge0Bar(i,i) = 0.0;
            hodge0BarInv(i,i) = 0.0;
        }
    }

    B = zeroFormWeakLaplacian.transpose() * hodge0Inv * hodge0Bar * hodge0Inv * zeroFormWeakLaplacian;
    // B = zeroFormWeakLaplacian.transpose() * hodge0BarInv * zeroFormWeakLaplacian;
    BBT = B + B.transpose();

    Beigen = B.toEigen();
    BBTeigen = BBT.toEigen();


    // Edge weights
    for(EdgePtr e : mesh->edges()) {
        // double eArea = faceAreas[e.halfedge().face()] / 3.0;
        double eArea = geometry->area(e.halfedge().face()) / 3.0;
        if(!e.isBoundary()) {
            // eArea += faceAreas[e.halfedge().twin().face()] / 3.0;
            eArea += geometry->area(e.halfedge().twin().face()) / 3.0;
        }

        edgeAreas.push_back(eArea);
        // edgeAreas.push_back(1.0);
    }
}


VertexData<Vector3> PositionSolver::solve(const EdgeData<double>& targetEdgeLengths, const VertexData<char>& vertexFixed) {

    using namespace Eigen;
    using namespace LBFGSpp;

    cout << endl << "Solving for vertex positions with specified edge lengths..." << endl;

    // Initialize with current positions
    VectorXd sol = VectorXd(3*nVert);
    for(VertexPtr v : mesh->vertices()) {
        size_t i = vInd[v];
        for(size_t j = 0; j < 3; j++) sol[3*i + j] = geometry->position(v)[j];
    }

    // Store target lengths
    targetLengths = VectorXd(nEdge);
    for(EdgePtr e : mesh->edges()) {
        targetLengths[eInd[e]] = targetEdgeLengths[e];
    }
    
    // Store fixed
    vCut.resize(nVert);
    for(VertexPtr v : mesh->vertices()) {
        vCut[vInd[v]] = vertexFixed[v];
        // vCut[vInd[v]] = false;
    }

    buildMatrices(targetEdgeLengths);

    // Solve with LBFGS
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1000;
    // param.max_iterations = 0; // unlimited
    param.linesearch = 2; // WOLFE

    // Create a solver
    LBFGSSolver<double> solver(param);

    // Solve
    double energy;
    int nIterations = solver.minimize(*this, sol, energy);

    // Print info about result
    cout << endl << "  LBFGS finished in " << nIterations << " iterations" << std::endl;
    cout << "  Final energy = " << energy << endl;

    // Store current positions
    VertexData<Vector3> newPositions(mesh);
    for(VertexPtr v : mesh->vertices()) {
        size_t i = vInd[v];
        for(size_t j = 0; j < 3; j++) newPositions[v][j] = sol[3*i + j];
    }
    
    std::vector<double> edgeErrors;
    double errorSum = 0;
    double errorWeightSum = 0;
    size_t iE = 0;
    for(EdgePtr e : mesh->edges()) {
        double targetLen = targetEdgeLengths[e];
        double actualLen = norm(newPositions[e.halfedge().vertex()] - newPositions[e.halfedge().twin().vertex()]);

        double err = std::abs(1.0 - actualLen / targetLen);
        edgeErrors.push_back(err);

        errorSum += edgeAreas[iE] * err;
        errorWeightSum += edgeAreas[iE];
        iE++;
    }

    double meanError = errorSum / errorWeightSum;
    cout << "Mean edge length deviation on surface = " << meanError << endl;

    // std::sort(edgeErrors.begin(), edgeErrors.end());
    // for(double e : edgeErrors) {
    //     cout << "Edge error: " << e << endl;
    // }
    // int ind = static_cast<int>(std::floor(edgeErrors.size() * .98));
    // cout << "98 percentile error = " << edgeErrors[ind] << endl;


    return newPositions;
}

double PositionSolver::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

    using namespace Eigen;

    grad = VectorXd::Zero(3*nVert);
    double energy = 0.0;

    for(size_t iE = 0; iE < edges.size(); iE++) {
        std::array<size_t, 2> e = edges[iE];
        double eW = edgeAreas[iE] * globalScale;
        // cout << eW << endl;

        Vector3 vA = Vector3{x[3*e[0]+0], x[3*e[0]+1], x[3*e[0]+2]};
        Vector3 vB = Vector3{x[3*e[1]+0], x[3*e[1]+1], x[3*e[1]+2]};

        Vector3 vAB = vB - vA;
        double l = norm(vAB);

        double targetL = targetLengths[iE];
        Vector3 vABg = 2.0 * unit(vAB) * (l - targetL);

        for(size_t j = 0; j < 3; j++) {
            grad[3*e[0]+j] += -vABg[j] * eW;
            grad[3*e[1]+j] +=  vABg[j] * eW;
        }

        energy += eW * (l - targetL) * (l - targetL);
    }
    double edgeLengthEnergy = energy;
    cout << "Edge length energy = " << edgeLengthEnergy << endl;


    // Bending term
    for(int j = 0; j < 3; j++) {

        // GC::DenseVector<double> posVec(nVert);
        // for(size_t i = 0; i < nVert; i++) {
        //     posVec(i) = x[3*i + j];
        // }
        
        Eigen::VectorXd posVec(nVert);
        for(size_t i = 0; i < nVert; i++) {
            posVec[i] = x[3*i + j];
        }

        // Energy
        // GC::DenseVector<double> Lf = zeroFormLaplacian * posVec;
        // double bTerm = bendingWeight * (Lf.transpose() * (hodge0Bar * Lf))(0,0);
        double bTerm = bendingWeight * (posVec.transpose() * Beigen * posVec)(0,0);
        cout << "Bending energy = " << bTerm << endl;
        energy += globalScale *bTerm;

        // Gradient
        // GC::DenseVector<double> gradVec = 2.0 * bendingWeight * hodge0Bar * Lf;
        Eigen::VectorXd gradVec = bendingWeight * BBTeigen * posVec;
        for(size_t i = 0; i < nVert; i++) {
            grad[3*i+j] += globalScale * gradVec[i];
        }
    }


    // Zero gradient at fixed vertices
    // for(size_t i = 0; i < nVert; i++) {
    //     if(vCut[i]) {
    //         for(size_t j = 0; j < 3; j++) {
    //             grad[3*i+j] = 0.0;
    //         }
    //     }
    // }

    cout << endl;

    return energy;
}