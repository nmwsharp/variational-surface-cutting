#include "developable_approximation.h"

#include "discrete_operators.h"


DevelopableApproximation::DevelopableApproximation(Geometry<Euclidean>* geometry_) :
    initGeometry(geometry_)
{
    mesh = initGeometry->getMesh();

    nVert = mesh->nVertices();
    vInd = mesh->getVertexIndices();

    // Construct an additional geometry object to mess with
    devGeom = new Geometry<Euclidean>(*mesh);
    for(VertexPtr v : mesh->vertices()) {
        (*devGeom)[v] = (*initGeometry)[v];
    } 

    initGeometry->getHalfedgeAngles(initAngles);
}


void DevelopableApproximation::solve(const VertexData<char>& vertexFixed) {

    cout << endl << "Solving for developable approximation..." << endl;

    // Initialize
    vFixed = vertexFixed;
    buildFilterMatrix();

    int iIter = 0;
    int maxIter = 100;
    double startingStepSize = 1.0;
    while(iIter < maxIter) {

        double currEnergy = computeEnergy();
        cout << endl << "Starting step " << iIter << " with energy " << currEnergy << endl;

        VertexData<Vector3> gradient = computeGradient();

        double gradSlope = 0;
        for(VertexPtr v : mesh->vertices()) {
            gradSlope += norm2(gradient[v]);
        }
        cout << "Gradient slope = " << gradSlope << endl;

        // Check for convergence
        if(gradSlope < 1e-5) {
            cout << "Converged!" << endl;
        }

        // Backtracking line search
        double stepSize = startingStepSize;
        double tauShortenFactor = 0.5;
        double cExpectedDecrease = 0.2;
        VertexData<Vector3> initPositions = *devGeom;
        int nDecreases = 0;
        while(true) {

            // Take a candidate step
            for(VertexPtr v : mesh->vertices()) {
                (*devGeom)[v] -= stepSize * gradient[v];
            }
            double newEnergy = computeEnergy();
    
            double energyDecrease = currEnergy - newEnergy; 
            double expectedDecrease = stepSize * gradSlope;
    
            cout << endl << "stepSize = " << stepSize << endl;
            cout << "decrease = " << energyDecrease << endl;
            cout << "expected decrease = " << expectedDecrease << endl;

            // Accept the step
            if(energyDecrease > cExpectedDecrease * expectedDecrease) {
                break;
            } 
            // Try a smaller step
            else {
                // Shorten the step size
                stepSize *= tauShortenFactor;
                nDecreases++;
    
                // Roll back the iteration
                for(VertexPtr v : mesh->vertices()) {
                    (*devGeom)[v] = initPositions[v];
                } 
            }
        }

        // Tweak the starting size so future searches go faster
        if(nDecreases == 0) {
            startingStepSize *= 2.0;
        } else if(nDecreases > 1) {
            startingStepSize /= 2.0;
        }

        cout << "Took step with size " << stepSize << endl;
        cout << "Total interior curvature = " << computeInteriorCurvature() << endl;

        iIter++;
    }

    if(iIter == maxIter) {
        cout << "Terminated on max iterations :(" << endl;
    }


    cout << "Optimization finished with energy " << computeEnergy() << endl;
}

VertexData<Vector3> DevelopableApproximation::computeGradient() {

    VertexData<Vector3> gradient(mesh, Vector3{0.0,0.0,0.0});

    // Useful data 
    HalfedgeData<double> currAngles;
    devGeom->getHalfedgeAngles(currAngles);
    VertexData<double> angleDefect;
    devGeom->getVertexAngleDefects(angleDefect);
    FaceData<Vector3> faceNormals;
    devGeom->getFaceNormals(faceNormals);

    // Add gradient terms with respect to the angle opposite this halfedge
    for(HalfedgePtr he : mesh->halfedges()) {

        VertexPtr vA = he.vertex();
        VertexPtr vB = he.next().vertex();
        VertexPtr vRoot = he.next().next().vertex();

        // Energy excludes fixed vertices
        if(vFixed[vRoot]) continue;

        double theta = currAngles[he];
        double defect = angleDefect[vRoot];
        Vector3 normal = faceNormals[he.face()];

        { // vA term
            Vector3 dThetadA = -devGeom->vector(he.next().next()).rotate_around(normal, PI/2.0);
            double r2 = norm2(dThetadA);
            dThetadA /= r2;
            gradient[vA] += -2 * defect * dThetadA;
        } 
        
        { // vB term
            Vector3 dThetadB = -devGeom->vector(he.next()).rotate_around(normal, PI/2.0);
            double r2 = norm2(dThetadB);
            dThetadB /= r2;
            gradient[vB] += -2 * defect * dThetadB;
        } 
        
        { // vRoot term
            Vector3 basisX = unit(devGeom->vector(he));
            Vector3 basisY = basisX.rotate_around(normal, PI/2.0);
            double L = norm(devGeom->vector(he));
            Vector3 vAR = -devGeom->vector(he.next().next());
            double pX = dot(basisX, vAR);
            double pY = dot(basisY, vAR);

            double gradPx = -pY / ((L-pX)*(L-pX) + pY*pY) + pY / (pX*pX + pY*pY);
            double gradPy = (-L + pX) / ((L-pX)*(L-pX) + pY*pY) - pX / (pX*pX + pY*pY);

            Vector3 dThetadR = basisX * gradPx + basisY * gradPy;
            gradient[vRoot] += -2 * defect * dThetadR;
        } 
    }

    // Zero gradient at fixed vertices
    for(VertexPtr v : mesh->vertices()) {
        if(vFixed[v]) {
            gradient[v] = Vector3{0.0, 0.0, 0.0};
        }
    }

    // Filter gradient
    filterGradient(gradient);

    return gradient;
}

double DevelopableApproximation::computeInteriorCurvature() {

    VertexData<double> angleDefect;
    devGeom->getVertexAngleDefects(angleDefect);

    double total = 0;
    for(VertexPtr v : mesh->vertices()) {

        // Energy excludes fixed vertices
        if(vFixed[v]) continue;

        total += std::abs(angleDefect[v]);
    }

    return total;
}

double DevelopableApproximation::computeEnergy() {

    HalfedgeData<double> currAngles;
    devGeom->getHalfedgeAngles(currAngles);

    VertexData<double> angleDefect;
    devGeom->getVertexAngleDefects(angleDefect);

    double energy = 0;
    for(VertexPtr v : mesh->vertices()) {

        // Energy excludes fixed vertices
        if(vFixed[v]) continue;

        energy += angleDefect[v] * angleDefect[v];
    }

    // for(HalfedgePtr he : mesh->halfedges()) {
    //     energy += 
    // }

    return energy;
}

void DevelopableApproximation::buildFilterMatrix() {


    GC::SparseMatrix<double> hodge1 = GC::buildHodge1<double>(initGeometry);
    GC::SparseMatrix<double> d0 = GC::buildDerivative0<double>(mesh);
    GC::SparseMatrix<double> d0T = d0.transpose();
    filterMatrix = d0T * hodge1 * d0;
    filterMatrix.shift(1e-4);

    // Boundary conditions
    GC::DenseVector<double> mask(nVert);
    for(VertexPtr v : mesh->vertices()) {
        size_t i = vInd[v];
        if(vFixed[v]) {
            mask(i) = 0.0;
        } else {
            mask(i) = 1.0 / filterMatrix(i,i);
        }
    }
    filterMatrix = mask.asDiagonal() * filterMatrix;
    for(VertexPtr v : mesh->vertices()) {
        size_t i = vInd[v];
        if(vFixed[v]) {
            filterMatrix(i,i) = 1.0;
        }
    }

}


void DevelopableApproximation::filterGradient(VertexData<Vector3>& gradient) {

    for(size_t j = 0; j < 3; j++) {

        // Copy in
        GC::DenseVector<double> initVals(nVert);
        for(VertexPtr v : mesh->vertices()) {
            initVals(vInd[v]) = gradient[v][j];
        }

        // Solve
        // TODO could be a Cholesky solve if we treat boundaries less lazily
        GC::DenseVector<double> filteredVals;
        GC::solve(filterMatrix, filteredVals, initVals);

        // Copy out
        for(VertexPtr v : mesh->vertices()) {
            gradient[v][j] = filteredVals(vInd[v]);
        }
    }

}