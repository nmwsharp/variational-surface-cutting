// Author: Rohan Sawhney
// Fall 2016

#include "bff.h"
#ifdef HAVE_SUITESPARSE
GC::LinearContext context;
#endif

BoundaryFirstFlattening::BoundaryFirstFlattening(Geometry<Euclidean> *geometry_):
Parameterization(geometry_),
performSetup(true),
wc(0),
iwc(0),
bwc(0),
indices(mesh)
{
    
}

void BoundaryFirstFlattening::setup(bool assignVertexIndices)
{
    if (assignVertexIndices) {
        // wc = mesh->nVertices();
        // bwc = mesh->nImaginaryHalfedges();
        // iwc = wc - bwc;
        
        // // Assign vertex indices
        // unsigned int ivc = 0;
        // unsigned int ivb = iwc;
        // for (VertexPtr v: mesh->vertices()) {
        //     unsigned int vIdx = onCutBoundary(v) ? ivb++ : ivc++;
        //     for (CornerPtr c: v.adjacentCorners()) indices[c] = vIdx;
        // }
        // Assign wedge indices
        assignWedgeIndices(indices, true);
        wc = wedgeCount;
        iwc = interiorWedgeCount;
        bwc = boundaryWedgeCount;
        
    } else {
        // Assign wedge indices
        assignWedgeIndices(indices, true);
        wc = wedgeCount;
        iwc = interiorWedgeCount;
        bwc = boundaryWedgeCount;
    }
    
    // Initialize variables
    u.resize(wc);
    v.resize(wc);
    
    // Compute solvers
    int minIwc = max(0, iwc-1);
    int maxIwc = min(iwc, wc-1);
    
    A = cotanMatrix<double>(geometry, indices, wc);
    for (int i = 0; i < wc; i++) A(i, i) += 1e-8;
    
    Aii = A.sub(0, minIwc, 0, minIwc);
    Aib = A.sub(0, minIwc, maxIwc, wc-1);
    Abb = A.sub(maxIwc, wc-1, maxIwc, wc-1);
 
#ifndef HAVE_SUITESPARSE
    solverLd.compute(Aii.toEigen());
    solverLn.compute(A.toEigen());
#endif
}

double BoundaryFirstFlattening::computeWedgeAngleSum(HalfedgePtr h) const
{
    double sum = 0.0;
    HalfedgePtr he = h.twin().next();
    do {
        sum += geometry->angle(he.next());
        
        if (he.edge().isCut()) break;
        he = he.twin().next();
    } while (he.isReal());
    
    return sum;
}

void BoundaryFirstFlattening::computeBoundaryScaleFactors(GC::DenseMatrix<double>& boundaryScaleFactors)
{
    for (HalfedgePtr h: mesh->cutBoundary()) {
        unsigned int j = indices[wedge(h)] - iwc;
        
        if (targetLengths.size() > 0) {
            // Compute scale factors from target lengths
            HalfedgePtr p = previous(h);
            double lij = targetLengths[wIndices[wedge(p)]];
            double ljk = targetLengths[wIndices[wedge(h)]];
            
            boundaryScaleFactors(j) = (lij*log(lij/geometry->length(p.edge())) +
                                       ljk*log(ljk/geometry->length(h.edge())))/(lij + ljk);
            
        } else {
            // Map boundary isometrically
            boundaryScaleFactors(j) = 0.0;
        }
    }
}

void invert2x2(GC::DenseMatrix<double>& m)
{
    double det = m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0);
    
    swap(m(0, 0), m(1, 1));
    m(0, 1) = -m(0, 1);
    m(1, 0) = -m(1, 0);
    m *= 1.0/det;
}

void BoundaryFirstFlattening::closeLengths(GC::DenseMatrix<double>& lengths,
                                           const GC::DenseMatrix<Complex>& tangents) const
{
    // Create a map from halfedges to indices that indexes into unique lengths
    unsigned int lengthVariables = 0;
    EdgeData<Bool> visited(mesh);
    HalfedgeData<unsigned int> hIndices(mesh);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        if (!visited[h.edge()]) {
            hIndices[h] = lengthVariables++;
            visited[h.edge()] = true;
            
        } else {
            hIndices[h] = hIndices[h.twin()];
        }
    }
    
    // Accumulate mass matrix and tangents
    GC::DenseMatrix<double> l(lengthVariables);
    GC::DenseMatrix<double> T(2, lengthVariables);
    GC::SparseMatrix<double> M(lengthVariables, lengthVariables);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        unsigned int j = indices[wedge(h)] - iwc;
        unsigned int jj = hIndices[h];
        
        l(jj) = lengths(j);
        M(jj, jj) += 1.0/geometry->length(h.edge());
        
        T(0, jj) += tangents(j).real();
        T(1, jj) += tangents(j).imag();
    }
    
    // Invert mass matrix
    GC::SparseMatrix<double> Minv(lengthVariables, lengthVariables);
    for (unsigned int i = 0; i < lengthVariables; i++) Minv(i, i) = 1.0/M(i, i);
    GC::DenseMatrix<double> m = T*(Minv*T.transpose());
    invert2x2(m);
    
    // Modify lengths to ensure gamma closes
    l -= Minv*(T.transpose()*m)*(T*l);

    // Copy modified lengths back
    for (HalfedgePtr h: mesh->cutBoundary()) {
        unsigned int j = indices[wedge(h)] - iwc;
        unsigned int jj = hIndices[h];
        
        lengths(j) = l(jj);
    }
}

void BoundaryFirstFlattening::computeCurveApproximation(GC::DenseMatrix<double>& gammaReal,
                                                        GC::DenseMatrix<double>& gammaImag,
                                                        GC::DenseMatrix<double>& lengths,
                                                        const GC::DenseMatrix<double>& curvatures) const
{
    // Compute tangents
    double phi = 0.0;
    GC::DenseMatrix<Complex> tangents(bwc);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        unsigned int j = indices[wedge(h)] - iwc;
        
        phi += curvatures(j);
        tangents(j) = Complex(cos(phi), sin(phi));
    }
    
    // Modify lengths to ensure gamma closes
    closeLengths(lengths, tangents);
    
    // Compute gamma
    Complex gamma(0.0, 0.0);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        unsigned int j = indices[wedge(h)] - iwc;
        
        gammaReal(j) = gamma.real();
        gammaImag(j) = gamma.imag();
        gamma += lengths(j)*tangents(j);
    }
}

void BoundaryFirstFlattening::solveDirichletLaplace(GC::DenseMatrix<double>& x,
                                                    GC::DenseMatrix<double>& rhs)
{
    // Backsolve
#ifdef HAVE_SUITESPARSE
    GC::DenseMatrix<double> xx;
    solvePositiveDefinite(Aii, xx, rhs);
    x.setSubMatrix(0, max(0, iwc-1), 0, 0, xx);
    
#else
    Eigen::Matrix<double, Eigen::Dynamic, 1> xx = solverLd.solve(rhs.toEigen());
    for (int i = 0; i < iwc; i++) x(i) = xx(i);
#endif
}

void BoundaryFirstFlattening::computeHarmonicExtension(const GC::DenseMatrix<double>& boundaryData)
{
    // Compute right hand side
    GC::DenseMatrix<double> rhs = -(Aib*boundaryData);
    
    // Solve
    solveDirichletLaplace(u, rhs);
    u.setSubMatrix(iwc, wc-1, 0, 0, boundaryData);
}

void BoundaryFirstFlattening::solveNeumannLaplace(GC::DenseMatrix<double>& x,
                                                  GC::DenseMatrix<double>& rhs)
{
    // Backsolve
#ifdef HAVE_SUITESPARSE
    solvePositiveDefinite(A, x, rhs);
    
#else
    Eigen::Matrix<double, Eigen::Dynamic, 1> xx = solverLn.solve(rhs.toEigen());
    for (int i = 0; i < wc; i++) x(i) = xx(i);
#endif
    x.removeMean();
}

void BoundaryFirstFlattening::computeHarmonicConjugate()
{
    // Compute neumann boundary data
    GC::DenseMatrix<double> rhs(wc);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        HalfedgePtr p = previous(h);
        CornerPtr w = wedge(h);
        unsigned int i = indices[wedge(p)];
        unsigned int j = indices[w];
        unsigned int k = indices[w.prev()];
        
        rhs(j) = -0.5*(u(k) - u(i));
    }
    
    // Solve
    solveNeumannLaplace(v, rhs);
}

void BoundaryFirstFlattening::computeTargetLengthFlattening()
{
    // Compute boundary scale factors
    GC::DenseMatrix<double> boundaryScaleFactors(bwc);
    computeBoundaryScaleFactors(boundaryScaleFactors);

    GC::DenseMatrix<double> K0(iwc);
    if (iwc > 0) {
        // Solve del(u) = K0 | bdy(u) = g
        for (VertexPtr v: mesh->vertices()) {
            if (!onCutBoundary(v)) {
                K0(indices[v.corner()]) = geometry->angleDefect(v);
            }
        }
        
        GC::DenseMatrix<double> rhs = -(K0 + Aib*boundaryScaleFactors);
        solveDirichletLaplace(u, rhs);
    }
    
    // Evaluate normal derivative dudn = bdy(-del(u))
    GC::DenseMatrix<double> dudn = -(Aib.transpose()*u.sub(0, max(0, iwc-1)) + Abb*boundaryScaleFactors);
    
    // Compute target curvatures and lengths
    GC::DenseMatrix<double> curvatures(bwc);
    GC::DenseMatrix<double> lengths(bwc);
    double k0Sum = 0;
    for (HalfedgePtr h: mesh->cutBoundary()) {
        CornerPtr w = wedge(h);
        unsigned int j = indices[w] - iwc;
        unsigned int k = indices[w.prev()] - iwc;
        
        // Compute curvature
        double k0 = M_PI - computeWedgeAngleSum(h);
        k0Sum += k0;
        curvatures(j) = k0 - dudn(j);
        
        // Compute target length
        double uj = boundaryScaleFactors(j);
        double uk = boundaryScaleFactors(k);
        lengths(j) = exp(0.5*(uj + uk))*geometry->length(h.edge());
    }

    // Compute gamma
    GC::DenseMatrix<double> gammaReal(bwc);
    GC::DenseMatrix<double> gammaImag(bwc);
    computeCurveApproximation(gammaReal, gammaImag, lengths, curvatures);
    
    // Compute harmonic extension
    computeHarmonicExtension(gammaReal);
    
    // Compute harmonic conjugate
    computeHarmonicConjugate();
}

void BoundaryFirstFlattening::extendBoundaryCurve(const GC::DenseMatrix<double>& gammeReal,
                                                  const GC::DenseMatrix<double>& gammaImag)
{
    // Extend imaginary part
    computeHarmonicExtension(gammaImag);
    
    // Swap
    swap(u, v);
    
    // Extend real part
    computeHarmonicExtension(gammeReal);
}

void BoundaryFirstFlattening::computeTargetAngleFlattening()
{
    // Solve del(u) = K0 | dudn = k0 - k
    GC::DenseMatrix<double> rhs(wc);
    if (iwc > 0) {
        GC::DenseMatrix<double> K0(iwc);
        for (VertexPtr v: mesh->vertices()) {
            if (!onCutBoundary(v)) K0(indices[v.corner()]) = geometry->angleDefect(v);
        }
        
        rhs.setSubMatrix(0, max(0, iwc-1), 0, 0, -K0);
    }
    
    GC::DenseMatrix<double> dudn(bwc);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        CornerPtr w = wedge(h);
        unsigned int j = indices[w] - iwc;
        
        // Evaluate dudn
        double k0 = M_PI - computeWedgeAngleSum(h);
        dudn(j) = k0 - targetCurvatures[wIndices[w]];
    }
    
    rhs.setSubMatrix(iwc, wc-1, 0, 0, -dudn);
    solveNeumannLaplace(u, rhs);
    
    // Compute target curvatures and lengths
    GC::DenseMatrix<double> curvatures(bwc);
    GC::DenseMatrix<double> lengths(bwc);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        CornerPtr w = wedge(h);
        unsigned int j = indices[w];
        unsigned int k = indices[w.prev()];
        
        // Set target curvature
        curvatures(j-iwc) = targetCurvatures[wIndices[w]];
        
        // Compute target length
        lengths(j-iwc) = exp(0.5*(u(j) + u(k)))*geometry->length(h.edge());
    }
    
    // Compute gamma
    GC::DenseMatrix<double> gammaReal(bwc);
    GC::DenseMatrix<double> gammaImag(bwc);
    computeCurveApproximation(gammaReal, gammaImag, lengths, curvatures);
    
    // Extend
    extendBoundaryCurve(gammaReal, gammaImag);
}

void BoundaryFirstFlattening::computeConeFlattening()
{
    // Solve del(u) = K0 - K | bdy(u) = 0
    GC::DenseMatrix<double> K0(iwc), K(iwc);
    for (VertexPtr v: mesh->vertices()) {
        unsigned int j = indices[v.corner()];
        if (!v.isBoundary()) {
            K0(j) = geometry->angleDefect(v);
            if (coneAngles.find(vIndices[v]) != coneAngles.end()) K(j) = coneAngles[vIndices[v]];
        }
    }
    
    u.zero();
    GC::DenseMatrix<double> rhs = -(K0 - K);
    solveDirichletLaplace(u, rhs);
    
    // Make local copies of scale factors and indices
    GC::DenseMatrix<double> uCopy = u;
    CornerData<size_t> indicesCopy = indices;
    
    // Perform setup
    setup();
    
    // Copy scale factors into new array
    for (VertexPtr v: mesh->vertices()) {
        for (CornerPtr c: v.adjacentCorners()) {
            u(indices[c]) = uCopy(indicesCopy[c]);
        }
    }
    
    // Evaluate normal derivative dudn = bdy(-del(u))
    GC::DenseMatrix<double> dudn = -(Aib.transpose()*u.sub(0, max(0, iwc-1)) + Abb*u.sub(iwc, wc-1));
    
    // Compute target curvatures and lengths
    GC::DenseMatrix<double> curvatures(bwc);
    GC::DenseMatrix<double> lengths(bwc);
    for (HalfedgePtr h: mesh->cutBoundary()) {
        CornerPtr w = wedge(h);
        unsigned int j = indices[w];
        unsigned int k = indices[w.prev()];
        
        // Compute curvature
        double k0 = M_PI - computeWedgeAngleSum(h);
        curvatures(j-iwc) = k0 - dudn(j-iwc);
        
        // Compute target length
        lengths(j-iwc) = exp(0.5*(u(j) + u(k)))*geometry->length(h.edge());
    }

    // Compute gamma
    GC::DenseMatrix<double> gammaReal(bwc);
    GC::DenseMatrix<double> gammaImag(bwc);
    computeCurveApproximation(gammaReal, gammaImag, lengths, curvatures);
    
    // Extend
    extendBoundaryCurve(gammaReal, gammaImag);
}

void BoundaryFirstFlattening::setUVs()
{
    for (VertexPtr vv: mesh->vertices()) {
        for (CornerPtr c: vv.adjacentCorners()) {
            Vector2& uv = uvs[c];
            unsigned int j = indices[c];
            
            uv[0] = -u(j); // Accounting for clockwise boundary traversal
            uv[1] = v(j);
        }
    }
    
    normalize();
}

bool BoundaryFirstFlattening::flatten()
{
    if (mesh->nBoundaryLoops() > 1) {
        cout << "BFF does not work with holes. Continuing anyway because this isn't a good test." << endl;
        // return false;
    }
    
    // Perform setup
    if (performSetup) {
        setup(coneAngles.size() > 0);
        performSetup = false;
    }
    
    if (coneAngles.size() > 0) {
        // Compute flattening with prescribed cone angles
        computeConeFlattening();
    
    } else if (targetCurvatures.size() > 0) {
        // Compute target angles
        computeTargetAngleFlattening();
    
    } else {
        // Compute target lengths or default flattening with natural boundary conditions
        computeTargetLengthFlattening();
    }
    
    // Set uv coords
    setUVs();
    
    return true;
}
