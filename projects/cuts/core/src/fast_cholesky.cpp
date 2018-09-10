#include "fast_cholesky.h"

#include <algorithm>
#include <numeric>


#include "dense_matrix.h"
#include "timing.h"

// Undef to disable inner loop checks
#define FASTCHOLESKY_CHECKS

// Phases
#define PHASE_CREATED (0)
#define PHASE_RESERVED (1)
#define PHASE_FACTORED (2)

using namespace GC;


FastCholesky::FastCholesky(size_t N_) :
    N(N_)
{
    cholmod_l_start(&common);
    // common.supernodal = CHOLMOD_SUPERNODAL;

    // TODO supernodal vs simplicial?
}


FastCholesky::FastCholesky(GC::SparseMatrix<double>& A) {

    N = A.nRows();

    cholmod_l_start(&common);
    // common.supernodal = CHOLMOD_SUPERNODAL;

    // Convert representation
    mat = A.to_cholmod();
    ownMatrixData = false;

    // Right now it would break if the uesr tried to modify, so  just factor immediately and call it a day
    phase = PHASE_RESERVED;
    factor(); 
}

void FastCholesky::reserveColumns(const std::vector<size_t>& maxColSizes) {

    if(phase != PHASE_CREATED) {
        throw std::runtime_error("Cannot reserve except immediately after creation.");
    }

    if(maxColSizes.size() != N) {
        throw std::runtime_error("colSizes must be a vector of size N");
    }

    // Count nonzeros
    Nnonzero = std::accumulate(maxColSizes.begin(), maxColSizes.end(), (size_t)0);

    // Allocate the matrix
    int sorted = false;
    int packed = false;
    int stype = 1; // upper triangular
    mat = cholmod_l_allocate_sparse(N, N, Nnonzero, sorted, packed, stype, CHOLMOD_REAL, &common);

    // Pull out useful pointers
    values = (double*) mat->x;
    rowIndices = (SuiteSparse_long*) mat->i;
    colStart = (SuiteSparse_long*) mat->p;
    colSize = (SuiteSparse_long*) mat->nz;

    // Fill out the column stars and identity locations, since we know them
    size_t elemCount = 0;
    for(size_t iCol = 0; iCol < N; iCol++) {

        // Initialize diagonal
        rowIndices[elemCount] = iCol;
        colSize[iCol] = 1;

        // Fill start locations
        colStart[iCol] = elemCount;
        elemCount += maxColSizes[iCol];
    }
    colStart[N] = elemCount;

    // Zero out values so we can build incrementally
    std::fill(values, values+Nnonzero, 0.0);

    phase++;
}

// Convention: columns are stored unsorted, except that the diagonal entry is always first.
// Recall that we only store the upper triangular
// TODO implicitly swap row/column for better access pattern, since symmetric?
void FastCholesky::addValue(size_t iRow_, size_t iCol_, double val, bool definitelyNew) {

    SuiteSparse_long iRow = iRow_;
    SuiteSparse_long iCol = iCol_;

    // Only use upper triangular
    if(iCol < iRow) return;

#ifdef FASTCHOLESKY_CHECKS
    if(phase != PHASE_RESERVED) {
        throw std::runtime_error("Cannot insert value except after reservation.");
    }
    if(iRow_ >= N) {
        throw std::runtime_error("Invalid row index of " + std::to_string(iRow) + " matrix size is " + std::to_string(N));
    }
    if(iCol_ >= N) {
        throw std::runtime_error("Invalid col index of " + std::to_string(iCol) + " matrix size is " + std::to_string(N));
    }
#endif // FASTCHOLESKY_CHECKS

    // Identity gets handled special
    if(iRow == iCol) {
        values[colStart[iCol]] += val;
        return;
    }

    // Look for an existing index for this value
    SuiteSparse_long location = -77; 
    if(!definitelyNew) {
        for(SuiteSparse_long i = colStart[iCol] + 1; i < colStart[iCol] + colSize[iCol]; i++) {
            if(rowIndices[i] == iRow) {
                location = i;
                break;
            }
        }
    }

    // If needed, create the new entry
    if(location == -77) {

        // New location
        location = colStart[iCol] + colSize[iCol];

#ifdef FASTCHOLESKY_CHECKS
    if(location == colStart[iCol+1]) {
        throw std::runtime_error("Attempted to add too many entries to a column.");
    }
#endif // FASTCHOLESKY_CHECKS

        rowIndices[location] = iRow;
        colSize[iCol]++;
    }

    // Add the value
    values[location] += val; 
}
        
void FastCholesky::shiftDiagonal(double shiftVal) {

    if(phase != PHASE_RESERVED) {
        throw std::runtime_error("Cannot modify except after reservation.");
    }
    
    for(size_t i = 0; i < N; i++) { 
        values[colStart[i]] += shiftVal;
    }
}
        
void FastCholesky::factor() {

    if(phase != PHASE_RESERVED) {
        throw std::runtime_error("Cannot factor except after reservation.");
    }

    START_TIMING(factor)

    factorization = cholmod_l_analyze(mat, &common);
    cholmod_l_factorize(mat, factorization, &common);

    auto factorTime = FINISH_TIMING(factor);
    cout << "[FastCholesky] Factorization took " << pretty_time(factorTime) << endl;
    
    phase++;    
}

DenseVector<double> FastCholesky::solve(const DenseVector<double>& rhs) {
    
    if(phase != PHASE_FACTORED) {
        throw std::runtime_error("Cannot solve except after factorization.");
    }

    START_TIMING(solve)

    // Convert rhs to suitesparse
    cholmod_dense* rhsC = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_REAL, &common);
    double* rhsPtr = (double*) rhsC->x;
    for(size_t i = 0; i < N; i++) {
        rhsPtr[i] = rhs(i);
    }
    
    // Solve
    cholmod_dense* xC = cholmod_l_solve(CHOLMOD_A, factorization, rhsC, &common);

    // Compute residual
    if(computeResidual) {

        cholmod_dense* Ax = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_REAL, &common);
        double alpha[] = {1.0, 1.0};
        double beta[] = {0.0, 0.0};
        cholmod_l_sdmult(mat, 0, alpha, beta, xC, Ax, &common);

        double resMax = 0;
        double rhsMax = 0;
        double* AxPtr = (double*) Ax->x;
        for(size_t i = 0; i < N; i++) {
            resMax = std::max(resMax, std::abs(AxPtr[i] - rhsPtr[i]));
            rhsMax = std::max(rhsMax, std::abs(rhsPtr[i]));
        }

        cout << "[FastCholesky] Relative residual = " << (resMax / rhsMax) << endl;
        cholmod_l_free_dense(&Ax, &common);
    }

    // Copy result
    DenseVector<double> x(N);
    double* xPtr = (double*) xC->x;
    for(size_t i = 0; i < N; i++) {
        x(i) = xPtr[i];
    }

    cholmod_l_free_dense(&rhsC, &common);
    cholmod_l_free_dense(&xC, &common);
    
    auto solveTime = FINISH_TIMING(solve);
    cout << "[FastCholesky] Solve took " << pretty_time(solveTime) << endl;

    return x;
}

FastCholesky::~FastCholesky() {

    if(mat && ownMatrixData) {
        cholmod_l_free_sparse(&mat, &common);
    }
    if(factorization) {
        cholmod_l_free_factor(&factorization, &common);
    }

    cholmod_l_finish(&common);
}

