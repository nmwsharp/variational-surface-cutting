#pragma once

#include <vector>

// Suitesparse includes
#include <cholmod.h>


#include "dense_matrix.h"
#include "sparse_matrix.h"

/*

    A fast Cholesky solver for sparse SPD systems using SuiteSparse.

    Assumptions: SPD. Has an element on every diagonal entry.

*/


class FastCholesky {

    public:

        FastCholesky(size_t N);
        FastCholesky(GC::SparseMatrix<double>& A);
        ~FastCholesky();

        // Building the matrix
        void reserveColumns(const std::vector<size_t>& colSizes);
        void addValue(size_t iRow, size_t iCol, double val, bool definitelyNew = false);
        void shiftDiagonal(double shiftVal);

        // Fun stuff
        FastCholesky setRowsToIdentiy(const std::vector<size_t> idRows); // returns a new matrix

        // Solving systems
        void factor();
        GC::DenseVector<double> solve(const GC::DenseVector<double>& rhs);
        bool computeResidual = true;


        // Hide copy and move constructors, we don't wanna mess with that
        FastCholesky(const FastCholesky&other) = delete;
        FastCholesky&operator=(const FastCholesky&other) = delete;
        FastCholesky(FastCholesky&&other) = delete;
        FastCholesky&operator=(FastCholesky&&other) = delete;


    private:

        // State
        size_t N;
        size_t Nnonzero = 0;
        int phase = 0; // {CREATED, RESERVED, FACTORED}
        cholmod_common common;

        // Data
        cholmod_sparse* mat = nullptr;
        cholmod_factor* factorization = nullptr;
        bool ownMatrixData = true; // avoid double-free

        // Convenience pointers to matrix components
        double* values;
        SuiteSparse_long* rowIndices = nullptr;
        SuiteSparse_long* colStart = nullptr;
        SuiteSparse_long* colSize = nullptr;
        

};