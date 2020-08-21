#pragma once

#include <vector>

#include "Matrix.h"

#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl.h"

// TODO(anton): Support lapack INFO return cases, currently I only allow completely successful calls.
#define LAPACK_CHECK(expr) { \
ASSERT(expr == 0); \
}

typedef LA::Matrix<Complex> ZMatrix;

namespace LAPACK {

    enum MatrixLayout {
        // Ref mkl_lapacke.h
        RowMajor = LAPACK_ROW_MAJOR,
        ColMajor = LAPACK_COL_MAJOR
    };

    struct EigenParameters {
        MatrixLayout matrixLayout = MatrixLayout::RowMajor;
        bool computeLeftEigenvectors = false;
        bool computeRightEigenvectors = true;
        u32 squareMatrixOrder = 1; // set to N for a N x N-matrix.
    };

    class Eigenproblems {
    public:
        Eigenproblems();

        void GeneralisedComplexSolver(EigenParameters params, const ZMatrix &A_Matrix, const ZMatrix &B_matrix,
                                      std::vector<Complex> &out_eigenvalues, ZMatrix &out_eigenvectors);

    private:

    };

}
