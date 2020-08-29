#pragma once

#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl.h"

// TODO(anton): Support lapack INFO return cases, currently I only allow completely successful calls.
#define LAPACK_CHECK(expr) { \
ASSERT(expr == 0); \
}

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

    struct LinearEqParameters {
        MatrixLayout matrixLayout = MatrixLayout::RowMajor;
        bool computeLeftEigenvectors = false;
        bool computeRightEigenvectors = true;
        u32 squareMatrixOrder = 1; // set to N for a N x N-matrix.
    };
}