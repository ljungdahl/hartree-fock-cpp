#pragma once

#include "MKL_include_and_define.h"

#include "linear_algebra_typedefs.h"


typedef LA::Matrix<Complex> ZMatrix;
typedef LA::Vector<Complex> ZVector;

namespace LAPACK {


    struct EigenParameters {
        Layout matrixLayout = Layout::RowMajor;
        bool computeLeftEigenvectors = false;
        bool computeRightEigenvectors = true;
        u32 squareMatrixOrder = 1; // set to N for a N x N-matrix.
    };

    class Eigenproblems {
    public:
        Eigenproblems();

        void GeneralisedComplexSolver(EigenParameters params, const ZMatrix &A_Matrix, const ZMatrix &B_matrix,
                                      ZVector &out_eigenvalues, ZMatrix &out_eigenvectors);
    private:

    };

}
