#pragma once

#define MKL_INTERFACE_LAYER LP64

#include <vector>

#include "Matrix.h"

#include "LapackCommon.h"


typedef LA::Matrix<Complex> ZMatrix;

namespace LAPACK {


    class Eigenproblems {
    public:
        Eigenproblems();

        void GeneralisedComplexSolver(EigenParameters params, const ZMatrix &A_Matrix, const ZMatrix &B_matrix,
                                      ZVector &out_eigenvalues, ZMatrix &out_eigenvectors);
    private:

    };

}
