#include "Eigenproblems.h"

LAPACK::Eigenproblems::Eigenproblems() {
}

void LAPACK::Eigenproblems::GeneralisedComplexSolver(LAPACK::EigenParameters params, const ZMatrix &A_Matrix,
                                                 const ZMatrix &B_Matrix, std::vector<Complex> &out_eigenvalues,
                                                 ZMatrix &out_eigenvectors) {
// This function wraps the zggev LAPACK-routine which solves the generalised eigenproblem
// Ac = EBc, where A,B are nonsymmetric n x n-matrices, and c an eigenvector.
//TODO(anton): Also implement left eigenvectors, currently only working for right eigenvectors as written above.

    u32 n = params.squareMatrixOrder;
    ASSERT(n*n == A_Matrix.TotalSize() && n*n == B_Matrix.TotalSize());
    u32 leadingDimension_A = n, leadingDimension_B = n;
    u32 leadingDimRightEigvec = n;
    u32 leadingDimLeftEigvec = 1; // NOTE(anton): currently only right eigvecs supported here.

    auto computeLeft = 'N';
    auto computeRight = params.computeRightEigenvectors ? 'V' : 'N';
    if (computeRight == 'N') {
        Logger::Warn("LAPACK::Eigenproblems::GeneralComplexSolver(): computeRightEigvecs == 'N'.");
    }

    LAPACK::MatrixLayout layout = params.matrixLayout;

    // Lapack destroys the input matrices, so we need to make local copies.
    // ZMatrix internal storage is a vector anyway and LAPACK routine just takes pointers.
    // Thus it is easy to use std::vectors for the local data and we can hope the compiler can do some
    // nice things for the copy operation.
    std::vector<Complex> temp_A; //
    temp_A.resize(n*n);
    A_Matrix.copyToVector(temp_A);

    std::vector<Complex> temp_B;
    temp_B.resize(n*n);
    B_Matrix.copyToVector(temp_B);

    // The generalised eigenvalue is a ratio E = alpha / beta, and what we get back from zggev are arrays
    // with the n values alpha and beta.
    std::vector<Complex> alpha;
    alpha.resize(n);

    std::vector<Complex> beta;
    beta.resize(n);

    std::vector<Complex> leftEigVecs;
    leftEigVecs.resize(1);

    std::vector<Complex> rightEigVecs;
    rightEigVecs.resize(n*n);

    LAPACK_CHECK(
            LAPACKE_zggev(layout, computeLeft, computeRight, n,
                    temp_A.data(), leadingDimension_A,
                    temp_B.data(), leadingDimension_B,
                    alpha.data(), beta.data(),
                    leftEigVecs.data(), leadingDimLeftEigvec,
                    rightEigVecs.data(), leadingDimRightEigvec)
            );


    // From Intel MKL documentation on ?ggev:
    // For complex flavors:
    //"The i-th component of the j-th eigenvector vj is stored in vr[(i - 1) + (j - 1)*ldvr] for column major layout
    // and in vr[(i - 1)*ldvr + (j - 1)] for row major layout."
    // NOTE(anton): I think they use one-indexing (Fortran) here, so we should instead say (for row major):
    // Eigenvector_j(i) = vr[i*ldvr + j]
    ASSERT(out_eigenvectors.numCols() == n && out_eigenvectors.numRows() == n);
    ASSERT(layout == LAPACK::MatrixLayout::RowMajor);
    u32 ldvr = leadingDimRightEigvec;
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            out_eigenvectors(i,j) = rightEigVecs[i*ldvr + j];
        }
    }


 Logger::Log("LAPACK::Eigenproblems::GeneralComplexSolver(): Call to LAPACKE_zggev finished successfully (INFO == 0)");

}

