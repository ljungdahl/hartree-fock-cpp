#include <cmath>
#include <algorithm> // std::sort, std::stable_sort
#include <numeric> // std::iota, just to fill a container with sequentally increasing values.

#include "Eigenproblems.h"

std::vector<u32> sortZVectorIndices(const ZVector &v) {
    std::vector<u32> indices(v.size());
    std::iota(indices.begin(), indices.end(), 0); // Fills indices with integers from 0 to v.size()-1

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort(indices.begin(), indices.end(),
                     [&v](u32 i, u32 j) {
                         Complex a = v[i];
                         Complex b = v[j];
                         return a.real() < b.real() || a.real() == b.real() && a.imag() < b.imag();
                     });
    return indices;
}

LAPACK::Eigenproblems::Eigenproblems() {
}

void LAPACK::Eigenproblems::GeneralisedComplexSolver(LAPACK::EigenParameters params, const ZMatrix &A_Matrix,
                                                     const ZMatrix &B_Matrix, ZVector &out_eigenvalues,
                                                     ZMatrix &out_eigenvectors) {
// This function wraps the zggev LAPACK-routine which solves the generalised eigenproblem
// Ac = EBc, where A,B are nonsymmetric n x n-matrices, and c an eigenvector.
//TODO(anton): Also implement left eigenvectors, currently only working for right eigenvectors as written above.

    u32 n = params.squareMatrixOrder;
    ASSERT(n * n == A_Matrix.totalSize() && n * n == B_Matrix.totalSize());
    u32 leadingDimension_A = n, leadingDimension_B = n;
    u32 leadingDimRightEigvec = n;
    u32 leadingDimLeftEigvec = 1; // NOTE(anton): currently only right eigvecs supported here.

    auto computeLeft = 'N';
    auto computeRight = params.computeRightEigenvectors ? 'V' : 'N';
    if (computeRight == 'N') {
        Logger::Warn("LAPACK::Eigenproblems::GeneralComplexSolver(): computeRightEigvecs == 'N'.");
    }

    LAPACK::Layout layout = params.matrixLayout;

    // Lapack destroys the input matrices, so we need to make local copies.
    // ZMatrix internal storage is a vector anyway and LAPACK routine just takes pointers.
    // Thus it is easy to use std::vectors for the local data and we can hope the compiler can do some
    // nice things for the copy operation.
    ZVector temp_A; //
    temp_A.resize(n * n);
    A_Matrix.copyToVector(temp_A.data());

    ZVector temp_B;
    temp_B.resize(n * n);
    B_Matrix.copyToVector(temp_B.data());

    // The generalised eigenvalue is a ratio E = alpha / beta, and what we get back from zggev are arrays
    // with the n values alpha and beta.
    ZVector alpha;
    alpha.resize(n);

    ZVector beta;
    beta.resize(n);

    ZVector leftEigVecs;
    leftEigVecs.resize(1);

    ZVector rightEigVecs;
    rightEigVecs.resize(n * n);

    LAPACK_CHECK(
            LAPACKE_zggev3(layout, computeLeft, computeRight, n,
                          temp_A.dataPtr(), leadingDimension_A,
                          temp_B.dataPtr(), leadingDimension_B,
                          alpha.dataPtr(), beta.dataPtr(),
                          leftEigVecs.dataPtr(), leadingDimLeftEigvec,
                          rightEigVecs.dataPtr(), leadingDimRightEigvec)
    );

    Logger::Log(
            "LAPACK::Eigenproblems::GeneralComplexSolver(): Call to LAPACKE_zggev finished successfully (INFO == 0)");


//    Logger::Trace("i alpha[i]:              beta[i]:       alpha[i]/beta[i]     ");
    ZVector unsorted_eigenvalues;
    for (int i = 0; i < n; i++) {
        Complex val;
        if (std::abs(beta[i]) < 1e-8) {
            val = alpha[i];
            Logger::Warn("LAPACK::Eigenproblems::GeneralComplexSolver(): \n"
                         "beta[%i] < 1e-8, eigenvalue %i uses only alpha component. (should be infinity)", i);
        }


        val = alpha[i] / beta[i];
        unsorted_eigenvalues.push_back(val);
//        Logger::Trace("%i        (%f,%f)            (%f, %f)          (%f,%f)",
//                i, alpha[i].real(), alpha[i].imag(),
//                beta[i].real(), beta[i].imag(),
//                val.real(), val.imag());
    }
    auto old_indices = sortZVectorIndices(unsorted_eigenvalues);

    ZVector eigenvalues;
    eigenvalues.resize(unsorted_eigenvalues.size(), Complex(0.0));
    ZVector sorted_alpha;
    ZVector sorted_beta;
    u32 i = 0;
    for (auto i_old : old_indices) {
        eigenvalues[i] = unsorted_eigenvalues[i_old];
        sorted_alpha.push_back(alpha[i_old]);
        sorted_beta.push_back(beta[i_old]);
//        eigenvalues.push_back(unsorted_eigenvalues[i]);
        i++;
    }

    if (eigenvalues.size() == out_eigenvalues.size()) {
        out_eigenvalues = eigenvalues;
    } else {
        out_eigenvalues.clear();
        out_eigenvalues.resize(eigenvalues.size(), Complex(0.0));
        out_eigenvalues = eigenvalues;
    }

    // From Intel MKL documentation on ?ggev:
    // For complex flavors:
    //"The i-th component of the j-th eigenvector vj is stored in vr[(i - 1) + (j - 1)*ldvr] for column major layout
    // and in vr[(i - 1)*ldvr + (j - 1)] for row major layout."
    // NOTE(anton): I think they use one-indexing (Fortran) here, so we should instead say (for row major):
    // Eigenvector_j(i) = vr[i*ldvr + j]. Note that we use the sorted indices here to get the eigenvector.
    ASSERT(out_eigenvectors.numCols() == n && out_eigenvectors.numRows() == n);
    ASSERT(layout == LAPACK::Layout::RowMajor);
    u32 ldvr = leadingDimRightEigvec;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u32 unsorted_index = old_indices[j];
            out_eigenvectors(i, j) = rightEigVecs[i * ldvr + unsorted_index];
        }
    }

}
