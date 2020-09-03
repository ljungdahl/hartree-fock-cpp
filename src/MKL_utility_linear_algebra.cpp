#include "MKL_utility_linear_algebra.h"

void LAPACK::InvertMatrixInPlace(Matrix &A) {

    u32 rows = A.numRows(), cols = A.numCols();
    //TODO(anton): Make this work for general types (or at least also real doubles)
    auto layout = LAPACK::Layout::RowMajor;
    Vector A_copy = Vector(rows * cols);
    A.copyToVector(A_copy.data()); // Is this needed?... we will overwrite A with its inverse eventually though.

    // zgetrf - performs LU factorisation. A_copy is overwritten by resulting LU matrix.
    u32 leadingDimensionA = cols; // must be at least max(1, cols) for row major layout
    if (layout == LAPACK::Layout::ColMajor) {
        leadingDimensionA = rows; // or at least max(1, rows) for column major layout.
    }
    // From MKL docs: "ipiv contains the pivot indices for 1 <= i <= min(rows, cols), row i was interchanged with row ipiv(i)"
    std::vector<lapack_int> ipiv;
    ipiv.resize(rows);


    u32 info;
    info = LAPACKE_zgetrf(layout, rows, cols,
                          A_copy.dataPtr(), leadingDimensionA,
                          ipiv.data());
    if (info != 0) {
        Logger::Error("LAPACKE_zgetrf returned info = %i. \n"
                      "From Lapack documentation: \n"
                      "If info = 0 the execution is successful.\n"
                      "If info= -i, parameter i had an illegal value. \n"
                      "If info = i, element u_ii is zero; \n "
                      "the factorisation is complete but U is exactly singular. Division by 0 will occur if U is used to solve lineqs.",
                      info);
        debugBreak();
    }

    // zgetri - Computes the inverse of an LU-factored general double complex-valued matrix
    u32 n = rows;

    info = LAPACKE_zgetri(layout, n, A_copy.dataPtr(), leadingDimensionA, ipiv.data());
    if (info != 0) {
        Logger::Error("LAPACKE_zgetri returned info = %i. \n"
                      "From Lapack documentation: \n"
                      "If info = 0 the execution is successful.\n"
                      "If info= -i, parameter i had an illegal value. \n"
                      "If info = i the i-th diagonal element of the factor U is zero, Uis singular, and the inversion could not be completed.",
                      info);
        debugBreak();
    }

    // Now A_copy is the inverted matrix.
    // Put it back into A. Vector type has a member data() that gives back _reference_ to its container (usually std::vector<T>&).
    // Compare with std::vector<T>.data(), which is a pointer to the underlying memory array. Vector instead has a .dataPtr() with this corresponding function.
    A.setFromVector(A_copy.data());
}

Matrix MatMul(const Matrix &A, const Matrix &B) {
    // This routine performs a matrix multiplication C = A * B;
    // The CBLAS call computes C = alpha*A*B+beta*C, but this wrapper uses alpha = 1 and beta = 0.

    // TODO(anton): Make this work for other types than complex. Also other types than row major.
    // also test if maybe zgemm3m is better and faster for most usecases?

    auto layout = static_cast<const CBLAS_LAYOUT>(CBLAS_LAYOUT::CblasRowMajor);
    auto transposeA = CBLAS_TRANSPOSE::CblasNoTrans;
    auto transposeB = transposeA;

    ASSERT(A.numCols() == B.numRows()); // we need to multiply m x k-matrix by a k x n-matrix to get a m x n-matrix.
    u32 k = A.numCols();
    u32 m = A.numRows();
    u32 n = B.numCols();

    Complex alpha = Complex(1.0, 0.0);
    Complex beta = Complex(0.0, 0.0);

    u32 leadingDimA = k;
    u32 leadingDimB = n;
    u32 leadingDimC = n;

    if (layout == CBLAS_LAYOUT::CblasRowMajor) {
        leadingDimA = m;
        leadingDimB = k;
        leadingDimC = m;
    }

    Matrix C = Matrix(m, n);

    cblas_zgemm(layout, transposeA, transposeB,
                m, n, k, &alpha,
                A.GetDataPtr(), leadingDimA,
                B.GetDataPtr(), leadingDimB,
                &beta,
                C.GetDataPtr(), leadingDimC);

    return C;
}
