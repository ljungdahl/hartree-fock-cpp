#include "logger.h"
#include "HartreeFock.h"
#include "MKL_utility_linear_algebra.h"
#include "Eigenproblems.h"

Atom::HartreeFock::HartreeFock() {
    m_eig = new LAPACK::Eigenproblems();
}

void Atom::HartreeFock::PerformInitialStep(Matrix& H, Matrix& B_inverse) {
    /* This function calculates the inital approximation to the orbital wavefunctions
     * in the Hartree-Fock approximation.
     * We start with the atomic Hamiltonian without any electron-electron interaction,
     * i.e. H = -0.5 dr^2 + 0.5 l(l+1)/r^2 - Z/r.
     *
     * The relevant equation using Bsplines is H'c = EBc, where
     * H' is H multiplied by a Bspline integral from the left. This is a
     * generalised eigenvalue problem, but we can do better. By preemptively inverting
     * the RHS Bspline matrix B -> B^-1, we can write a standard eigenvalue problem.
     * B^-1H'c = Ec, and we can solve for the Eigenvalues E and the
     * eigenvectors (Bspline coefficients) c.
     */

    Matrix LHS = MatMul(B_inverse, H); // LHS = B^-1 * H.

    //TODO(anton): HERE COMES ZGEEVSOLVE!
    u32 rows = LHS.numRows();
    u32 cols = LHS.numCols();

    m_eigenvectors.resize(rows, cols);
    m_eigenvalues.resize(rows);

    LAPACK::EigenParameters params;
    params.squareMatrixOrder = rows;

    m_eig->NonsymmetricComplexSolver(params, LHS, m_eigenvalues, m_eigenvectors);
//    TestZGEEV();
}

void Atom::HartreeFock::SelfConsistentSolution() {

}
void Atom::HartreeFock::TestZGEEV() {
/* Locals */
/* Parameters */
constexpr u32 N = 4;
constexpr u32 LDA = N;
constexpr u32 LDVL = N;
constexpr u32 LDVR = N;

    MKL_INT n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
    /* Local arrays */
    MKL_Complex16 w[N], vl[LDVL*N], vr[LDVR*N];
    MKL_Complex16 a[LDA*N] = {
        {-3.84,  2.25}, {-8.94, -4.75}, { 8.95, -6.53}, {-9.87,  4.82},
        {-0.66,  0.83}, {-4.40, -3.82}, {-3.50, -4.26}, {-3.15,  7.36},
        {-3.99, -4.73}, {-5.88, -6.60}, {-3.36, -0.40}, {-0.75,  5.23},
        { 7.74,  4.18}, { 3.66, -7.53}, { 2.58,  3.60}, { 4.59,  5.41}
    };
    /* Executable statements */
    printf( "LAPACKE_zgeev (row-major, high-level) Example Program Results\n" );
    /* Solve eigenproblem */
    info = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'V', 'V', n, a, lda, w, vl,
                          ldvl, vr, ldvr );
    Logger::Trace("info %i", info);
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
}
