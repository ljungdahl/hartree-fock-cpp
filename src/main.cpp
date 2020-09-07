/* NOTE(anton):
 * This is a test program that solves the Hartree-Fock equations for some nice test problem, probably Neon?
 * The numerical method is using Bsplines, and everything is written around that. The aim is to make some generalisations
 * so that different numerical methods might be used in the future.
 */
#include "typedefs.h"
#include "logger.h"
#include "custom_asserts.h"

#include "Parameters.h"
#include "Vector.h"
#include "Grid.h"
#include "Bsplines.h"
#include "GaussLegendre.h"
#include "HartreeFock.h"
#include "temporary_tests.h"

typedef LA::Vector<Complex> Vector;

// Forward declaration of local functions. Definitions after main().
namespace Atom {
    void SetupHandBMatrices(Matrix &H, Matrix &B, Atom::Bsplines &Bsplines,
                            Atom::SystemParameters params, GaussLegendre::Integration &GLI);
}
void TestGLI(GaussLegendre::Integration &GLI);

int main(int argc, char *argv[]) {

    // Here we get atomic system details, ie what Z, which orbitals are to be calculated, etc.
    auto systemParameters = Atom::GetAtomicSystemParameters();
    systemParameters.Log();

    auto gridParams = Atom::GetGridParameters();
    gridParams.Log();

    // Representation of a radial grid with a linear spacing of dr = (end-start)/numGridPts;
    Atom::Grid Grid = Atom::Grid(gridParams.numberOfGridPoints,
                                 gridParams.gridStart, gridParams.gridEnd);

    auto bsplineParams = Atom::GetBsplineParameters();
    bsplineParams.Log();

    // In principle a valid minimal situation would be numGridPoints == numPhysicalKnotPoints,
    // but for all interesting applications we should have the case that numGridPoints >> numKnotPoints.
    ASSERT(gridParams.numberOfGridPoints > bsplineParams.numberOfKnotPoints);
    Atom::Bsplines Bsplines = Atom::Bsplines(bsplineParams.numberOfKnotPoints, bsplineParams.bsplineOrder);

    Bsplines.setupKnotPoints(Grid.getGridPoints(), Atom::knotSequenceType::Linear);

    // See this function for how boundary conditions rationalise what bsplines are
    // used in the solution.
    Atom::SetupUsedBsplineIndicesFromBoundaryConditions(Bsplines);

    GaussLegendre::Integration GaussLegendreIntegration = GaussLegendre::Integration();

    // We start with the atomic Hamiltonian without any electron-electron interaction,
    // i.e. H = -0.5 dr^2 + 0.5 l(l+1)/r^2 - Z/r.
    //
    // The relevant equation using Bsplines is H'c = EBc, where
    // H' is the Hamiltonian H multiplied by a Bspline integral from the left. This is a
    // generalised eigenvalue problem, but we can do better. By preemptively inverting
    // the RHS Bspline matrix B -> B^-1, we can write a standard eigenvalue problem.
    // B^-1H'c = Ec, and we can solve for the Eigenvalues E and the
    // eigenvectors (Bspline coefficients) c.

    Matrix B_inverse = Matrix();
    Matrix H = Matrix();

    {
        Matrix B = Matrix(); // A temporary B matrix that will be inverted and saved into B_inverse.

        // This resizes H and B and sets up the matrix elements using gauss legendre integration of Bsplines.
        Atom::SetupHandBMatrices(H, B, Bsplines, systemParameters, GaussLegendreIntegration);
        LAPACK::InvertMatrixInPlace(B);
        B_inverse.resize(B.numRows(), B.numCols());
        B_inverse.setFromMatrix(B);
    }

//    TestGLI(GaussLegendreIntegration);

    // TODO(anton): We need to make an array of H_matrices, one for each subshell.
    // Or perhaps just make a "setup l(l+1)-term" routine and calculate on the fly to save memory?
    // Currently only He 1s^2, so we don't have a l(l+1)-term ( s -> l = 0).

    // TODO(anton): The HF instantiation should take parameters so that
    // arrays with all relevant data can be allocated and initialised.
    Atom::HartreeFock hf = Atom::HartreeFock();

    // The starting point for the HF solver will be just the independent particle
    // model without any electron-electron interaction.
    hf.PerformInitialStep(H, B_inverse);

//    hf.SelfConsistentSolution();

    return 0;
}

void TestGLI(GaussLegendre::Integration &GLI) {

    auto integrate_interval = [&](Complex a, Complex b) {
        u32 numIntegrationPoints = 6;
        auto abscissae = GLI.getShiftedAbscissae(a, b, numIntegrationPoints);
        ASSERT(abscissae.size() == numIntegrationPoints);
        auto prefactor = GLI.b_minus_a_half(a, b);
        if (std::abs(prefactor) < 1e-8) {
            return Complex(0.0);
        }
        auto pWeights = GLI.getPointerToZWeights(numIntegrationPoints);

        Complex term = Complex(0.0);
        for (int m = 0; m < numIntegrationPoints; ++m) {
            auto r = abscissae[m];
            auto r2 = r*r;
            term += prefactor * pWeights[m] * r2;
        }

        return term;
    };

    Complex integrate_x2 = Complex(0.0);
    integrate_x2 += integrate_interval(0.0, 0.5);
    integrate_x2 += integrate_interval(0.5, 1.0);

    Logger::Trace("integrate_x2 = (%f, %f), "
                  "The integral of x^2 from zero to one is 1/3.", integrate_x2.real(), integrate_x2.imag());

}

void Atom::SetupHandBMatrices(Matrix &H, Matrix &B, Atom::Bsplines &Bsplines,
        Atom::SystemParameters params, GaussLegendre::Integration &GLI) {
    Logger::Trace("entered setupHandBMatrices");
    // The matrix dimension n (for square matrices of size n x n) are determined by
    // the number of bsplines used, thus we get the dimension from the number of elements in the
    // usedBsplineIndices array.
    u32 n = Bsplines.m_usedBsplineIndices.size();

    B.resize(n,n);
    H.resize(n,n);

    auto calculate_H_matrix_element = [&](u32 i, u32 j, f64 l, f64 Z) {
        // Since we skip the first (index 0) and last (index numBsplines-1) Bsplines to
        // account for the boundary conditions, the indices into the matrix are not the same as the Bspline indices.
        i32 bsplineIndex_j = Bsplines.m_usedBsplineIndices[j];
        i32 bsplineIndex_i = Bsplines.m_usedBsplineIndices[i];

        if (std::abs(bsplineIndex_i - bsplineIndex_j) > Bsplines.m_order) {
            return Complex(0.0);
        }

        Complex a = Bsplines.m_knotPoints[bsplineIndex_i];
        Complex b = Bsplines.m_knotPoints[bsplineIndex_i + Bsplines.m_order];

        u32 numIntegrationPoints = Bsplines.m_order;
        auto abscissae = GLI.getShiftedAbscissae(a, b, numIntegrationPoints);
        ASSERT(abscissae.size() == numIntegrationPoints);
        auto prefactor = GLI.b_minus_a_half(a, b);
        if (std::abs(prefactor) < 1e-8) {
            return Complex(0.0);
        }
        auto pWeights = GLI.getPointerToZWeights(numIntegrationPoints);

        auto term1 = Complex(0.0, 0.0);
        auto term2 = Complex(0.0, 0.0);
        auto term3 = Complex(0.0, 0.0);
        // We have three terms here for H B_i (r)
        // Integrals are perfomed using gaussian quadrature.
        for (int m = 0; m < numIntegrationPoints; ++m) {
            // First \int dr B_j (-0.5*d_r^2) B_i. We can show that this in fact is equal to
            // 0.5\int dr dB_j/dr dB_i/dr, without 0.5 factor.
            Complex dB_i = Bsplines.GetDerivativeAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex dB_j = Bsplines.GetDerivativeAtCoordinate(abscissae[m], bsplineIndex_j);
            auto f = dB_j * dB_i;
            term1 += prefactor * pWeights[m] * f;

            Complex B_i = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex B_j = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_j);
            if (l > 0.0) {
                // term2: \int B_j l(l+1)/r^2 B_i dr, without l(l+1)
                Complex r2 = abscissae[m] * abscissae[m];
                f = B_j * B_i / r2;
                term2 += prefactor * pWeights[m] * f;
            }
            // term3: \int B_j Z/r B_i dr, without Z factor.
            Complex r = abscissae[m];
            f = B_j * B_i / r;
            term3 += prefactor * pWeights[m] * f;
        }

        // H_ij = \int B_j (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r) B_i
        f64 lfactor = l * (l + 1.0);
        Complex return_value = 0.5 * term1;
//        Complex return_value = 0.5 * Complex(0.0);
        if (lfactor > 0.0) {
            return_value += (0.5 * lfactor) * term2;
        }
        return_value += (-1.0 * Z) * term3;
        return return_value;
    };

    auto calculate_B_matrix_element = [&](u32 i, u32 j) {
        u32 bsplineIndex_j = Bsplines.m_usedBsplineIndices[j];
        u32 bsplineIndex_i = Bsplines.m_usedBsplineIndices[i];

        if (std::abs((f64) bsplineIndex_i - (f64) bsplineIndex_j) >= (f64) Bsplines.m_order) {
            return Complex(0.0);
        }

        Complex a = Bsplines.m_knotPoints[bsplineIndex_i];
        Complex b = Bsplines.m_knotPoints[bsplineIndex_i + Bsplines.m_order];

        u32 points = Bsplines.m_order;
        auto abscissae = GLI.getShiftedAbscissae(a, b, points);
        ASSERT(abscissae.size() == points);
        auto prefactor = GLI.b_minus_a_half(a, b);
        auto pWeights = GLI.getPointerToZWeights(points);

        auto term1 = Complex(0.0, 0.0);
        for (int m = 0; m < points; ++m) {
            Complex B_i = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex B_j = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_j);
            auto f = B_j * B_i;
            term1 += prefactor * pWeights[m] * f;
        }

        Complex return_value = term1;
        return return_value;
    };

    // TODO(anton): We're testing on He 1s^2 right now. Need to make this work with general
    // subshells later!!!
    f64 l_f64 = (f64) params.shells[0].l;
    f64 Z_f64 = (f64) params.Z;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H(i, j) = calculate_H_matrix_element(i, j, l_f64, Z_f64);
            B(i, j) = calculate_B_matrix_element(i, j);
        }
    }

}