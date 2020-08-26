#include <cstdio>
#include <cmath>
#include <string>

#define MKL_INTERFACE_LAYER LP64
#define MKL_THREADING_LAYER INTEL

#include "typedefs.h"
#include "custom_asserts.h"

#include "FileIO.h"
#include "Grid.h"
#include "Bsplines.h"
#include "Matrix.h"
#include "GaussLegendre.h"
#include "Eigenproblems.h"


typedef LA::Matrix<Complex> ZMatrix;

namespace Atom {
    enum AngularMomentumQuantumNumber {
        s = 0,
        p = 1,
        d = 2,
        f = 3,
        g = 4,
    };

}

// fwd declare
void setupGeneralisedEVPMatrices(ZMatrix &H, ZMatrix &B,
                                 Atom::Bsplines &Bsplines, GaussLegendre::Integration &GLI,
                                 std::vector<u32> bsplineIndices, u32 l, u32 Z);

void writeKnotPointsToFile(const ZVector &knotPts) {
    DVector output;
    for (auto point : knotPts) {
        output.push_back(point.real());
    }

    FileIO::writeColDataToFile(output, "../dat/knotpoints.dat");
}

void writeGridToFile(const ZVector &grid) {
    DVector output;
    for (auto point : grid) {
        output.push_back(point.real());
    }

    FileIO::writeColDataToFile(output, "../dat/grid.dat");
}

int main() {
    constexpr f64 eVperHartree = 27.211396641308;

    // Atomic units. Let's do a 10 Bohr radii sized grid.
    constexpr u32 numGridPoints = 1000;
    constexpr f64 gridStart = 0.0, gridEnd = 60.0;
    Atom::Grid Grid = Atom::Grid(numGridPoints, gridStart, gridEnd);
    writeGridToFile(Grid.getGridPoints());

    constexpr u32 bsplineOrder = 6;
    constexpr u32 numKnotPoints = 80;
    ASSERT(numGridPoints > numKnotPoints);
    Atom::Bsplines Bsplines = Atom::Bsplines(numKnotPoints, bsplineOrder);

    Bsplines.setupKnotPoints(Grid.getGridPoints(), Atom::knotSequenceType::Linear);
    writeKnotPointsToFile(Bsplines.m_knotPoints);

    GaussLegendre::Integration GaussLegendreIntegration = GaussLegendre::Integration();

    // Boundary conditions are \Psi(r = 0) = 0, and \Psi(r = \infty) = 0,
    // so we don't use the first and last Bsplines.
    // This gives numKnotPoints-bsplineOrder-2 Bsplines in total.
    // We make an array with the relevant bspline indices for ease of access.
    std::vector<u32> bsplineIndices;
    // Note that I am using zero indexing for Bspline, ie the first Bspline has index 0,
    // and the last Bspline has index numBsplines-1.
    for (u32 i = 1; i <= numKnotPoints - bsplineOrder - 2; i++) {
        bsplineIndices.push_back(i);
    }

    // We use Hartree atomic units. hbar = 1, electron charge e = 1, bohr radius a_0 = 1, electron mass m_e = 1.
    // We also have that 4\pi \epsilon_0 = 1,
    // since the Bohr radius a_0 = (4\pi epsilon_0 hbar^2)/(m_e e^2) = 1, and hbar = m_e = e = 1 in atomic units.
    // Then the radial Hamiltonian H' is
    // H' = (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r).
    constexpr u32 Z = 1; // Hydrogen atom.

    u32 l = Atom::AngularMomentumQuantumNumber::s; // l = 0,1,2,3... <-> s,p,d,f,...

    // Homogeneous boundaries at r = 0, and r = \infty, so we remove first and last bspline.
    u32 matrixDimension = bsplineIndices.size();
    Logger::Trace("number of Bsplines included in calculation (matrixDimension): %i", matrixDimension);
    ZMatrix H = ZMatrix(matrixDimension, matrixDimension);
    H.setToZero();
    ZMatrix B = ZMatrix(matrixDimension, matrixDimension);
    B.setToZero();

    // Create and initialise generalised eigenvalue problem solution (gevp) resources.
    ZVector gevp_Eigenvalues;
    gevp_Eigenvalues.resize(matrixDimension, Complex(0.0));
    // Eigenvectors v_j = \sum c_i B_i are expressed in their BSpline coefficients c_i
    ZMatrix gevp_Eigenvectors = ZMatrix(matrixDimension, matrixDimension);

    // Here we setup the B matrix where each element B_ij is an integral \int dr B_j B_i.
    // We also setup the zeroth step in the iteration for the LHS, which is an integral
    // H_ij = \int B_j (-0.5* d^2/dr^2+0.5(l(l+1))/r^2 - Z/r) B_k dr.
    // Subsequently we will add the potential term V_ij = \int B_j V_ee B_i dr to the H-matrix,
    // but for the zeroth step we take it to be zero.
    setupGeneralisedEVPMatrices(H, B, Bsplines, GaussLegendreIntegration, bsplineIndices, l, Z);

    LAPACK::EigenParameters eigenParameters;
    eigenParameters.computeRightEigenvectors = true;
    eigenParameters.matrixLayout = LAPACK::MatrixLayout::RowMajor;
    eigenParameters.squareMatrixOrder = matrixDimension;

    LAPACK::Eigenproblems eigenSolvers = LAPACK::Eigenproblems();
    eigenSolvers.GeneralisedComplexSolver(eigenParameters, H, B, gevp_Eigenvalues, gevp_Eigenvectors);

//    for (int i = 0; i < gevp_Eigenvectors.numRows(); i++) {
//        printf("%i ", i);
//        for (int j = 0; j < gevp_Eigenvectors.numCols(); j++) {
//            printf(" %i (%f, %f) ", j, gevp_Eigenvectors(i,j).real(), gevp_Eigenvectors(i,j).imag());
//        }
//        printf("\n");
//    }
//

Logger::Trace(" ");
Logger::Trace("Negative eigenvalues:");
    for (int i = 0; i < gevp_Eigenvalues.size(); i++) {
        auto val = gevp_Eigenvalues[i];
        if (val.real() < 0.0 && i >= 20) {
            if ( i == 20) {
                auto val2 = gevp_Eigenvalues[i - 1];
                Logger::Trace("Nonphysical before first physical? %i, (%f, %f) Hartree", i - 1, val2.real(),
                              val2.imag());
            }
            Logger::Trace("%i (%f, %f) Hartree,    (%f, %f) eV", i,
                          val.real(), val.imag(),
                          val.real() * eVperHartree, val.imag() * eVperHartree);
        }
    }

    return 0;
}

void setupGeneralisedEVPMatrices(ZMatrix &H, ZMatrix &B, Atom::Bsplines &Bsplines,
                                 GaussLegendre::Integration &GLI, std::vector<u32> bsplineIndices,
                                 u32 l, u32 Z) {

    u32 matrixDimension = H.numRows();

    auto calculate_H_matrix_element = [&](u32 i, u32 j, f64 l, f64 Z) {
        // Since we skip the first (index 0) and last (index numBsplines-1) Bsplines to
        // account for the boundary conditions, the indices into the matrix are not the same as the Bspline indices.
        u32 bsplineIndex_j = bsplineIndices[j];
        u32 bsplineIndex_i = bsplineIndices[i];

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
            // 0.5\int dr dB_j/dr dB_i/dr.
            Complex dB_i = Bsplines.GetDerivativeAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex dB_j = Bsplines.GetDerivativeAtCoordinate(abscissae[m], bsplineIndex_j);
            auto f = dB_j * dB_i;
            term1 += prefactor * pWeights[m] * f;

            Complex B_i = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex B_j = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_j);
            if (l > 0.0) {
                // term2: \int B_j l(l+1)/r^2 B_i dr
                Complex r2 = abscissae[m] * abscissae[m];
                f = B_j * B_i / r2;
                term2 += prefactor * pWeights[m] * f;
            }
            // term3: \int B_j Z/r B_i dr
//            B_i = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_i);
//            B_j = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_j);
            Complex r = abscissae[m];
            f = B_j * B_i / r;
            term3 += prefactor * pWeights[m] * f;
        }

        // H_ij = \int B_j (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r) B_i
        f64 lfactor = l*(l+1.0);
        Complex return_value = 0.5 * term1;
        if (lfactor > 0.0) {
            return_value += (0.5 * lfactor) * term2;
        }
        return_value += (-1.0 * Z ) * term3;
        return return_value;
    };

    auto calculate_B_matrix_element = [&](u32 i, u32 j) {
        u32 bsplineIndex_j = bsplineIndices[j];
        u32 bsplineIndex_i = bsplineIndices[i];

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

    f64 l_double = (f64) l;
    f64 Z_double = (f64) Z;

    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            H(i, j) = calculate_H_matrix_element(i, j, l_double, Z_double);
            B(i, j) = calculate_B_matrix_element(i, j);
        }
    }


}