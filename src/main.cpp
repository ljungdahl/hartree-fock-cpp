#include <iostream>
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

namespace Atom {
enum AngularMomentumQuantumNumber {
    s = 0,
    p = 1,
    d = 2,
    f = 3,
    g = 4,
};
}

typedef LA::Matrix<Complex> ZMatrix;

void writeKnotPointsToFile(const std::vector<Complex> &knotPts) {
    std::vector<f64> output;
    for (auto point : knotPts) {
        output.push_back(point.real());
    }

    FileIO::writeColDataToFile(output, "../dat/knotpoints.dat");
}

void writeGridToFile(const std::vector<Complex> &grid) {
    std::vector<f64> output;
    for (auto point : grid) {
        output.push_back(point.real());
    }

    FileIO::writeColDataToFile(output, "../dat/grid.dat");
}

int main() {
    constexpr f64 eVperHartree = 27.211396641308;

    // Atomic units. Let's do a 10 Bohr radii sized grid.
    constexpr u32 numGridPoints = 2000;
    constexpr f64 gridStart = 0.0, gridEnd = 80.0;
    Atom::Grid Grid = Atom::Grid(numGridPoints, gridStart, gridEnd);
    writeGridToFile(Grid.getGridPoints());

    constexpr u32 bsplineOrder = 6;
    constexpr u32 numKnotPoints = 100;
    Atom::Bsplines Bsplines = Atom::Bsplines(numKnotPoints, bsplineOrder);

    Bsplines.setupKnotPoints(Grid.getGridPoints(), Atom::knotSequenceType::Linear);
    writeKnotPointsToFile(Bsplines.m_knotPoints);

    GaussLegendre::Integration GaussLegendreIntegration = GaussLegendre::Integration();

    // We are looking to solve the Eigenvalue problem for the hydrogen atom,
    // H'\Psi = E\Psi. If we say that \Psi(r) = \sum c_i B_i(r) we can transform the problem
    // into an eigenvalue equation that solves for the coefficient vector c = ( c_0,...,c_n-1 ).
    // This equation reads Hc = EBc, where c is as above, and H and B are matrices. We get to this point by
    // multiplying H\sum c_i B_i = E\sum c_i B_i by B_j from the left, and then integrating over all space.
    // This means that an element H_ij of H' is
    // H_ij = \int dr B_j(r) H' B_i(r),
    // and an element of B is
    // B_ij = \int dr B_j(r)B_i(r).

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
    Logger::Trace("bsplineIndices.size(): %i, last element: %i",
                  bsplineIndices.size(), bsplineIndices[bsplineIndices.size() - 1]);

    // We use Hartree atomic units. hbar = 1, electron charge e = 1, bohr radius a_0 = 1, electron mass m_e = 1.
    // We also have that 4\pi \epsilon_0 = 1,
    // since the Bohr radius a_0 = (4\pi epsilon_0 hbar^2)/(m_e e^2) = 1, and hbar = m_e = e = 1 in atomic units.
    // Then the radial Hamiltonian H' is
    // H' = (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r).
    constexpr u32 Z = 1; // Hydrogen atom.

    u32 l = Atom::AngularMomentumQuantumNumber::d; // l = 0,1,2,3... <-> s,p,d,f,...

    // Homogeneous boundaries at r = 0, and r = infty, so we remove first and last bspline.
    u32 matrixDimension = bsplineIndices.size();
    Logger::Trace("matrixDimension: %i", matrixDimension);
    ZMatrix H = ZMatrix(matrixDimension, matrixDimension);
    H.setToZero();

    auto calculate_H_matrix_element = [&](u32 i, u32 j, f64 l, f64 Z) {
        u32 bsplineIndex_j = bsplineIndices[j];
        u32 bsplineIndex_i = bsplineIndices[i];

        Complex a = Bsplines.m_knotPoints[bsplineIndex_i];
        Complex b = Bsplines.m_knotPoints[bsplineIndex_i + Bsplines.m_order];

        u32 points = Bsplines.m_order;
        auto abscissae = GaussLegendreIntegration.getShiftedAbscissae(a, b, points);
        ASSERT(abscissae.size() == points);
        auto prefactor = GaussLegendreIntegration.b_minus_a_half(a, b);
        auto pWeights = GaussLegendreIntegration.getPointerToZWeights(points);

        // We have three terms here for H B_i (r)
        // First \int dr B_j (-0.5*d_r^2) B_i. We can show that this in fact is equal to
        // 0.5\int dr dB_j/dr dB_i/dr.
        // Integrals are perfomed using gaussian quadrature.
        auto term1 = Complex(0.0, 0.0);
        for (int m = 0; m < points; ++m) {
            Complex dB_i = Bsplines.GetDerivativeAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex dB_j = Bsplines.GetDerivativeAtCoordinate(abscissae[m], bsplineIndex_j);
            auto f = dB_j * dB_i;
            term1 += prefactor * pWeights[m] * f;
        }

        // term2: \int B_j l(l+1)/r^2 B_i dr
        auto term2 = Complex(0.0, 0.0);
        for (int m = 0; m < points; ++m) {
            Complex B_i = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex B_j = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_j);
            Complex r2 = abscissae[m] * abscissae[m];
            auto f = B_j * B_i / r2;
            term2 += prefactor * pWeights[m] * f;
        }

        // term3: \int B_j Z/r
        auto term3 = Complex(0.0, 0.0);
        for (int m = 0; m < points; ++m) {
            Complex B_i = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_i);
            Complex B_j = Bsplines.GetBsplineAtCoordinate(abscissae[m], bsplineIndex_j);
            Complex r = abscissae[m];
            auto f = B_j * B_i / r;
            term3 += prefactor * pWeights[m] * f;
        }

        Complex return_value = 0.5 * term1 + (0.5 * l * (l + 1)) * term2 - Z * term3;
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
        auto abscissae = GaussLegendreIntegration.getShiftedAbscissae(a, b, points);
        ASSERT(abscissae.size() == points);
        auto prefactor = GaussLegendreIntegration.b_minus_a_half(a, b);
        auto pWeights = GaussLegendreIntegration.getPointerToZWeights(points);

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

    ZMatrix B = ZMatrix(matrixDimension, matrixDimension);
    B.setToZero();

    for (int i = 0; i < matrixDimension; i++) {
        for (int j = 0; j < matrixDimension; j++) {
            H(i, j) = calculate_H_matrix_element(i, j, l, Z);
            B(i, j) = calculate_B_matrix_element(i, j);
        }
    }

    LAPACK::EigenParameters eigenParams;
    eigenParams.matrixLayout = LAPACK::RowMajor;
    eigenParams.computeRightEigenvectors = true;
    eigenParams.squareMatrixOrder = matrixDimension;

    std::vector<Complex> eigenvalues;
    eigenvalues.resize(matrixDimension);
    ZMatrix eigenvectorCoefficients = ZMatrix(matrixDimension, matrixDimension);

    LAPACK::Eigenproblems eigen = LAPACK::Eigenproblems();
    eigen.GeneralisedComplexSolver(eigenParams, H, B, eigenvalues, eigenvectorCoefficients);

    // Get the bound states, ie the states that have lowest energy for ground dstate ~0.5 Hartree
    // and still negative energy  (implying bound state).
    // TODO(anton): Understand where all of the other negative energies are coming from?
    std::vector<Complex> selectedEigenvalues;
    std::vector<u32> selectedIndices;
    u32 valIdx = 0;

    f64 lowestEnergies[3] = {/*1s:*/-14.0, /*2p:*/-4.0, /*3d:*/-1.6};

    for (auto val : eigenvalues) {
        auto lowerLimit = lowestEnergies[l];
        if (val.real() > lowerLimit && val.real() < 0.0) {
            selectedEigenvalues.push_back(val);
//            Logger::Trace("output eigenvalue %i: (%E, %E) Hartree, (%E, %E) eV",
//                          valIdx, val.real(), val.imag(), val.real() * eVperHartree, val.imag() * eVperHartree);
            selectedIndices.push_back(valIdx);
        }
        valIdx++;
    }

    std::vector<std::vector<Complex>> selected_eigenvectors;
    for (auto index : selectedIndices) {
        std::vector<Complex> eigenvector;
        for (int i = 0; i < eigenvectorCoefficients.numRows(); i++) {
            eigenvector.push_back(eigenvectorCoefficients(i, index));
        }
        selected_eigenvectors.push_back(eigenvector);
    }

    u32 vec_index = 0;
    for (auto &vec : selected_eigenvectors) {

        std::string filename = "../dat/eigenvector_function_" + std::to_string(vec_index) + ".dat";
        vec_index++;

        std::vector<Complex> outputFunction;

        for (auto r : Grid.getGridPoints()) {
            Complex value_at_r = Complex(0.0);
            u32 bspl_index = 0;

            for (auto coefficient : vec) {
                value_at_r += coefficient * Bsplines.GetBsplineAtCoordinate(r, bsplineIndices[bspl_index]);
                bspl_index += 1;
            }
            outputFunction.push_back(value_at_r);
        }

        FileIO::writeComplexVectorToFile(outputFunction, filename);
    }

    constexpr f64 eVperkHz = 1000.0 * 4.1356 * 1e-15;

    constexpr f64 Hydrogen_n1_l0_NIST_level_kHz = -3288086857127.6; // kHz
    Logger::Trace("Hydrogen 1s (NIST): %f", Hydrogen_n1_l0_NIST_level_kHz * eVperkHz);

    constexpr f64 Hydrogen_n2_l0_NIST_level_kHz = -822025443940.5; // kHz
    Logger::Trace("Hydrogen 2s (NIST): %f", Hydrogen_n2_l0_NIST_level_kHz * eVperkHz);

    constexpr f64 Hydrogen_n2_l1_NIST_level_kHz = -822015532742.95; // kHz
    Logger::Trace("Hydrogen 2p (NIST): %f", Hydrogen_n2_l1_NIST_level_kHz * eVperkHz);


    constexpr f64 Hydrogen_n3_l2_NIST_level_kHz = -365339565239.0; // kHz
    Logger::Trace("Hydrogen 3d (NIST): %f", Hydrogen_n3_l2_NIST_level_kHz * eVperkHz);

    return 0;
}