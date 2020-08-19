#include <cstdio>
#include <iostream>

#include "mkl.h"

#include "typedefs.h"
#include "custom_asserts.h"

#include "FileIO.h"
#include "Grid.h"
#include "Bsplines.h"
#include "Matrix.h"
#include "GaussLegendre.h"

void writeToFile(Atom::Bsplines *Bsplines, Atom::Grid *Grid) {
    std::vector<std::vector<f64>> knotPointsOutput;

    for (auto knotPt : Bsplines->m_knotPoints) {
        std::vector<f64> row;
        row.push_back(knotPt.real());
        knotPointsOutput.push_back(row);
    }

    FileIO::writeRowColDataToFile(knotPointsOutput, "../knotpoints.dat");

    std::vector<std::vector<f64>> outputData;
    for (auto val : Grid->getGridPoints()) {
        std::vector<f64> rowVector;

        rowVector.push_back(val.real());

        for (u32 bsplineIndex = 0; bsplineIndex < Bsplines->m_numBsplines; bsplineIndex++) {
            auto result = Bsplines->GetBsplineAtCoordinate(val, bsplineIndex);
            rowVector.push_back(result.real());
        }

        outputData.push_back(rowVector);
    }

    FileIO::writeRowColDataToFile(outputData, "../bsplines.dat");
}



void testLapack() {
    Logger::Trace("Entered testLapack()");

    Logger::Trace("Testing CBLAS zgemv y = A*x");

    Complex rhs[3] = {Complex(1.0), Complex(2.0), Complex(3.0)};
    Complex A[3][3];
    for (int i = 0; i < 1; i++) {
        A[i][i] = Complex(2.0, 0.0);
    }
    A[0][1] = Complex(1.0, 0.0);
    A[1][1] = Complex(2.0, 0.0);
    A[2][2] = Complex(3.0, 1.0);

    CBLAS_LAYOUT Layout = CblasRowMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;

    MKL_INT m = 3, n = m;
    Complex alpha = Complex(1.0, 0.0);
    Complex beta = Complex(0.0, 0.0);
    MKL_INT lda = n, incy = 1, incx = 1;

    Complex result[3] = {Complex(1.0)};

    Complex aa[2] = {Complex(1.0, 2.0), Complex(3.0, 2.0)};
    Complex bb[2] = {Complex(3.0, 1.0), Complex(4.0, 2.0)};
    Complex c;

    Logger::Trace("testing zdotu: c = [(1.0, 2.0), (3.0, 2.0)] . [(3.0, 1.0), (4.0, 2.0)]:");
    cblas_zdotu_sub(2, aa, 1, bb, 1, &c);
    Logger::Trace("c = (%f, %f)", c.real(), c.imag());

    Logger::Trace("Testing cblas_zgemv!");

    cblas_zgemv(Layout, trans, m, n, &alpha, A, lda, rhs, incx, &beta, result, incy);
    Logger::Trace("zgemv done");
    for (int i = 0; i < 3; i++) {

        Logger::Trace("rhs[%i]: (%f, %f), result[%i]: (%f, %f)",
                      i, rhs[i].real(), rhs[i].imag(), i, result[i].real(), result[i].imag()
        );

    }
}

void MKL_Test() {
    Logger::Trace("Entered MKL_Test()");

    MKLVersion ver;
    int len = 198;
    char buf[198];

    MKL_Get_Version_String(buf, len);
    printf("\n%s\n", buf);
    printf("\n");

    MKL_Get_Version(&ver);
    printf("Major version:          %d\n", ver.MajorVersion);
    printf("Minor version:          %d\n", ver.MinorVersion);
    printf("Update version:         %d\n", ver.UpdateVersion);
    printf("Product status:         %s\n", ver.ProductStatus);
    printf("Build:                  %s\n", ver.Build);
    printf("Processor optimization: %s\n", ver.Processor);
    printf("================================================================\n");
    printf("\n");

    double a[2] = {1.0, 2.0};

    double result = cblas_ddot(2, a, 1, a, 1);
    Logger::Trace("result %f", result);
    testLapack();
}

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

int main() {
//    MKL_Test();

    // Atomic units. Let's do a 10 Bohr radii sized grid.
    constexpr u32 numGridPoints = 1000;
    constexpr f64 gridStart = 0.0, gridEnd = 10.0;
    Atom::Grid Grid = Atom::Grid(numGridPoints, gridStart, gridEnd);

    constexpr u32 bsplineOrder = 6;
    constexpr u32 numKnotPoints = 30;
    Atom::Bsplines Bsplines = Atom::Bsplines(numKnotPoints, bsplineOrder);
    // Linear knotpoint sequence;
    Bsplines.setupKnotPoints(Grid.getGridPoints());

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

    // We use Hartree atomic units. hbar = 1, electron charge e = 1, bohr radius a_0 = 1, electron mass m_e = 1.
    // We also have that 4\pi \epsilon_0 = 1,
    // since the Bohr radius a_0 = (4\pi epsilon_0 hbar^2)/(m_e e^2) = 1, and hbar = m_e = e = 1 in atomic units.
    // Then the radial Hamiltonian H' is
    // H' = (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r).
    constexpr u32 Z = 1; // Hydrogen atom.

    u32 l = Atom::AngularMomentumQuantumNumber::s; // l = 0,1,2,3... <-> s,p,d,f,...

    u32 matrixDimension = bsplineIndices.size();

    ZMatrix H = ZMatrix(matrixDimension, matrixDimension);
    H.setToZero();

    auto testGLQuadrature = [&](u32 points) {
        // Test Guass Legendre quadrature by integrating \int_0^1 x^2 = 1/3.
        auto a = Complex(0.0, 0.0);
        auto b = Complex(1.0, 0.0);

        auto abscissae = GaussLegendreIntegration.getShiftedAbscissae(a, b, points);
        ASSERT(abscissae.size() == points);
        auto prefactor = GaussLegendreIntegration.b_minus_a_half(a, b);
        auto pWeights = GaussLegendreIntegration.getPointerToZWeights(points);

        auto result = Complex(0.0, 0.0);
        for(int i = 0; i < points; ++i) {
            auto f = abscissae[i]*abscissae[i]; // Function is x^2.
            result += prefactor*pWeights[i]*f;
        }

        Logger::Trace("Result of GL integration of x^2 from a = 0, to b = 1 (should be 1/3 = 0.333333):");
        Logger::Trace("(%f, %f)", result.real(), result.imag());
    };

    testGLQuadrature(bsplineOrder);

    auto calculate_H_matrix_element = [&](u32 i, u32 j, u32 l) {
        // We have three terms here for H B_i (r)
        // First \int dr B_j (-0.5*d_r^2) B_i. We can show that this in fact is equal to
        // \int dr dB_j/dr dB_i/dr.
        // Integrals are perfomed using gaussian quadrature.

        Complex return_value = Complex(0.0);
        return return_value;
    };

    ZMatrix B = ZMatrix(matrixDimension, matrixDimension);
    B.setToZero();




    return 0;
}