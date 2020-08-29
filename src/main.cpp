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

static constexpr f64 g_OneOverFourPi = 0.07957747155;
static constexpr f64 g_FourPi = 12.56637061436;
typedef LA::Matrix<Complex> ZMatrix;

namespace Atom {
    enum AngularMomentumQuantumNumber {
        s = 0,
        p = 1,
        d = 2,
        f = 3,
        g = 4,
    };

    struct SubShell {
        u32 n; // PrincipalQuantumNumber
        u32 l; // AngularMomentumQuantumNumber
        u32 NumberOfElectrons;
        f64 radialFunctionIntegratedOverAllSpace = 1.0;
    };

    struct SystemParameters {
        std::string Name;
        u32 Z;
        std::vector<Atom::SubShell> shells;
        u32 TotalOccupationNumber = 0;
    };

}


// fwd declare
void setupPoissonLHS(Atom::Bsplines &Bsplines, ZMatrix &rhs, GaussLegendre::Integration &GLI);

void setupGeneralisedEVPMatrices(ZMatrix &H, ZMatrix &B,
                                 Atom::Bsplines &Bsplines, GaussLegendre::Integration &GLI,
                                 std::vector<u32> bsplineIndices, u32 l, u32 Z);

Complex ChargeDensityAtCoordinate(Complex r, Atom::Bsplines &Bsplines,
                                  const ZMatrix &coefficientMatrix, const std::vector<u32> &physicalIndices,
                                  const Atom::SystemParameters &params);

Complex GetRadialFunctionAtCoordinate(Complex r, Atom::Bsplines &Bsplines, u32 functionIndex,
                                      const ZMatrix &coefficients);

void CalculateRadialFunctionNorms(Atom::Bsplines &Bsplines,
                                  const ZMatrix &coefficientMatrix, const std::vector<u32> &physicalIndices,
                                  Atom::SystemParameters &params, GaussLegendre::Integration &GLI);

void IntegrateChargeDensity(Atom::Bsplines &Bsplines, ZMatrix &coeffs, std::vector<u32> &physicalIndices,
                            Atom::SystemParameters AtomParameters, GaussLegendre::Integration &GLI);

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

void writeBsplineBasisFunctionToFile(const ZVector &basisFunctionOfR, u32 bsplineIndex, u32 derivativeOrder) {
    DVector output;
    for (auto point : basisFunctionOfR) {
        output.push_back(point.real());
    }
    std::string filename;
    if (derivativeOrder == 2) {
        filename = "../dat/dB2_" + std::to_string(bsplineIndex) + ".dat"; // dB2
    } else if (derivativeOrder == 1) { // dB
        filename = "../dat/dB_" + std::to_string(bsplineIndex) + ".dat"; // dB2
    } else { // B
        filename = "../dat/B_" + std::to_string(bsplineIndex) + ".dat"; // dB2
    }
    ASSERT(!filename.empty());

    FileIO::writeColDataToFile(output, filename);
}

void TestBsplineDerivatives(Atom::Bsplines &Bsplines, Atom::Grid &Grid) {
    // We test k = 4 Bsplines, regular, first derivative, and second derivative.


    std::vector<u32> orders = {0, 1, 2};

    u32 bsplineIndex = Bsplines.m_usedBsplineIndices[Bsplines.m_usedBsplineIndices.size()-1];
    for (auto order : orders) {
        u32 derivOrder = order;
        ZVector basisFunction;
        for (auto r : Grid.getGridPoints()) {
            Complex val = Complex(0.0);
            if (derivOrder == 0) {
                val = Bsplines.GetBsplineAtCoordinate(r, bsplineIndex);
            }
            if (derivOrder == 1) {
                val = Bsplines.GetDerivativeAtCoordinate(r, bsplineIndex);
//                Logger::Trace("first derivative for bspl #%i: (%f, %f)", bsplineIndex, val.real(), val.imag());
            }
            if (derivOrder == 2) {
                val = Bsplines.GetSecondDerivativeAtCoordinate(r, bsplineIndex);
//                Logger::Trace("second derivative for bspl #%i: (%f, %f)", bsplineIndex, val.real(), val.imag());
            }
            basisFunction.push_back(val);
        }
        writeBsplineBasisFunctionToFile(basisFunction, bsplineIndex, derivOrder);
    }
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
    Bsplines.SetBoundaryConditionBsplineIndices(bsplineIndices);

    // We use Hartree atomic units. hbar = 1, electron charge e = 1, bohr radius a_0 = 1, electron mass m_e = 1.
    // We also have that 4\pi \epsilon_0 = 1,
    // since the Bohr radius a_0 = (4\pi epsilon_0 hbar^2)/(m_e e^2) = 1, and hbar = m_e = e = 1 in atomic units.
    // Then the radial, hydrogen-like, Hamiltonian H' is
    // H' = (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r).
    constexpr u32 Z = 2; // Helium
    u32 l = Atom::AngularMomentumQuantumNumber::s; // l = 0,1,2,3... <-> s,p,d,f,...

    Atom::SubShell _1s;
    _1s.n = 1;
    _1s.l = l;
    _1s.NumberOfElectrons = 2;

    Atom::SystemParameters AtomParameters;
    AtomParameters.shells.push_back(_1s);
    AtomParameters.Z = Z;
    AtomParameters.Name = "Helium";
    for (auto shell : AtomParameters.shells) {
        AtomParameters.TotalOccupationNumber += shell.NumberOfElectrons;
    }

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

    // We get more solutions than are physical bound solutions, for some reason it seems to be always after the 20th
    // solution? At least for enough Bsplines. We will save the physical solution eigenvectors and use them to
    // compute the radial charge density.
    std::vector<u32> physicalIndices;
    Logger::Trace(" ");
    Logger::Trace("Negative eigenvalues:");
    for (int i = 0; i < gevp_Eigenvalues.size(); i++) {

        auto val = gevp_Eigenvalues[i];

        if (val.real() < 0.0 && i >= 20) {

            if (i == 20) {
                auto val2 = gevp_Eigenvalues[i - 1];
                Logger::Trace("Nonphysical before first physical? %i, (%f, %f) Hartree", i - 1, val2.real(),
                              val2.imag());
            }

            Logger::Trace("%i (%f, %f) Hartree,    (%f, %f) eV", i,
                          val.real(), val.imag(),
                          val.real() * eVperHartree, val.imag() * eVperHartree);
            physicalIndices.push_back(i);
        }

    }

    // Now we calculate radial charge density (averaged over angular parts)
    // This we can show is
    // rho(r) = 1/4pi sum_j^occupied orbitals e * N_j (P_njlj(r)/r)^2
    // where N_j is occupation number for subshell n_jl_j.
    // If we do this correctly we should have 4pi \int rho(r) r^2 dr = N_occ.

    // Normalize radial functions. We save the resulting normSquardes in the relevant SubShell in AtomParameters.shells.
    // We need to use this when calculating the charge density with properly normalised radial functions!!
    CalculateRadialFunctionNorms(Bsplines, gevp_Eigenvectors, physicalIndices, AtomParameters,
                                 GaussLegendreIntegration);

    // Integrate charge density. This is for sanity checking.
    IntegrateChargeDensity(Bsplines, gevp_Eigenvectors, physicalIndices, AtomParameters, GaussLegendreIntegration);

    // Now we solve the poisson equation for direct electron-electron potential due to this charge density
    // Nabla^2 V_dir = 1/r del^2 phi = -4\pi rho, where have used V_dir = phi(r)/r.
    // The equation we solve then is del^2 phi(r) = -r4\pi rho(r).

    // Here we use collocation method with a different bspline order, but with a lot more knotpts to get
    // a good representation of the potential.
    constexpr u32 poissonBsplineOrder = 4;
    constexpr u32 poissonNumKnotPoints = 16; // So many knotpoints for this?
    ASSERT(numGridPoints > numKnotPoints);
    Atom::Bsplines poissonBsplines = Atom::Bsplines(poissonNumKnotPoints, poissonBsplineOrder);
    poissonBsplines.setupKnotPoints(Grid.getGridPoints(), Atom::knotSequenceType::Linear);
    writeKnotPointsToFile(poissonBsplines.m_knotPoints);



    // For the poisson equation we say that at the origin everyting is zero, so we skip the first bspline;
    // the one with index zero.
    // We also say that the derivative of the solution (phi) is zero at the last knotpoint. The physical situation
    // is that the potential is approximately constant "far out", so the derivative vanishes.
    // A difference from the r = 0 zero condition, and the r -> infty zero condition from the eigenvalue problem,
    // is that we actually want this condition on the second derivative to be present in the solution. SO for this
    // endpoint boundary condition we can't just remove the coefficient frome the calculation, it is the
    // rhs that should be zero here.
    // We will be using all bsplines except the first one (index 0)
    std::vector<u32> poissonBsplineIndices;
    // Note that I am using zero indexing for Bspline, ie the first Bspline has index 0,
    // and the last Bspline has index numBsplines-1.
    for (u32 i = 1; i < poissonBsplines.numberOfBsplines(); i++) {
        poissonBsplineIndices.push_back(i);
    }
    poissonBsplines.SetBoundaryConditionBsplineIndices(poissonBsplineIndices);
    TestBsplineDerivatives(poissonBsplines, Grid);

    // We have N-k unknowns (the c_n coeffiecients, one for each Bspline), but only really N-k-2 equations.
    // If we use that c_0 = 0 we have N-k-1 unknowns, but still not enough equations. We can add an equation then for
    // the last condition, so that the system size that we plug into LAPACK solver routines is N-k-1.
    u32 poissonDimension = poissonBsplines.m_numKnotPoints - poissonBsplines.m_order - 1;
    ZMatrix poissonLHS = ZMatrix(poissonDimension, poissonDimension);
    setupPoissonLHS(poissonBsplines, poissonLHS, GaussLegendreIntegration);
    Logger::Trace("Dimension for poisson equation LHS matrix: %i", poissonDimension);
    for (u32 i = 0; i < poissonDimension; i++) {
        for (u32 j = 0; j < poissonDimension; j++) {
            printf("%f ", poissonLHS(i, j).real());
        }
        printf("\n");
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
        i32 bsplineIndex_j = bsplineIndices[j];
        i32 bsplineIndex_i = bsplineIndices[i];

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
        f64 lfactor = l * (l + 1.0);
        Complex return_value = 0.5 * term1;
        if (lfactor > 0.0) {
            return_value += (0.5 * lfactor) * term2;
        }
        return_value += (-1.0 * Z) * term3;
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

Complex ChargeDensityAtCoordinate(Complex r, Atom::Bsplines &Bsplines,
                                  const ZMatrix &coefficientMatrix, const std::vector<u32> &physicalIndices,
                                  const Atom::SystemParameters &params) {

    auto prefactor = g_OneOverFourPi; // 1/4pi

    Complex rho_at_r = Complex(0.0);
    if (std::abs(r) < 1e-8) {
        return rho_at_r; // Return zero if we're close to the origin.
    }

    // Loop over sub shell n_j l_j
    u32 j = 0;
    for (auto &shell : params.shells) {
        u32 N_j = shell.NumberOfElectrons; // Occupation number
        // Get the index into eigenvector coefficient matrix that correspond to coefficients
        // for the P_njlj (this shell's) radial function
        u32 coefficientIndex = physicalIndices[j];
        auto P_nl = GetRadialFunctionAtCoordinate(r, Bsplines,
                                                  coefficientIndex,
                                                  coefficientMatrix); // Radial function at coordinate r
        auto P_normSquared = shell.radialFunctionIntegratedOverAllSpace;
        P_nl = P_nl / sqrt(P_normSquared);
        rho_at_r += prefactor * N_j * P_nl * P_nl / (r * r);
        ++j;
    }
    return rho_at_r;
};

Complex GetRadialFunctionAtCoordinate(Complex r, Atom::Bsplines &Bsplines, u32 functionIndex,
                                      const ZMatrix &coefficients) {

    auto return_value = Complex(0.0);
    // Loop over all bspline coefficients and create superposition
    // P(r) = sum_i c_i B_i(r)
    ASSERT(coefficients.numRows() == Bsplines.m_usedBsplineIndices.size());
    for (u32 i = 0; i < coefficients.numRows(); ++i) {
        u32 bsplineIndex = Bsplines.m_usedBsplineIndices[i];
        auto B_i = Bsplines.GetBsplineAtCoordinate(r, bsplineIndex);
        auto c_i = coefficients(i, functionIndex);
        return_value += c_i * B_i;
    }

    return return_value;
};

void CalculateRadialFunctionNorms(Atom::Bsplines &Bsplines,
                                  const ZMatrix &coefficientMatrix, const std::vector<u32> &physicalIndices,
                                  Atom::SystemParameters &params, GaussLegendre::Integration &GLI) {


    // Loop over sub shell n_j l_j
    u32 j = 0;
    for (auto &shell : params.shells) {
        // Get the index into eigenvector coefficient matrix that correspond to coefficients
        // for the P_njlj (this shell's) radial function
        u32 coefficientIndex = physicalIndices[j];
        Complex norm = Complex(0.0);
        for (u32 i = 0; i < Bsplines.m_numBsplines - 1; i++) {
            // Loop here is so that the last point b will be at i+1, and at the last physical point;
            // Shift so we don't integrate in ghost points
            u32 knotIndex = i + Bsplines.m_order - 1;

            Complex a = Bsplines.m_knotPoints[knotIndex];
            Complex b = Bsplines.m_knotPoints[knotIndex + 1];

            u32 numIntegrationPoints = Bsplines.m_order;
            auto abscissae = GLI.getShiftedAbscissae(a, b, numIntegrationPoints);
            ASSERT(abscissae.size() == numIntegrationPoints);
            auto prefactor = GLI.b_minus_a_half(a, b);

            if (std::abs(prefactor) < 1e-8) {
                continue;
            }

            auto pWeights = GLI.getPointerToZWeights(numIntegrationPoints);

            auto term1 = Complex(0.0, 0.0);
            for (int m = 0; m < numIntegrationPoints; ++m) {
                auto r = abscissae[m];
                auto f = GetRadialFunctionAtCoordinate(r, Bsplines,
                                                       coefficientIndex,
                                                       coefficientMatrix); // Radial function at coordinate r
                term1 += prefactor * pWeights[m] * f * f;
            }

            norm += term1;
        }
        shell.radialFunctionIntegratedOverAllSpace = std::abs(norm);
        ++j;
    }
};

void IntegrateChargeDensity(Atom::Bsplines &Bsplines, ZMatrix &coeffs, std::vector<u32> &physicalIndices,
                            Atom::SystemParameters AtomParameters, GaussLegendre::Integration &GLI) {

    auto ChargeDensityIntegratedOverAllSpace = Complex(0.0);
    for (u32 i = 0; i < Bsplines.m_numBsplines - 1; i++) {
        // Loop here is so that the last point b will be at i+1, and at the last physical point;
        // Shift so we don't integrate in ghost points
        u32 knotIndex = i + Bsplines.m_order - 1;

        Complex a = Bsplines.m_knotPoints[knotIndex];
        Complex b = Bsplines.m_knotPoints[knotIndex + 1];

        u32 numIntegrationPoints = Bsplines.m_order;
        auto abscissae = GLI.getShiftedAbscissae(a, b, numIntegrationPoints);
        ASSERT(abscissae.size() == numIntegrationPoints);
        auto prefactor = GLI.b_minus_a_half(a, b);
        if (std::abs(prefactor) < 1e-8) {
            continue;
        }
        auto pWeights = GLI.getPointerToZWeights(numIntegrationPoints);

        auto term1 = Complex(0.0, 0.0);
        for (int m = 0; m < numIntegrationPoints; ++m) {
            auto r = abscissae[m];
            auto r2 = r * r;
            auto f = r2 * ChargeDensityAtCoordinate(r, Bsplines, coeffs, physicalIndices, AtomParameters);
            term1 += prefactor * pWeights[m] * f;
        }
        ChargeDensityIntegratedOverAllSpace += g_FourPi * term1;
    }

    Logger::Trace("4pi integral over all r rho(r)*r^2 dr = (%f, %f). For %s this should equal %i",
                  ChargeDensityIntegratedOverAllSpace.real(), ChargeDensityIntegratedOverAllSpace.imag(),
                  AtomParameters.Name.c_str(), AtomParameters.TotalOccupationNumber);
}

void setupPoissonLHS(Atom::Bsplines &Bsplines, ZMatrix &lhs, GaussLegendre::Integration &GLI) {
    u32 numCols = lhs.numCols();
    u32 numRows = lhs.numRows();
    u32 k = Bsplines.m_order;

    // First we loop over all rows that can be fully "banded", this excludes the first and the last row.
    ASSERT(Bsplines.m_usedBsplineIndices.size() > 0);
    for (u32 i = 1; i < numRows-1; i++) {
        u32 knotPointIndex = i + k - 1; // The first physical point has index 3 in this array, for k = 4.
        // But in this loop we are actually in the second physical knotpoint.
        auto r = Bsplines.m_knotPoints[knotPointIndex];
        u32 bsplineIndex = Bsplines.m_usedBsplineIndices[i-1]; // we take care of the matrix indexing relative to Bspline numbering here.
        lhs(i, i-1) = Bsplines.GetSecondDerivativeAtCoordinate(r, bsplineIndex);//bsplineIndex;
        lhs(i, i) = Bsplines.GetSecondDerivativeAtCoordinate(r, bsplineIndex + 1);//bsplineIndex + 1;
        lhs(i, i+1) = Bsplines.GetSecondDerivativeAtCoordinate(r, bsplineIndex + 2);
    }
    // Then we fix the first row:
    auto r_first = Bsplines.m_knotPoints[k-1];
    lhs(0,0) = Bsplines.GetSecondDerivativeAtCoordinate(r_first, 1);
    lhs(0,1) = Bsplines.GetSecondDerivativeAtCoordinate(r_first, 2);

    // also the last row. Note here how we are using the _first_ derivative here.
    auto r_last = Bsplines.m_knotPoints[Bsplines.m_numKnotPoints -(k-1)-1];
    u32 lastBsplineIndex = Bsplines.m_usedBsplineIndices[Bsplines.m_usedBsplineIndices.size()-1];
    Logger::Trace("lastBsplineIndex: %i",lastBsplineIndex);
    Logger::Trace("numBsplines (poisson):, %i", Bsplines.numberOfBsplines());
    lhs(numRows-1, numCols-2) = Bsplines.GetDerivativeAtCoordinate(r_last, lastBsplineIndex-1);//Bsplines.GetDerivativeAtCoordinate(r_last, lastBsplineIndex-1);
    lhs(numRows-1, numCols-1) = Bsplines.GetDerivativeAtCoordinate(r_last, lastBsplineIndex);//Bsplines.GetDerivativeAtCoordinate(r_last, lastBsplineIndex);
}
