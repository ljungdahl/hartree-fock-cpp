#include "Parameters.h"

// TODO(anton): These functions should read parameters from input files.

Atom::GridParameters Atom::GetGridParameters() {

    constexpr u32 numGridPoints = 1000;
    constexpr f64 gridStart = 0.0, gridEnd = 20.0;

    GridParameters params;
    params.numberOfGridPoints = numGridPoints;
    params.gridStart = gridStart;
    params.gridEnd = gridEnd;
    params.isInitialised = true;

    return params;
}

Atom::BsplineParameters Atom::GetBsplineParameters() {

    constexpr u32 bsplineOrder = 6;
    constexpr u32 numKnotPoints = 60;

    u32 numberOfBsplines = numKnotPoints - bsplineOrder;
    ASSERT(numberOfBsplines >= bsplineOrder);

    BsplineParameters params;
    params.bsplineOrder = bsplineOrder;
    params.numberOfKnotPoints = numKnotPoints;
    params.isInitialised = true;

    return params;
}

Atom::SystemParameters Atom::GetAtomicSystemParameters() {
    // We use Hartree atomic units. hbar = 1, electron charge e = 1, bohr radius a_0 = 1, electron mass m_e = 1.
    // We also have that 4\pi \epsilon_0 = 1,
    // since the Bohr radius a_0 = (4\pi epsilon_0 hbar^2)/(m_e e^2) = 1, and hbar = m_e = e = 1 in atomic units.
    // Then the radial, hydrogen-like, Hamiltonian H' is
    // H' = (-0.5*d_r^2 + 0.5*l(l+1)/r^2 - Z/r).
    constexpr u32 Z = 1; // HYdrogen

    Atom::SubShell _1s;
    _1s.n = 1;
    _1s.l = Atom::AngularMomentumQuantumNumber::s; // l = 0,1,2,3... <-> s,p,d,f,...;
    _1s.numberOfElectrons = 1;

    Atom::SystemParameters AtomParameters;
    AtomParameters.shells.push_back(_1s);
    AtomParameters.Z = Z;
    AtomParameters.name = "Hydrogen";
    for (auto shell : AtomParameters.shells) {
        AtomParameters.totalOccupationNumber += shell.numberOfElectrons;
    }

    AtomParameters.isInitialised = true;

    return AtomParameters;
}