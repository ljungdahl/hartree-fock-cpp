#include "logger.h"
#include "HartreeFock.h"

Atom::HartreeFock::HartreeFock() {

}

void Atom::HartreeFock::PerformInitialStep() {
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




}

void Atom::HartreeFock::SelfConsistentSolution() {

}
