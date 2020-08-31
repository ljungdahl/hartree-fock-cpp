/* NOTE(anton):
 * This is a test program that solves the Hartree-Fock equations for some nice test problem, probably Neon?
 * The numerical method is using Bsplines, and everything is written around that. The aim is to make some generalisations
 * so that different numerical methods might be used in the future.
 */

#include "typedefs.h"
#include "logger.h"
#include "custom_asserts.h"

#include "Vector.h"
#include "Grid.h"
#include "Bsplines.h"
#include "HartreeFock.h"

typedef LA::Vector<Complex> Vector;


int main(int argc, char* argv[]) {
// Atomic units. Let's do a 10 Bohr radii sized grid.
    constexpr u32 numGridPoints = 1000;
    constexpr f64 gridStart = 0.0, gridEnd = 10.0;
    Atom::Grid Grid = Atom::Grid(numGridPoints, gridStart, gridEnd);

    constexpr u32 bsplineOrder = 6;
    constexpr u32 numKnotPoints = 80;
    ASSERT(numGridPoints > numKnotPoints);
    Atom::Bsplines Bsplines = Atom::Bsplines(numKnotPoints, bsplineOrder);

    Bsplines.setupKnotPoints(Grid.getGridPoints(), Atom::knotSequenceType::Linear);

    // TODO(anton): The hf instantiation should take parameters so that
    // arrays with all relevant data can be allocated and initialised.
    Atom::HartreeFock hf = Atom::HartreeFock();


    // The starting point for the HF solver will be just the independent particle
    // model without any electron-electron interaction.
    hf.PerformInitialStep();

    hf.SelfConsistentSolution();

    return 0;
}