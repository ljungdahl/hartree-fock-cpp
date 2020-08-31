/* NOTE(anton):
 * This is a test program that solves the Hartree-Fock equations for some nice test problem, probably Neon?
 * The numerical method is using Bsplines, and everything is written around that. The aim is to make some generalisations
 * so that different numerical methods might be used in the future.
 */

#include "typedefs.h"
#include "logger.h"
#include "Vector.h"

typedef LA::Vector<Complex> Vector;


int main(int argc, char* argv[]) {
    Vector vec = Vector(4);

    u32 i = 0;
    for (auto& el : vec.data()) {
        el = Complex(i);
        i++;
    }

    for (auto el : vec.data()) {
        Logger::Trace("(%f, %f)", el.real(), el.imag());
    }

    return 0;
}