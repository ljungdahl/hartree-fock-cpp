#pragma once

#include "typedefs.h"
#include "linear_algebra_typedefs.h" // Matrix, Vector

namespace Atom {

class HartreeFock {
public:
    HartreeFock();

    void PerformInitialStep(Matrix &H, Matrix &B_inverse);
    void SelfConsistentSolution();

private:

};

}


