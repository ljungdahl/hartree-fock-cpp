#pragma once

#include "typedefs.h"
#include "linear_algebra_typedefs.h" // Matrix, Vector

namespace LAPACK {
class Eigenproblems; // fwd declare.
}
namespace Atom {

class HartreeFock {
public:
    HartreeFock();

    void PerformInitialStep(Matrix &H, Matrix &B_inverse);
    void SelfConsistentSolution();

public:

    Matrix m_eigenvectors;
    Vector m_eigenvalues;

private:
    LAPACK::Eigenproblems* m_eig;

    void TestZGEEV();
};

}


