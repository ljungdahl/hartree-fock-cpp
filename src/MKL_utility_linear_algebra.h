#pragma once

/*
 * This header and the corresponding .cpp-file contains functions that uses LAPACK routines
 * not specifically tied to a certain type of solver.
 * For example we take care of matrix inversion here. This is used for the Hartree-Fock equations
 * to solve a standard eigenvalue problem instead of a generalised eigenvalue problem. Ie we
 * can use ?geev instead of ?ggev (?ggev performs matrix inversion anyway, so we only do it once).
 */
#include "MKL_include_and_define.h"
#include "linear_algebra_typedefs.h"

namespace LAPACK {

void InvertMatrixInPlace(Matrix &A);

}

Matrix MatMul(const Matrix &A, const Matrix &B);