#pragma once
#include <complex>

#define MKL_INTERFACE_LAYER LP64

#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl.h"

// TODO(anton): Support lapack INFO return cases, currently I only allow completely successful calls.
#define LAPACK_CHECK(expr) { \
ASSERT(expr == 0); \
}

namespace LAPACK {
    enum Layout {
        // Ref mkl_lapacke.h
        RowMajor = LAPACK_ROW_MAJOR,
        ColMajor = LAPACK_COL_MAJOR
    };
}