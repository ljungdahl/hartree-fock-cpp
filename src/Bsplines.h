#pragma once

#include "logger.h"
#include "typedefs.h"
#include "ComplexType.h"

namespace Bsplines {
    void bsplvb_Complex(const Complex *knotPoints, const u32 bsplineOrder,
            const Complex x, const u32 left,
            Complex* outputArray);

    Complex GetValueAtCoordinate(const Complex x_value, const Complex *knotPoints,
                                 const u32 numKnotPoints, const u32 bsplineOrder,
                                 const u32 bsplineIndex, const u32 numberOfBsplines);
}
