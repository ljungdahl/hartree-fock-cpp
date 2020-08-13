#include <cmath>
#include "custom_asserts.h"
#include "Bsplines.h"

void Bsplines::bsplvb_Complex(const Complex *knotPoints, const u32 bsplineOrder,
                              const Complex x, const u32 left,
                              Complex *outputArray) {
    Complex saved;
    Complex term;

    // Assume index == 1, ie starting order is j = 0;
    outputArray[0] = Complex(1.0, 0.0);

    u32 k = bsplineOrder;
    auto *t = knotPoints;
    auto *dR = new Complex[k];
    auto *dL = new Complex[k];

    for (u32 j = 0; j < k - 1; ++j) {
        dR[j] = t[left + j + 1] - x;
        dL[j] = x - t[left - j];

//        Logger::Trace("dR[%i]: %f, dL[%i]: %f", j, dR[j].Re, j, dL[j].Re);
        saved = Complex(0.0);

        for (int i = 0; i <= j; i++) {
            term = outputArray[i] / (dR[i] + dL[j - i]);
//            Logger::Trace("term: %f", term.Re);
            outputArray[i] = saved + dR[i] * term;
            saved = dL[j - i] * term;
        }

//        Logger::Trace("saved: %f", saved.Re);
        outputArray[j + 1] = saved;
    }

    delete[] dR;
    delete[] dL;
}

Complex Bsplines::GetValueAtCoordinate(const Complex x_value, const Complex *knotPoints,
                                       const u32 numKnotPoints, const u32 bsplineOrder,
                                       const u32 bsplineIndex, const u32 numberOfBsplines) {
    f64 difference_tolerance = 1e-8;

    if (x_value.Re < knotPoints[0].Re) {
//        Logger::Trace("x_val is less than the first grid point. Return Complex(0.0, 0.0)");
        return Complex(0.0, 0.0);
    }

    if (x_value.Re > knotPoints[numKnotPoints - 1].Re) {
//        Logger::Trace("x_val is greater than the last grid point. Return Complex(0.0, 0.0)");
        return Complex(0.0, 0.0);
    }

//    if (bsplineIndex == 0) {
//        // We're looking at the first bspline
//        if (abs(x_value.Re - knotPoints[0].Re) < difference_tolerance) {
////            Logger::Trace("Difference too small between x_val and last knot point. Return Complex(1.0, 0.0)");
//            return Complex(1.0, 0.0);
//        }
//    }

    if (bsplineIndex == numberOfBsplines - 1) {
        // We're looking at the last bspline
        if (abs(x_value.Re - knotPoints[numKnotPoints - 1].Re) < difference_tolerance) {
//            Logger::Trace("Difference too small between x_val and last knot point. Return Complex(1.0, 0.0)");
            return Complex(1.0, 0.0);
        }
    }

    u32 left_knotPoint_index = 0;
    for (int t = 0; t < numKnotPoints; t++) {
        if (x_value.Re >= knotPoints[t].Re) {

            left_knotPoint_index = t;
        }
    }

    Logger::Trace("left_knotPoint_index %i", left_knotPoint_index);
    u32 accessIndex = bsplineIndex - left_knotPoint_index + bsplineOrder;
    if (
            accessIndex < 0 || accessIndex > bsplineOrder-1) {
//        Logger::Trace("Index condition not fullfilled for:");
//        Logger::Trace("indexCondition: %i, bsplineIndex: %i, left_knotPoint_index: %i, bsplineOrder: %i",
//                      indexCondition, bsplineIndex, left_knotPoint_index, bsplineOrder);
        return Complex(0.0, 0.0);
    }

//    Logger::Trace("Calculating Bspline value for");
//    Logger::Trace("bsplineIndex: %i, left_knotPoint_index: %i, bsplineOrder: %i",
//                  bsplineIndex, left_knotPoint_index, bsplineOrder);

    auto *Sp = new Complex[bsplineOrder]; // TODO(anton): Explain this.
    Bsplines::bsplvb_Complex(knotPoints, bsplineOrder, x_value, left_knotPoint_index, Sp);

    Complex ret = Sp[accessIndex];
    delete[] Sp;
    return ret;
}