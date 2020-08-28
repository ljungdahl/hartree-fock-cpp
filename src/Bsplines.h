#pragma once
#include <vector>

#include "logger.h"
#include "typedefs.h"

namespace Atom {

    enum knotSequenceType {
        Linear = 0,
        firstPointsCloserThenLinear = 1
    };

    class Bsplines {
    public:
        u32 m_numKnotPoints;
        u32 m_order;
        u32 m_numBsplines;
        std::vector<Complex> m_knotPoints;
        std::vector<u32> m_usedBsplineIndices;
    public:
        Bsplines(u32 numKnotPoints_, u32 bsplineOrder_);
        void setupKnotPoints(const std::vector<Complex> &gridPoints, knotSequenceType sequenceType);
        void DebugLogKnotSequence();
        u32 numberOfBsplines();
        Complex GetBsplineAtCoordinate(Complex coordinate, u32 bsplineIndex);
        Complex GetBsplineAtCoordinate(Complex coordinate, u32 bsplineIndex, u32 order);
        Complex GetDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex);
        Complex GetDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex, u32 order);
        Complex GetSecondDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex);
        Complex GetSecondDerivativeAtCoordinate(Complex coordinate, u32 bsplineIndex, u32 order);
        void SetBoundaryConditionBsplineIndices(std::vector<u32>& indices);

    private:
        std::vector<Complex> m_bsplvb_dL;
        std::vector<Complex> m_bsplvb_dR;
        std::vector<Complex> m_Sp;
    private:
        void bsplvb_Complex(Complex coordinate, u32 left, u32 k_order); // For derivative calculations
    };

}
















//namespace Bsplines {
//    void bsplvb_Complex(const Complex *knotPoints, u32 bsplineOrder,
//                        Complex x, u32 left,
//                        Complex *outputArray);
//
//    Complex GetValueAtCoordinate(Complex x_value, const Complex *knotPoints,
//                                 u32 numKnotPoints, u32 bsplineOrder,
//                                 u32 bsplineIndex,  u32 numberOfBsplines);
//}
