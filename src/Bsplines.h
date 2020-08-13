#pragma once
#include <vector>

#include "logger.h"
#include "typedefs.h"

namespace Atom {

    class Bsplines {
    public:
        u32 m_numKnotPoints;
        u32 m_order;
        u32 m_numBsplines;
        std::vector<Complex> m_knotPoints;
    public:
        Bsplines(u32 numKnotPoints_, u32 bsplineOrder_);
        void setupKnotPoints(const std::vector<Complex> &gridPoints);
        void LogKnotSequence();
        u32 numberOfBsplines();
        Complex GetBsplineAtCoordinate(Complex coordinate, u32 bsplineIndex);

    private:
        std::vector<Complex> m_bsplvb_dL;
        std::vector<Complex> m_bsplvb_dR;
        std::vector<Complex> m_Sp;
    private:
        void bsplvb_Complex(Complex coordinate, u32 left);
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
