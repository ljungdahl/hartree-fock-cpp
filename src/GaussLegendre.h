#pragma once
#include <vector>
#include <array>

#include "typedefs.h"

namespace GaussLegendre {
  class Integration {
  public:
      Integration();
      const Complex* getPointerToZAbscissae(u32 points);
      const Complex* getPointerToZWeights(u32 points);
      std::vector<Complex> getShiftedAbscissae(Complex a, Complex b, u32 points);
      Complex b_minus_a_half(Complex a, Complex b);
  private:
      std::array<Complex, 6> m_ZAbscissaeSixPoints;
      std::array<Complex, 6> m_ZWeightsSixPoints;
      std::array<Complex, 4> m_ZAbscissaeFourPoints;
      std::array<Complex, 4> m_ZWeightsFourPoints;
  };
}