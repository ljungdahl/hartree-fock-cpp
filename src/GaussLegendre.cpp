//
// Created by anton on 2020-08-19.
//

#include "logger.h"
#include "GaussLegendre.h"

GaussLegendre::Integration::Integration() {
// We hardcode constants here for a few number of integration points using
// the webpage https://www.efunda.com/math/num_integration/findgausslegendre.cfm.

// TODO(anton): Write my own function to generate the weights and abscissae given a number of integration points.


// k = 4
    m_ZAbscissaeFourPoints = {
            Complex(-0.3399810435848563, 0.0),
            Complex(0.3399810435848563, 0.0),
            Complex(-0.8611363115940526, 0.0),
            Complex(0.8611363115940526, 0.0)
    };

    m_ZWeightsFourPoints = {
            Complex(0.6521451548625461, 0.0),
            Complex(0.6521451548625461, 0.0),
            Complex(0.3478548451374538, 0.0),
            Complex(0.3478548451374538, 0.0)
    };

    // k = 6
    m_ZAbscissaeSixPoints = {
            Complex(0.6612093864662645, 0.0),
            Complex(-0.6612093864662645, 0.0),
            Complex(-0.2386191860831969, 0.0),
            Complex(0.2386191860831969, 0.0),
            Complex(-0.932469514203152, 0.0),
            Complex(0.932469514203152, 0.0)
    };

    m_ZWeightsSixPoints = {
            Complex(0.3607615730481386, 0.0),
            Complex(0.3607615730481386, 0.0),
            Complex(0.4679139345726910, 0.0),
            Complex(0.4679139345726910, 0.0),
            Complex(0.1713244923791704, 0.0),
            Complex(0.1713244923791704, 0.0)
    };
}

const Complex *GaussLegendre::Integration::getPointerToZAbscissae(u32 points) {
    if (points == 6) {
        return m_ZAbscissaeSixPoints.data();
    }

    if (points == 4) {
        return m_ZAbscissaeFourPoints.data();
    }

    Logger::Warn("getPointerToZAbscissae(k = %i): returned nullptr!", points);
    return nullptr;
}

const Complex *GaussLegendre::Integration::getPointerToZWeights(u32 points) {
    if (points == 6) {
        return m_ZWeightsSixPoints.data();
    }

    if (points == 4) {
        return m_ZWeightsFourPoints.data();
    }
    Logger::Warn("getPointerToZWeights(k = %i): returned nullptr!", points);
    return nullptr;
}

Complex GaussLegendre::Integration::b_minus_a_half(Complex a, Complex b) {
    Complex return_value = 0.5 * (b - a);
    return return_value;
}

std::vector<Complex> GaussLegendre::Integration::getShiftedAbscissae(Complex a, Complex b, u32 points) {
    std::vector<Complex> shiftedAbscissae;
    Complex b_minus_a_half = 0.5 * (b - a);
    Complex b_plus_a_half = 0.5 * (b + a);

    if (points == 6) {
        for (auto abscissa : m_ZAbscissaeSixPoints) {
            shiftedAbscissae.push_back(abscissa * b_minus_a_half + b_plus_a_half);
        }
        return shiftedAbscissae;
    }

    if (points == 4) {
        for (auto abscissa : m_ZAbscissaeFourPoints) {
            shiftedAbscissae.push_back(abscissa * b_minus_a_half + b_plus_a_half);
        }
        return shiftedAbscissae;
    }

    Logger::Warn("getShiftedAbscissae(a = (%f,%f), b = (%f,%f), k = %i) returned with zero size vector.",
                 a.real(), a.imag(), b.real(), b.imag(), points);
    return shiftedAbscissae;
}
