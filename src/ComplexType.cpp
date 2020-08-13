#include "typedefs.h"
#include "ComplexType.h"

Complex::Complex() {
    Re = 0.0;
    Im = 0.0;
}

Complex::Complex(Real a) {
    Re = a;
    Im = a;
}

Complex::Complex(Real a, Real b) {
    Re = a;
    Im = b;
}

Complex operator+(const Complex &lhs, const Complex &rhs) {
    Complex result = Complex(lhs.Re + rhs.Re, lhs.Im + rhs.Im);
    return result;
}

Complex operator-(const Complex &lhs, const Complex &rhs) {
    Complex result = Complex(lhs.Re - rhs.Re, lhs.Im - rhs.Im);
    return result;
}

Complex operator*(const Complex &lhs, const Complex &rhs) {
    Real ac_minus_bd = lhs.Re * rhs.Re - lhs.Im * rhs.Im;
    Real ad_plus_bc = lhs.Re * rhs.Im + lhs.Im*rhs.Re;
    Complex result = Complex(ac_minus_bd, ad_plus_bc);
    return result;
}

Complex operator/(const Complex &lhs, const Complex &rhs) {
    Real ac_plus_bd = lhs.Re * rhs.Re + lhs.Im * rhs.Im;
    Real bc_minus_ad = lhs.Im * rhs.Re - lhs.Re * rhs.Im;
    Real c2 = rhs.Re * rhs.Re;
    Real d2 = rhs.Im * rhs.Im;
    Complex result = Complex(ac_plus_bd / (c2 + d2), bc_minus_ad / (c2 + d2));
    return result;
}
