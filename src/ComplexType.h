#pragma once
#include "typedefs.h"

struct Complex {
    Real Re;
    Real Im;
    Complex();
    Complex(Real a);
    Complex(Real a, Real b);
};


Complex operator+(const Complex &lhs, const Complex &rhs);

Complex operator-(const Complex &lhs, const Complex &rhs);

Complex operator*(const Complex &lhs, const Complex &rhs);

Complex operator/(const Complex &lhs, const Complex &rhs);