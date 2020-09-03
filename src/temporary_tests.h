#pragma once

#include "MKL_utility_linear_algebra.h"

void testMatrixInversion(const Matrix& B) {

//    u32 n = 2, m = 2;
//    Matrix B = Matrix(m, n);
//
//    B(0, 0) = Complex(5.0);
//    B(0, 1) = Complex(-2.0);
//    B(1, 0) = Complex(-2.0);
//    B(1, 1) = Complex(1.0);

//    for (u32 i = 0; i < B.numRows(); i++) {
////        for (u32 j = 0; j < B.numCols(); j++) {
////            B(i,j) = Complex((f64)i, (f64)j);
////        }
//        B(i, i) = Complex((f64) i);
//    }

//    B(1, 2) = Complex(3.0);
//
//    for (u32 i = 0; i < B.numRows(); i++) {
//        for (u32 j = 0; j < B.numCols(); j++) {
//            auto val = B(i, j);
//            printf("(%f, %f) ", val.real(), val.imag());
//        }
//        printf("\n");
//    }

    Matrix B_inverse = Matrix(B.numRows(), B.numCols());
//    B_inverse(0,0) = Complex(1.0);
//    B_inverse(0,1) = Complex(2.0);
//    B_inverse(1,0) = Complex(2.0);
//    B_inverse(1,1) = Complex(5.0);

    Logger::Trace("calling setFromMatrix");
    B_inverse.setFromMatrix(B);
//    Logger::Trace("print B_inverse after copy from B:");
//    for (u32 i = 0; i < B_inverse.numRows(); i++) {
//        for (u32 j = 0; j < B_inverse.numCols(); j++) {
//            auto val = B_inverse(i, j);
//            printf("(%f, %f) ", val.real(), val.imag());
//        }
//        printf("\n");
//    }

    Logger::Trace("Calling InvertInPlace");
    LAPACK::InvertMatrixInPlace(B_inverse);
//    for (u32 i = 0; i < B_inverse.numRows(); i++) {
//        for (u32 j = 0; j < B_inverse.numCols(); j++) {
//            auto val = B_inverse(i, j);
//            printf("(%f, %f) ", val.real(), val.imag());
//        }
//        printf("\n");
//    }

    Logger::Trace("Calling MatMul");
    Matrix I = MatMul(B_inverse, B);

    Logger::Trace("B^-1 * B: ");
    for (u32 i = 0; i < I.numRows(); i++) {
        for (u32 j = 0; j < I.numCols(); j++) {
            auto val = I(i, j);
//            printf("(%f, %f)", val.real(), val.imag());
            printf("%f ", val.real());
        }
        printf("\n");
    }

    Matrix I2 = MatMul(B, B_inverse);
    Logger::Trace("B * B^-1: ");
    for (u32 i = 0; i < I2.numRows(); i++) {
        for (u32 j = 0; j < I2.numCols(); j++) {
            auto val = I2(i, j);
//            printf("(%f, %f)", val.real(), val.imag());
            printf("%f ", val.real());
        }
        printf("\n");
    }
}