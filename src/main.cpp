#include <cstdio>
#include <iostream>

#include "typedefs.h"

#include "mkl.h"

#include "FileIO.h"
#include "Grid.h"
#include "Bsplines.h"

void writeToFile(Atom::Bsplines *Bsplines, Atom::Grid *Grid) {
    std::vector<std::vector<f64>> knotPointsOutput;

    for (auto knotPt : Bsplines->m_knotPoints) {
        std::vector<f64> row;
        row.push_back(knotPt.real());
        knotPointsOutput.push_back(row);
    }

    FileIO::writeRowColDataToFile(knotPointsOutput, "../knotpoints.dat");

    std::vector<std::vector<f64>> outputData;
    for (auto val : Grid->getGridPoints()) {
        std::vector<f64> rowVector;

        rowVector.push_back(val.real());

        for (u32 bsplineIndex = 0; bsplineIndex < Bsplines->m_numBsplines; bsplineIndex++) {
            auto result = Bsplines->GetBsplineAtCoordinate(val, bsplineIndex);
            rowVector.push_back(result.real());
        }

        outputData.push_back(rowVector);
    }

    FileIO::writeRowColDataToFile(outputData, "../bsplines.dat");
}

//void testLapack() {
//    Logger::Trace("Entered testLapack()");
//
//    Logger::Trace("Testing CBLAS zgemv y = A*x");
//
//    MKL_Complex16 rhs[3] = { Complex(1.0), Complex(2.0), Complex(3.0) };
//    MKL_Complex16 A[3][3];
//    for (int i = 0; i < 3; i++) {
//        A[i][i] = Complex(1.0, 0.0);
//    }
//
//    //A[0][0] = Complex(2.0);
//    //A[1][0] = Complex(1.0, 0.0);
//    //A[2][0] = Complex(1.0, 0.0);
//
//    CBLAS_LAYOUT Layout = CblasRowMajor;
//    CBLAS_TRANSPOSE trans = CblasNoTrans;
//
//    MKL_INT m = 3, n = m;
//    MKL_Complex16 alpha = Complex(1.0, 0.0);
//    void* beta = nullptr;
//    MKL_INT lda = n, incy = 1, incx = 1;
//
//    MKL_Complex16 result[3] = { Complex(1.0) };
//
//    Logger::Trace("Before lapack call");
//
//
//    cblas_zgemv(Layout, trans, m, n, &alpha, &A, lda, &rhs, incx, beta, &result, incy);
//    for (int i = 0; i < 3; i++) {
//
//        Logger::Trace("rhs[%i]: (%f, %f), result[%i]: (%f, %f)",
//                      i, rhs[i].real(), rhs[i].imag(), i, result[i].real(), result[i].imag()
//        );
//
//    }
//}

int main() {

    Logger::Trace("Entered main");

    MKLVersion ver;
    int len=198;
    char buf[198];

    MKL_Get_Version_String(buf, len);
    printf("\n%s\n",buf);
    printf("\n");

    MKL_Get_Version(&ver);
    printf("Major version:          %d\n",ver.MajorVersion);
    printf("Minor version:          %d\n",ver.MinorVersion);
    printf("Update version:         %d\n",ver.UpdateVersion);
    printf("Product status:         %s\n",ver.ProductStatus);
    printf("Build:                  %s\n",ver.Build);
    printf("Processor optimization: %s\n",ver.Processor);
    printf("================================================================\n");
    printf("\n");

//    Complex a[2] = { Complex(1.0, 0.0), Complex(2.0, 0.0) };
//    Complex b[2] = { Complex(1.0, 0.0), Complex(1.0, 0.0) };
//    Complex dotc;
    double a[2] = { 1.0, 2.0 };

    double result = cblas_ddot(2, a, 1, a, 1);
    Logger::Trace("result %f", result);
//    testLapack();

    return 0;
}