#include <cstdio>
#include <iostream>

#include "Grid.h"
#include "Bsplines.h"

template<class T>
void writeDataToFile(std::vector<std::vector<T>> out, std::string fileName) {
    FILE *pFile;
    pFile = fopen(fileName.c_str(), "w");

    u32 numRows = out.size();
    u32 numCols = out[0].size();

    for (u32 i = 0; i < numRows; i++ ) {

        for (u32 j = 0; j < numCols; j++) {
            T value = out[i][j];
            fprintf(pFile, "%E ", value);
        }

        fprintf(pFile, "\n");
    }

    fclose(pFile);
    Logger::Trace("wrote %s to disk.", fileName.c_str());
}

int main() {

    constexpr u32 NUMBER_OF_KNOTPOINTS = 20;
    constexpr u32 BSPLINE_ORDER = 4;
    constexpr u32 GRID_SIZE = 1000;
    constexpr f64 GRID_START = 0.0, GRID_END = 1.0;

    auto Grid = Atom::Grid(GRID_SIZE, GRID_START, GRID_END);
    auto Bsplines = Atom::Bsplines(NUMBER_OF_KNOTPOINTS, BSPLINE_ORDER);

    Bsplines.setupKnotPoints(Grid.getGridPoints());

    std::vector<std::vector<f64>> outputData;
    for (auto val : Grid.getGridPoints()) {
        std::vector<f64> rowVector;
        rowVector.push_back(val.real());
        for (u32 bsplineIndex = 0; bsplineIndex < Bsplines.m_numBsplines; bsplineIndex++) {
            auto result = Bsplines.GetBsplineAtCoordinate(val, bsplineIndex);
            rowVector.push_back(result.real());
        }
        outputData.push_back(rowVector);
    }

    writeDataToFile(outputData, "../bsplines.dat");
    return 0;
}
//
////    Logger::Trace("entered main");
//    Grid radialGrid;
//
//    // Linspace grid
//    {
//        for (int i = 0; i < radialGrid.size; i++) {
//            f64 r = i * radialGrid.stepSize;
//            radialGrid.values.push_back(r);
//        }
//    }
////
////    Logger::Trace("%f", radialGrid.values[0]);
////    Logger::Trace("%f", radialGrid.values[radialGrid.values.size()-1]);
////    Logger::Trace("%i", radialGrid.values.size());
//
//    std::vector<Complex> knotPoints;
//    u32 bsplineOrder = 4;
//    u32 numberOfGhostPointsInEachEnd = (bsplineOrder - 1);
//
//    u32 numberOfPhysicalKnotPoints = 11;
//    u32 numberOfKnotPoints = numberOfPhysicalKnotPoints + 2 * numberOfGhostPointsInEachEnd;
//    u32 numberOfBsplines = numberOfKnotPoints - bsplineOrder;
//    Logger::Trace("Number of Bsplines: %i", numberOfBsplines);
//
//    knotPoints.resize(numberOfKnotPoints, Complex(0.0, 0.0));
//
//    // linear knot sequence
//    {
//        // set ghost knot points
//        for (int i = 0; i < numberOfGhostPointsInEachEnd; ++i) {
//            knotPoints[i] = Complex(radialGrid.start, 0.0);
//            u32 endIndex = (knotPoints.size() - i) - 1;
//            knotPoints[endIndex] = Complex(radialGrid.end, 0.0);
//        }
//
//        // set real knot points
//        u32 offset = numberOfGhostPointsInEachEnd;
//        u32 stride = (u32) (1000 / numberOfPhysicalKnotPoints);
//        for (int i = 0; i < numberOfPhysicalKnotPoints; ++i) {
//            u32 knotPointIndex = i + offset;
//            u32 gridIndex = i * stride;
//            ASSERT(gridIndex <= radialGrid.size - 1);
//            Complex knotPointValue = Complex(radialGrid.values[gridIndex], 0.0);
//            knotPoints[knotPointIndex] = knotPointValue;
////            Logger::Trace("%i %i %f", knotPointIndex, gridIndex, knotPointValue);
//        }
//    }
//
//    Logger::Trace("knotpoints[0].Re: %f", knotPoints[0].Re);
//    Logger::Trace("radialGrid.values[0]: %f", radialGrid.values[0]);
//
//    std::vector<std::vector<f64>> outputData;
////    Logger::Trace("radialGrid.values.size(): %i", radialGrid.values.size());
//    outputData.push_back(radialGrid.values);
//    // NOTE(anton): Bspline indices from zero. The first Bspline has index 0, and the last has index numBsplines-1.
//    for (u32 bsplineIndex = 0; bsplineIndex < numberOfBsplines; bsplineIndex++) {
//        std::vector<f64> bsplineValues;
//        for (auto val : radialGrid.values) {
//            Complex x_value = Complex(val, 0.0);
//
//            auto result = Bsplines::GetValueAtCoordinate(x_value, knotPoints.data(), knotPoints.size(), bsplineOrder,
//                                                         bsplineIndex, numberOfBsplines);
////        Logger::Trace("%i, result.Re: %f at x: %f", idx, result.Re, x_value.Re);
//            bsplineValues.push_back(result.Re);
//        }
//        outputData.push_back(bsplineValues);
//    }
//
//    Logger::Trace("Done with bspline calc.");
//    std::string dataFileName = "../bsplines.dat";
//    writeDataToFile(outputData, dataFileName);
//    return 0;
//}
//
//
//
