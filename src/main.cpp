#include <cstdio>
#include <vector>
#include <string>
#include <array>

#include "custom_asserts.h"
#include "logger.h"
#include "ComplexType.h"
#include "Bsplines.h"

struct Grid {
    u32 size = 1000;
    f64 start = 0.0, end = 1.0;
    f64 stepSize = (end - start) / (f64) (size - 1);
    std::vector<f64> values;
};

template<class T>
void writeDataToFile(std::vector<std::vector<T>> out, std::string fileName) {
    FILE *pFile;
    pFile = fopen(fileName.c_str(), "w");

    u32 numRows = out[0].size();
    u32 numCols = out.size();

    for (u32 i = 0; i < numRows; i++ ) {
        for (u32 j = 0; j < numCols; j++) {
            T value = out[j][i];
            fprintf(pFile, "%E ", value);
        }
        fprintf(pFile, "\n");
    }

    fclose(pFile);
    Logger::Trace("wrote %s to disk.", fileName.c_str());
}

int main() {

//    Logger::Trace("entered main");
    Grid radialGrid;

    // Linspace grid
    {
        for (int i = 0; i < radialGrid.size; i++) {
            f64 r = i * radialGrid.stepSize;
            radialGrid.values.push_back(r);
        }
    }
//
//    Logger::Trace("%f", radialGrid.values[0]);
//    Logger::Trace("%f", radialGrid.values[radialGrid.values.size()-1]);
//    Logger::Trace("%i", radialGrid.values.size());

    std::vector<Complex> knotPoints;
    u32 bsplineOrder = 4;
    u32 numberOfGhostPointsInEachEnd = (bsplineOrder - 1);

    u32 numberOfPhysicalKnotPoints = 11;
    u32 numberOfKnotPoints = numberOfPhysicalKnotPoints + 2 * numberOfGhostPointsInEachEnd;
    u32 numberOfBsplines = numberOfKnotPoints - bsplineOrder;
    Logger::Trace("Number of Bsplines: %i", numberOfBsplines);

    knotPoints.resize(numberOfKnotPoints, Complex(0.0, 0.0));

    // linear knot sequence
    {
        // set ghost knot points
        for (int i = 0; i < numberOfGhostPointsInEachEnd; ++i) {
            knotPoints[i] = Complex(radialGrid.start, 0.0);
            u32 endIndex = (knotPoints.size() - i) - 1;
            knotPoints[endIndex] = Complex(radialGrid.end, 0.0);
        }

        // set real knot points
        u32 offset = numberOfGhostPointsInEachEnd;
        u32 stride = (u32) (1000 / numberOfPhysicalKnotPoints);
        for (int i = 0; i < numberOfPhysicalKnotPoints; ++i) {
            u32 knotPointIndex = i + offset;
            u32 gridIndex = i * stride;
            ASSERT(gridIndex <= radialGrid.size - 1);
            Complex knotPointValue = Complex(radialGrid.values[gridIndex], 0.0);
            knotPoints[knotPointIndex] = knotPointValue;
//            Logger::Trace("%i %i %f", knotPointIndex, gridIndex, knotPointValue);
        }
    }

    Logger::Trace("knotpoints[0].Re: %f", knotPoints[0].Re);
    Logger::Trace("radialGrid.values[0]: %f", radialGrid.values[0]);

    std::vector<std::vector<f64>> outputData;
//    Logger::Trace("radialGrid.values.size(): %i", radialGrid.values.size());
    outputData.push_back(radialGrid.values);
    // NOTE(anton): Bspline indices from zero. The first Bspline has index 0, and the last has index numBsplines-1.
    for (u32 bsplineIndex = 0; bsplineIndex < numberOfBsplines; bsplineIndex++) {
        std::vector<f64> bsplineValues;
        for (auto val : radialGrid.values) {
            Complex x_value = Complex(val, 0.0);

            auto result = Bsplines::GetValueAtCoordinate(x_value, knotPoints.data(), knotPoints.size(), bsplineOrder,
                                                         bsplineIndex, numberOfBsplines);
//        Logger::Trace("%i, result.Re: %f at x: %f", idx, result.Re, x_value.Re);
            bsplineValues.push_back(result.Re);
        }
        outputData.push_back(bsplineValues);
    }

    Logger::Trace("Done with bspline calc.");
    std::string dataFileName = "../bsplines.dat";
    writeDataToFile(outputData, dataFileName);
    return 0;
}



