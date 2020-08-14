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
    std::vector<std::vector<f64>> knotPointsOutput;
    for(auto knotPt : Bsplines.m_knotPoints) {
        std::vector<f64> row;
        row.push_back(knotPt.real());
        knotPointsOutput.push_back(row);
    }
    writeDataToFile(knotPointsOutput, "../knotpoints.dat");

    std::vector<std::vector<f64>> outputData;
    for (auto val : Grid.getGridPoints()) {
        std::vector<f64> rowVector;
        rowVector.push_back(val.real());
        for (u32 bsplineIndex = 0; bsplineIndex < Bsplines.m_numBsplines; bsplineIndex++) {
            auto result = Bsplines.GetBsplineAtCoordinate(val, bsplineIndex);
//            auto result = Bsplines.GetBsplineAtCoordinate(val, /*bsplineIndex*/0);
            rowVector.push_back(result.real());
        }
        outputData.push_back(rowVector);
    }

    writeDataToFile(outputData, "../bsplines.dat");
    return 0;
}