#include "logger.h"
#include "FileIO.h"

void FileIO::writeComplexVectorToFile(std::vector<Complex> out, std::string fileName) {
    FILE *pFile;
    pFile = fopen(fileName.c_str(), "w");

    u32 numCols = out.size();

    for (u32 j = 0; j < numCols; j++) {
        auto value = out[j];
        fprintf(pFile, "%E %E\n", value.real(), value.imag());
    }
    fclose(pFile);
    Logger::Log("wrote %s to disk.", fileName.c_str());
}

void FileIO::writeRowColDataToFile(std::vector<std::vector<f64>> out, std::string fileName) {
    FILE *pFile;
    pFile = fopen(fileName.c_str(), "w");

    u32 numRows = out.size();
    u32 numCols = out[0].size();

    for (u32 i = 0; i < numRows; i++) {

        for (u32 j = 0; j < numCols; j++) {
            f64 value = out[i][j];
            fprintf(pFile, "%E ", value);
        }
        fprintf(pFile, "\n");
    }
    fclose(pFile);
    Logger::Log("wrote %s to disk.", fileName.c_str());
}

void FileIO::writeColDataToFile(std::vector<f64> out, std::string fileName) {
    FILE *pFile;
    pFile = fopen(fileName.c_str(), "w");

    u32 numCols = out.size();

    for (u32 j = 0; j < numCols; j++) {
        f64 value = out[j];
        fprintf(pFile, "%E\n", value);
    }
    fclose(pFile);
    Logger::Log("wrote %s to disk.", fileName.c_str());
}