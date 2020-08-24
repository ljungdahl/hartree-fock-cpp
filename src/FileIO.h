#pragma once
#include <vector>
#include <string>

#include "typedefs.h"

namespace FileIO {
void writeComplexVectorToFile(std::vector<Complex> out, std::string fileName);

void writeRowColDataToFile(std::vector<std::vector<f64>> out, std::string fileName);

void writeColDataToFile(std::vector<f64> out, std::string fileName);
}