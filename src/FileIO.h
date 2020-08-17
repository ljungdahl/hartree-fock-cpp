#pragma once
#include <vector>
#include <string>

#include "typedefs.h"

namespace FileIO {
    void writeRowColDataToFile(std::vector<std::vector<f64>> out, std::string fileName);
}