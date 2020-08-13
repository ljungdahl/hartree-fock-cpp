#pragma once
#include <vector>
#include <string>

#include "logger.h"
#include "typedefs.h"

namespace Atom {

    class Grid {
    public:

    public:
        Grid(u32 size_, u32 start_, u32 end_);
        void Log();
        const std::vector<Complex> &getGridPoints();
    private:
        u32 m_size;
        Complex m_start;
        Complex m_end;
        std::vector<Complex> m_gridPoints;
        std::string m_logPrefix = "[GRID]: ";
    private:
        void SetLinspace();
    };

}

