#pragma once
#include <string>

#include "logger.h"
#include "typedefs.h"

namespace Atom {

    class Grid {
    public:

    public:
        Grid(u32 size_, u32 start_, u32 end_);
        void Log();
        const ZVector &getGridPoints() const;
    private:
        u32 m_size;
        Complex m_start;
        Complex m_end;
        ZVector m_gridPoints;
        std::string m_logPrefix = "[GRID]: ";
    private:
        void SetLinspace();
    };

}

