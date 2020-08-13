#include "Grid.h"

Atom::Grid::Grid(u32 size_, u32 start_, u32 end_)
        : m_size(size_),
          m_start(start_),
          m_end(end_)
          {

    m_gridPoints.resize(m_size, Complex(0.0));
    SetLinspace();
}

void Atom::Grid::SetLinspace() {
    f64 start = m_start.real(), end = m_end.real();
    f64 stepSize = (end-start) / (f64) (m_size-1);

    u32 pointIndex = 0;
    for (auto& point : m_gridPoints) {
        f64 r = pointIndex * stepSize;
        point = r;
        ++pointIndex;
    }
}

void Atom::Grid::Log() {
    u32 increment = 0;
    std::string logString = m_logPrefix+"gridPoint[%i]: %E";

    for (auto point : m_gridPoints) {
        Logger::Log(logString.c_str(), increment, point.real());
        increment += 1;
    }
}

const std::vector<Complex>& Atom::Grid::getGridPoints() {
    return m_gridPoints;
}