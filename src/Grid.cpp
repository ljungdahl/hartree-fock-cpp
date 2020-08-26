#include "Grid.h"

Atom::Grid::Grid(u32 size_, u32 start_, u32 end_)
        : m_size(size_),
          m_start(start_),
          m_end(end_) {

    m_gridPoints.resize(m_size, Complex(0.0));
    SetLinspace();
}

void Atom::Grid::SetLinspace() {
    f64 start = m_start.real(), end = m_end.real();
    f64 stepSize = (end - start) / (f64) (m_size - 1);

    for (u32 i = 0; i < m_gridPoints.size(); i++) {
        f64 r = i * stepSize;
        m_gridPoints[i] = Complex(r,0.0);
    }
}

void Atom::Grid::Log() {
    u32 increment = 0;
    std::string logString = m_logPrefix + "gridPoint[%i]: %E";

    for (auto point : m_gridPoints) {
        Logger::Log(logString.c_str(), increment, point.real());
        increment += 1;
    }
}

const ZVector &Atom::Grid::getGridPoints() const {
    return m_gridPoints;
}