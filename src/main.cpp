#include <vector>

#include "custom_asserts.h"
#include "logger.h"
#include "typedefs.h"

struct Grid {
    u32 size = 1000;
    f64 start = 0.0, end = 1.0;
    f64 stepSize = (end - start) / (f64) (size - 1);
    std::vector<f64> values;
};

Complex get_coefficient_deBoor(const void *knotSequence, const u32 bsplineOrder, const f64 coordinate) {
    u32 k = bsplineOrder;
    Complex *t = (Complex *) knotSequence;
    Complex d[bsplineOrder];


}

int main() {

    Grid radialGrid;

    // Linspace grid
    {
        for (int i = 0; i < radialGrid.size; i++) {
            f64 r = i * radialGrid.stepSize;
            radialGrid.values.push_back(r);
        }
    }
//    Logger::Trace("%f", radialGrid.values[0]);
//    Logger::Trace("%f", radialGrid.values[radialGrid.values.size()-1]);
//    Logger::Trace("%i", radialGrid.values.size());

    std::vector<Complex> knotPoints;
    u32 bsplineOrder = 4;
    u32 numberOfGhostPointsInEachEnd = (bsplineOrder - 1);
    u32 numberOfBsplines = 13;
    u32 numberOfKnotPoints = numberOfBsplines + bsplineOrder;
    u32 numberOfPhysicalKnotPoints = numberOfKnotPoints - 2 * numberOfGhostPointsInEachEnd;

    knotPoints.resize(numberOfKnotPoints, Complex{0.0, 0.0});

    // linear knot sequence
    {
        // set ghost knot points
        for (int i = 0; i < numberOfGhostPointsInEachEnd; ++i) {
            knotPoints[i] = Complex{radialGrid.start, 0.0};
            u32 endIndex = (knotPoints.size() - i) - 1;
            knotPoints[endIndex] = Complex{radialGrid.end, 0.0};
        }

        // set real knot points
        u32 offset = numberOfGhostPointsInEachEnd;
        u32 stride = (u32) (1000 / numberOfPhysicalKnotPoints);
        for (int i = 0; i < numberOfPhysicalKnotPoints; ++i) {
            u32 knotPointIndex = i + offset;
            u32 gridIndex = i * stride;
            ASSERT(gridIndex <= radialGrid.size - 1);
            Complex knotPointValue = Complex{radialGrid.values[gridIndex], 0.0};
            knotPoints[knotPointIndex] = knotPointValue;
            //Logger::Trace("%i %i %f", knotPointIndex, gridIndex, knotPointValue);
        }
    }

    Complex S_at_x = get_coefficient_deBoor(knotPoints.data(), bsplineOrder, 0.5)
    return 0;
}



