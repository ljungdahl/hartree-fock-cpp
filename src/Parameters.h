#pragma once
#include <string>

#include "typedefs.h"
#include "logger.h"
#include "custom_asserts.h"

namespace Atom {
    enum AngularMomentumQuantumNumber {
        s = 0,
        p = 1,
        d = 2,
        f = 3,
        g = 4,
    };

    struct SubShell {
        u32 n; // PrincipalQuantumNumber
        u32 l; // AngularMomentumQuantumNumber
        u32 numberOfElectrons;
        f64 radialFunctionIntegratedOverAllSpace = 1.0;
    };

    struct SystemParameters {
        bool isInitialised = false;
        std::string name;
        u32 Z;
        std::vector<Atom::SubShell> shells;
        u32 totalOccupationNumber = 0;
        void Log() {
            ASSERT(isInitialised);

            Logger::Log("Atomic system: %s", name.c_str());
            Logger::Log("Z = %i", Z);
            Logger::Log("Shells: ");
            for(auto shell : shells) {
                Logger::Log("(n, l) = (%i, %i), occupied by %i electrons", shell.n, shell.l, shell.numberOfElectrons);
            }
            printf("\n");
        }
    };

    struct GridParameters {
        bool isInitialised = false;
        u32 numberOfGridPoints;
        f64 gridStart;
        f64 gridEnd;
        void Log() {
            ASSERT(isInitialised);
            Logger::Log("GridParameters: ");
            Logger::Log("Number of grid points: %i", numberOfGridPoints);
            Logger::Log("Grid start coordinate: %f", gridStart);
            Logger::Log("Grid end coordinate: %f \n", gridEnd);
        }
    };

    struct BsplineParameters {
        bool isInitialised = false;
        u32 bsplineOrder;
        u32 numberOfKnotPoints;
        void Log() {
            ASSERT(isInitialised);
            Logger::Log("BsplineParameters: ");
            Logger::Log("Bspline order: %i", bsplineOrder);
            Logger::Log("Number of knot points: %i \n", numberOfKnotPoints);
        }
    };

    GridParameters GetGridParameters();
    BsplineParameters GetBsplineParameters();
    SystemParameters GetAtomicSystemParameters();


}
