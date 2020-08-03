#include "custom_asserts.h"
#include "logger.h"

int main() {

    bool testAssert = true;
    Logger::Trace("Testing Logger: %i", 1);

    ASSERT(testAssert);
    return 0;
}


