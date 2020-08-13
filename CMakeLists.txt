cmake_minimum_required(VERSION 3.16)
project(hartree_fock_cpp)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(hartree_fock_cpp
        src/logger.cpp
        src/Bsplines.cpp
        src/main.cpp
        src/Grid.cpp)
