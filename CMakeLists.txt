cmake_minimum_required(VERSION 3.16)
project(hartree_fock_cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(hartree_fock_cpp
        src/main.cpp
        src/logger.cpp
        src/logger.h src/custom_asserts.h)