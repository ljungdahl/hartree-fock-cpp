cmake_minimum_required(VERSION 3.16)
project(hartree_fock_cpp)

set(CMAKE_CXX_STANDARD 14)

if (WIN32)
    message("Win32 MKL 2020 Include!")
    set(MKL_DIR
            "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020/windows/mkl")
    set(MKL_LIB_DIR ${MKL_DIR}/lib/intel64_win)

    set(MKL_COMPILER_DIR
            "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020/windows/compiler")

    set(MKL_COMPILER_REDIST_DIR
            "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020/windows/redist/intel64_win/compiler")

endif (WIN32)

set(MKL_INCLUDE_DIRS ${MKL_DIR}/include)

include_directories(${MKL_INCLUDE_DIRS})
include_directories(src)

add_executable(${PROJECT_NAME}
        src/logger.cpp
        src/Bsplines.cpp
        src/main.cpp
        src/Grid.cpp
        src/FileIO.cpp
        src/testing_mkl.cpp)


# Static linking MKL libs
set(MKL_INTERFACE_LAYER ${MKL_LIB_DIR}/mkl_intel_lp64.lib)
set(MKL_THREADING_LAYER ${MKL_LIB_DIR}/mkl_intel_thread.lib)
set(MKL_COMPUTATIONAL_LAYER ${MKL_LIB_DIR}/mkl_core.lib)

# This shit is to link the run time library and copy the dll to the same location as the .exe binary
# see https://stackoverflow.com/a/59594095
set(MKL_COMPILER_LIB ${MKL_COMPILER_DIR}/lib/intel64_win/libiomp5md.lib)
add_library(mkl_libiomp5md SHARED IMPORTED GLOBAL)
set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_IMPLIB_RELEASE ${MKL_COMPILER_LIB})
set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_LOCATION_RELEASE ${MKL_COMPILER_REDIST_DIR}/libiomp5md.dll)
set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_IMPLIB_DEBUG ${MKL_COMPILER_LIB})
set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_LOCATION_DEBUG ${MKL_COMPILER_REDIST_DIR}/libiomp5md.dll)

#Combine all libs
set(MKL_LIBRARIES ${MKL_INTERFACE_LAYER} ${MKL_THREADING_LAYER} ${MKL_COMPUTATIONAL_LAYER})

target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES} mkl_libiomp5md)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:mkl_libiomp5md> $<TARGET_FILE_DIR:${PROJECT_NAME}>
        )