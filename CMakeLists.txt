cmake_minimum_required(VERSION 3.16)
project(hartree_fock_cpp)

set(CMAKE_CXX_STANDARD 14)

if (WIN32)
    message("Win32 MKL 2020 Include!")
    set(MKL_WIN_PATH
            "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020/windows")

    set(MKL_DIR
            ${MKL_WIN_PATH}/mkl)

    set(MKL_LIB_DIR
            ${MKL_DIR}/lib/intel64_win)

    set(MKL_COMPILER_DIR
            ${MKL_WIN_PATH}/compiler/lib/intel64_win)

    set(MKL_COMPILER_REDIST_DIR
            ${MKL_WIN_PATH}/redist/intel64_win/compiler)

endif (WIN32)

# NOTE(anton): I am using arch linux and the intel-mkl/intel-mkl-static packages.
# The install might look different if you've downloaded from Intel's webpage.
if (UNIX)
    message("Unix MKL 2020 Include!")

        set(MKL_DIR
                /opt/intel/mkl)
#    set(MKL_DIR
#            /home/anton/lib/mkl/mkl)

    set(MKL_LIB_DIR ${MKL_DIR}/lib/intel64)

    set(MKL_COMPILER_DIR
            ${MKL_LIB_DIR})

    set(MKL_COMPILER_REDIST_DIR
            ${MKL_LIB_DIR})
endif (UNIX)


set(MKL_INCLUDE_DIRS ${MKL_DIR}/include)

include_directories(${MKL_INCLUDE_DIRS})
include_directories(src)

add_executable(${PROJECT_NAME}
        src/main.cpp
        src/logger.cpp
        src/Bsplines.cpp
        src/Grid.cpp
        src/FileIO.cpp
        src/Matrix.cpp
        src/GaussLegendre.cpp
        src/Eigenproblems.cpp)

if (UNIX)
    # Single Dynamic Link library used on my Arch Linu machine!

    # We get these options from the Intel link-line tool.
    set(GCC_COVERAGE_COMPILE_FLAGS "-m64")
    set(GCC_COVERAGE_LINK_FLAGS "-L${MKL_DIR}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl")

    # Append our custom flags
    set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif (UNIX)

if (WIN32)
    # Static linking MKL libs
    set(MKL_INTERFACE_LAYER ${MKL_LIB_DIR}/mkl_intel_lp64.lib)
    set(MKL_THREADING_LAYER ${MKL_LIB_DIR}/mkl_intel_thread.lib)
    set(MKL_COMPUTATIONAL_LAYER ${MKL_LIB_DIR}/mkl_core.lib)

    # This shit is to link the run time library and copy the dll to the same location as the .exe binary
    # see https://stackoverflow.com/a/59594095
    set(MKL_COMPILER_LIB ${MKL_COMPILER_DIR}/libiomp5md.lib)
    add_library(mkl_libiomp5md SHARED IMPORTED GLOBAL)
    set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_IMPLIB_RELEASE ${MKL_COMPILER_LIB})
    set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_LOCATION_RELEASE ${MKL_COMPILER_REDIST_DIR}/libiomp5md.dll)
    set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_IMPLIB_DEBUG ${MKL_COMPILER_LIB})
    set_property(TARGET mkl_libiomp5md PROPERTY IMPORTED_LOCATION_DEBUG ${MKL_COMPILER_REDIST_DIR}/libiomp5md.dll)

    #Combine all libs
    set(MKL_LIBRARIES ${MKL_INTERFACE_LAYER} ${MKL_THREADING_LAYER} ${MKL_COMPUTATIONAL_LAYER})

    target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES} mkl_libiomp5md)

    # This commands copies the .dll to same folder as the built executable.
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:mkl_libiomp5md> $<TARGET_FILE_DIR:${PROJECT_NAME}>)
endif (WIN32)

