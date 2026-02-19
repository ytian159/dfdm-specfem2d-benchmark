cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(mat_test LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})

find_package(OpenBLAS REQUIRED) 
if(OpenBLAS_FOUND)
    message(STATUS "OpenBLAS found at: ${OpenBLAS_INCLUDE_DIRS}")
endif()
add_subdirectory(matrix)
include_directories(matrix/include ${OpenBLAS_INCLUDE_DIRS})

option(DFDM_MATRIX_BUILD_TESTS "Build DFDM matrix tests" OFF)

if(DFDM_MATRIX_BUILD_TESTS)
    add_executable(mat_test tests/test.cpp)
    target_link_libraries(mat_test matrix_lib)
endif()