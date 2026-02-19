cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(matrixlib LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})

message(STATUS "Building DFDM Matrix libs...")
message(STATUS "Building DFDM Matrix Lib from directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "variable CMAKE_MODULE_PATH is set to: ${CMAKE_MODULE_PATH}")
message(STATUS "variable CMAKE_PREFIX_PATH is set to: ${CMAKE_PREFIX_PATH}")

find_package(OpenBLAS REQUIRED) 

file(GLOB SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if(OpenBLAS_FOUND)
    message(STATUS "OpenBLAS found at: ${OpenBLAS_INCLUDE_DIRS}")
    include_directories(${OpenBLAS_INCLUDE_DIRS}) 
    add_library(matrix_lib SHARED ${SOURCE_FILES})
    target_link_libraries(matrix_lib ${OpenBLAS_LIBRARIES})
else()
    message(FATAL_ERROR "OpenBlas not found, failed to build Matrix libs")
endif()


