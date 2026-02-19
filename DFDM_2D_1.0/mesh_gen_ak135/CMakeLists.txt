cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mesh_gen LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
# option(BUILD_TESTS "Building Tests" ON)

message(STATUS "Currently processing in directory: ${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})
# docker container openblas install: /home/mymac/seismo/install/lib/cmake
find_package(OpenBLAS REQUIRED) 
# find_package(MPI REQUIRED)
add_subdirectory(matrix)
# file(COPY tests DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "Currently binary directory: ${CMAKE_CURRENT_BINARY_DIR}")

if(OpenBLAS_FOUND)
    message(STATUS "Parent CMake found OpenBlas, will be included in DFDM")
else()
    message(FATAL_ERROR "openblas not found")
endif()

file(GLOB SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)




include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/toml11 ${CMAKE_CURRENT_SOURCE_DIR}/matrix/matrix/include ${OpenBLAS_INCLUDE_DIRS})
#  ${CMAKE_CURRENT_SOURCE_DIR}/toml11 
add_executable(mesh_gen ${SOURCE_FILES})
target_link_libraries(mesh_gen matrix_lib)

# if(BUILD_TESTS)
#     set (passRegex ".test Passed")
#     set (failRegex ".test Failed!")
#     message(STATUS "Adding DFDM tests")
#     include(CTest)
#     add_test(full_simulation python3 tests/test_run.py)
#     set_property (TEST full_simulation PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#     set_property (TEST full_simulation PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")


#     add_test(full_simulation_mpi python3 tests/test_run_mpi.py)
#     set_property (TEST full_simulation_mpi PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#     set_property (TEST full_simulation_mpi PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")
# endif()