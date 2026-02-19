cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(dfdm LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)



if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type (Debug, Release, etc.)" FORCE)
endif()
option(BUILD_TESTS "Building Tests" OFF)

message(STATUS "CMAKE BUILD TYPE set to ${CMAKE_BUILD_TYPE}")
message(STATUS "BUILD_TESTS set to ${BUILD_TESTS}")
message(STATUS "Currently processing in directory: ${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{OPEN_BLAS_CMAKE_PATH})
# docker container openblas install: /home/my_mac/blas_test/blas_install/lib/cmake
find_package(OpenBLAS REQUIRED) 
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
add_subdirectory(spdlog)
add_subdirectory(matrix/matrix)
file(COPY tests DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "Currently binary directory: ${CMAKE_CURRENT_BINARY_DIR}")

if(OpenBLAS_FOUND)
    message(STATUS "Parent CMake found OpenBlas, will be included in DFDM")
else()
    message(FATAL_ERROR "openblas not found")
endif()

file(GLOB SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/toml11 ${CMAKE_CURRENT_SOURCE_DIR}/CLI11/include/CLI ${CMAKE_CURRENT_SOURCE_DIR}/matrix/matrix/include ${CMAKE_CURRENT_SOURCE_DIR}/spdlog ${OpenBLAS_INCLUDE_DIRS})
add_executable(dfdm ${SOURCE_FILES})

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # target_compile_options(dfdm PRIVATE -fsanitize=address -g)
    # target_link_options(dfdm PRIVATE -fsanitize=address)
endif()

# target_link_libraries(dfdm matrix_lib spdlog MPI::MPI_CXX)
target_link_libraries(dfdm matrix_lib spdlog MPI::MPI_CXX OpenMP::OpenMP_CXX)
if(BUILD_TESTS STREQUAL "ON")
    target_compile_definitions(dfdm PUBLIC ENABLE_TEST)

    set (passRegex ".test Passed")
    set (failRegex ".test Failed!")
    message(STATUS "Adding DFDM tests")
    include(CTest)

    add_test(bspline_matrix python3 tests/bspline.py)
    set_property (TEST bspline_matrix PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
    set_property (TEST bspline_matrix PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")
    set_tests_properties(bspline_matrix PROPERTIES LABELS "operators")

    add_test(full_simulation_mpi python3 tests/test_run_mpi.py)
    set_property (TEST full_simulation_mpi PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
    set_property (TEST full_simulation_mpi PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")
    set_tests_properties(full_simulation_mpi PROPERTIES LABELS "full_sim")
    
    add_test(full_simulation_serial python3 tests/test_run.py)
    set_property (TEST full_simulation_serial PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
    set_property (TEST full_simulation_serial PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")
    set_tests_properties(full_simulation_serial PROPERTIES LABELS "full_sim")

    add_test(src_recv_loc python3 tests/test_source_location.py)
    set_property (TEST src_recv_loc PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
    set_property (TEST src_recv_loc PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")
    set_tests_properties(src_recv_loc PROPERTIES LABELS "src_recv_loc")
endif()