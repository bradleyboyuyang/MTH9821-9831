cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project (basic_transformation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

set(EIGEN_ROOT_DIR "/usr/local/include/eigen/3.4.0_1/include/eigen3/")

include_directories(
    ${EIGEN_ROOT_DIR}
)

add_executable (basic_transformation main.cpp regression.cpp)
