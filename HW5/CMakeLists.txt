cmake_minimum_required(VERSION 3.24)
project(problem2_correct)

set(CMAKE_CXX_STANDARD 14)

include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/")

add_executable(
        regression.h
        regression.cpp
        mc_regression.h
        mc_regression.cpp
        american_pricers.cpp main.cpp)
