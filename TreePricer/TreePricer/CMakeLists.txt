cmake_minimum_required(VERSION 3.13)  # CMake version check
project(TreePricer)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 20)            # Enable c++14 standard

# use C++11
# setting warning compiler flags
add_compile_options(-Wall -std=c++20)

# Add main.cpp file of the project root directory as a source file
set(SOURCE_FILES AmericanPut.cpp           BTree.cpp                 BTreeDiscreteDividend.cpp EuropeanOption.cpp        EuropeanPut.cpp           TTree.cpp                 TreeIV.cpp                main.cpp)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(simple_example ${SOURCE_FILES})
