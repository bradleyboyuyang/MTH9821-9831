cmake_minimum_required(VERSION 3.12)
project(hw1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add your source files
set(SOURCE_FILES main.cpp)

# Create an executable
add_executable(main ${SOURCE_FILES})

