cmake_minimum_required(VERSION 3.1.2)
project(HW2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find the Boost library
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

# compile writer and reader separately
add_executable(writer q1-writer.cpp)
target_link_libraries(writer ${Boost_LIBRARIES})

add_executable(reader q1-reader.cpp)
target_link_libraries(reader ${Boost_LIBRARIES})

# compile test cases for products
add_executable(test_product test_products.cpp)