cmake_minimum_required(VERSION 3.15)
project(Homework2)

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)

aux_source_directory(. DIR_SRCS)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(Homework2 ${DIR_SRCS})
target_link_libraries(Homework2 ${OpenCV_LIBS} Eigen3::Eigen)