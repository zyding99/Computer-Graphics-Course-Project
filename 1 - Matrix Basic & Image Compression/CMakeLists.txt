cmake_minimum_required(VERSION 3.15)
project(Homework1)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)
add_definitions(-DOPENCV_DISABLE_EIGEN_TENSOR_SUPPORT)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(Homework1 main.cpp)
target_link_libraries(Homework1 ${OpenCV_LIBS})