cmake_minimum_required(VERSION 3.3)
project(RayTracer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(RayTracer ${SOURCE_FILES})