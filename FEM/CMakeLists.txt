cmake_minimum_required(VERSION 3.26)
project(FEM)

set(CMAKE_CXX_STANDARD 17)

add_executable(FEM main.cpp
        Displaying.cpp
        Calculations.cpp)
