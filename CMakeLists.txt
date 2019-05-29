cmake_minimum_required(VERSION 3.13)
project(GAHW)

set(CMAKE_CXX_STANDARD 14)
add_definitions(-DFLT_EVAL_METHOD=0)
add_definitions(-D_USE_MATH_DEFINES)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

include_directories(Include)

add_executable(
        GAHW Source/main.cpp
        Source/TestGenerator.cpp
        Source/Descriptor.cpp
        Source/GA.cpp
        Include/GA.h
        Source/NaiveGA.cpp
        Include/NaiveGA.h
        Source/dubins.c
        Include/dubins.h
        Source/SchGA.cpp 
        Include/SchGA.h
        Source/SearchEngine.cpp
        )
