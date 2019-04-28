cmake_minimum_required(VERSION 3.13)
project(GAHW)

set(CMAKE_CXX_STANDARD 14)

add_executable(GAHW main.cpp GA.cpp GA.h NaiveGA.cpp NaiveGA.h)