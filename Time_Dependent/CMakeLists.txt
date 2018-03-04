cmake_minimum_required(VERSION 2.8)

PROJECT(Time_Dependent)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(Time_Dependent MACOSX_BUNDLE Time_Dependent.cpp )

target_link_libraries(Time_Dependent ${VTK_LIBRARIES})
