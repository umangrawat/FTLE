cmake_minimum_required(VERSION 2.8)

PROJECT(Velocity_field)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(Velocity_field MACOSX_BUNDLE Velocity_field.cpp )

target_link_libraries(Velocity_field ${VTK_LIBRARIES})
