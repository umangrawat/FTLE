cmake_minimum_required(VERSION 2.8)

PROJECT(lcs)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(lcs MACOSX_BUNDLE lcs.cpp )

target_link_libraries(lcs ${VTK_LIBRARIES})
