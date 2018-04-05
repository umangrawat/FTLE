cmake_minimum_required(VERSION 3.5)
project(ColorMap)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE RELAESE)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../libs/linalg)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../catch)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Eigen)

if (NOT ParaView_BINARY_DIR)
    find_package(ParaView REQUIRED)
    include(${PARAVIEW_USE_FILE})
endif()

include(ParaViewPlugins)

ADD_PARAVIEW_PLUGIN(ColorMap "1.0"
        SERVER_MANAGER_XML ColorMap.xml
        SERVER_MANAGER_SOURCES  vtkColorMap.cpp
        SOURCES vtkColorMap.h vtkColorMap.cpp
        )
