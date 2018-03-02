# Install script for directory: /export/home/urawat/Desktop/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/AdolcForward"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/AlignedVector3"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/ArpackSupport"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/AutoDiff"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/BVH"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/EulerAngles"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/FFT"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/IterativeSolvers"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/KroneckerProduct"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/MatrixFunctions"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/MoreVectorization"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/MPRealSupport"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/NonLinearOptimization"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/NumericalDiff"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/OpenGLSupport"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/Polynomials"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/Skyline"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/SparseExtra"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/SpecialFunctions"
    "/export/home/urawat/Desktop/eigen/unsupported/Eigen/Splines"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/export/home/urawat/Desktop/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/export/home/urawat/Desktop/eigen/build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

