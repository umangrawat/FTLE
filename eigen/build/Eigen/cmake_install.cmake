# Install script for directory: /export/home/urawat/Desktop/eigen/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/export/home/urawat/Desktop/eigen/Eigen/SPQRSupport"
    "/export/home/urawat/Desktop/eigen/Eigen/SparseCholesky"
    "/export/home/urawat/Desktop/eigen/Eigen/PardisoSupport"
    "/export/home/urawat/Desktop/eigen/Eigen/Geometry"
    "/export/home/urawat/Desktop/eigen/Eigen/StdDeque"
    "/export/home/urawat/Desktop/eigen/Eigen/Householder"
    "/export/home/urawat/Desktop/eigen/Eigen/LU"
    "/export/home/urawat/Desktop/eigen/Eigen/SparseQR"
    "/export/home/urawat/Desktop/eigen/Eigen/StdVector"
    "/export/home/urawat/Desktop/eigen/Eigen/Sparse"
    "/export/home/urawat/Desktop/eigen/Eigen/Cholesky"
    "/export/home/urawat/Desktop/eigen/Eigen/IterativeLinearSolvers"
    "/export/home/urawat/Desktop/eigen/Eigen/Eigenvalues"
    "/export/home/urawat/Desktop/eigen/Eigen/SparseCore"
    "/export/home/urawat/Desktop/eigen/Eigen/Dense"
    "/export/home/urawat/Desktop/eigen/Eigen/OrderingMethods"
    "/export/home/urawat/Desktop/eigen/Eigen/QtAlignedMalloc"
    "/export/home/urawat/Desktop/eigen/Eigen/StdList"
    "/export/home/urawat/Desktop/eigen/Eigen/Core"
    "/export/home/urawat/Desktop/eigen/Eigen/QR"
    "/export/home/urawat/Desktop/eigen/Eigen/SuperLUSupport"
    "/export/home/urawat/Desktop/eigen/Eigen/Jacobi"
    "/export/home/urawat/Desktop/eigen/Eigen/SparseLU"
    "/export/home/urawat/Desktop/eigen/Eigen/PaStiXSupport"
    "/export/home/urawat/Desktop/eigen/Eigen/Eigen"
    "/export/home/urawat/Desktop/eigen/Eigen/MetisSupport"
    "/export/home/urawat/Desktop/eigen/Eigen/SVD"
    "/export/home/urawat/Desktop/eigen/Eigen/CholmodSupport"
    "/export/home/urawat/Desktop/eigen/Eigen/UmfPackSupport"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/export/home/urawat/Desktop/eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

