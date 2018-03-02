#ifndef __CUDAINTEGRATOR_H_
#define __CUDAINTEGRATOR_H_

#include "linalg.h"


/// REMEBER TO UNDECLARE
//TODO remember to define/undefine for tests
//#undef UNITTEST
//#undef CPUEXEC
//#define CPUEXEC

#ifndef UNITTEST
#include <cuda_runtime.h>
#include <cublas.h>
#include "device_functions.hpp"
#endif

#include <vector>


#ifdef UNITTEST
    #define CPUEXEC
    #define MAKEDEVICE
    #define MAKEGLOBAL
    #warning "Swap __device__ defintion"
    #include <iostream>
#elif !defined CPUEXEC
    #define MAKEDEVICE __device__
    #define MAKEGLOBAL __global__

#endif


#if defined CPUEXEC && !defined UNITTEST
    #define MAKEDEVICE
    #define MAKEGLOBAL
    #warning "Swap __device__ defintion"
    #include <iostream>
//#elif defined UNITTEST
//    #define MAKEDEVICE __device__
//    #define MAKEGLOBAL __global__
#endif

/*
struct GpuInfo{
    int dimensionSeedGrid[3];
    int dimensionInputGrid[3];
    double originInput[3];
    double originSource[3];
    int numberOfPointsSeedGrid;
    int numberOfPointsInputGrid;
    double spacingSeedGrid[3];
    double spacingInputGrid[3];
    int integrationDirection;
    double stepSize;
    double integrationTime;
    int maxNumberOfSteps;
    double stagnationThreshold;
};

*/
////changed 3 to 2
struct GpuInfo{
    int dimensionSeedGrid[2];
    int dimensionInputGrid[2];
    double originInput[2];
    double originSource[2];
    int numberOfPointsSeedGrid;
    int numberOfPointsInputGrid;
    double spacingSeedGrid[2];
    double spacingInputGrid[2];
    int integrationDirection;
    double stepSize;
    double integrationTime;
    int maxNumberOfSteps;
    double stagnationThreshold;
};
/*
class CudaIntegrator
{

private:


    double* integratePointRK4(vec3 startLocation);

public:
    CudaIntegrator(int dimensionsSeedGrid_[3], int dimensionsInputGrid_[3], double origin_[3], double originSource_[3],
                   int numberOfPointsSeedGrid_, int numberOfPointsInputGrid_,
                   double* spacingSeedGrid_, double* spacingInputGrid_,
                   int integrationDirection_, double stepSize_,
                   double intTime, double stagnationThreshold_, int maxNumberOfSteps_)
            : numPointsSeedGrid(numberOfPointsSeedGrid_),
              numPointsInputGrid(numberOfPointsInputGrid_),
              spacingSeedGrid(spacingSeedGrid_),
              spacingInputGrid(spacingInputGrid_),
              integrationDirection(integrationDirection_),
              stepSize(stepSize_), integrationTime(intTime),
              stagnationThreshold(stagnationThreshold_),
              maxNumberOfSteps(maxNumberOfSteps_)
            {
                this->dimensionsSeedGrid[0] = dimensionsSeedGrid_[0];
                this->dimensionsSeedGrid[1] = dimensionsSeedGrid_[1];
                this->dimensionsSeedGrid[2] = dimensionsSeedGrid_[2];
                this->dimensionsInputGrid[0] = dimensionsInputGrid_[0];
                this->dimensionsInputGrid[1] = dimensionsInputGrid_[1];
                this->dimensionsInputGrid[2] = dimensionsInputGrid_[2];
                this->origin[0] = origin_[0];
                this->origin[1] = origin_[1];
                this->origin[2] = origin_[2];
                this->originSource[0] = originSource_[0];
                this->originSource[1] = originSource_[1];
                this->originSource[2] = originSource_[2];
                this->info.numberOfPointsSeedGrid = numberOfPointsSeedGrid_;
                this->info.numberOfPointsInputGrid = numberOfPointsInputGrid_;
                this->info.maxNumberOfSteps = maxNumberOfSteps_;

            }
    ~CudaIntegrator();


    double origin[3];
    double originSource[3];
    int integrationDirection;
    double stepSize;
    double integrationTime;
    double stagnationThreshold;
    int dimensionsSeedGrid[3];
    int dimensionsInputGrid[3];
    int numPointsSeedGrid;
    int numPointsInputGrid;
    double* spacingSeedGrid;
    double* spacingInputGrid;
    int maxNumberOfSteps;

    GpuInfo info;




    void trilinInterpolator(vec3 a,vec3 b, vec3 c, vec3 d, vec3 e, vec3 f, vec3 g, vec3 h, double spacingX, double spacingY, double spacingZ, vec3 location, vec3& output);
    bool cellLocator(vec3& output,vec3 location);



    double* integrate(float* inputArray, double* flowMap);

    bool pointInterpolator(vec3 currentLocation, vec3& vecNextLocation);





};
*/

////changed 3D to 2D
class CudaIntegrator
{

private:


    double* integratePointRK4(vec2 startLocation);

public:
    CudaIntegrator(int dimensionsSeedGrid_[2], int dimensionsInputGrid_[2], double origin_[2], double originSource_[2],
                   int numberOfPointsSeedGrid_, int numberOfPointsInputGrid_,
                   double* spacingSeedGrid_, double* spacingInputGrid_,
                   int integrationDirection_, double stepSize_,
                   double intTime, double stagnationThreshold_, int maxNumberOfSteps_)
            : numPointsSeedGrid(numberOfPointsSeedGrid_),
              numPointsInputGrid(numberOfPointsInputGrid_),
              spacingSeedGrid(spacingSeedGrid_),
              spacingInputGrid(spacingInputGrid_),
              integrationDirection(integrationDirection_),
              stepSize(stepSize_), integrationTime(intTime),
              stagnationThreshold(stagnationThreshold_),
              maxNumberOfSteps(maxNumberOfSteps_)
            {
                this->dimensionsSeedGrid[0] = dimensionsSeedGrid_[0];
                this->dimensionsSeedGrid[1] = dimensionsSeedGrid_[1];
                ////this->dimensionsSeedGrid[2] = dimensionsSeedGrid_[2];
                this->dimensionsInputGrid[0] = dimensionsInputGrid_[0];
                this->dimensionsInputGrid[1] = dimensionsInputGrid_[1];
                ////this->dimensionsInputGrid[2] = dimensionsInputGrid_[2];
                this->origin[0] = origin_[0];
                this->origin[1] = origin_[1];
                ////this->origin[2] = origin_[2];
                this->originSource[0] = originSource_[0];
                this->originSource[1] = originSource_[1];
                ////this->originSource[2] = originSource_[2];
                this->info.numberOfPointsSeedGrid = numberOfPointsSeedGrid_;
                this->info.numberOfPointsInputGrid = numberOfPointsInputGrid_;
                this->info.maxNumberOfSteps = maxNumberOfSteps_;

            }
    ~CudaIntegrator();


    double origin[2];
    double originSource[2];
    int integrationDirection;
    double stepSize;
    double integrationTime;
    double stagnationThreshold;
    int dimensionsSeedGrid[2];
    int dimensionsInputGrid[2];
    int numPointsSeedGrid;
    int numPointsInputGrid;
    double* spacingSeedGrid;
    double* spacingInputGrid;
    int maxNumberOfSteps;

    GpuInfo info;



////changed 3 to 2, removed e,f,g,h and spacingZ and changed the name of trilin to bilin
    void bilinInterpolator(vec2 a,vec2 b, vec2 c, vec2 d, double spacingX, double spacingY, vec2 location, vec2& output);
    bool cellLocator(vec2& output,vec2 location);



    double* integrate(float* inputArray, double* flowMap);

    bool pointInterpolator(vec2 currentLocation, vec2& vecNextLocation);





};

/*
MAKEDEVICE int getIndexInputGrid(int x , int y , int z, int* dim);
MAKEDEVICE bool interpolate(vec3 location, vec3 dataVec, float *inputGrid, double *origin,
                           int *dimInputGrid, double *spacingInputGrid, int intDirection);
MAKEDEVICE void getDataPoint(vec3 dataVec, int xPos, int yPos, int zPos,
                             float* inputGrid, int* dimInputGrid);
MAKEDEVICE double integratePoint(vec3 &location, float *inputGrid, double *origin, int *dimInputGrid, double *spacingInputGrid,
                                 int integrationDirection, double stepize, double &stagnation);
*/

////removed z component and changed 3 to 2
MAKEDEVICE int getIndexInputGrid(int x , int y , int* dim);
MAKEDEVICE bool interpolate(vec2 location, vec2 dataVec, float *inputGrid, double *origin,
                           int *dimInputGrid, double *spacingInputGrid, int intDirection);
MAKEDEVICE void getDataPoint(vec2 dataVec, int xPos, int yPos,
                             float* inputGrid, int* dimInputGrid);
MAKEDEVICE double integratePoint(vec2 &location, float *inputGrid, double *origin, int *dimInputGrid, double *spacingInputGrid,
                                 int integrationDirection, double stepize, double &stagnation);



#define CUDA_OUTPUT_STEPS
#endif
