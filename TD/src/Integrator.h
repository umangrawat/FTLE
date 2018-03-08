#ifndef __Integrator_h
#define __Integrator_h

#include "linalg.h"

#include "vtkPolyData.h"

#include "vtkImageData.h"

#include <iostream>
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "cudaIntegrator.h"
#include <vtkDataArray.h>
#include "vtkPoints.h"
#include <vector>

using namespace std;


const double LIMIT_DOUBLE = 1000* std::numeric_limits<double>::epsilon();

/*
class Integrator
{


private:


    vtkSmartPointer<vtkPoints> integratePointRK4(vec3 startLocation);

public:
    Integrator();
    Integrator(int intDirection, double stepsize, double intTime, double stagThresh,
               vtkSmartPointer<vtkImageData> inputGrid, vtkSmartPointer<vtkImageData> source);
    ~Integrator();



    int integrationDirection;
    double stepSize;
    double integrationTime;
    int maxNumberOfSteps;
    double integrationLength;
    double stagnationThreshold;
    vtkSmartPointer<vtkImageData> input;
    vtkSmartPointer<vtkImageData> source;

    double origin[3];
    double originSource[3];
    bool calcStreamlines;



    void trilinInterpolator(vec3 a,vec3 b, vec3 c, vec3 d, vec3 e, vec3 f, vec3 g, vec3 h, double spacingX, double spacingY, double spacingZ, vec3 location, vec3& output);
    bool cellLocator(vec3& output,vec3 location);



    vtkSmartPointer<vtkPolyData> integrateRK4();
    void integrateRK4GPUWrapper();
    void setOrigin(double ori[3]);
    void setOriginSource(double ori[3]);

    bool pointInterpolator(vec3 currentLocation, vec3& vecNextLocation);


};
*/

////changed 3D to 2D
class Integrator
{


private:


    vtkSmartPointer<vtkPoints> integratePointRK4(vec2 startLocation);

public:
    Integrator();
    Integrator(int intDirection, double stepsize, double starttime, double endtime, double stagThresh,
               vtkSmartPointer<vtkImageData> inputGrid, vtkSmartPointer<vtkImageData> source,
               vtkSmartPointer<vtkImageData> nextinputGrid, int blockNum,  int initialnum, int finalnum, double* time, std::vector<vector<double>> previousPts);
    ~Integrator();

    int integrationDirection;
    double stepSize;
    double integrationTime;
    int maxNumberOfSteps;
    double integrationLength;
    double stagnationThreshold;
    vtkSmartPointer<vtkImageData> input;
    vtkSmartPointer<vtkImageData> source;
    vtkSmartPointer<vtkImageData> nextinput;

    std::vector<vector<double>> NewEndPoints;
    std::vector<vector<double>> PathlineEndPoints;

    double startTime;               ////
    double endTime;                 ////
    int initialNum;
    int finalNum;
    double origin[2];
    double originSource[2];
    bool calcStreamlines;
    int iter;
    double* timestep;


////changed 3 to 2 and removed e,f,g,h and spacingZ and changed the name from tri to bi
    void bilinInterpolator(vec2 a,vec2 b, vec2 c, vec2 d, double spacingX, double spacingY, vec2 location, vec2& output);
    bool cellLocator(vec2& output, vec2& nextoutput, vec2 location);
    vtkSmartPointer<vtkPolyData> integrateRK4();
    void integrateRK4GPUWrapper();
    void setOrigin(double ori[2]);
    void setOriginSource(double ori[2]);

    bool pointInterpolator(vec2 currentLocation, vec2& vecNextLocation, vec2& vecNextSliceLocation);


};

#endif
