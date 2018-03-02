#include "cudaIntegrator.h"
#include "linalgCUDA.h"

/*
MAKEDEVICE int getIndexInputGrid(int x , int y , int z, int* dim)
{
    /// return index, 3 times since values
    int index = (z*dim[1]*dim[0]+ (y)*dim[0]+x)*3;
    return index;
}
*/

////changed 3D to 2D
MAKEDEVICE int getIndexInputGrid(int x , int y, int* dim)
{
    /// return index, 2 times since values
    int index = ((y)*dim[0]+x)*2;
    return index;
}

/*
MAKEDEVICE void getDataPoint(vec3 dataVec, int xPos, int yPos, int zPos,
                             float* inputGrid, int* dimInputGrid)
{
    int index = getIndexInputGrid(xPos,yPos,zPos,dimInputGrid);
    dataVec[0] = inputGrid[index+0];
    dataVec[1] = inputGrid[index+1];
    dataVec[2] = inputGrid[index+2];
}
*/

////changed 3 to 2
MAKEDEVICE void getDataPoint(vec2 dataVec, int xPos, int yPos,
                             float* inputGrid, int* dimInputGrid)
{
    int index = getIndexInputGrid(xPos,yPos,dimInputGrid);
    dataVec[0] = inputGrid[index+0];
    dataVec[1] = inputGrid[index+1];
    ////dataVec[2] = inputGrid[index+2];
}

/*
MAKEDEVICE bool interpolate(vec3 location, vec3 dataVec, float *inputGrid, double *origin,
            int *dimInputGrid, double *spacingInputGrid, int intDirection) {
    vec3 loc;
    d_vec3copy(location,loc);
    //
    //loc[0] -= origin[0];
    //loc[1] -= origin[1];
    //loc[2] -= origin[2];
    //TODO Check local index calculation
    int xPos = (int) (loc[0] / spacingInputGrid[0]);
    int yPos = (int) (loc[1] / spacingInputGrid[1]);
    int zPos = (int) (loc[2] / spacingInputGrid[2]);

    double xLocal = loc[0] - double(xPos)*spacingInputGrid[0];
    double yLocal = loc[1] - double(yPos)*spacingInputGrid[1];
    double zLocal = loc[2] - double(zPos)*spacingInputGrid[2];
    xLocal /= spacingInputGrid[0];
    yLocal /= spacingInputGrid[1];
    zLocal /= spacingInputGrid[2];


    if(xPos > dimInputGrid[0] || yPos > dimInputGrid[1]
       || zPos > dimInputGrid[2]  || xPos < 0 || yPos < 0 || zPos <0) {
        return false;
    }


    vec3 a,b,c,d,e,f,g,h;
    getDataPoint(a,xPos,yPos,zPos,inputGrid, dimInputGrid);
    getDataPoint(b,xPos+1,yPos,zPos,inputGrid, dimInputGrid);
    getDataPoint(c,xPos,yPos+1,zPos,inputGrid, dimInputGrid);
    getDataPoint(d,xPos+1,yPos+1,zPos,inputGrid, dimInputGrid);
    getDataPoint(e,xPos,yPos,zPos+1,inputGrid, dimInputGrid);
    getDataPoint(f,xPos+1,yPos,zPos+1,inputGrid, dimInputGrid);
    getDataPoint(g,xPos,yPos+1,zPos+1,inputGrid, dimInputGrid);
    getDataPoint(h,xPos+1,yPos+1,zPos+1,inputGrid, dimInputGrid);


    d_vec3trilint(a,b,c,d,e,f,g,h, xLocal, yLocal, zLocal, dataVec);
    if(intDirection == -1)
        d_vec3scal(dataVec, -1., dataVec);
    return true;
}
*/

////changed 3 to 2
MAKEDEVICE bool interpolate(vec2 location, vec2 dataVec, float *inputGrid, double *origin,
            int *dimInputGrid, double *spacingInputGrid, int intDirection) {
    vec2 loc;
    d_vec2copy(location,loc);
    //
    //loc[0] -= origin[0];
    //loc[1] -= origin[1];
    //loc[2] -= origin[2];
    //TODO Check local index calculation
    int xPos = (int) (loc[0] / spacingInputGrid[0]);
    int yPos = (int) (loc[1] / spacingInputGrid[1]);
    ////int zPos = (int) (loc[2] / spacingInputGrid[2]);

    double xLocal = loc[0] - double(xPos)*spacingInputGrid[0];
    double yLocal = loc[1] - double(yPos)*spacingInputGrid[1];
    ////double zLocal = loc[2] - double(zPos)*spacingInputGrid[2];
    xLocal /= spacingInputGrid[0];
    yLocal /= spacingInputGrid[1];
    ////zLocal /= spacingInputGrid[2];


    if(xPos > dimInputGrid[0] || yPos > dimInputGrid[1]
       || xPos < 0 || yPos < 0) {
        return false;
    }

////removed e,f,g,h
    vec2 a,b,c,d;
    getDataPoint(a,xPos,yPos,inputGrid, dimInputGrid);
    getDataPoint(b,xPos+1,yPos,inputGrid, dimInputGrid);
    getDataPoint(c,xPos,yPos+1,inputGrid, dimInputGrid);
    getDataPoint(d,xPos+1,yPos+1,inputGrid, dimInputGrid);
    ////getDataPoint(e,xPos,yPos,zPos+1,inputGrid, dimInputGrid);
    ////getDataPoint(f,xPos+1,yPos,zPos+1,inputGrid, dimInputGrid);
    ////getDataPoint(g,xPos,yPos+1,zPos+1,inputGrid, dimInputGrid);
    ////getDataPoint(h,xPos+1,yPos+1,zPos+1,inputGrid, dimInputGrid);

    ////changed trilint to bilint
    d_vec2bilint(a,b,c,d, xLocal, yLocal, dataVec);
    if(intDirection == -1)
        d_vec2scal(dataVec, -1., dataVec);
    return true;
}

/*
MAKEDEVICE double integratePoint(vec3 &location, float *inputGrid, double *origin, int *dimInputGrid, double *spacingInputGrid,
               int integrationDirection, double stepSize, double &stagnation) {


    double intTime = 1;
    vec3 dataVec;
    d_vec3set(dataVec,0,0,0);
    stagnation = 0.;
    ///


    if(false) {

        if (!interpolate(location, dataVec, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection)) {
            stagnation = 0.;
            return intTime;
        }
        ///Simple Euler step

        intTime = stepSize / d_vec3mag(dataVec);

        stagnation = d_vec3mag(dataVec);


        d_vec3nrm(dataVec, dataVec);
        d_vec3scal(dataVec, stepSize, dataVec);
        //d_vec3add(location, dataVec, location);
        d_vec3copy(location,dataVec);
        stagnation = 0.;
        return intTime;
    }
    ///
    /// Runge Kutta integration scheme






    double currentIntTime = 0;
    double currentIntLength = 0;
    vec3 currentLocation;
    vec3 vecNextLocation;
    vec3 currentLocationSave;
    vec3 k1,k2,k3,k4;
    vec3 k1norm, k2norm, k3norm, k4norm;
    vec3 k1scal, k2scal, k3scal, k4scal;
    vec3 step;
    vec3 temp;
    d_vec3copy(location,currentLocation);
    double interStepWeights[3] = {.5,.5,1.0};
    double finalStepWeights[4] = {1.0/6.0, 1.0/3.0,1.0/3.0,1.0/6.0};


    d_vec3copy(currentLocation, currentLocationSave);
    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }

    d_vec3copy(vecNextLocation,k1);
    d_vec3nrm(k1,k1norm);
    d_vec3scal(k1, stepSize,k1scal);

    d_vec3copy(k1scal,temp);
    d_vec3scal(temp, interStepWeights[0],temp);
    d_vec3copy(currentLocationSave, currentLocation);
    d_vec3add(currentLocation,temp,currentLocation);

    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }
    d_vec3copy(vecNextLocation,k2);
    d_vec3nrm(k2,k2norm);
    d_vec3scal(k2,stepSize,k2scal);

    d_vec3copy(k2scal,temp);
    d_vec3scal(temp,interStepWeights[1],temp);
    d_vec3copy(currentLocationSave, currentLocation);
    d_vec3add(currentLocation,temp,currentLocation);

    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }

    d_vec3copy(vecNextLocation,k3);
    d_vec3nrm(k3,k3norm);
    d_vec3scal(k3,stepSize,k3scal);

    d_vec3copy(k3scal,temp);
    d_vec3scal(temp,interStepWeights[2],temp);
    d_vec3copy(currentLocationSave, currentLocation);
    d_vec3add(currentLocation,temp,currentLocation);

    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }

    d_vec3copy(vecNextLocation,k4);
    d_vec3nrm(k4,k4norm);
    d_vec3scal(k4,stepSize,k4scal);

    d_vec3scal(k1norm, finalStepWeights[0],temp);
    d_vec3copy(temp, step);
    d_vec3scal(k2norm, finalStepWeights[1], temp);
    d_vec3add(step, temp,step);
    d_vec3scal(k3norm, finalStepWeights[2], temp);
    d_vec3add(step, temp,step);
    d_vec3scal(k4norm, finalStepWeights[3],temp);
    d_vec3add(step,temp,step);
    d_vec3nrm(step,step);
    d_vec3scal(step, stepSize, step);
    currentIntLength += d_vec3mag(step);

    vec3 vel;
    d_vec3copy(k1,vel);
    d_vec3scal(k1,finalStepWeights[0],vel);
    d_vec3scal(k2,finalStepWeights[1],temp);
    d_vec3add(vel,temp,vel);
    d_vec3scal(k3,finalStepWeights[2],temp);
    d_vec3add(vel,temp,vel);
    d_vec3scal(k4,finalStepWeights[3],temp);
    d_vec3add(vel,temp,vel);
    currentIntTime += stepSize/d_vec3mag(vel);

    d_vec3copy(currentLocationSave, currentLocation);
    d_vec3add(currentLocation,step,currentLocation);

    stagnation = d_vec3mag(vel);

    ///End of Runge Kutta

    intTime = stepSize/d_vec3mag(vel);


    //stagnation = d_vec3mag(step);


    //d_vec3nrm(vel,vel);
    //d_vec3scal(vel, stepSize, vel);
    d_vec3add(location, step, location);
    return intTime;

}
*/


////changed 3 to 2
MAKEDEVICE double integratePoint(vec2 &location, float *inputGrid, double *origin, int *dimInputGrid, double *spacingInputGrid,
               int integrationDirection, double stepSize, double &stagnation) {


    double intTime = 1;
    vec2 dataVec;
    d_vec2set(dataVec,0,0);
    stagnation = 0.;
    ///


    if(false) {

        if (!interpolate(location, dataVec, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection)) {
            stagnation = 0.;
            return intTime;
        }
        ///Simple Euler step

        intTime = stepSize / d_vec2mag(dataVec);

        stagnation = d_vec2mag(dataVec);


        d_vec2nrm(dataVec, dataVec);
        d_vec2scal(dataVec, stepSize, dataVec);
        //d_vec3add(location, dataVec, location);
        d_vec2copy(location,dataVec);
        stagnation = 0.;
        return intTime;
    }
    ///
    /// Runge Kutta integration scheme


    double currentIntTime = 0;
    double currentIntLength = 0;
    vec2 currentLocation;
    vec2 vecNextLocation;
    vec2 currentLocationSave;
    vec2 k1,k2,k3,k4;
    vec2 k1norm, k2norm, k3norm, k4norm;
    vec2 k1scal, k2scal, k3scal, k4scal;
    vec2 step;
    vec2 temp;
    d_vec2copy(location,currentLocation);
    double interStepWeights[3] = {.5,.5,1.0};
    double finalStepWeights[4] = {1.0/6.0, 1.0/3.0,1.0/3.0,1.0/6.0};


    d_vec2copy(currentLocation, currentLocationSave);
    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }

    d_vec2copy(vecNextLocation,k1);
    d_vec2nrm(k1,k1norm);
    d_vec2scal(k1, stepSize,k1scal);

    d_vec2copy(k1scal,temp);
    d_vec2scal(temp, interStepWeights[0],temp);
    d_vec2copy(currentLocationSave, currentLocation);
    d_vec2add(currentLocation,temp,currentLocation);

    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }
    d_vec2copy(vecNextLocation,k2);
    d_vec2nrm(k2,k2norm);
    d_vec2scal(k2,stepSize,k2scal);

    d_vec2copy(k2scal,temp);
    d_vec2scal(temp,interStepWeights[1],temp);
    d_vec2copy(currentLocationSave, currentLocation);
    d_vec2add(currentLocation,temp,currentLocation);

    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }

    d_vec2copy(vecNextLocation,k3);
    d_vec2nrm(k3,k3norm);
    d_vec2scal(k3,stepSize,k3scal);

    d_vec2copy(k3scal,temp);
    d_vec2scal(temp,interStepWeights[2],temp);
    d_vec2copy(currentLocationSave, currentLocation);
    d_vec2add(currentLocation,temp,currentLocation);

    if(!interpolate(currentLocation, vecNextLocation, inputGrid, origin, dimInputGrid, spacingInputGrid, integrationDirection))
    {

        stagnation = 0.;
        return intTime;
    }

    d_vec2copy(vecNextLocation,k4);
    d_vec2nrm(k4,k4norm);
    d_vec2scal(k4,stepSize,k4scal);

    d_vec2scal(k1norm, finalStepWeights[0],temp);
    d_vec2copy(temp, step);
    d_vec2scal(k2norm, finalStepWeights[1], temp);
    d_vec2add(step, temp,step);
    d_vec2scal(k3norm, finalStepWeights[2], temp);
    d_vec2add(step, temp,step);
    d_vec2scal(k4norm, finalStepWeights[3],temp);
    d_vec2add(step,temp,step);
    d_vec2nrm(step,step);
    d_vec2scal(step, stepSize, step);
    currentIntLength += d_vec2mag(step);

    vec2 vel;
    d_vec2copy(k1,vel);
    d_vec2scal(k1,finalStepWeights[0],vel);
    d_vec2scal(k2,finalStepWeights[1],temp);
    d_vec2add(vel,temp,vel);
    d_vec2scal(k3,finalStepWeights[2],temp);
    d_vec2add(vel,temp,vel);
    d_vec2scal(k4,finalStepWeights[3],temp);
    d_vec2add(vel,temp,vel);
    currentIntTime += stepSize/d_vec2mag(vel);

    d_vec2copy(currentLocationSave, currentLocation);
    d_vec2add(currentLocation,step,currentLocation);

    stagnation = d_vec2mag(vel);

    ///End of Runge Kutta

    intTime = stepSize/d_vec2mag(vel);


    //stagnation = d_vec3mag(step);


    //d_vec3nrm(vel,vel);
    //d_vec3scal(vel, stepSize, vel);
    d_vec2add(location, step, location);
    return intTime;

}

/*
#if !defined UNITTEST && !defined CPUEXEC
MAKEDEVICE int getGlobalIdx() {
    int blockId = blockIdx.x
                  + blockIdx.y * gridDim.x
                  + gridDim.x * gridDim.y * blockIdx.z;
    int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
                   + (threadIdx.z * (blockDim.x * blockDim.y))
                   + (threadIdx.y * blockDim.x)
                   + threadIdx.x;
    return threadId;
}
#endif
*/

////removed z component
#if !defined UNITTEST && !defined CPUEXEC
MAKEDEVICE int getGlobalIdx() {
    int blockId = blockIdx.x
                  + blockIdx.y * gridDim.x;
    int threadId = (threadIdx.y * blockDim.x)
                   + threadIdx.x;
    return threadId;
}
#endif


/// fy refers to the current function value
/// fx and fz to the i-1 and i+1 function value
MAKEDEVICE double d_centralDiff(double fx, double fz, double dist)
{

    return (fz-fx)/(2*dist);
}
MAKEDEVICE double d_forwardDiff(double fy, double fz, double dist)
{
    return (fz-fy)/dist;
}
///Changed  fx and fy because of wrong function call further dowm
MAKEDEVICE double d_backwardDiff(double fx, double fy, double dist)
{
    return (fy-fx)/dist;

}

/*
#if !defined UNITTEST && !defined CPUEXEC
MAKEGLOBAL void computeFTLE(double* ftleField, double* flowMap, int* dimSeedGrid,
                            double* spacingSeedGrid, double intTime, double d_LIMIT_DOUBLE)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;


    int z = index / (dimSeedGrid[0] * dimSeedGrid[1]);
    int y = (index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[0];
    int x = (index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[0];
    //if(index < dimSeedGrid[0]*dimSeedGrid[1]*dimSeedGrid[2]) {
     //   int x = blockIdx.x * blockDim.x + threadIdx.x;//(index / (dimSeedGrid[0] * dimSeedGrid[1]));
     //   int y = blockIdx.y * blockDim.y + threadIdx.y;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[1];
     //   int z = blockIdx.z * blockDim.z + threadIdx.z;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[1];
    if (x < dimSeedGrid[0] && y < dimSeedGrid[1] && z < dimSeedGrid[2] &&
            x >=0 && y >=0 && z>=0){
        vec3 t1, t2;
        int id;
        vec3 du, dv, dw;


        if (z <= 0) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y, z + 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[2] = d_forwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[2] = d_forwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[2] = d_forwardDiff(t1[2], t2[2], spacingSeedGrid[2]);

        } else if (z >= dimSeedGrid[2] - 1) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y, z - 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[2] = d_backwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[2] = d_backwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[2] = d_backwardDiff(t1[2], t2[2], spacingSeedGrid[2]);
        } else {
            id = getIndexInputGrid(x, y, z - 1, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y, z + 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[2] = d_centralDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[2] = d_centralDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[2] = d_centralDiff(t1[2], t2[2], spacingSeedGrid[2]);

        }


        /// y Component
        if (y <= 0) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y + 1, z, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[1] = d_forwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[1] = d_forwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[1] = d_forwardDiff(t1[2], t2[2], spacingSeedGrid[2]);

        } else if (y >= dimSeedGrid[1] - 1) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y - 1, z, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[1] = d_backwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[1] = d_backwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[1] = d_backwardDiff(t1[2], t2[2], spacingSeedGrid[2]);

        } else {
            id = getIndexInputGrid(x, y - 1, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y + 1, z, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[1] = d_centralDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[1] = d_centralDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[1] = d_centralDiff(t1[2], t2[2], spacingSeedGrid[2]);
        }


        ///x Component
        if (x <= 0) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x + 1, y, z, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[0] = d_forwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[0] = d_forwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[0] = d_forwardDiff(t1[2], t2[2], spacingSeedGrid[2]);
        } else if (x >= dimSeedGrid[0] - 1) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x - 1, y, z, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[0] = d_backwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[0] = d_backwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[0] = d_backwardDiff(t1[2], t2[2], spacingSeedGrid[2]);


        } else {
            id = getIndexInputGrid(x - 1, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x + 1, y, z, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[0] = d_centralDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[0] = d_centralDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[0] = d_centralDiff(t1[2], t2[2], spacingSeedGrid[2]);

        }

        mat3 jacobi;
        d_mat3setrows(jacobi, du, dv, dw);
        mat3 jacobiT;
        mat3 cgTensor;
        d_mat3trp(jacobi, jacobiT);
        d_mat3mul(jacobiT, jacobi, cgTensor);
        double eMax = 0;
        vec3 eigenV;
        int realEigen = d_mat3eigenvalues(cgTensor, eigenV);
        eMax = fmax(fmax(eigenV[0], eigenV[1]), eigenV[2]);

        id = getIndexInputGrid(x, y, z, dimSeedGrid) / 3;
        if (realEigen != 3 || fabs(eMax) < d_LIMIT_DOUBLE) {
            // std::cerr<<"Eigenvalues of Cauchy Green Tensor are not real"<<endl;
            //cout<<realEigen<< " eigenvalues"<< eigenV[0]<< " " <<eigenV[1]<< " "<< eigenV[2]<<" "<<LIMIT_DOUBLE<<endl;
            eMax = 0;
            ftleField[id] = 0;
            return;
        }


        eMax = 1.0 / fabs(intTime) * log(sqrt(eMax));
        ftleField[id] = eMax;

    }
}
#endif
*/

////removing z components
#if !defined UNITTEST && !defined CPUEXEC
MAKEGLOBAL void computeFTLE(double* ftleField, double* flowMap, int* dimSeedGrid,
                            double* spacingSeedGrid, double intTime, double d_LIMIT_DOUBLE)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;


    ////int z = index / (dimSeedGrid[0] * dimSeedGrid[1]);
    int y = (index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[0];
    int x = (index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[0];
    //if(index < dimSeedGrid[0]*dimSeedGrid[1]*dimSeedGrid[2]) {
     //   int x = blockIdx.x * blockDim.x + threadIdx.x;//(index / (dimSeedGrid[0] * dimSeedGrid[1]));
     //   int y = blockIdx.y * blockDim.y + threadIdx.y;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[1];
     //   int z = blockIdx.z * blockDim.z + threadIdx.z;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[1];
    if (x < dimSeedGrid[0] && y < dimSeedGrid[1] &&
            x >=0 && y >=0){
        vec2 t1, t2;
        int id;
        vec2 du, dv, dw;

/*
        if (z <= 0) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y, z + 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[2] = d_forwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[2] = d_forwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[2] = d_forwardDiff(t1[2], t2[2], spacingSeedGrid[2]);

        } else if (z >= dimSeedGrid[2] - 1) {
            id = getIndexInputGrid(x, y, z, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y, z - 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[2] = d_backwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[2] = d_backwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[2] = d_backwardDiff(t1[2], t2[2], spacingSeedGrid[2]);
        } else {
            id = getIndexInputGrid(x, y, z - 1, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y, z + 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            t2[2] = (double)flowMap[id + 2];

            du[2] = d_centralDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[2] = d_centralDiff(t1[1], t2[1], spacingSeedGrid[1]);
            dw[2] = d_centralDiff(t1[2], t2[2], spacingSeedGrid[2]);

        }
*/
	////removed z component
        /// y Component
        if (y <= 0) {
            id = getIndexInputGrid(x, y, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            ////t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y + 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            ////t2[2] = (double)flowMap[id + 2];

            du[1] = d_forwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[1] = d_forwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            ////dw[1] = d_forwardDiff(t1[2], t2[2], spacingSeedGrid[2]);

        } else if (y >= dimSeedGrid[1] - 1) {
            id = getIndexInputGrid(x, y, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            ////t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y - 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            ////t2[2] = (double)flowMap[id + 2];

            du[1] = d_backwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[1] = d_backwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            ////dw[1] = d_backwardDiff(t1[2], t2[2], spacingSeedGrid[2]);

        } else {
            id = getIndexInputGrid(x, y - 1, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            ////t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x, y + 1, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            ////t2[2] = (double)flowMap[id + 2];

            du[1] = d_centralDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[1] = d_centralDiff(t1[1], t2[1], spacingSeedGrid[1]);
            ////dw[1] = d_centralDiff(t1[2], t2[2], spacingSeedGrid[2]);
        }

	////removed z component
        ///x Component
        if (x <= 0) {
            id = getIndexInputGrid(x, y, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            ////t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x + 1, y, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            ////t2[2] = (double)flowMap[id + 2];

            du[0] = d_forwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[0] = d_forwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            ////dw[0] = d_forwardDiff(t1[2], t2[2], spacingSeedGrid[2]);
        } else if (x >= dimSeedGrid[0] - 1) {
            id = getIndexInputGrid(x, y, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            ////t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x - 1, y, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            ////t2[2] = (double)flowMap[id + 2];

            du[0] = d_backwardDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[0] = d_backwardDiff(t1[1], t2[1], spacingSeedGrid[1]);
            ////dw[0] = d_backwardDiff(t1[2], t2[2], spacingSeedGrid[2]);


        } else {
            id = getIndexInputGrid(x - 1, y, dimSeedGrid);
            t1[0] = (double)flowMap[id + 0];
            t1[1] = (double)flowMap[id + 1];
            ////t1[2] = (double)flowMap[id + 2];
            id = getIndexInputGrid(x + 1, y, dimSeedGrid);
            t2[0] = (double)flowMap[id + 0];
            t2[1] = (double)flowMap[id + 1];
            ////t2[2] = (double)flowMap[id + 2];

            du[0] = d_centralDiff(t1[0], t2[0], spacingSeedGrid[0]);
            dv[0] = d_centralDiff(t1[1], t2[1], spacingSeedGrid[1]);
            ////dw[0] = d_centralDiff(t1[2], t2[2], spacingSeedGrid[2]);

        }

        mat2 jacobi;
        d_mat2setrows(jacobi, du, dv);
        mat2 jacobiT;
        mat2 cgTensor;
        d_mat2trp(jacobi, jacobiT);
        d_mat2mul(jacobiT, jacobi, cgTensor);
        double eMax = 0;
        vec2 eigenV;
        int realEigen = d_mat2eigenvalues(cgTensor, eigenV);
	////changed fmax(fmax()) to fmax()
        eMax = fmax(eigenV[0], eigenV[1]);

        id = getIndexInputGrid(x, y, dimSeedGrid) / 2;
        if (realEigen != 2 || fabs(eMax) < d_LIMIT_DOUBLE) {
            // std::cerr<<"Eigenvalues of Cauchy Green Tensor are not real"<<endl;
            //cout<<realEigen<< " eigenvalues"<< eigenV[0]<< " " <<eigenV[1]<< " "<< eigenV[2]<<" "<<LIMIT_DOUBLE<<endl;
            eMax = 0;
            ftleField[id] = 0;
            return;
        }


        eMax = 1.0 / fabs(intTime) * log(sqrt(eMax));
        ftleField[id] = eMax;

    }
}
#endif

/*
MAKEGLOBAL void
calcFtleField(float *inputArray, double*flowMap, double *ftleField, int *dimInputGrid, int *dimSeedGrid,
              double *spacingInputGrid,
              double *spacingSeedGrid, double *origin, double* originSource, int integrationDirection, double stepSize, double intTime,
              double stagnationThresh, double d_LIMIT_DOUBLE, int maxNumberOfSteps, int CPUindex = 0) {
    //int index = getGlobalIdx();
    ///Calc global thresh index
    int index;
#ifdef CPUEXEC
    index = CPUindex;
#else
    index = blockDim.x * blockIdx.x + threadIdx.x;
#endif

    int zIndex = index / (dimSeedGrid[0] * dimSeedGrid[1]);
    //blockIdx.x * blockDim.x + threadIdx.x;//(index / (dimSeedGrid[0] * dimSeedGrid[1]));
    int yIndex = (index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[0];
    //blockIdx.y * blockDim.y + threadIdx.y;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[1];
    int xIndex = (index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[0];
    //blockIdx.z * blockDim.z + threadIdx.z;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[1];


    if (xIndex < dimSeedGrid[0] && yIndex < dimSeedGrid[1] && zIndex < dimSeedGrid[2]&&
            xIndex >=0 && yIndex >=0 && zIndex|>=0)

    {
        vec3 location;
        double xLocation = double(xIndex) *  spacingSeedGrid[0] + originSource[0];//blockIdx.x * blockDim.x + threadIdx.x * spacingSeedGrid[0];
        double yLocation = double(yIndex) *  spacingSeedGrid[1] + originSource[1];//blockIdx.y * blockDim.y + threadIdx.y * spacingSeedGrid[1];
        double zLocation = double(zIndex) *  spacingSeedGrid[2] + originSource[2];//blockIdx.z * blockDim.z + threadIdx.z * spacingSeedGrid[2];


        d_vec3set(location,xLocation,yLocation,zLocation);

        double currentIntTime = 0.;
        double currentStagnation = 0;
        int maxSteps = maxNumberOfSteps;
        int currentSteps = 0;

        /// Integrate Streamline


        while (currentIntTime <= intTime) {
            currentIntTime += integratePoint(location, inputArray, origin, dimInputGrid, spacingInputGrid,
                                             integrationDirection, stepSize, currentStagnation);
            ///Check for stagnation of integration
            if (currentStagnation <= stagnationThresh)
                break;
            ++currentSteps;
            if(currentSteps >= maxSteps)
                break;
        }
        ///Set flow map (x,y,z)
        flowMap[3 * index] = location[0];
        flowMap[3 * index + 1] = location[1];
        flowMap[3 * index + 2] = location[2];
        //ftleField[index] = d_vec3mag(location);
        //TODO add finite differences when sync issue is fixed
        //computeFTLE(ftleField, flowMap , dimSeedGrid , spacingSeedGrid , intTime, d_LIMIT_DOUBLE);
        //d_vec3print(location, " location: ");
        //printf("SeedGrid: %d %d %d \n", dimSeedGrid[0],dimSeedGrid[1],dimSeedGrid[2]);
        //printf("Index: %d %d %d %d %f \n",xIndex, yIndex,zIndex, index, ftleField[index]);





    }

    return;
}
*/

MAKEGLOBAL void
calcFtleField(float *inputArray, double*flowMap, double *ftleField, int *dimInputGrid, int *dimSeedGrid,
              double *spacingInputGrid,
              double *spacingSeedGrid, double *origin, double* originSource, int integrationDirection, double stepSize, double intTime,
              double stagnationThresh, double d_LIMIT_DOUBLE, int maxNumberOfSteps, int CPUindex = 0) {
    //int index = getGlobalIdx();
    ///Calc global thresh index
    int index;
#ifdef CPUEXEC
    index = CPUindex;
#else
    index = blockDim.x * blockIdx.x + threadIdx.x;
#endif

////removed z component
    ////int zIndex = index / (dimSeedGrid[0] * dimSeedGrid[1]);
    //blockIdx.x * blockDim.x + threadIdx.x;//(index / (dimSeedGrid[0] * dimSeedGrid[1]));
    int yIndex = (index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[0];
    //blockIdx.y * blockDim.y + threadIdx.y;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) / dimSeedGrid[1];
    int xIndex = (index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[0];
    //blockIdx.z * blockDim.z + threadIdx.z;//(index % (dimSeedGrid[0] * dimSeedGrid[1])) % dimSeedGrid[1];

////removed z component
    if (xIndex < dimSeedGrid[0] && yIndex < dimSeedGrid[1] &&
            xIndex >=0 && yIndex >=0)


    {
        vec2 location;
        double xLocation = double(xIndex) *  spacingSeedGrid[0] + originSource[0];//blockIdx.x * blockDim.x + threadIdx.x * spacingSeedGrid[0];
        double yLocation = double(yIndex) *  spacingSeedGrid[1] + originSource[1];//blockIdx.y * blockDim.y + threadIdx.y * spacingSeedGrid[1];
        ////double zLocation = double(zIndex) *  spacingSeedGrid[2] + originSource[2];//blockIdx.z * blockDim.z + threadIdx.z * spacingSeedGrid[2];


        d_vec2set(location,xLocation,yLocation);

        double currentIntTime = 0.;
        double currentStagnation = 0;
        int maxSteps = maxNumberOfSteps;
        int currentSteps = 0;

        /// Integrate Streamline


        while (currentIntTime <= intTime) {
            currentIntTime += integratePoint(location, inputArray, origin, dimInputGrid, spacingInputGrid,
                                             integrationDirection, stepSize, currentStagnation);
            ///Check for stagnation of integration
            if (currentStagnation <= stagnationThresh)
                break;
            ++currentSteps;
            if(currentSteps >= maxSteps)
                break;
        }

////changed 3 to 2
        ///Set flow map (x,y)
        flowMap[2 * index] = location[0];
        flowMap[2 * index + 1] = location[1];
        ////flowMap[3 * index + 2] = location[2];
        //ftleField[index] = d_vec3mag(location);
        //TODO add finite differences when sync issue is fixed
        //computeFTLE(ftleField, flowMap , dimSeedGrid , spacingSeedGrid , intTime, d_LIMIT_DOUBLE);
        //d_vec3print(location, " location: ");
        //printf("SeedGrid: %d %d %d \n", dimSeedGrid[0],dimSeedGrid[1],dimSeedGrid[2]);
        //printf("Index: %d %d %d %d %f \n",xIndex, yIndex,zIndex, index, ftleField[index]);





    }

    return;
}



CudaIntegrator::~CudaIntegrator() {}

/*
void CudaIntegrator::
inInterpolator(vec3 a, vec3 b, vec3 c, vec3 d, vec3 e, vec3 f, vec3 g, vec3 h, double spacingX,
                                        double spacingY, double spacingZ, vec3 location, vec3 &output) {


}
*/
////removed component e,f,g,h, spacingZ and changed 3 to 2 and changed name to bilin
void CudaIntegrator::bilinInterpolator(vec2 a, vec2 b, vec2 c, vec2 d, double spacingX,
                                        double spacingY, vec2 location, vec2 &output) {


}

/*
double *CudaIntegrator::integrate(float *inputArray, double* flowMap) {
    //double* flowMap;
    double *d_flowMap;
    double *ftleField;
    double *d_ftleField;
    float *d_inputArray;
    int *d_dimSeedGrid;
    int *d_dimInputGrid;
    double *d_spacingSeedGrid;
    double *d_spacingInputGrid;
    double *d_origin;
    double *d_originSource;
    int sizeFlowMap = sizeof(double) * this->numPointsSeedGrid * 3;
    int sizeFtleField = sizeof(double) * this->numPointsSeedGrid;
    int sizeInputArray = sizeof(float) * this->numPointsInputGrid * 3;
    int sizeDim = sizeof(int) * 3;
    int sizeSpacing = sizeof(double) * 3;


    ///Alloc host and device memory
    //flowMap = (double*) malloc(sizeFlowMap);
    ftleField = (double *) malloc(sizeFtleField);
*/

////changed 3 to 2
double *CudaIntegrator::integrate(float *inputArray, double* flowMap) {
    //double* flowMap;
    double *d_flowMap;
    double *ftleField;
    double *d_ftleField;
    float *d_inputArray;
    int *d_dimSeedGrid;
    int *d_dimInputGrid;
    double *d_spacingSeedGrid;
    double *d_spacingInputGrid;
    double *d_origin;
    double *d_originSource;
    int sizeFlowMap = sizeof(double) * this->numPointsSeedGrid * 2;
    int sizeFtleField = sizeof(double) * this->numPointsSeedGrid;
    int sizeInputArray = sizeof(float) * this->numPointsInputGrid * 2;
    int sizeDim = sizeof(int) * 2;
    int sizeSpacing = sizeof(double) * 2;


    ///Alloc host and device memory
    //flowMap = (double*) malloc(sizeFlowMap);
    ftleField = (double *) malloc(sizeFtleField);

/*
#ifndef UNITTEST
    cudaMalloc((void **) &d_flowMap, sizeFlowMap);
    cudaMalloc((void **) &d_ftleField, sizeFtleField);
    cudaMalloc((void **) &d_inputArray, sizeInputArray);
    cudaMalloc((void **) &d_dimSeedGrid, sizeDim);
    cudaMalloc((void **) &d_dimInputGrid, sizeDim);
    cudaMalloc((void **) &d_spacingSeedGrid, sizeSpacing);
    cudaMalloc((void **) &d_spacingInputGrid, sizeSpacing);
    cudaMalloc((void **) &d_origin, sizeSpacing);
    cudaMalloc((void **) &d_originSource, sizeSpacing);

    ///Set sizes for block and threads
    dim3 gridDim = dim3(this->dimensionsSeedGrid[0] /8 +1 , this->dimensionsSeedGrid[1] / 8 +1 , this->dimensionsSeedGrid[2] /8 +1);
    dim3 blockDim = dim3(8, 8, 8);

    int M = this->dimensionsSeedGrid[0] * this->dimensionsSeedGrid[1] *
            this->dimensionsSeedGrid[2];
    int N = 256;

    //for(int i=0 ; i<numPointsSeedGrid; ++i)
      //  printf("Location.x %i: %i \n",i ,ftleField[i]);

#endif
*/

#ifndef UNITTEST
    cudaMalloc((void **) &d_flowMap, sizeFlowMap);
    cudaMalloc((void **) &d_ftleField, sizeFtleField);
    cudaMalloc((void **) &d_inputArray, sizeInputArray);
    cudaMalloc((void **) &d_dimSeedGrid, sizeDim);
    cudaMalloc((void **) &d_dimInputGrid, sizeDim);
    cudaMalloc((void **) &d_spacingSeedGrid, sizeSpacing);
    cudaMalloc((void **) &d_spacingInputGrid, sizeSpacing);
    cudaMalloc((void **) &d_origin, sizeSpacing);
    cudaMalloc((void **) &d_originSource, sizeSpacing);

////changed 3 to 2
    ///Set sizes for block and threads
    dim3 gridDim = dim3(this->dimensionsSeedGrid[0] /8 +1 , this->dimensionsSeedGrid[1] /8 +1);
    dim3 blockDim = dim3(8, 8);

    ////changed 256 to 32
    int M = this->dimensionsSeedGrid[0] * this->dimensionsSeedGrid[1];
    int N = 32;

    //for(int i=0 ; i<numPointsSeedGrid; ++i)
      //  printf("Location.x %i: %i \n",i ,ftleField[i]);

#endif



#ifndef CPUEXEC
    ///Copy data to device
    printf("Start GPU Kernel\n");
    cudaMemcpy(d_flowMap, flowMap, sizeFlowMap, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ftleField, ftleField, sizeFtleField, cudaMemcpyHostToDevice);
    cudaMemcpy(d_inputArray, inputArray, sizeInputArray, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dimSeedGrid, this->dimensionsSeedGrid, sizeDim, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dimInputGrid, this->dimensionsInputGrid, sizeDim, cudaMemcpyHostToDevice);
    cudaMemcpy(d_spacingSeedGrid, this->spacingSeedGrid, sizeSpacing, cudaMemcpyHostToDevice);
    cudaMemcpy(d_spacingInputGrid, this->spacingInputGrid, sizeSpacing, cudaMemcpyHostToDevice);
    cudaMemcpy(d_origin, this->origin, sizeSpacing, cudaMemcpyHostToDevice);
    cudaMemcpy(d_originSource, this->originSource, sizeSpacing, cudaMemcpyHostToDevice);




    calcFtleField <<< M / N +1, N >>>
                                     (d_inputArray, d_flowMap, d_ftleField, d_dimInputGrid, d_dimSeedGrid,
                                             d_spacingInputGrid, d_spacingSeedGrid,
                                             d_origin, d_originSource, this->integrationDirection,
                                             this->stepSize, this->integrationTime,
                                             this->stagnationThreshold, LIMIT_DOUBLE, this->maxNumberOfSteps);


    cudaDeviceSynchronize();

    computeFTLE<<< M / N +1, N >>>
            (d_ftleField, d_flowMap , d_dimSeedGrid , d_spacingSeedGrid , this->integrationTime, LIMIT_DOUBLE);

    ///Copy data from device to host
    cudaDeviceSynchronize();
    cudaMemcpy(ftleField, d_ftleField, sizeFtleField, cudaMemcpyDeviceToHost);
    cudaMemcpy(flowMap, d_flowMap, sizeFlowMap, cudaMemcpyDeviceToHost);



    //for(int i=0 ; i<numPointsSeedGrid; ++i)
      //  printf("Location.x %i: %f \n",i ,ftleField[i]);


    cudaError_t error = cudaGetLastError();
    if(error!=cudaSuccess)
    {
        fprintf(stderr,"ERROR: %s\n", cudaGetErrorString(error) );
        exit(-1);
    }
    //free(flowMap);
    //free(ftleField);
    cudaFree(d_flowMap);
    cudaFree(d_ftleField);
    cudaFree(d_inputArray);
    cudaFree(d_dimSeedGrid);
    cudaFree(d_dimInputGrid);
    cudaFree(d_spacingSeedGrid);
    cudaFree(d_spacingInputGrid);
    cudaFree(d_origin);
#else
    for(int i = 0; i < this->numPointsSeedGrid; ++i) {
        calcFtleField(inputArray, flowMap, ftleField, dimensionsInputGrid, dimensionsSeedGrid,
                      spacingInputGrid, spacingSeedGrid,
                      origin, originSource, this->integrationDirection,
                      this->stepSize, this->integrationTime,
                      this->stagnationThreshold, LIMIT_DOUBLE, this->maxNumberOfSteps, i);


       // computeFTLE(d_ftleField, d_flowMap, d_dimSeedGrid, d_spacingSeedGrid, this->integrationTime, LIMIT_DOUBLE);
    }
#endif

    return ftleField;

}

