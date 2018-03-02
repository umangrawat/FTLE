#define CATCH_CONFIG_MAIN


#include <iostream>
#include "catch.h"
//#include "vtkFTLE.h"
#include "Integrator.h"
#include "cudaIntegrator.h"





unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}
/*
TEST_CASE("GetInputIndex", "[getInputIndex]"){
    int dimension[3] = {10,10,10};
    CHECK(getIndexInputGrid(0,0,0,&dimension[0]) == 0);
    CHECK(getIndexInputGrid(9,9,9,&dimension[0]) == 999 *3);

    int dimension1[3] = {10,18,6};
    CHECK(getIndexInputGrid(0,0,0,&dimension1[0]) == 0);
    CHECK(getIndexInputGrid(3,8,2,&dimension1[0]) == (2*18*10+8*10+3)*3);
    CHECK(getIndexInputGrid(9,17,5,&dimension1[0]) == (10*18*6-1) *3);

}
*/

////changed 3D to 2D and changed 999 to 99
TEST_CASE("GetInputIndex", "[getInputIndex]"){
    int dimension[2] = {10,10};
    CHECK(getIndexInputGrid(0,0,&dimension[0]) == 0);
    CHECK(getIndexInputGrid(9,9,&dimension[0]) == 99 *2);

    int dimension1[2] = {10,18};
    CHECK(getIndexInputGrid(0,0,&dimension1[0]) == 0);
    CHECK(getIndexInputGrid(3,8,&dimension1[0]) == (8*10+3)*2);
    CHECK(getIndexInputGrid(9,17,&dimension1[0]) == (10*18-1) *2);

}

/*
TEST_CASE("Test Integration", "[Test integration]"){

    vec3 loc1,loc2;
    vec3set(loc1, 0,0,0);
    vec3set(loc2, 3.9,1.6,10.9);
    vec3 dataVec;

    double originInput[3]= {0,0,0};
    double originSeed[3] = {3,1,6};
    double spacingInput[3] = {1,1,1};
    double spacingSeed[3] = {9,9,9};
    int dimInputGrid[3] = {401,501,21};
    int dimSeedGrid[3] = {2,2,2};
    double stepSize = .1;
    double stagnation = LIMIT_DOUBLE;
    float* inputGrid = (float*) malloc(sizeof(float) *
                                              dimInputGrid[0]*dimInputGrid[1]*dimInputGrid[2]*3);

    for(int i = 0; i <dimInputGrid[0]*dimInputGrid[1]*dimInputGrid[2]*3;i+=3 ){
        inputGrid[i+0] = i ;
        inputGrid[i+1] = 2.*i;
        inputGrid[i+2] = 3.*i;
    }
    double intTime = integratePoint(loc2,&inputGrid[0], &originInput[0], &dimInputGrid[0],
                                    &spacingInput[0],1,stepSize,stagnation);

}
*/

///changed 3D to 2D
TEST_CASE("Test Integration", "[Test integration]"){

    vec2 loc1,loc2;
    vec2set(loc1, 0,0);
    vec2set(loc2, 3.9,1.6);
    vec2 dataVec;

    double originInput[2]= {0,0};
    double originSeed[2] = {3,1};
    double spacingInput[2] = {1,1};
    double spacingSeed[2] = {9,9};
    int dimInputGrid[2] = {401,501};
    int dimSeedGrid[2] = {2,2};
    double stepSize = .1;
    double stagnation = LIMIT_DOUBLE;
    float* inputGrid = (float*) malloc(sizeof(float) *
                                              dimInputGrid[0]*dimInputGrid[1]*2);

    for(int i = 0; i <dimInputGrid[0]*dimInputGrid[1]*2;i+=2 ){
        inputGrid[i+0] = i ;
        inputGrid[i+1] = 2.*i;
        ////inputGrid[i+2] = 3.*i;
    }
    double intTime = integratePoint(loc2,&inputGrid[0], &originInput[0], &dimInputGrid[0],
                                    &spacingInput[0],1,stepSize,stagnation);

}


/*
TEST_CASE("Interpolat/Locate Cell/Get input data", "[InterPolate]")
{
    vec3 location, location1;
    vec3set(location, 3.5,7.5,3.5);
    vec3set(location1, 12,9,9);
    vec3 locationLocal;
    vec3set(locationLocal,.5,.5,.5);
    vec3 dataVec;
    int dimInputGrid[3] = {10,10,10};
    float inputGrid[dimInputGrid[0]*dimInputGrid[1]*dimInputGrid[2]*3];
    double origin[3] = {0.,0.,0.};
    double spacingInputGrid[3] = {1.,1.,1.};
    int direction = 0;
    for(int i = 0; i<dimInputGrid[0]*dimInputGrid[1]*dimInputGrid[2]*3;i+=3 )
    {
        inputGrid[i+0] = i ;
        inputGrid[i+1] = 2.*i;
        inputGrid[i+2] = 3.*i;
    }

    CHECK(interpolate(location1,dataVec,inputGrid,origin,
                      dimInputGrid,spacingInputGrid,direction)==false);

    int xPos = 3;
    int yPos = 7;
    int zPos = 3;

    getDataPoint(dataVec,xPos,yPos,zPos,inputGrid,dimInputGrid);
    vec3 testVec;

    vec3set(testVec,inputGrid[getIndexInputGrid(xPos,yPos,zPos,dimInputGrid)+0],
            inputGrid[getIndexInputGrid(xPos,yPos,zPos,dimInputGrid)+1],
            inputGrid[getIndexInputGrid(xPos,yPos,zPos,dimInputGrid)+2]);

    CHECK(vec3eq(dataVec,testVec));
    vec3set(dataVec,0,0,0);
    interpolate(location,dataVec,inputGrid,origin,dimInputGrid,spacingInputGrid,direction);

    vec3 a,b,c,d,e,f,g,h;
    getDataPoint(a,xPos,yPos,zPos,inputGrid, dimInputGrid);
    getDataPoint(b,xPos+1,yPos,zPos,inputGrid, dimInputGrid);
    getDataPoint(c,xPos,yPos+1,zPos,inputGrid, dimInputGrid);
    getDataPoint(d,xPos+1,yPos+1,zPos,inputGrid, dimInputGrid);
    getDataPoint(e,xPos,yPos,zPos+1,inputGrid, dimInputGrid);
    getDataPoint(f,xPos+1,yPos,zPos+1,inputGrid, dimInputGrid);
    getDataPoint(g,xPos,yPos+1,zPos+1,inputGrid, dimInputGrid);
    getDataPoint(h,xPos+1,yPos+1,zPos+1,inputGrid, dimInputGrid);

    double xLocal = .5;
    double yLocal = .5;
    double zLocal = .5;
    vec3trilint(a,b,c,d,e,f,g,h, xLocal, yLocal, zLocal, testVec);
    REQUIRE(vec3eq(dataVec,testVec));

    double stepSize = 100;
    double stagnation  = -10.;
    vec3 loc;
    vec3copy(location,loc);
    double intTime = 1;
    double currentIntTime = 0;
    while(currentIntTime < intTime) {
        currentIntTime += integratePoint(location, inputGrid, origin,
                                  dimInputGrid, spacingInputGrid,
                                  direction, stepSize, stagnation);
           
    }
    CHECK(-10. != stagnation);

    for(int i = 0; i<dimInputGrid[0]*dimInputGrid[1]*dimInputGrid[2]*3;i+=3 )
    {
        inputGrid[i+0] = LIMIT_DOUBLE/10.;
        inputGrid[i+1] = LIMIT_DOUBLE/10.;
        inputGrid[i+2] = LIMIT_DOUBLE/10.;
    }

    while(currentIntTime < intTime) {
        intTime += integratePoint(location, inputGrid, origin,
                                  dimInputGrid, spacingInputGrid,
                                  direction, stepSize, stagnation);
        if(stagnation <= LIMIT_DOUBLE)
            break;
    }
    CHECK(stagnation <= LIMIT_DOUBLE );

    origin[0] = 100.;
    origin[1] = 100.;
    origin[2] = 100.;
    intTime += integratePoint(location, inputGrid, origin,
                              dimInputGrid, spacingInputGrid,
                              direction, stepSize, stagnation);

    CHECK(stagnation == Approx(0.));
    //CHECK(intTime == stepSize/vec3mag(testVec));

    vec3nrm(testVec,testVec);
    vec3scal(testVec,stepSize,testVec);
    vec3add(loc,testVec,loc);
    //CHECK(vec3eq(loc,location));




}
*/


////changed 3D to 2D
TEST_CASE("Interpolat/Locate Cell/Get input data", "[InterPolate]")
{
    vec2 location, location1;
    vec2set(location, 3.5,7.5);
    vec2set(location1, 12,9);
    vec2 locationLocal;
    vec2set(locationLocal,.5,.5);
    vec2 dataVec;
    int dimInputGrid[2] = {10,10};
    float inputGrid[dimInputGrid[0]*dimInputGrid[1]*2];
    double origin[2] = {0.,0.};
    double spacingInputGrid[2] = {1.,1.};
    int direction = 0;
    for(int i = 0; i<dimInputGrid[0]*dimInputGrid[1]*2;i+=2 )
    {
        inputGrid[i+0] = i ;
        inputGrid[i+1] = 2.*i;
        ////inputGrid[i+2] = 3.*i;
    }

    CHECK(interpolate(location1,dataVec,inputGrid,origin,
                      dimInputGrid,spacingInputGrid,direction)==false);

    int xPos = 3;
    int yPos = 7;
    ////int zPos = 3;

    ////removed zpos
    getDataPoint(dataVec,xPos,yPos,inputGrid,dimInputGrid);
    vec2 testVec;

    vec2set(testVec,inputGrid[getIndexInputGrid(xPos,yPos,dimInputGrid)+0],
            inputGrid[getIndexInputGrid(xPos,yPos,dimInputGrid)+1]);

    CHECK(vec2eq(dataVec,testVec));
    vec2set(dataVec,0,0);
    interpolate(location,dataVec,inputGrid,origin,dimInputGrid,spacingInputGrid,direction);

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

    double xLocal = .5;
    double yLocal = .5;
    ////double zLocal = .5;
    ////changed tri to bi
    vec2bilint(a,b,c,d, xLocal, yLocal, testVec);
    REQUIRE(vec2eq(dataVec,testVec));

    double stepSize = 100;
    double stagnation  = -10.;
    vec2 loc;
    vec2copy(location,loc);
    double intTime = 1;
    double currentIntTime = 0;
    while(currentIntTime < intTime) {
        currentIntTime += integratePoint(location, inputGrid, origin,
                                  dimInputGrid, spacingInputGrid,
                                  direction, stepSize, stagnation);

    }
    CHECK(-10. != stagnation);

    for(int i = 0; i<dimInputGrid[0]*dimInputGrid[1]*2;i+=2 )
    {
        inputGrid[i+0] = LIMIT_DOUBLE/10.;
        inputGrid[i+1] = LIMIT_DOUBLE/10.;
        ////inputGrid[i+2] = LIMIT_DOUBLE/10.;
    }

    while(currentIntTime < intTime) {
        intTime += integratePoint(location, inputGrid, origin,
                                  dimInputGrid, spacingInputGrid,
                                  direction, stepSize, stagnation);
        if(stagnation <= LIMIT_DOUBLE)
            break;
    }
    CHECK(stagnation <= LIMIT_DOUBLE );

    origin[0] = 100.;
    origin[1] = 100.;
    ////origin[2] = 100.;
    intTime += integratePoint(location, inputGrid, origin,
                              dimInputGrid, spacingInputGrid,
                              direction, stepSize, stagnation);

    CHECK(stagnation == Approx(0.));
    //CHECK(intTime == stepSize/vec3mag(testVec));

    vec2nrm(testVec,testVec);
    vec2scal(testVec,stepSize,testVec);
    vec2add(loc,testVec,loc);
    //CHECK(vec3eq(loc,location));




}
