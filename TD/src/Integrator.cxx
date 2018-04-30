#include "Integrator.h"
#include "vtkCellArray.h"
#include "vtkCell.h"

#include <vtkDataArray.h>
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include <math.h>


//int iter1 = 0;
//int iter2 = 0;
//int iter3 = 0;
int iter1,iter2,iter3,iter4, iter5, iter6 = 0;



////changed 3D to 2D
int getIndex1(int x,int y, int* dataDims)
{
    return (y*dataDims[0]+x);
}

////changed 3D to 2D
void vec2print(vec2 a)
{
    std::cout<< a[0] << " " << a[1]<< " ";
}


Integrator::Integrator()

{

}


Integrator::Integrator(int intDirection, double stepsize, double starttime, double endtime, double stagThresh,
                       vtkSmartPointer<vtkImageData> inputGrid, vtkSmartPointer<vtkImageData> sourcePoints, vtkSmartPointer<vtkImageData> nextinputGrid,
                        int blockNum,  int initialnum,  int finalnum, double* time, std::vector<vector<double>> previousPts)
{

    this->integrationDirection = intDirection;
    this->stepSize = stepsize;
    this->integrationTime = (starttime + endtime) * intDirection;
    this->stagnationThreshold = stagThresh;
    this->input = inputGrid;
    this->source = sourcePoints;
    this->origin[0] = 0;
    this->origin[1] = 0;
    this->nextinput = nextinputGrid;
    this->iter = blockNum;
    //this->currentIntLength = 0;
    this->startTime = starttime;
    this->endTime = endtime;
    this->initialNum = initialnum;
    this->finalNum = finalnum;
    this->timestep = time;


    if (blockNum != this->initialNum) {
        this->NewEndPoints = previousPts;
    }

    if(!this->input || ! this->source)
    {
        std::cerr<<" Input or source not set correctly"<< std::endl;
    }
}


Integrator::~Integrator() {}

void Integrator::setOrigin(double *ori) {
    this->origin[0] = ori[0];
    this->origin[1] = ori[1];
}


void Integrator::setOriginSource(double *ori) {
    this->originSource[0] = ori[0];
    this->originSource[1] = ori[1];
}


////changed 3D to 2D, removed components e,f,g,h (from cube to square), change the name from tri to bi
void Integrator::bilinInterpolator(vec2 a,vec2 b, vec2 c, vec2 d,
                                   double spacingX, double spacingY,
                                    vec2 location, vec2& output)
{

    ////changed r,s,t to r,s
    double r,s;
    double t1 = location[0];
    double t2 = location[1];
    ////double t3 = location[2];


    r= t1/(spacingX);
    s= t2/(spacingY);
    ////t= t3/(spacingZ);


    if( fabs(r)< LIMIT_DOUBLE)
        r = 0;
    if( fabs(s)< LIMIT_DOUBLE)
        s = 0;
    /*
    if( fabs(t)< LIMIT_DOUBLE)
    {
        t = 0;
    }
    */

    //// 3D to 2D
    vec2bilint(a,b,c,d, r, s, output);

}

/*
bool Integrator::cellLocator(vec3& output, vec3 location)
{

    vec3 loc;
    vec3copy(location, loc);
    //std::cout<<"Cell Locator"<<std::endl;
    //vec3print(location);

    double spacing[3];
    this->input->GetSpacing(spacing);
    int dim[3];
    this->input->GetDimensions(dim);

    ///Translate to origin of inputGrid to calculate the correct cell index
   // loc[0] -= origin[0];
   // loc[1] -= origin[1];
   // loc[2] -= origin[2];

    int xIndex = 0;
    xIndex = static_cast<int>(loc[0]/spacing[0]);
    int yIndex = 0;
    yIndex = static_cast<int>(loc[1]/spacing[1]);
    int zIndex = 0;
    zIndex = static_cast<int>(loc[2]/spacing[2]);

    if(xIndex > dim[0]-2 || yIndex > dim[1]-2 || zIndex > dim[2]-2
       || xIndex <0 || yIndex < 0 || zIndex <0)
    {
        return false;

    }

    //std::cout<< xIndex << " "<< yIndex <<" "<< zIndex<<std::endl;
    vec3 a,b,c,d,e,f,g,h;
    double p[3];
    vtkIdType id;
    vtkSmartPointer<vtkDataArray> data = this->input->GetPointData()->GetArray(0);
    if(!data)
    {
        std::cerr<<"error getting data array"<<std::endl;
        return false;
    }

    int index = getIndex1(xIndex, yIndex,zIndex, dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(a,p[0],p[1],p[2]);


    index = getIndex1(xIndex+1, yIndex,zIndex,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(b,p[0],p[1],p[2]);

    index = getIndex1(xIndex, yIndex+1,zIndex,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(c,p[0],p[1],p[2]);

    index = getIndex1(xIndex+1, yIndex+1,zIndex,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(d,p[0],p[1],p[2]);

    index = getIndex1(xIndex, yIndex,zIndex+1,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(e,p[0],p[1],p[2]);

    index = getIndex1(xIndex+1, yIndex,zIndex+1,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(f,p[0],p[1],p[2]);

    index = getIndex1(xIndex, yIndex+1,zIndex+1,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(g,p[0],p[1],p[2]);

    index = getIndex1(xIndex+1, yIndex+1,zIndex+1,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    vec3set(h,p[0],p[1],p[2]);

    ///Translate location to local cell coordinates
    loc[0] -= double(xIndex) * spacing[0];
    loc[1] -= double(yIndex) * spacing[1];
    loc[2] -= double(zIndex) * spacing[2];

    this->trilinInterpolator(a,b,c,d,e,f,g,h,spacing[0],spacing[1],spacing[2],loc,output);



    return true;

}
*/

////changed 3D to 2D
bool Integrator::cellLocator(vec2& output, vec2& nextoutput, vec2 location)
{

    ////changed 3 to 2
    vec2 loc;
    vec2copy(location, loc);
    //std::cout<<"Cell Locator"<<std::endl;
    //vec3print(location);

    ////changed 3 to 2
    double spacing[2];                        ////
    this->input->GetSpacing(spacing);
    int dim[2];                               ////
    this->input->GetDimensions(dim);


    double nextspacing[2];                        ////
    this->nextinput->GetSpacing(nextspacing);
    int nextdim[2];                               ////
    this->nextinput->GetDimensions(nextdim);

    ///Translate to origin of inputGrid to calculate the correct cell index
   // loc[0] -= origin[0];
   // loc[1] -= origin[1];
   // loc[2] -= origin[2];

    int xIndex = 0;
    xIndex = static_cast<int>(loc[0]/spacing[0]);
    int yIndex = 0;
    yIndex = static_cast<int>(loc[1]/spacing[1]);
    ////int zIndex = 0;
    ////zIndex = static_cast<int>(loc[2]/spacing[2]);

    int nextxIndex = 0;
    nextxIndex = static_cast<int>(loc[0]/nextspacing[0]);
    int nextyIndex = 0;
    nextyIndex = static_cast<int>(loc[1]/nextspacing[1]);

    ////changed 3 to 2 and removed e,f,g,h
    vec2 a,b,c,d;
    vec2 nexta,nextb,nextc,nextd;
    double p[2];
    double nextp[2];
    vtkIdType id;
    vtkSmartPointer<vtkDataArray> data = this->input->GetPointData()->GetArray(0);
    vtkSmartPointer<vtkDataArray> nextdata = this->nextinput->GetPointData()->GetArray(0);
    if(!data)
    {
        std::cerr<<"error getting data array"<<std::endl;
        return false;
    }

    if(xIndex > dim[0]-1 || yIndex > dim[1] -1
       || xIndex < 0 || yIndex < 0)
    {
        return false;
    }

    else if (xIndex == dim[0] - 1 && yIndex != dim[1] - 1)
    {
        int index = getIndex1(xIndex - 1, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(a, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nexta, nextp[0], nextp[1]);

        ////removed zIndex
        index = getIndex1(xIndex, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(b, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextb, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex -1, yIndex + 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(c, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextc, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex, yIndex + 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(d, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextd, nextp[0], nextp[1]);
    }

    else if ( xIndex != dim[0] - 1 && yIndex == dim[1] - 1) {

        int index = getIndex1(xIndex - 1, yIndex - 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(a, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nexta, nextp[0], nextp[1]);

        ////removed zIndex
        index = getIndex1(xIndex + 1, yIndex - 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(b, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextb, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(c, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextc, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex + 1, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(d, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextd, nextp[0], nextp[1]);

    }

    else if (xIndex == dim[0] - 1 && yIndex == dim[1] - 1)
    {
        int index = getIndex1(xIndex - 1, yIndex - 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(a, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nexta, nextp[0], nextp[1]);

        ////removed zIndex
        index = getIndex1(xIndex, yIndex - 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(b, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextb, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex - 1, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(c, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextc, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(d, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextd, nextp[0], nextp[1]);
    }
    else {

        ////removed zIndex
        int index = getIndex1(xIndex, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(a, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nexta, nextp[0], nextp[1]);

        ////removed zIndex
        index = getIndex1(xIndex + 1, yIndex, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(b, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextb, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex, yIndex + 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(c, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextc, nextp[0], nextp[1]);


        ////removed zIndex
        index = getIndex1(xIndex + 1, yIndex + 1, dim);
        id = vtkIdType(index);
        this->input->GetPoint(id, p);
        data->GetTuple(id, p);
        ////changed 3 to 2 and removed p[2]
        vec2set(d, p[0], p[1]);

        this->nextinput->GetPoint(id, nextp);
        nextdata->GetTuple(id, nextp);
        vec2set(nextd, nextp[0], nextp[1]);
    }

    ///Translate location to local cell coordinates
    loc[0] -= double(xIndex) * spacing[0];
    loc[1] -= double(yIndex) * spacing[1];
    //loc[2] -= double(zIndex) * spacing[2];

    ////removed e,f,g,h and spacing[2] and changed the name from tri to bi
    this->bilinInterpolator(a,b,c,d,spacing[0],spacing[1],loc,output);
    this->bilinInterpolator(nexta, nextb, nextc, nextd, nextspacing[0], nextspacing[1], loc, nextoutput);


    return true;

}


////changed 3 to 2
bool Integrator::pointInterpolator(vec2 currentLocation, vec2& vecNextLocation, vec2& vecNextSliceLocation)
{
    vec2 location;
    vec2copy(currentLocation,location);
    bool result = this->cellLocator(vecNextLocation, vecNextSliceLocation, location);
    if(this->integrationDirection < 0)
    {
        vecNextLocation[0] *=-1;
        vecNextLocation[1] *=-1;
        ////vecNextLocation[2] *=-1;
        vecNextSliceLocation[0] *=-1;
        vecNextSliceLocation[1] *=-1;
    }

    return result;

}


vtkSmartPointer<vtkPolyData> Integrator::integrateRK4()                                             ////changed from vtkPolyData
{

    vtkSmartPointer<vtkCellArray> outputLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();


    if(!this->source || ! this->input)
    {
        std::cerr<<"Source or Input not set correctly"<<std::endl;
        return output;
    }


    //// if it's first iteration, then take the starting points as defined in the grid
    if( this->iter == this->initialNum ) {

        int numOfSeedPoints = this->source->GetNumberOfPoints();
        std::cout << "num of seed points" << numOfSeedPoints << std::endl;

        vtkSmartPointer<vtkPoints> line;

        int currentProgress = 0;

        for (int i = 0; i < numOfSeedPoints; i++)
        {
            if (int((double(i) / double(numOfSeedPoints)) * 100.0) >= currentProgress) {
                std::cout << " Progress: " << int((double(i) / double(numOfSeedPoints)) * 100.0)
                          << "%\r";
                std::cout.flush();
                ++currentProgress;
            }
            ////changed 3 to 2
            vec2 startLocation;
            double p[2];
            double point[3];
            vtkIdType currentId = vtkIdType(i);

            this->source->GetPoint(currentId, p);

            std::vector<double> startLoc = {p[0], p[1]};
            this->StartPoint.push_back(startLoc);

            ////changed 3 to 2 and removed p[2]
            vec2set(startLocation, p[0], p[1]);

            line = this->integratePointRK4(startLocation);
            //std::cout<< "Number of Points in stream line "<< i<< " " <<line->GetNumberOfPoints()<<std::endl;

            outputLines->InsertNextCell(line->GetNumberOfPoints());

            for (int j = 0; j < line->GetNumberOfPoints(); j++) {
                vtkIdType id = vtkIdType(j);
                line->GetPoint(id, point);
                outputPoints->InsertNextPoint(point);
                vtkIdType nextPoint = outputPoints->InsertNextPoint(point);
                outputLines->InsertCellPoint(nextPoint);

            }
        }
    }

    //// After the first iteration, change the input as end points of last iteration
    else
    {

        int numOfSeedPoints = NewEndPoints.size();
        std::cout << "num of seed points" << numOfSeedPoints << std::endl;

        vtkSmartPointer<vtkPoints> line;

        int currentProgress = 0;

        for (int i = 0; i < numOfSeedPoints; i++)
        {
            if (int((double(i) / double(numOfSeedPoints)) * 100.0) >= currentProgress) {
                std::cout << " Progress: " << int((double(i) / double(numOfSeedPoints)) * 100.0)
                          << "%\r";
                std::cout.flush();
                ++currentProgress;
            }
            ////changed 3 to 2
            vec2 startLocation;
            std::vector<double> p = NewEndPoints.at(i);

            double point[3];

            vec2set(startLocation, p[0], p[1]);

            line = this->integratePointRK4(startLocation);

            outputLines->InsertNextCell(line->GetNumberOfPoints());

            for (int j = 0; j < line->GetNumberOfPoints(); j++) {
                vtkIdType id = vtkIdType(j);
                line->GetPoint(id, point);
                outputPoints->InsertNextPoint(point);
                vtkIdType nextPoint = outputPoints->InsertNextPoint(point);
                outputLines->InsertCellPoint(nextPoint);
            }
        }

    }

    std::cout << "iter 1 " <<iter1 << " iter 2 " << iter2 << " iter 3 " << iter3 << " iter 4 " << iter4<< " iter 5 " << iter5 <<" iter 6 " << iter6<<std::endl;
    vtkIdType numPoints = outputPoints->GetNumberOfPoints();
    std::cout<< "number of points " << numPoints <<std::endl;

    vtkIdType numCells = outputLines->GetNumberOfCells();
    std::cout<< "number of cells in intergrator " << numCells <<std::endl;

    output->SetPoints(outputPoints);
    output->SetLines(outputLines);

    return output;

}

////changed 3 to 2
vtkSmartPointer<vtkPoints> Integrator::integratePointRK4(vec2 startLocation)
{

    vec2 currentLocation;
    vec2 velCurrLocation;
    vec2 currentLocationSave;
    vec2 traverseDist;
    vec2 temp;
    vec2copy(startLocation,currentLocation);
    vec2 velNextLocation;
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    double intFac;
    double currentIntTime = 0.;
    double currentIntLength = this->timestep[this->iter] * this->integrationDirection;                             /// start from the particular iteration time step
    double timeDiff = (this->timestep[this->iter + 1 * this->integrationDirection ] - this->timestep[this->iter]) * this->integrationDirection;                ////                /// Time gap between 2 consecutive intervals

    vec2 k1,k2,k3,k4;
    vec2 k1norm, k2norm, k3norm, k4norm;
    vec2 k1scal, k2scal, k3scal, k4scal;
    vec2 step;
    //vec2 k1, k1scal;

    double interStepWeights[3] = {.5,.5,1.0};
    double finalStepWeights[4] = {1.0/6.0, 1.0/3.0,1.0/3.0,1.0/6.0};


    ////setting initial currenttime only for the first iteration
    if (this->iter == this->initialNum)
    {
        currentIntTime = (this->startTime - this->timestep[this->iter]) * this->integrationDirection;
        currentIntLength = this->startTime * this->integrationDirection;
    }


    ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
    double point[3] = {currentLocation[0] + this->origin[0],
                       currentLocation[1] + this->origin[1],0.};

    vtkIdType id = outputPoints->InsertNextPoint(point);
    id = outputPoints->InsertNextPoint(point);

    bool pointValid = true;


    /// for every iteration run between the last time step and for the last iteration stop at the integration time

    while ( currentIntTime < timeDiff && currentIntLength < this->integrationTime )
    {
        pointValid = true;
        vec2copy(currentLocation, currentLocationSave);

        if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {

            iter1 += 1;
            std::cout <<"outside loop " <<  currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << std::endl;
            //pointValid = false;
            //break;
        }


        int iterator = 0;

        //// very first Starting point, first step, need to calculate the new timestep
        if (this->iter == this->initialNum && currentIntTime == (this->startTime - this->timestep[this->iter]) * this->integrationDirection )
        {

            int count = 0;
            double starttimevalue = fabs(this->startTime);
            starttimevalue = starttimevalue - static_cast<int>(starttimevalue);

            while (starttimevalue != 0)
            {
                starttimevalue = starttimevalue * 10;
                count += 1;
                starttimevalue = starttimevalue - static_cast<int>(starttimevalue);
            }

            double nearestValue = floor(this->startTime * pow(10.0, count - 1)) / pow(10.0, count - 1) + stepSize;
            double firstTimestep = nearestValue - this->startTime;

            if (firstTimestep < LIMIT_DOUBLE)
            {
                firstTimestep = this->stepSize;
            }

            intFac = currentIntTime / timeDiff;

            vec2 vel;
            vec2 nextvel;
            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k1);
            //vec2nrm(k1,k1norm);
            vec2scal(k1,firstTimestep,k1scal);

            vec2copy(k1scal,temp);
            vec2scal(temp, interStepWeights[0],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);


            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter2 += 1;
                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);

            intFac = (currentIntTime + (firstTimestep / 2.) )/ timeDiff;

            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k2);
            //vec2nrm(k2,k2norm);
            vec2scal(k2,firstTimestep,k2scal);

            vec2copy(k2scal,temp);
            vec2scal(temp, interStepWeights[1],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter3+=1;
                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k3);
            //vec2nrm(k3,k3norm);
            vec2scal(k3,firstTimestep,k3scal);

            vec2copy(k3scal,temp);
            vec2scal(temp, interStepWeights[2],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter4+=1;

                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            intFac = (currentIntTime + firstTimestep) / timeDiff;                                        ///
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k4);
            //vec2nrm(k4,k4norm);
            vec2scal(k4,firstTimestep,k4scal);


            vec2copy(k1scal,vel);
            vec2scal(k1scal,finalStepWeights[0],vel);
            vec2scal(k2scal,finalStepWeights[1],temp);
            vec2add(vel,temp,vel);
            vec2scal(k3scal,finalStepWeights[2],temp);
            vec2add(vel,temp,vel);
            vec2scal(k4scal,finalStepWeights[3],temp);
            vec2add(vel,temp,vel);
            currentIntTime += firstTimestep;
            currentIntLength += firstTimestep;

            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, vel, currentLocation);
            //vec2add(currentLocation, step, currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter5+=1;
                std::cout <<"iter5 first loop " <<  currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << std::endl;
                // vec3print(currentLocation);
                // std::cout<< "k"<<std::endl;
                pointValid = false;
                break;
            }

        }

            //// for the last point, need to again calculate the new timestep
        else if ( this->integrationTime - currentIntLength < this->stepSize )
        {

            double lastTimestep = this->integrationTime - currentIntLength;    //// new step size

            intFac = currentIntTime / timeDiff;

            vec2 vel;
            vec2 nextvel;
            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k1);
            //vec2nrm(k1,k1norm);
            vec2scal(k1,lastTimestep,k1scal);

            vec2copy(k1scal,temp);
            vec2scal(temp, interStepWeights[0],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);


            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter2 += 1;
                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);

            intFac = (currentIntTime + lastTimestep / 2. )/ timeDiff;

            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k2);
            //vec2nrm(k2,k2norm);
            vec2scal(k2,lastTimestep,k2scal);

            vec2copy(k2scal,temp);
            vec2scal(temp, interStepWeights[1],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter3+=1;
                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k3);
            //vec2nrm(k3,k3norm);
            vec2scal(k3,lastTimestep,k3scal);

            vec2copy(k3scal,temp);
            vec2scal(temp, interStepWeights[2],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter4+=1;

                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            intFac = (currentIntTime + lastTimestep) / timeDiff;                                        ///
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k4);
            //vec2nrm(k4,k4norm);
            vec2scal(k4,lastTimestep,k4scal);


            vec2copy(k1scal,vel);
            vec2scal(k1scal,finalStepWeights[0],vel);
            vec2scal(k2scal,finalStepWeights[1],temp);
            vec2add(vel,temp,vel);
            vec2scal(k3scal,finalStepWeights[2],temp);
            vec2add(vel,temp,vel);
            vec2scal(k4scal,finalStepWeights[3],temp);
            vec2add(vel,temp,vel);
            currentIntTime += lastTimestep;
            currentIntLength += lastTimestep;

            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, vel, currentLocation);
            //vec2add(currentLocation, step, currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter5+=1;
                std::cout <<"iter5 second loop " <<  currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << std::endl;
                // vec3print(currentLocation);
                // std::cout<< "k"<<std::endl;
                pointValid = false;
                break;
            }

        }
            //// implementing RK4
            //// rest points

        else {

            intFac = currentIntTime / timeDiff;

            vec2 vel;
            vec2 nextvel;
            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k1);
            //vec2nrm(k1,k1norm);
            vec2scal(k1,this->stepSize,k1scal);

            vec2copy(k1scal,temp);
            vec2scal(temp, interStepWeights[0],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);


            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter2 += 1;
                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);

            intFac = (currentIntTime + this->stepSize / 2. )/ timeDiff;

            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k2);
            //vec2nrm(k2,k2norm);
            vec2scal(k2,this->stepSize,k2scal);

            vec2copy(k2scal,temp);
            vec2scal(temp, interStepWeights[1],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter3+=1;
                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k3);
            //vec2nrm(k3,k3norm);
            vec2scal(k3,this->stepSize,k3scal);

            vec2copy(k3scal,temp);
            vec2scal(temp, interStepWeights[2],temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation,temp,currentLocation);

            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter4+=1;

                //vec3print(currentLocation);
                // std::cout<< "k2"<<std::endl;
                pointValid = false;
                break;
            }

            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            intFac = (currentIntTime + this->stepSize) / timeDiff;                                        ///
            vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
            vec2copy(temp, vel);

            vec2copy(vel, k4);
            //vec2nrm(k4,k4norm);
            vec2scal(k4,this->stepSize,k4scal);


            vec2copy(k1scal,vel);
            vec2scal(k1scal,finalStepWeights[0],vel);
            vec2scal(k2scal,finalStepWeights[1],temp);
            vec2add(vel,temp,vel);
            vec2scal(k3scal,finalStepWeights[2],temp);
            vec2add(vel,temp,vel);
            vec2scal(k4scal,finalStepWeights[3],temp);
            vec2add(vel,temp,vel);
            currentIntTime += this->stepSize;
            currentIntLength += this->stepSize;

            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, vel, currentLocation);
            //vec2add(currentLocation, step, currentLocation);


            if(!this->pointInterpolator(currentLocation,velCurrLocation, velNextLocation))
            {
                iter5+=1;
                //std::cout << "third loop " << currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << " " << this->timestep[this->iter] <<std::endl;

                //std::cout<< "step "<< step[0] <<" " << step[1] <<std::endl;
                //std::cout<< "temp "<< temp[0] <<" " << temp[1] <<std::endl;

                // vec3print(currentLocation);
                // std::cout<< "k"<<std::endl;
                pointValid = false;
                break;
            }

            //iter3 += 1;
        }

        if (currentIntTime >= timeDiff || currentIntLength >= this->integrationTime)
        {
            //std::cout<< currentLocation[0] << " " << currentLocation[1] << std::endl;
            std::vector<double> currloc = {currentLocation[0], currentLocation[1]};
            this->PathlineEndPoints.push_back(currloc);
        }


        if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {

            iter6 += 1;
            //mstd::cout << "outside loop " << currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << " " << this->timestep[this->iter] <<std::endl;
            pointValid = false;
            break;
        }

        ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
        /// Add each point of the integration to the streamline if flag is set
        if (this->calcStreamlines) {
            double point2[2] = {currentLocation[0] + this->origin[0],
                                currentLocation[1] + this->origin[1]};
            vtkIdType id = outputPoints->InsertNextPoint(point2);
        }


    }


    ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
    ///If flag is not set only set end location of integration


    if(pointValid) {
        double point2[3] = {currentLocation[0] + this->origin[0] ,
                            currentLocation[1] + this->origin[1], 0.};
        id = outputPoints->InsertNextPoint(point2);
    }

    return outputPoints;


}

/*
////changed 3 to 2
vtkSmartPointer<vtkPoints> Integrator::integratePointRK4(vec2 startLocation)
{

    vec2 currentLocation;
    vec2 velCurrLocation;
    vec2 currentLocationSave;
    vec2 traverseDist;
    vec2 temp;
    vec2copy(startLocation,currentLocation);
    vec2 velNextLocation;
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    double intFac;
    double currentIntTime = 0.;
    double currentIntLength = this->timestep[this->iter] * this->integrationDirection;                             /// start from the particular iteration time step
    double timeDiff = (this->timestep[this->iter + 1 * this->integrationDirection ] - this->timestep[this->iter]) * this->integrationDirection;                ////                /// Time gap between 2 consecutive intervals


    ////setting initial currenttime only for the first iteration
    if (this->iter == this->initialNum)
    {
        currentIntTime = (this->startTime - this->timestep[this->iter]) * this->integrationDirection;
        currentIntLength = this->startTime * this->integrationDirection;
    }


    ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
    double point[2] = {currentLocation[0] + this->origin[0],
                        currentLocation[1] + this->origin[1]};

    vtkIdType id = outputPoints->InsertNextPoint(point);
    id = outputPoints->InsertNextPoint(point);

    bool pointValid = true;


    /// for every iteration run between the last time step and for the last iteration stop at the integration time

    while ( currentIntTime < timeDiff && currentIntLength < this->integrationTime )
    {
        pointValid = true;
        vec2copy(currentLocation, currentLocationSave);

        if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {

            //iter2 += 1;
            //std::cout << currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << std::endl;
            pointValid = false;
            break;
        }

        intFac = currentIntTime / timeDiff;

        vec2 vel;
        vec2 nextvel;
        vec2copy(velCurrLocation, vel);
        vec2copy(velNextLocation, nextvel);
        vec2lint(vel, nextvel, intFac, temp);                                                   ///// edited
        vec2copy(temp, vel);

        int iterator = 0;

        //// very first Starting point, first step, need to calculate the new timestep
        if (this->iter == this->initialNum && currentIntTime == (this->startTime - this->timestep[this->iter]) * this->integrationDirection )
        {
            int count = 0;
            double starttimevalue = fabs(this->startTime);
            starttimevalue = starttimevalue - static_cast<int>(starttimevalue);

            while (starttimevalue != 0)
            {
                starttimevalue = starttimevalue * 10;
                count += 1;
                starttimevalue = starttimevalue - static_cast<int>(starttimevalue);
            }

            double nearestValue = floor(this->startTime * pow(10.0, count - 1)) / pow(10.0, count - 1) + stepSize;
            double firstTimestep = nearestValue - this->startTime;

            if (firstTimestep < LIMIT_DOUBLE)
            {
                firstTimestep = this->stepSize;
            }

            vec2scal(vel, firstTimestep, traverseDist);

            vec2copy(traverseDist, temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, temp, currentLocation);

            currentIntTime += firstTimestep;
            currentIntLength += firstTimestep;

            //iter1 += 1;

        }

            //// for the last point, need to again calculate the new timestep
        else if ( this->integrationTime - currentIntLength < this->stepSize )
        {
            double lastTimestep = this->integrationTime - currentIntLength;    //// new step size


            vec2scal(vel, lastTimestep, traverseDist);

            vec2copy(traverseDist, temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, temp, currentLocation);

            currentIntTime += lastTimestep;
            currentIntLength += lastTimestep;


        }

        //// rest points
        else {

            vec2scal(vel, stepSize, traverseDist);

            vec2copy(traverseDist, temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, temp, currentLocation);

            currentIntTime += this->stepSize;
            currentIntLength += this->stepSize;

            //iter3 += 1;
        }

        if (currentIntTime >= timeDiff || currentIntLength >= this->integrationTime)
        {
            std::vector<double> currloc = {currentLocation[0], currentLocation[1]};
            this->PathlineEndPoints.push_back(currloc);
        }

        if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {

            //iter2 += 1;
            //std::cout << currentLocation[0] << " " << currentLocation[1] << " " << startLocation[0] << " " << startLocation[1] << " " << this->timestep[this->iter] <<std::endl;
            pointValid = false;
            break;
        }

        ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
        /// Add each point of the integration to the streamline if flag is set
        if (this->calcStreamlines) {
            double point2[2] = {currentLocation[0] + this->origin[0],
                                    currentLocation[1] + this->origin[1]};
            vtkIdType id = outputPoints->InsertNextPoint(point2);
        }


    }


    ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
    ///If flag is not set only set end location of integration


    if(pointValid) {
        double point2[2] = {currentLocation[0] + this->origin[0] ,
                            currentLocation[1] + this->origin[1]};
        id = outputPoints->InsertNextPoint(point2);
    }

    return outputPoints;


}
*/

void Integrator::integrateRK4GPUWrapper()
{
    ///TODO check if this is still needed
    //this->originSource[0] -= this->origin[0];
    //this->originSource[1] -= this->origin[1];
    //this->originSource[2] -= this->origin[2];

    CudaIntegrator gpuIntegrator = CudaIntegrator(this->source->GetDimensions(),
                                                  this->input->GetDimensions(), this->origin,
                                                  this->originSource,
                                                  this->source->GetNumberOfPoints(),
                                                  this->input->GetNumberOfPoints(),
                                                  this->source->GetSpacing(),
                                                  this->input->GetSpacing(),
                                                  this->integrationDirection,
                                                  this->stepSize, this->integrationTime,
                                                  this->stagnationThreshold,
                                                  this->maxNumberOfSteps);


    auto floatInput = vtkSmartPointer<vtkFloatArray>::New();
    ////changed 3 to 2
    floatInput->SetNumberOfComponents(2);
    floatInput->SetNumberOfTuples(this->input->GetPointData()->GetArray(0)->GetNumberOfTuples());
    for(unsigned int index = 0; index < floatInput->GetNumberOfTuples();++index)
    {
        floatInput->SetTuple(index,this->input->GetPointData()->GetArray(0)->GetTuple(index));
    }
    float* inputArray = (float*)floatInput->GetVoidPointer(0);
    //int sizeFlowMap = sizeof(float) * this->source->GetNumberOfPoints()* 3;
    ////changed 3 to 2
    int sizeFlowMap = sizeof(double) * this->source->GetNumberOfPoints()* 2;
#ifndef UNITTEST
    //float* flowArray = (float*) malloc(sizeFlowMap);
    double* flowArray = (double*) malloc(sizeFlowMap);
    double* ftleArray = gpuIntegrator.integrate(inputArray,flowArray);
#else
    double* ftleArray;
    double* flowArray;
#endif

    //auto flowMap = vtkSmartPointer<vtkFloatArray>::New();
    auto flowMap = vtkSmartPointer<vtkDoubleArray>::New();
    ////changed 3 to 2
    flowMap->SetNumberOfComponents(2);
    //Black Pointer Voodoo Magic
#ifndef CPUEXEC
    auto ftleField = vtkSmartPointer<vtkDoubleArray>::New();
    ftleField->SetVoidArray(ftleArray, this->source->GetNumberOfPoints(),0);
    ftleField->SetName("FTLE");
    this->source->GetPointData()->AddArray(ftleField);
#endif
    ////changed 3 to 2
    flowMap->SetVoidArray(flowArray, this->source->GetNumberOfPoints()*2,0);
    //ftleField->Print(std::cout);
    flowMap->SetName("FlowMap");
    this->source->GetPointData()->AddArray(flowMap);

    //return ftleField;

}

