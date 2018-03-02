#include "Integrator.h"
#include "vtkCellArray.h"
#include "vtkCell.h"

#include <vtkDataArray.h>
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"


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

/*
Integrator::Integrator(int intDirection, double stepsize, double intTime, double stagThresh, vtkSmartPointer<vtkImageData> inputGrid, vtkSmartPointer<vtkImageData> sourcePoints)
{

    this->integrationDirection =  intDirection;
    this->stepSize = stepsize;
    this->integrationTime = intTime;
    this->integrationLength = intTime;
    this->stagnationThreshold = stagThresh;
    this->input = inputGrid;
    this->source = sourcePoints;
    this->origin[0] = 0;
    this->origin[1] = 0;
    this->origin[2] = 0;

    if(!this->input || ! this->source)
    {
        std::cerr<<" Input or source not set correctly"<< std::endl;

    }

}
*/

Integrator::Integrator(int intDirection, double stepsize, double starttime, double endtime, double stagThresh,
                       vtkSmartPointer<vtkImageData> inputGrid, vtkSmartPointer<vtkImageData> sourcePoints, vtkSmartPointer<vtkImageData> nextinputGrid,
                       unsigned int blockNum, unsigned int initialnum, unsigned int finalnum, double* time, std::vector<vector<double>> previousPts)
{

    this->integrationDirection = intDirection;
    this->stepSize = stepsize;
    //this->integrationTime = endtime - starttime;
    //this->integrationLength = endtime - starttime;
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

    if (blockNum != initialnum) {
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

    if(xIndex > dim[0]-2 || yIndex > dim[1]-2
       || xIndex < 0 || yIndex < 0)
    {
        return false;
    }


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

    ////removed zIndex
    int index = getIndex1(xIndex, yIndex, dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    ////changed 3 to 2 and removed p[2]
    vec2set(a,p[0],p[1]);

    this->nextinput->GetPoint(id,nextp);
    nextdata->GetTuple(id,nextp);
    vec2set(nexta,nextp[0],nextp[1]);

    ////removed zIndex
    index = getIndex1(xIndex+1, yIndex,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    ////changed 3 to 2 and removed p[2]
    vec2set(b,p[0],p[1]);

    this->nextinput->GetPoint(id,nextp);
    nextdata->GetTuple(id,nextp);
    vec2set(nextb,nextp[0],nextp[1]);


    ////removed zIndex
    index = getIndex1(xIndex, yIndex+1,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    ////changed 3 to 2 and removed p[2]
    vec2set(c,p[0],p[1]);

    this->nextinput->GetPoint(id,nextp);
    nextdata->GetTuple(id,nextp);
    vec2set(nextc,nextp[0],nextp[1]);


    ////removed zIndex
    index = getIndex1(xIndex+1, yIndex+1,dim);
    id = vtkIdType(index);
    this->input->GetPoint(id,p);
    data->GetTuple(id,p);
    ////changed 3 to 2 and removed p[2]
    vec2set(d,p[0],p[1]);

    this->nextinput->GetPoint(id,nextp);
    nextdata->GetTuple(id,nextp);
    vec2set(nextd,nextp[0],nextp[1]);


    ///Translate location to local cell coordinates
    loc[0] -= double(xIndex) * spacing[0];
    loc[1] -= double(yIndex) * spacing[1];
    //loc[2] -= double(zIndex) * spacing[2];

    ////removed e,f,g,h and spacing[2] and changed the name from tri to bi
    this->bilinInterpolator(a,b,c,d,spacing[0],spacing[1],loc,output);
    this->bilinInterpolator(nexta, nextb, nextc, nextd, nextspacing[0], nextspacing[1], loc, nextoutput);


    return true;

}

/*
bool Integrator::pointInterpolator(vec3 currentLocation, vec3& vecNextLocation)
{
    vec3 location;
    vec3copy(currentLocation,location);
    bool result = this->cellLocator(vecNextLocation, location);
    if(this->integrationDirection < 0)
    {
        vecNextLocation[0] *=-1;
        vecNextLocation[1] *=-1;
        vecNextLocation[2] *=-1;
    }

    return result;

}
*/
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
    }

    return result;

}

/*
vtkSmartPointer<vtkPolyData> Integrator::integrateRK4()
{


    vtkSmartPointer<vtkCellArray> outputLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
    if(!this->source || ! this->input)
    {
        std::cerr<<"Source or Input not set correctly"<<std::endl;
        return output;
    }

    int numOfSeedPoints = this->source->GetNumberOfPoints();
    vtkSmartPointer<vtkPoints> line;



    int currentProgress = 0;
    for(int i = 0; i < numOfSeedPoints; i++)
    {
        if(int((double(i)/double(numOfSeedPoints))  * 100.0) >= currentProgress)
        {
            std::cout << " Progress: " << int((double(i)/double(numOfSeedPoints))  * 100.0)
                      << "%\r";
            std::cout.flush();
            ++currentProgress;
        }
        vec3 startLocation;
        double p[3];
        double point[3];
        vtkIdType currentId = vtkIdType(i);
        this->source->GetPoint(currentId,p);
        vec3set(startLocation,p[0],p[1],p[2]);


        line =  this->integratePointRK4(startLocation);
        //std::cout<< "Number of Points in stream line "<< i<< " " <<line->GetNumberOfPoints()<<std::endl;
        outputLines->InsertNextCell(line->GetNumberOfPoints());
        for(int j = 0; j < line->GetNumberOfPoints(); j++)
        {
            vtkIdType id = vtkIdType(j);
            line->GetPoint(id,point);
            vtkIdType nextPoint = outputPoints->InsertNextPoint(point);
            outputLines->InsertCellPoint(nextPoint);


        }
    }

    output->SetPoints(outputPoints);
    //std::cout<<" \n"<<outputPoints->GetNumberOfPoints()<<std::endl;
    //if(outputPoints->GetNumberOfPoints() > 1)
    output->SetLines(outputLines);


    return output;

}

*/

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

    std::cout << "iterator" << (this->iter) <<std::endl;

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
            double point[2];
            vtkIdType currentId = vtkIdType(i);

            this->source->GetPoint(currentId, p);

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

        std::cout << "iteration "<<  this->iter + 1 << std::endl;

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

            double point[2];

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

    vtkIdType numPoints = outputPoints->GetNumberOfPoints();
    std::cout<< "number of points " << numPoints <<std::endl;

    vtkIdType numCells = outputLines->GetNumberOfCells();
    std::cout<< "number of cells in intergrator " << numCells <<std::endl;

    output->SetPoints(outputPoints);
    output->SetLines(outputLines);

    return output;

}

/*
vtkSmartPointer<vtkPoints> Integrator::integratePointRK4(vec3 startLocation)
{
    double currentIntTime = 0;
    double currentIntLength = 0;
    vec3 currentLocation;
    vec3 nextLocation;
    vec3 vecNextLocation;
    vec3 currentLocationSave;
    vec3 k1,k2,k3,k4;
    vec3 k1norm, k2norm, k3norm, k4norm;
    vec3 k1scal, k2scal, k3scal, k4scal;
    vec3 step;
    vec3 temp;
    vec3copy(startLocation,currentLocation);
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();

    double interStepWeights[3] = {.5,.5,1.0};
    double finalStepWeights[4] = {1.0/6.0, 1.0/3.0,1.0/3.0,1.0/6.0};

    // std::cout<<"Start Location"<<std::endl;
    // vec3print(startLocation);

    double point[3] = {currentLocation[0] + this->origin[0] ,
                        currentLocation[1] + this->origin[1],
                        currentLocation[2] + this->origin[2]};
    vtkIdType id = outputPoints->InsertNextPoint(point);
     id = outputPoints->InsertNextPoint(point);

    bool pointValid = true;
    while(currentIntTime < this->integrationTime)
        //while(currentIntLength < this->integrationLength)
    {
        pointValid = true;
        vec3copy(currentLocation, currentLocationSave);
        if(!this->pointInterpolator(currentLocation, vecNextLocation))
        {
            //vec3print(currentLocation);
            //std::cout<< "k1"<<std::endl;
            pointValid = false;
            break;
        }

        vec3copy(vecNextLocation,k1);
        vec3nrm(k1,k1norm);
        vec3scal(k1, stepSize,k1scal);

        vec3copy(k1scal,temp);
        vec3scal(temp, interStepWeights[0],temp);
        vec3copy(currentLocationSave, currentLocation);
        vec3add(currentLocation,temp,currentLocation);

        if(!this->pointInterpolator(currentLocation,vecNextLocation))
        {
            //vec3print(currentLocation);
            // std::cout<< "k2"<<std::endl;
            pointValid = false;
            break;pointInterpolator(currentLocation, vecNextLocation, vecNextSliceLocation))
        }
        vec3copy(vecNextLocation,k2);
        vec3nrm(k2,k2norm);
        vec3scal(k2,stepSize,k2scal);

        vec3copy(k2scal,temp);
        vec3scal(temp,interStepWeights[1],temp);
        vec3copy(currentLocationSave, currentLocation);
        vec3add(currentLocation,temp,currentLocation);

        if(!this->pointInterpolator(currentLocation,vecNextLocation))
        {
            // vec3print(currentLocation);
            // std::cout<< "k3"<<std::endl;
            pointValid = false;
            break;
        }

        vec3copy(vecNextLocation,k3);
        vec3nrm(k3,k3norm);
        vec3scal(k3,stepSize,k3scal);

        vec3copy(k3scal,temp);
        vec3scal(temp,interStepWeights[2],temp);
        vec3copy(currentLocationSave, currentLocation);
        vec3add(currentLocation,temp,currentLocation);

        if(!this->pointInterpolator(currentLocation,vecNextLocation))
        {
            // vec3print(currentLocation);
            // std::cout<< "k4"<<std::endl;
            pointValid = false;
            break;
        }

        vec3copy(vecNextLocation,k4);
        vec3nrm(k4,k4norm);
        vec3scal(k4,stepSize,k4scal);

        vec3scal(k1norm, finalStepWeights[0],temp);
        vec3copy(temp, step);
        vec3scal(k2norm, finalStepWeights[1], temp);
        vec3add(step, temp,step);
        vec3scal(k3norm, finalStepWeights[2], temp);
        vec3add(step, temp,step);
        vec3scal(k4norm, finalStepWeights[3],temp);
        vec3add(step,temp,step);
        vec3nrm(step,step);
        vec3scal(step, stepSize, step);
        currentIntLength += vec3mag(step);

        vec3 vel;
        vec3copy(k1,vel);
        vec3scal(k1,finalStepWeights[0],vel);
        vec3scal(k2,finalStepWeights[1],temp);
        vec3add(vel,temp,vel);
        vec3scal(k3,finalStepWeights[2],temp);
        vec3add(vel,temp,vel);
        vec3scal(k4,finalStepWeights[3],temp);
        vec3add(vel,temp,vel);
        currentIntTime += stepSize/vec3mag(vel);

        vec3copy(currentLocationSave, currentLocation);

        //vec3copy(k1,step);
        //vec3nrm(step,step);
        //vec3scal(step,stepSize,step);
        // std::cout<<"Step"<<std::endl;
        // vec3print(currentLocation);
        vec3add(currentLocation,step,currentLocation);
        // vec3print(currentLocation);

        if(!this->pointInterpolator(currentLocation,vecNextLocation))
        {
            // vec3print(currentLocation);
            // std::cout<< "k"<<std::endl;
            pointValid = false;
            break;
        }
        if(vec3mag(vel)<this->stagnationThreshold)
        {
            //std::cout<< "Stagnation break "<<this->stagnationThreshold<<std::endl;
            //vec3print(step);
            break;
        }
        // vec3print(step);

        //std::cout<<  "next Point inserted "<< outputPoints->GetNumberOfPoints()<<endl;
        /// Add each point of the integration to the streamline if flag is set
        if(this->calcStreamlines) {
            double point2[3] = {currentLocation[0] + this->origin[0] ,
                                currentLocation[1] + this->origin[1],
                                currentLocation[2] + this->origin[2]};
            vtkIdType id = outputPoints->InsertNextPoint(point2);
        }


        //std::cout<< " currentIntTime "<< currentIntTime<< " this->integrationTime "<< this->integrationTime << std::endl;
    }
    ///If flag is not set only set end location of integration
    if(pointValid) {
        double point2[3] = {currentLocation[0] + this->origin[0] ,
                            currentLocation[1] + this->origin[1],
                            currentLocation[2] + this->origin[2]};
        id = outputPoints->InsertNextPoint(point2);
    }


    //vec3print(currentLocation);

    //std::cout<<"new line Numer of points "<< outputPoints->GetNumberOfPoints()<< std::endl;
    return outputPoints;
}
*/

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
    double currentIntTime = 0;


    if (this->iter == this->initialNum)
    {
        currentIntTime = this->startTime - static_cast<int>(this->startTime);
    }


    ////changed 3 to 2 and removed currentLocation[2] + this->origin[2]
    double point[2] = {currentLocation[0] + this->origin[0],
                        currentLocation[1] + this->origin[1]};

    vtkIdType id = outputPoints->InsertNextPoint(point);
    id = outputPoints->InsertNextPoint(point);

    bool pointValid = true;

    if ( this->iter == this->finalNum - 1)
    {
        while (currentIntTime <= this->endTime)
        {

            pointValid = true;

            vec2copy(currentLocation, currentLocationSave);

            if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {
                pointValid = false;
                break;
            }

            //// changed 3 to 2

            vec2 vel;
            vec2 nextvel;
            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, currentIntTime, temp);
            vec2copy(temp, vel);

            vec2scal(vel, stepSize, traverseDist);

            vec2copy(traverseDist, temp);
            vec2copy(currentLocationSave, currentLocation);
            vec2add(currentLocation, temp, currentLocation);

            currentIntTime += stepSize;

            if ( this->endTime - currentIntTime < stepSize )     ///points at irregular time step
            {
                double interpolator = this->endTime - currentIntTime;    //// new step size
                double offset = this->endTime - static_cast<int>(this->endTime);   //// offset from last timestep for velocity calculations

                if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {
                    pointValid = false;
                    break;
                }

                vec2copy(velCurrLocation, vel);
                vec2copy(velNextLocation, nextvel);
                vec2lint(vel, nextvel, offset , temp);
                vec2copy(temp, vel);

                vec2scal(vel, interpolator, traverseDist);

                vec2copy(traverseDist, temp);
                vec2copy(currentLocationSave, currentLocation);
                vec2add(currentLocation, temp, currentLocation);

                std::vector<double> currloc = {currentLocation[0], currentLocation[1]};
                this->PathlineEndPoints.push_back(currloc);

                currentIntTime += stepSize;
            }

            if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {
                pointValid = false;
                break;
            }

            if (vec2mag(vel) < this->stagnationThreshold) {
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


    }

    else
    {
        while (currentIntTime <= this->timestep[this->iter + 1])
        {
            pointValid = true;
            ////changed 3 to 2
            vec2copy(currentLocation, currentLocationSave);

            if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {
                pointValid = false;
                break;
            }

            vec2 vel;
            vec2 nextvel;
            vec2copy(velCurrLocation, vel);
            vec2copy(velNextLocation, nextvel);
            vec2lint(vel, nextvel, currentIntTime, temp);
            vec2copy(temp, vel);

            //// very first Starting point, at first iteration
            if (this->iter == this->initialNum && currentIntTime == this->startTime - static_cast<int>(this->startTime) )
            {
                auto nearestValue = static_cast<double>(ceilf(static_cast<float>(currentIntTime) * 10)/10);
                double step = nearestValue - currentIntTime;

                vec2scal(vel, step, traverseDist);
                currentIntTime += step;

                vec2copy(traverseDist, temp);
                vec2copy(currentLocationSave, currentLocation);
                vec2add(currentLocation, temp, currentLocation);

                std::cout<<"current time" << currentIntTime <<std::endl;

            }

            //// rest points
            else {
                vec2scal(vel, stepSize, traverseDist);
                currentIntTime += stepSize;

                vec2copy(traverseDist, temp);
                vec2copy(currentLocationSave, currentLocation);
                vec2add(currentLocation, temp, currentLocation);
            }


            if (currentIntTime > this->timestep[this->iter + 1])
            {
                std::vector<double> currloc = {currentLocation[0], currentLocation[1]};
                this->PathlineEndPoints.push_back(currloc);
            }

            if (!this->pointInterpolator(currentLocation, velCurrLocation, velNextLocation)) {
                pointValid = false;
                break;
            }

            if (vec2mag(vel) < this->stagnationThreshold) {
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

/*
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
    floatInput->SetNumberOfComponents(3);
    floatInput->SetNumberOfTuples(this->input->GetPointData()->GetArray(0)->GetNumberOfTuples());
    for(unsigned int index = 0; index < floatInput->GetNumberOfTuples();++index)
    {
        floatInput->SetTuple(index,this->input->GetPointData()->GetArray(0)->GetTuple(index));
    }
    float* inputArray = (float*)floatInput->GetVoidPointer(0);
    //int sizeFlowMap = sizeof(float) * this->source->GetNumberOfPoints()* 3;
    int sizeFlowMap = sizeof(double) * this->source->GetNumberOfPoints()* 3;
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
    flowMap->SetNumberOfComponents(3);
    //Black Pointer Voodoo Magic
#ifndef CPUEXEC
    auto ftleField = vtkSmartPointer<vtkDoubleArray>::New();
    ftleField->SetVoidArray(ftleArray, this->source->GetNumberOfPoints(),0);
    ftleField->SetName("FTLE");
    this->source->GetPointData()->AddArray(ftleField);
#endif
    flowMap->SetVoidArray(flowArray, this->source->GetNumberOfPoints()*3,0);
    //ftleField->Print(std::cout);
    flowMap->SetName("FlowMap");
    this->source->GetPointData()->AddArray(flowMap);

    //return ftleField;

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