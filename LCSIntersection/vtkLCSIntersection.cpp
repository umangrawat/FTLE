#include <vtkVersion.h>
#include "vtkLCSIntersection.h"
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vtkDataSetAlgorithm.h>
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include <vtkDataSet.h>
#include "vtkPointData.h"

#include <vtkPolyData.h>
#include <vtkImageData.h>

#include "chrono"

#include "linalg.h"

#include <vtkGradientFilter.h>
#include <vtkPoints.h>
#include <vtkDataObject.h>

#include "vtkMath.h"
#include "vtkArrayCalculator.h"
#include "vtkNew.h"
#include "vtkMultiBlockDataSet.h"
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPolyDataConnectivityFilter.h>

#include <vtkXMLImageDataWriter.h>
#include <vtkFieldData.h>
#include <stdio.h>
#include <vtkStructuredGrid.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkArrayDataWriter.h>
#include <vtkAssignAttribute.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

#include <vtkPolyDataAlgorithm.h>

#include <Eigenvalues>
#include <Core>
#include <Dense>

#include "vtkPCAStatistics.h"
#include "vtkTable.h"

#include <random>
#include <vector>
#include <math.h>
#include <limits>
#include <iomanip>

#include <vtkPassArrays.h>
#include <vtkAppendPolyData.h>
#include <vtkXMLPolyDataReader.h>
//#include "vtkIntersectionPolyDataFilter.h"


using namespace std;
using namespace std::chrono;

const double LIMIT_DOUBLE = 1000* std::numeric_limits<double>::epsilon();

int getIndex(int x, int y, int* resolution)
{
    return ( y  * resolution[1]   + x );
}

vtkStandardNewMacro(vtkLCSIntersection);

//-----------------------------------------------------------------------------
vtkLCSIntersection::vtkLCSIntersection()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);                                                        ////changing to 2
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
vtkLCSIntersection::~vtkLCSIntersection()
{

}


//----------------------------------------------------------------------------
int vtkLCSIntersection::FillInputPortInformation( int port, vtkInformation* info )
{

    if (port == 0) {
        //info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}


//----------------------------------------------------------------------------
int vtkLCSIntersection::FillOutputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 )
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }

    return 0;
}

//----------------------------------------------------------------------------

int vtkLCSIntersection::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);


    int resolution[6];
    resolution[0] = 0;
    resolution[1] = static_cast<int>(this->cellsNumber[0]);
    resolution[2] = 0;
    resolution[3] = static_cast<int>(this->cellsNumber[1]);
    resolution[4] = 0;
    resolution[5] = static_cast<int>(this->cellsNumber[2]);

    return 1;

}

//----------------------------------------------------------------------------


int vtkLCSIntersection::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    //// Get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkSmartPointer<vtkPolyData> input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));


    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    double* bounds = input->GetBounds();
    //std::cout<<bounds[0] << " " <<bounds[1]<< " "<<bounds[2] << " " <<bounds[3]<<" "<<bounds[4] << " " <<bounds[5];
    int number = input->GetNumberOfCells();

    int resolution[6];
    resolution[0] = 0;
    resolution[1] = static_cast<int>(this->cellsNumber[0]);
    resolution[2] = 0;
    resolution[3] = static_cast<int>(this->cellsNumber[1]);
    resolution[4] = 0;
    resolution[5] = static_cast<int>(this->cellsNumber[2]);


    double spacing[3] ={0.};
    spacing[0] = bounds[1] / (static_cast<double>(resolution[1]));
    spacing[1] = bounds[3] / (static_cast<double>(resolution[3]));
    spacing[2] = bounds[5] / (static_cast<double>(resolution[5]));

    //std::cout<<spacing[0] << " " <<spacing[1]<< " "<<spacing[2] <<std::endl;
    //std::cout<<resolution[1] << " " <<resolution[3]<< " "<<resolution[5] <<std::endl;

    std::cout<< number <<std::endl;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points = input->GetPoints();

    vtkSmartPointer<vtkDoubleArray> LocationIndexValues = vtkSmartPointer<vtkDoubleArray>::New();
    LocationIndexValues->SetNumberOfComponents(6);

    int numpoints = points->GetNumberOfPoints();
    //std::cout<<"points" << numpoints <<std::endl;

    int n = numpoints / 2;

    int arr0[n] = {0};
    int arr1[n] = {0};

    std::vector<pair<int, int>> ids;

    for (int i = 0; i < points->GetNumberOfPoints(); i+=2)
    //for (int i = 0; i < points->GetNumberOfPoints(); i++)
    {
        double point[3] = {0.};
        vtkIdType id = vtkIdType(i);
        points->GetPoint(id, point);

        double point1[3] = {0.};
        vtkIdType id1 = vtkIdType(i+1);
        points->GetPoint(id1, point1);

        //if ( i < 100)
        //std::cout<< id << " "<< point[0] <<" "<<point[1] <<std::endl;
        //std::cout<< id1 << " "<< point1[0] <<" "<<point1[1] <<std::endl;


        vtkIdType idloc = vtkIdType(i / 2);

        int index = cellLocator(point,bounds, spacing, resolution);
        LocationIndexValues->InsertTuple6(idloc, index, point[0], point[1], point1[0], point1[1], 0.);
        arr0[i/2] = index;
        arr1[i/2] = idloc;
    }


    for (int m = 0; m < n; m++) {

        ids.push_back( make_pair (arr0[m], arr1[m]));
    }

    std::sort(ids.begin(), ids.end());

    int numTuples = LocationIndexValues->GetNumberOfTuples();
    std::cout<< numTuples<<std::endl;

    std::cout<<"vec size" <<ids.size()<<std::endl;
    int increment = 1;

    int tracker = 0;

    //for (int j = 0; j < ids.size(); j+= increment)
    while (ids.size() > 0)
    {
        increment = 1;
        vtkIdType idx = vtkIdType(ids[0].second);
        double point1[6] ={0.};
        LocationIndexValues->GetTuple(idx, point1);

        ///std::cout <<"Sorting " << ids[0].first << " " << ids[0].second <<" "<<  tuplevals[1] << " " <<tuplevals[2] << " " <<tuplevals[3] << " " <<tuplevals[4] << std::endl;

                                                                                     ///A0                      B0                  A1                      B1
        double A0 = point1[1];
        double B0 = point1[2];
        double A1 = point1[3];
        double B1 = point1[4];

        for ( int k = 1; k < ids.size(); k++)
        {
            if (ids[0].first == ids[k].first)
            {
                vtkIdType idx1 = vtkIdType(ids[k].second);
                double point2[6] ={0.};
                LocationIndexValues->GetTuple(idx1, point2);

                ///std::cout <<"Sorting " << ids[0].first << " " << ids[0].second <<" "<<  tuplevals1[1] << " " <<tuplevals1[2] << " " <<tuplevals1[3] << " " <<tuplevals1[4] << std::endl;

                                                                                    ///A2                      B2                  A3                      B3
                double A2 = point2[1];
                double B2 = point2[2];
                double A3 = point2[3];
                double B3 = point2[4];

                ///check if the two points intersect

                if (((A2-A0)*(B1-B0) - (B2-B0)*(A1-A0)) * ((A3-A0)*(B1-B0) - (B3-B0)*(A1-A0)) < 0
                    &&
                    ((A0-A2)*(B3-B2) - (B0-B2)*(A3-A2)) * ((A1-A2)*(B3-B2) - (B1-B2)*(A3-A2)) < 0)
                {
                    tracker+=1;
                }

                increment+=1;
            }
        }

        ids.erase(ids.begin());

    }

    //std::cout<< "repeated values " << repeatedIndices.size()<<std::endl;
    std::cout<<"intersection number "<< tracker <<std::endl;


    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines = input->GetLines();

    //int numlines = lines->GetNumberOfCells();
    //std::cout<<"numlines "<<numlines <<std::endl;

    //// Getting the scalar values in an array
    //vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
    //field = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("FTLE");
    return 1;

}


void vtkLCSIntersection::PrintSelf(ostream &os, vtkIndent indent)
{
}

//// Dot Product

double DotProd(double x[2], double y[2])
{
    return x[0]*y[0] + x[1]*y[1];
}

double distance(std::vector<double> x, std::vector<double> y)
{
    double xdirec = (y[0] - x[0]) * (y[0] - x[0]);
    double ydirec = (y[1] - x[1]) * (y[1] - x[1]);
    double zdirec = (y[2] - x[2]) * (y[2] - x[2]);

    return std::sqrt( xdirec + ydirec + zdirec );
}

int vtkLCSIntersection:: cellLocator( double Location[3], double bounds[6], double spacing[3], int* resolution )
{

    int xIndex, yIndex = 0;


    ///points which are on or more than x-boundary
    if ( Location[0] >= bounds[1] && Location[1] < bounds[3] )
    {
        xIndex = static_cast<int>( Location[0] / spacing[0] ) - 1;
        yIndex = static_cast<int>( Location[1] / spacing[1] );
    }

        ///points which are on or more than y-boundary
    else if ( Location[0] < bounds[1] && Location[1] >= bounds[3] )
    {
        xIndex = static_cast<int>( Location[0] / spacing[0] );
        yIndex = static_cast<int>( Location[1] / spacing[1] ) - 1;
    }

        ////both x and y boundaries
    else if ( Location[0] >= bounds[1] && Location[1] >= bounds[3] )
    {
        xIndex = static_cast<int>( Location[0] / spacing[0] ) - 1;
        yIndex = static_cast<int>( Location[1] / spacing[1] ) - 1;
    }

        ////rest of the points
    else {
        xIndex = static_cast<int>( Location[0] / spacing[0] );
        yIndex = static_cast<int>( Location[1] / spacing[1] );
    }

    int index = getIndex(xIndex, yIndex, resolution );

    return index;
}