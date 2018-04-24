#include <vtkVersion.h>
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
#include "vtkRidgeLine.h"

#include <vtkGradientFilter.h>
#include <vtkPoints.h>
#include <vtkDataObject.h>

#include "vtkMath.h"
#include "vtkArrayCalculator.h"
#include "vtkNew.h"

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

#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>

#include <vtkArrayDataWriter.h>
#include <vtkAssignAttribute.h>

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


using namespace std;
using namespace std::chrono;

const double LIMIT_DOUBLE = 1000* std::numeric_limits<double>::epsilon();

int getIndex(int z,int y, int x, int* dataDims);
double centralDiff(double fx, double fz, double dist);
double forwardDiff(double fy, double fz, double dist);
double backwardDiff(double fx, double fy, double dist);
double DotProd(double x[2], double y[2]);
double distance(std::vector<double> x, std::vector<double> y);
//void Interpolation (vtkImageData* oldinput, vtkDoubleArray* oldfield, int* olddims, int* dims,  int ResolutionFactor, vtkStructuredGrid* input);
void Gradient (vtkImageData* input, vtkDoubleArray* field, int* dims, vtkDoubleArray* gradient);
void Hessian ( vtkImageData* input, vtkDoubleArray* gradient, int* dims, vtkDoubleArray* hessian );
void MinEigenVector (vtkDoubleArray* hessian, int* dims, vtkDoubleArray* MinEigenvecs);
void PCA( vtkDoubleArray* MinEigenVecs ,int* dims, vtkDoubleArray* pcaMaxEigenVec );
void FlipMinEigenVecs(vtkDoubleArray* pcaMaxEigenVec, vtkDoubleArray* MinEigenvecs, int* dims);
void Ridges( vtkImageData* input, vtkDoubleArray* gradient, vtkDoubleArray* hessian, vtkDoubleArray* MinEigenvecs, int* dims, vtkPoints* RidgePoints );
bool RidgeEigenValue (double IntHessian[4]);

double limit = 0.;


vtkStandardNewMacro(vtkRidgeLine);

//-----------------------------------------------------------------------------
vtkRidgeLine::vtkRidgeLine()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);                                                        ////changing to 2
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
vtkRidgeLine::~vtkRidgeLine()
{

}


//----------------------------------------------------------------------------
int vtkRidgeLine::FillInputPortInformation( int port, vtkInformation* info )
{

    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
        return 1;
    }
    return 0;
}


//----------------------------------------------------------------------------
int vtkRidgeLine::FillOutputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 )
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }

    return 0;
}
//----------------------------------------------------------------------------


/*
int vtkRidgeLine::RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector)
{
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);



    int extentSeedGrid[6];
    extentSeedGrid[0] = 0;
    extentSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0])-1;
    extentSeedGrid[2] = 0;
    extentSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1])-1;
    extentSeedGrid[4] = 0;
    extentSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2])-1;

    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extentSeedGrid, 6);

    return 1;
}
*/
//----------------------------------------------------------------------------

int vtkRidgeLine::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    return 1;

}

//----------------------------------------------------------------------------


int vtkRidgeLine::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    //// Get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


    limit = this->LimitEigenValue;

    //// Input dimensions
    int* dims = input->GetDimensions();

    //// Getting the scalar values in an array
    vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
    //field = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetArray(0));
    field = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("FTLE"));

    std::cout<<"Starts" <<std::endl;
    //// Gradient
    vtkSmartPointer<vtkDoubleArray> gradient = vtkSmartPointer<vtkDoubleArray>::New();
    gradient->SetNumberOfComponents(2);
    Gradient(input, field, dims, gradient);


    //// Hessian
    vtkSmartPointer<vtkDoubleArray> hessian = vtkSmartPointer<vtkDoubleArray>::New();
    hessian->SetNumberOfComponents(4);
    Hessian(input, gradient, dims, hessian);


    //// Eigenvalues and Eigenvectors and get the minimum eigenvalue and corresponding eigenvector
    vtkSmartPointer<vtkDoubleArray> MinEigenvecs = vtkSmartPointer<vtkDoubleArray>::New();
    MinEigenvecs->SetNumberOfComponents(2);
    MinEigenVector(hessian, dims, MinEigenvecs);


    //// PCA and get the maximum eigenvector corresponding to max. eigenvalue
    vtkSmartPointer<vtkDoubleArray> pcaMaxEigenVec = vtkSmartPointer<vtkDoubleArray>::New();
    pcaMaxEigenVec->SetNumberOfComponents(2);
    PCA(MinEigenvecs, dims, pcaMaxEigenVec);



    //// Dot product between PCA max eigenvector and minimum eigenvector of each node and flipping if opposite
    FlipMinEigenVecs(pcaMaxEigenVec, MinEigenvecs, dims);


    //// Ridge Points using lookup table
    vtkSmartPointer<vtkPoints> RidgePoints = vtkSmartPointer <vtkPoints>::New();
    Ridges(input, gradient, hessian, MinEigenvecs, dims, RidgePoints);

    vtkIdType numRidgePoints = RidgePoints->GetNumberOfPoints();

    std::cout<< "Points calculated"<<std::endl;

    //// Create a line

    vtkSmartPointer<vtkCellArray> outputLines = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();

    for ( int k = 0; k < numRidgePoints - 1; k+= 2)
    {
        outputLines->InsertNextCell(2);

        for (int j = 0; j < 2; j++) {

            double point[3] = {0.};

            vtkIdType id = vtkIdType(k + j);
            RidgePoints->GetPoint(id, point);

            //outputPoints->InsertNextPoint(point);
            vtkIdType nextPoint = outputPoints->InsertNextPoint(point);
            outputLines->InsertCellPoint(nextPoint);
        }
    }

    std::cout<< "Creating lines"<<std::endl;

    vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();

    linesPolyData->SetPoints(outputPoints);
    linesPolyData->SetLines(outputLines);

    //// Copying the lines output array
    output->ShallowCopy(linesPolyData);



    /////////////////// Get the Length of ridge lines ///////////////////////////////

    //// storing all the cellpoints of ridges
    std::vector<vector<double>> cellpoints;
    std::vector<double> ridgeLinelength;

    ////storing all the ridge points in a vector
    for( int pointIndex = 0; pointIndex < numRidgePoints; pointIndex++ )
    {
        double cellpoint[3] = {0.};
        vtkIdType id = vtkIdType(pointIndex);
        RidgePoints->GetPoint(id, cellpoint);
        std::vector<double> point = {cellpoint[0], cellpoint[1], cellpoint[2]};

        cellpoints.push_back(point);
    }

    while (cellpoints.size() > 0)
    {
        int iter = 0;
        double ridgeLength = 0.;
        vtkSmartPointer<vtkPoints> linepts = vtkSmartPointer<vtkPoints>::New();

        std::vector<double> firstpoint = cellpoints.at( 0 );
        std::vector<double> nextpoint = cellpoints.at( 1 );

        ridgeLength += distance(firstpoint, nextpoint);                                                            /// length between first two points

        std::vector<int> storeid;

        for (int j = 2; j < cellpoints.size() - 1; j+=2)
        {
            std::vector<double> nextcellpoint = cellpoints.at(j);

            if ((std::fabs(nextpoint[0]) - std::fabs(nextcellpoint[0])) < LIMIT_DOUBLE
                && (std::fabs(nextpoint[1]) - std::fabs(nextcellpoint[1])) < LIMIT_DOUBLE
                && (std::fabs(nextpoint[2]) - std::fabs(nextcellpoint[2])) < LIMIT_DOUBLE)
            {
                std::vector<double> nextcellSecondpoint = cellpoints.at(j + 1);

                if (iter == 0)
                {
                    linepts->InsertNextPoint(firstpoint[0], firstpoint[1], firstpoint[2]);
                }

                linepts->InsertNextPoint(nextpoint[0], nextpoint[1], nextpoint[2]);
                linepts->InsertNextPoint(nextcellpoint[0], nextcellpoint[1], nextcellpoint[2]);

                ridgeLength += distance(nextcellpoint, nextcellSecondpoint);                                        /// length in the next cell and hence forth

                nextpoint = nextcellSecondpoint;

                iter += 1;
                storeid.push_back( j );
                storeid.push_back( j+1 );

            }

        }

        //// storing the length in the vector
        ridgeLinelength.push_back(ridgeLength);

        int idPoints = storeid.size();

        for (int l = 0; l < storeid.size(); l++)
        {
            int indexvalue = storeid.at(l);
            std::vector<vector<double>>::iterator indexIt = cellpoints.begin() + indexvalue - l;
            cellpoints.erase(indexIt);
        }

        cellpoints.erase( cellpoints.begin(), cellpoints.begin() + 2 );
        storeid.erase(storeid.begin(), storeid.end());

    }

    vtkSmartPointer<vtkDoubleArray> RidgeLineLengthArray = vtkSmartPointer<vtkDoubleArray>::New();
    RidgeLineLengthArray->SetNumberOfTuples(numRidgePoints);
    RidgeLineLengthArray->SetName("RidgeLineLengthArray");

    for ( int null = 0; null < numRidgePoints; null ++)
    {
        vtkIdType idNull = vtkIdType(null);
        RidgeLineLengthArray->InsertTuple1( idNull, 0. );
    }

    double totalLength = 0.;
    for (int len  = 0; len < ridgeLinelength.size(); len++)
    {
        vtkIdType idLen = vtkIdType(len);
        double length = 0.;
        length = ridgeLinelength.at(len);
        totalLength += length;

        RidgeLineLengthArray->SetTuple1(idLen, length);

        std::cout <<"Length of ridge " << len << " " << ridgeLinelength.at(len) <<std::endl;
    }

    std::cout <<"Total length " << totalLength << std::endl;
    std::cout<< "numpoints" << numRidgePoints << std::endl;

    //output->GetPointData()->AddArray(RidgeLineLengthArray);

    return 1;
}


void vtkRidgeLine::PrintSelf(ostream &os, vtkIndent indent)
{
}


/// Get the index of a 3 dim rectilinear grid

int getIndex(int z,int y, int x, int* dataDims)
{
    return (z*dataDims[1]*dataDims[0]+ y*dataDims[0]+x);
}

//// Central Differences

double centralDiff(double fx, double fz, double dist)
{
    return (fz-fx)/dist;
}


//// Forward Differences

double forwardDiff(double fy, double fz, double dist)
{
    return (fz-fy)/dist;
}


//// Backward Differnces

double backwardDiff(double fx, double fy, double dist)
{
    return (fx-fy)/dist;

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


///Gradient Calculation

void Gradient(vtkImageData* input, vtkDoubleArray* field, int* dims, vtkDoubleArray* gradient)
{
    double spaceX, spaceY = 0.;

    double du[2] = {0.};

    for ( int k = 0; k < 1; k++)
    {
        for (int j = 0; j < dims[1]; j++)                                              ////
        {
            for (int i = 0; i < dims[0]; i++)                                          ////
            {

                int index_gradient = getIndex(k, j, i, dims);
                vtkIdType id1 = vtkIdType(index_gradient);
                vtkIdType id2 = vtkIdType(index_gradient);
                double t1[1], t2[1] = {0.};
                double pt1[3], pt2[3] = {0.};


                //// y component
                if (j == 0)
                {

                    id1 = getIndex (k, j, i, dims);
                    field->GetTuple( id1, t1 );
                    input->GetPoint( id1, pt1);

                    id2 = getIndex ( k, j+1, i, dims);
                    field->GetTuple( id2, t2 );
                    input->GetPoint( id2, pt2);

                    spaceY = pt2[1] - pt1[1];

                    du[1] = forwardDiff(t1[0],t2[0], spaceY);

                }

                else if ( j == dims[1] - 1 )                                                ////
                {

                    id1 = getIndex(k, j, i, dims);
                    field->GetTuple(id1,t1);
                    input->GetPoint( id1, pt1);

                    id2 = getIndex(k, j-1, i, dims);
                    field->GetTuple(id2,t2);
                    input->GetPoint( id2, pt2);

                    spaceY = (pt2[1] - pt1[1]) * ( - 1 );

                    du[1] = backwardDiff( t1[0],t2[0], spaceY );

                }

                else
                {

                    id1 = getIndex(k, j-1, i, dims);
                    field->GetTuple(id1,t1);
                    input->GetPoint( id1, pt1);

                    id2 = getIndex(k, j+1, i, dims);
                    field->GetTuple(id2,t2);
                    input->GetPoint( id2, pt2);

                    spaceY = pt2[1] - pt1[1];

                    id1 = getIndex (k, j, i, dims);

                    du[1] = centralDiff(t1[0],t2[0], spaceY);

                }



                //// x Component

                if( i == 0 )
                {

                    id1 = getIndex(k,j,i,dims);
                    field->GetTuple(id1,t1);
                    input->GetPoint(id1,pt1);

                    id2 = getIndex(k,j,i+1,dims);
                    field->GetTuple(id2,t2);
                    input->GetPoint(id2,pt2);

                    spaceX = pt2[0] - pt1[0];

                    du[0] = forwardDiff(t1[0],t2[0], spaceX);
                }

                else if ( i == dims[0] - 1 )
                {

                    id1 = getIndex(k,j,i,dims);
                    field->GetTuple(id1,t1);
                    input->GetPoint(id1,pt1);

                    id2 = getIndex(k,j,i-1,dims);
                    field->GetTuple(id2,t2);
                    input->GetPoint(id2,pt2);

                    spaceX = (pt2[0] - pt1[0]) * ( -1 );

                    du[0] = backwardDiff(t1[0],t2[0], spaceX);

                }

                else
                {

                    id1 = getIndex(k,j,i-1,dims);
                    field->GetTuple(id1,t1);
                    input->GetPoint(id1,pt1);

                    id2 = getIndex(k,j,i+1,dims);
                    field->GetTuple(id2,t2);
                    input->GetPoint(id2,pt2);

                    spaceX = pt2[0] - pt1[0];

                    id1 = getIndex (k, j, i, dims);

                    du[0] = centralDiff(t1[0],t2[0], spaceX);

                }

                gradient->InsertTuple2 (id1, du[0], du[1] );


            }

        }

    }

    std::cout<<"gradient calculation finished" <<std::endl;

    return;

}



//// Hessian Calculation

void Hessian ( vtkImageData* input, vtkDoubleArray* gradient, int* dims, vtkDoubleArray* hessian )
{

    double dxdx, dxdy, dydx, dydy = 0.;

    for ( int k = 0; k < 1; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {

                int index_hessian = getIndex(k, j, i, dims);
                vtkIdType id1_hess = vtkIdType(index_hessian);
                vtkIdType id2_hess = vtkIdType(index_hessian);
                double s1[2], s2[2];
                double point1[3], point2[3] = {0.};

                double spaceX, spaceY = 0.;

                ////y component near boundary needs forward diff
                if (j == 0)
                {
                    id1_hess = getIndex(k, j, i, dims);
                    gradient->GetTuple(id1_hess, s1);
                    input->GetPoint(id1_hess,point1);

                    id2_hess = getIndex(k, j+1, i, dims);
                    gradient->GetTuple(id2_hess, s2);
                    input->GetPoint(id2_hess,point2);

                    spaceY = point2[1] - point1[1];

                    dydx = forwardDiff(s1[0], s2[0], spaceY);
                    dydy = forwardDiff(s1[1], s2[1], spaceY);

                }

                    ////y component near boundary needs backward diff
                else if (j == dims[1] - 1 )
                {
                    id1_hess = getIndex(k, j, i, dims);
                    gradient->GetTuple(id1_hess, s1);
                    input->GetPoint(id1_hess,point1);

                    id2_hess = getIndex(k, j-1, i, dims);
                    gradient->GetTuple(id2_hess, s2);
                    input->GetPoint(id2_hess,point2);

                    spaceY = (point2[1] - point1[1]) * (-1);

                    dydx = backwardDiff(s1[0], s2[0], spaceY);
                    dydy = backwardDiff(s1[1], s2[1], spaceY);
                }

                    ////rest central differences

                else
                {
                    id1_hess = getIndex(k, j-1, i, dims);
                    gradient->GetTuple(id1_hess, s1);
                    input->GetPoint(id1_hess,point1);

                    id2_hess = getIndex(k, j+1, i, dims);
                    gradient->GetTuple(id2_hess, s2);
                    input->GetPoint(id2_hess,point2);

                    id1_hess = getIndex (k, j, i, dims);

                    spaceY = point2[1] - point1[1];

                    dydx = centralDiff(s1[0], s2[0], spaceY);
                    dydy = centralDiff(s1[1], s2[1], spaceY);
                }



                if (i == 0)
                {
                    id1_hess = getIndex(k, j, i, dims);
                    gradient->GetTuple(id1_hess, s1);
                    input->GetPoint(id1_hess,point1);

                    id2_hess = getIndex(k, j, i+1, dims);
                    gradient->GetTuple(id2_hess, s2);
                    input->GetPoint(id2_hess,point2);

                    spaceX = point2[0] - point1[0];

                    dxdx = forwardDiff(s1[0], s2[0], spaceX);
                    dxdy = forwardDiff(s1[1], s2[1], spaceX);

                }

                else if (i == dims[0] - 1)
                {
                    id1_hess = getIndex(k, j, i, dims);
                    gradient->GetTuple(id1_hess, s1);
                    input->GetPoint(id1_hess,point1);

                    id2_hess = getIndex(k, j, i-1, dims);
                    gradient->GetTuple(id2_hess, s2);
                    input->GetPoint(id2_hess,point2);

                    spaceX = (point2[0] - point1[0]) * ( -1 );

                    dxdx = backwardDiff(s1[0], s2[0], spaceX);
                    dxdy = backwardDiff(s1[1], s2[1], spaceX);

                }

                else
                {
                    id1_hess = getIndex(k, j, i-1, dims);
                    gradient->GetTuple(id1_hess, s1);
                    input->GetPoint(id1_hess,point1);

                    id2_hess = getIndex(k, j, i+1, dims);
                    gradient->GetTuple(id2_hess, s2);
                    input->GetPoint(id2_hess,point2);

                    id1_hess = getIndex (k, j, i, dims);

                    spaceX = point2[0] - point1[0];

                    dxdx = centralDiff(s1[0], s2[0], spaceX);
                    dxdy = centralDiff(s1[1], s2[1], spaceX);

                }

                ////storing the hessian components
                hessian->InsertTuple4(id1_hess,dxdx, dxdy, dydx, dydy);
            }

        }

    }

    std::cout<<"Hessian calculation finished" <<std::endl;
    return;

}



//// Calculate eigenvalues and eigenvectors and get the minimum eigenvalue and corresponding eigenvector

void MinEigenVector (vtkDoubleArray* hessian, int* dims, vtkDoubleArray* MinEigenvecs)
{

    for ( int k = 0; k < 1; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int index_eigen = getIndex(k, j, i, dims);
                vtkIdType id_eigen = vtkIdType(index_eigen);
                double hessian_components[4];

                ////get hessian tuple
                hessian->GetTuple(id_eigen, hessian_components);


                ////storing hessian elements in a matrix
                Eigen::Matrix2d hessianmatrix(2,2);
                hessianmatrix(0, 0) = hessian_components[0];
                hessianmatrix(0, 1) = hessian_components[1];
                hessianmatrix(1, 0) = hessian_components[2];
                hessianmatrix(1, 1) = hessian_components[3];


                ////calculating eigenvectors and eigenvalues using Eigen
                Eigen::EigenSolver<Eigen::MatrixXd> es;
                es.compute(hessianmatrix, true);


                ////Get minimum eigenvalue and the corresponding eigenvector
                Eigen::MatrixXd ValueReal(2,1);
                Eigen::MatrixXd VectorReal(2,2);

                ValueReal = es.eigenvalues().real();     ///to get the real part of Evalues and Evectors
                VectorReal = es.eigenvectors().real();


                double lambda[2];
                double RealEV[2];

                lambda[0] = ValueReal(0,0);
                lambda[1] = ValueReal(1,0);


                //std::cout<< lambda[0] << ", " << lambda[1] << std::endl;

                if (lambda[0] < lambda[1] )
                {
                    RealEV[0] = VectorReal(0,0);
                    RealEV[1] = VectorReal(1,0);
                    MinEigenvecs->InsertTuple2(id_eigen, RealEV[0], RealEV[1]);
                }

                else if (lambda[1] < lambda[0])
                {
                    RealEV[0] = VectorReal(0,1);
                    RealEV[1] = VectorReal(1,1);
                    MinEigenvecs->InsertTuple2(id_eigen, RealEV[0], RealEV[1]);
                }

            }
        }
    }

    std::cout<<"Eigenvalue calculation finished" <<std::endl;

    return;
}

////changed from here


//// Principal Component Analysis, PCA

void PCA ( vtkDoubleArray* MinEigenvecs , int* dims, vtkDoubleArray* pcaMaxEigenVec )
{

    vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
    vtkSmartPointer<vtkDoubleArray> pcaEigenvecs = vtkSmartPointer<vtkDoubleArray>::New();


    for ( int k = 0; k < 1; k++)  ////doing the iteration till x-1 and y-1 to avoid boundary problems
    {
        for (int j = 0; j < dims[1] - 1 ; j++)
        {
            std::cout<<"PCA Progress "<< int(double (j+1)/double (dims[1])* 100 ) <<"%\r";
            for (int i = 0; i < dims[0] - 1; i++) {
                int nodeCnt = 4;
                double evec1[2], evec2[2], evec3[2], evec4[2] = {0.};
                //double XComp[8], YComp[8] = {0.};
                double X[2][2 * nodeCnt] = {0.};

                int pcaIndex = getIndex(k, j, i, dims);
                vtkIdType pcaId1 = vtkIdType(pcaIndex);
                vtkIdType pcaId2 = vtkIdType(pcaIndex);
                vtkIdType pcaId3 = vtkIdType(pcaIndex);
                vtkIdType pcaId4 = vtkIdType(pcaIndex);


                ////Get all the eigenvectors of a cell
                pcaId1 = getIndex(k, j, i, dims);
                MinEigenvecs->GetTuple(pcaId1, evec1);

                pcaId2 = getIndex(k, j, i + 1, dims);
                MinEigenvecs->GetTuple(pcaId2, evec2);

                pcaId3 = getIndex(k, j + 1, i, dims);
                MinEigenvecs->GetTuple(pcaId3, evec3);

                pcaId4 = getIndex(k, j + 1, i + 1, dims);
                MinEigenvecs->GetTuple(pcaId4, evec4);



                ////store X & Y component differently and negative of each

                X[0][0] = evec1[0];
                X[1][0] = evec1[1];
                X[0][1] = -1 * evec1[0];
                X[1][1] = -1 * evec1[1];
                X[0][2] = evec2[0];
                X[1][2] = evec2[1];
                X[0][3] = -1 * evec2[0];
                X[1][3] = -1 * evec2[1];
                X[0][4] = evec3[0];
                X[1][4] = evec3[1];
                X[0][5] = -1 * evec3[0];
                X[1][5] = -1 * evec3[1];
                X[0][6] = evec4[0];
                X[1][6] = evec4[1];
                X[0][7] = -1 * evec4[0];
                X[1][7] = -1 * evec4[1];

                mat2 C;
                {
                    for (int j = 0; j < 2; j++) {
                        for (int i = 0; i < 2; i++) {

                            double sum = 0.0;
                            for (int k = 0; k < 2 * nodeCnt; k++) {
                                sum += X[i][k] * X[j][k];
                            }

                            C[i][j] = sum / (2 * nodeCnt - 1);
                        }
                    }
                }

                // compute eigenvalues and eigenvectors
                vec2 eigenvalues;
                vec2 eigenvectors[2];
                {
                    // force C to be symmetric (added 2007-08-15, untested)
                    mat2symm(C, C);

                    // eigenvalues
                    bool allReal = (mat2eigenvalues(C, eigenvalues) == 2);

                    if (!allReal) {
                        //printf("got complex eigenvalues: %g, %g, %g, returning zero\n", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
                        //mat3dump(C, stdout);

                        return;
                    }
                    // eigenvectors
                    mat2realEigenvector(C, eigenvalues[0], eigenvectors[0]);
                    mat2realEigenvector(C, eigenvalues[1], eigenvectors[1]);
                }

#if 0
                // get largest eigenvalue
                int maxEVIdx;
                {
                        if (eigenvalues[0] > eigenvalues[1])
                        {
                                maxEVIdx = 0;
                         }
                        else {
                            maxEVIdx = 1;
                        }

                }
#else
                // sort eigenvalues in descending order
                int evalDescIndices[2];
                {
                    if (eigenvalues[0] > eigenvalues[1]) {

                        evalDescIndices[0] = 0;
                    } else {
                        evalDescIndices[0] = 1;
                    }


                    int remainingIndices[2];
                    switch (evalDescIndices[0]) {
                        case 0:
                            remainingIndices[0] = 1;
                            remainingIndices[1] = 2;
                            break;
                        case 1:
                            remainingIndices[0] = 0;
                            remainingIndices[1] = 2;
                            break;
                        case 2:
                            remainingIndices[0] = 0;
                            remainingIndices[1] = 1;
                            break;
                    }

                    if (eigenvalues[remainingIndices[0]] > eigenvalues[remainingIndices[1]]) {
                        evalDescIndices[1] = remainingIndices[0];
                        //evalDescIndices[2] = remainingIndices[1];
                    } else {
                        evalDescIndices[1] = remainingIndices[1];
                        //evalDescIndices[2] = remainingIndices[0];
                    }

                    /*
                    if (eigenValuesDesc) {
                        eigenValuesDesc[0] = eigenvalues[evalDescIndices[0]];
                        eigenValuesDesc[1] = eigenvalues[evalDescIndices[1]];
                        //eigenValuesDesc[2] = eigenvalues[evalDescIndices[2]];
                    }*/
                }
#endif

            vec2 evMax;
            {
                vec2copy(eigenvectors[evalDescIndices[0]], evMax);
            }

            pcaMaxEigenVec->InsertTuple2( pcaId1, evMax[0], evMax[1] );

                /*
                    ////store X & Y component differently and negative of each

                    XComp[0] = evec1[0];
                    XComp[1] = -1 * evec1[0];
                    XComp[2] = evec2[0];
                    XComp[3] = -1 * evec2[0];
                    XComp[4] = evec3[0];
                    XComp[5] = -1 * evec3[0];
                    XComp[6] = evec4[0];
                    XComp[7] = -1 * evec4[0];

                    YComp[0] = evec1[1];
                    YComp[1] = -1 * evec1[1];
                    YComp[2] = evec2[1];
                    YComp[3] = -1 * evec2[1];
                    YComp[4] = evec3[1];
                    YComp[5] = -1 * evec3[1];
                    YComp[6] = evec4[1];
                    YComp[7] = -1 * evec4[1];

                    vtkSmartPointer<vtkDoubleArray> xArray = vtkSmartPointer<vtkDoubleArray>::New();
                    xArray->SetNumberOfComponents(1);
                    xArray->SetName("x");

                    vtkSmartPointer<vtkDoubleArray> yArray = vtkSmartPointer<vtkDoubleArray>::New();
                    yArray->SetNumberOfComponents(1);
                    yArray->SetName("y");


                    for ( int l = 0; l < 8; l++)
                    {

                        vtkIdType idx_e;
                        idx_e = l;

                        xArray->InsertTuple1(idx_e, XComp[l]);
                        yArray->InsertTuple1(idx_e, YComp[l]);
                    }


                    vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();

                    datasetTable->AddColumn(xArray);
                    datasetTable->AddColumn(yArray);


                    #if VTK_MAJOR_VERSION <= 5
                    pcaStatistics->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );
                    #else
                    pcaStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );
                    #endif

                    pcaStatistics->SetColumnStatus("x", 1 );
                    pcaStatistics->SetColumnStatus("y", 1 );

                    pcaStatistics->RequestSelectedColumns();
                    pcaStatistics->SetDeriveOption(true);
                    pcaStatistics->Update();



                    ///////// Eigenvalues ////////////
                    vtkSmartPointer<vtkDoubleArray> eval = vtkSmartPointer<vtkDoubleArray>::New();
                    pcaStatistics->GetEigenvalues(eval);

                    ///////// Eigenvectors ////////////

                    pcaStatistics->GetEigenvectors(pcaEigenvecs);


                    for(vtkIdType m = 0; m < pcaEigenvecs->GetNumberOfTuples(); m++)
                    {
                        //double* evector = new double[pcaEigenvecs->GetNumberOfComponents()];
                        double evector[2] = {0.};
                        pcaEigenvecs->GetTuple( m, evector );

                        ////store the maximum eigenvector for each cell, which is the first one
                        if (m == 0)
                        {
                            pcaMaxEigenVec->InsertTuple2( pcaId1, evector[0], evector[1] );
                        }

                    }
                     */
            }
        }
    }

    std::cout<<"PCA calculation finished" <<std::endl;
    return;
}


//// Dot product between PCA max eigenvector and minimum eigenvector of each node and flipping opposite if <0

void FlipMinEigenVecs(vtkDoubleArray* pcaMaxEigenVec, vtkDoubleArray* MinEigenvecs, int* dims)
{

    for ( int k = 0; k < 1; k++)
    {
        for (int j = 0; j < dims[1] - 1; j++)
        {
            std::cout<<"Flip Progress "<< int(double (j+1)/double (dims[1])* 100 ) <<"%\r";

            for (int i = 0; i < dims[0] - 1; i++)
            {
                double minEvec1[2], minEvec2[2], minEvec3[2], minEvec4[2], pcaMaxEvec[2] = {0.};

                int idME = getIndex(k, j, i, dims);
                vtkIdType idMinEvec1 = vtkIdType(idME);
                vtkIdType idMinEvec2 = vtkIdType(idME);
                vtkIdType idMinEvec3 = vtkIdType(idME);
                vtkIdType idMinEvec4 = vtkIdType(idME);

                //// x,y

                idMinEvec1 = getIndex(k, j, i, dims);
                MinEigenvecs->GetTuple(idMinEvec1, minEvec1);
                pcaMaxEigenVec->GetTuple(idMinEvec1, pcaMaxEvec);

                double DotProdEigens1 = DotProd(minEvec1, pcaMaxEvec);

                //std::cout<< DotProd (minEvec1, pcaMaxEvec) <<std::endl;

                if (DotProdEigens1 < 0)
                {
                    minEvec1[0] = minEvec1[0] * (- 1.);
                    minEvec1[1] = minEvec1[1] * (- 1.);

                    MinEigenvecs->SetTuple2(idMinEvec1, minEvec1[0], minEvec1[1]);

                    //std::cout<< DotProd (minEvec1, pcaMaxEvec) <<std::endl;

                }

                //// x + 1, y

                idMinEvec2 = getIndex(k, j, i+1, dims);
                MinEigenvecs->GetTuple(idMinEvec2, minEvec2);

                double DotProdEigens2 = DotProd(minEvec2, pcaMaxEvec);

                //std::cout<< DotProd (minEvec2, pcaMaxEvec) <<std::endl;

                if (DotProdEigens2 < 0)
                {
                    minEvec2[0] = minEvec2[0] * (- 1.);
                    minEvec2[1] = minEvec2[1] * (- 1.);

                    MinEigenvecs->SetTuple2(idMinEvec2, minEvec2[0], minEvec2[1]);
                    //std::cout<< DotProd (minEvec2, pcaMaxEvec) <<std::endl;

                }

                //// x, y + 1

                idMinEvec3 = getIndex(k, j+1, i, dims);
                MinEigenvecs->GetTuple(idMinEvec3, minEvec3);

                double DotProdEigens3 = DotProd(minEvec3, pcaMaxEvec);

                //std::cout<< DotProd (minEvec3, pcaMaxEvec) <<std::endl;

                if (DotProdEigens3 < 0)
                {
                    minEvec3[0] = minEvec3[0] * (- 1.);
                    minEvec3[1] = minEvec3[1] * (- 1.);

                    MinEigenvecs->SetTuple2(idMinEvec3, minEvec3[0], minEvec3[1]);
                    //std::cout<< DotProd (minEvec3, pcaMaxEvec) <<std::endl;

                }

                //// x + 1, y + 1

                idMinEvec4 = getIndex(k, j+1, i+1, dims);
                MinEigenvecs->GetTuple(idMinEvec4, minEvec4);

                double DotProdEigens4 = DotProd(minEvec4, pcaMaxEvec);

                //std::cout<< DotProd (minEvec4, pcaMaxEvec) <<std::endl;

                if (DotProdEigens4 < 0)
                {
                    minEvec4[0] = minEvec4[0] * (- 1.);
                    minEvec4[1] = minEvec4[1] * (- 1.);

                    MinEigenvecs->SetTuple2(idMinEvec4, minEvec4[0], minEvec4[1]);
                    //std::cout<< DotProd (minEvec4, pcaMaxEvec) <<std::endl;

                }

                //std::cout<< DotProdEigens1 << " " << DotProdEigens2 << " "  << " " << DotProdEigens3 << " "<<DotProdEigens4<<std::endl;

            }
        }
    }

    std::cout<<"DotProd calculation finished" <<std::endl;
    return;
}


//// Ridge Points calculation

void Ridges( vtkImageData* input, vtkDoubleArray* gradient, vtkDoubleArray* hessian, vtkDoubleArray* MinEigenvecs, int* dims, vtkPoints* RidgePoints)
{

    int id_ridge = 0;
    vtkIdType X_Id = vtkIdType(id_ridge);
    X_Id = id_ridge;
    vtkIdType Y_Id = vtkIdType(id_ridge);
    Y_Id = id_ridge + 1;

    int number = 0;

    //// Dot product of each gradient with updated minimum eigenvectors

    for (int k = 0; k < 1; k++)
    {
        for (int j = 0; j < dims[1] - 1; j++)
        {
            std::cout<<"Ridge Pt Progress " << int(double (j+1)/double (dims[1])* 100 ) <<"%\r";
            for (int i = 0; i < dims[0] - 1; i++)
            {
                double grad1[2], grad2[2], grad3[2], grad4[2], hess1[4], hess2[4], hess3[4], hess4[4], IntHess1[4], IntHess2[4], IntHess3[4], IntHess4[4],
                        minEV1[2], minEV2[2], minEV3[2], minEV4[2], point1[3], point2[3], point3[3], point4[3] = {0.};

                int idxx = getIndex(k, j, i, dims);
                vtkIdType ideg1 = vtkIdType(idxx);
                vtkIdType ideg2 = vtkIdType(idxx);
                vtkIdType ideg3 = vtkIdType(idxx);
                vtkIdType ideg4 = vtkIdType(idxx);


                ////Get gradient, minimum eigenvector and coordinate values at each node


                //// x, y

                ideg1 = getIndex(k, j, i, dims);
                input->GetPoint(ideg1, point1);
                gradient->GetTuple(ideg1, grad1);
                hessian->GetTuple(ideg1, hess1);
                MinEigenvecs->GetTuple(ideg1, minEV1);

                double DotProd1 = DotProd(minEV1, grad1);

                //std::cout << point1[0] << " " << point1[1] << " " << point1[2] << std::endl;


                //// x + 1, y

                ideg2 = getIndex(k, j, i + 1, dims);
                input->GetPoint(ideg2, point2);
                gradient->GetTuple(ideg2, grad2);
                hessian->GetTuple(ideg2, hess2);
                MinEigenvecs->GetTuple(ideg2, minEV2);

                double DotProd2 = DotProd(minEV2, grad2);

                //std::cout << DotProd2 << std::endl;
                //std::cout << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;


                //// x, y + 1
                ideg3 = getIndex(k, j + 1, i, dims);
                input->GetPoint(ideg3, point3);
                gradient->GetTuple(ideg3, grad3);
                hessian->GetTuple(ideg3, hess3);
                MinEigenvecs->GetTuple(ideg3, minEV3);

                double DotProd3 = DotProd(minEV3, grad3);

                //std::cout << DotProd3 << std::endl;
                //std::cout << point3[0] << " " << point3[1] << " " << point3[2] << std::endl;


                //// x + 1, y + 1
                ideg4 = getIndex(k, j + 1, i + 1, dims);
                input->GetPoint(ideg4, point4);
                gradient->GetTuple(ideg4, grad4);
                hessian->GetTuple(ideg4, hess4);
                MinEigenvecs->GetTuple(ideg4, minEV4);


                double DotProd4 = DotProd(minEV4, grad4);



                ///////// Lookup Table cases //////////


                double phi, ro, phi2, ro2 = 0.;
                double RidgeX[3], RidgeY[3], RidgeX2[3], RidgeY2[3] = {0.};

            //// Case 1: ----, No ridges

                if (DotProd1 < 0 && DotProd2 < 0 && DotProd3 < 0 && DotProd4 < 0)
                {
                    continue;
                }


            //// Case 2: +---

                if (DotProd1 > 0 && DotProd2 < 0 && DotProd3 < 0 && DotProd4 < 0)
                {

                    //// Calculate x coordinate of ridge line on bottom edge

                    phi = DotProd1 / (DotProd1 - DotProd2);                      /// lamda = f_l / (f_l - f_r);
                    RidgeX[0] = ( 1 - phi ) * point1[0] + phi * point2[0];       /// x = (1 - lambda) * x_l + lambda * x_r
                    RidgeX[1] = point1[1];
                    RidgeX[2] = 0.;

                    /// Calculate y coordinate of ridge line on left edge

                    ro = DotProd1 / (DotProd1 - DotProd3);
                    RidgeY[0] = point1[0];
                    RidgeY[1] = ( 1 - ro ) * point1[1] + ro * point3[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                        IntHess2[l] = (1 - ro) * hess1[l] + ro * hess3[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout << "Case 2: " << "X line " << RidgeX[0] << " " << RidgeX[1] << " " << RidgeX[2]
                                      << std::endl;
                            std::cout << "Case 2: " << "Y line " << RidgeY[0] << " " << RidgeY[1] << " " << RidgeY[2]
                                      << std::endl;
                        }*/
                    }




                }


            //// Case 3: -+--

                else if (DotProd1 < 0 && DotProd2 > 0 && DotProd3 < 0 && DotProd4 < 0)
                {

                    //// X on bottom edge

                    phi = DotProd1 / (DotProd1 - DotProd2);
                    RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                    RidgeX[1] = point2[1];
                    RidgeX[2] = 0.;

                    //// Y on Right edge

                    ro = DotProd2 / (DotProd2 - DotProd4);
                    RidgeY[0] = point2[0];
                    RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                    RidgeY[2] = 0.;


                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                        IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {

                            std::cout<<"Case 3: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 3: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }

                }


            //// Case 4: --+-

                else if (DotProd1 < 0 && DotProd2 < 0 && DotProd3 > 0 && DotProd4 < 0)
                {
                    //// X on top edge

                    phi = DotProd3 / (DotProd3 - DotProd4);
                    RidgeX[0] = ( 1 - phi ) * point3[0] + phi * point4[0];
                    RidgeX[1] = point3[1];
                    RidgeX[2] = 0.;


                    //// Y on left edge

                    ro = DotProd1 / (DotProd1 - DotProd3);
                    RidgeY[0] = point3[0];
                    RidgeY[1] = ( 1 - ro ) * point1[1] + ro * point3[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess3[l] + phi * hess4[l];
                        IntHess2[l] = (1 - ro) * hess1[l] + ro * hess3[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout<<"Case 4: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 4: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }


                }


            //// Case 5: ---+

                else if (DotProd1 < 0 && DotProd2 < 0 && DotProd3 < 0 && DotProd4 > 0)
                {

                    //// X on top edge

                    phi = DotProd3 / (DotProd3 - DotProd4);
                    RidgeX[0] = ( 1 - phi ) * point3[0] + phi * point4[0];
                    RidgeX[1] = point4[1];
                    RidgeX[2] = 0.;


                    //// Y on right edge
                    ro = DotProd2 / (DotProd2 - DotProd4);
                    RidgeY[0] = point4[0];
                    RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess3[l] + phi * hess4[l];
                        IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout<<"Case 5: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 5: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }

                }


            //// Case 6: ++--

                else if (DotProd1 > 0 && DotProd2 > 0 && DotProd3 < 0 && DotProd4 < 0)
                {

                    //// Y on left edge

                    phi = DotProd1 / ( DotProd1 - DotProd3 );
                    RidgeX[0] = point1[0];
                    RidgeX[1] = (1 - phi) * point1[1] + phi * point3[1];
                    RidgeX[2] = 0.;

                    //// Y on right edge

                    ro = DotProd2/( DotProd2 - DotProd4 );
                    RidgeY[0] = point2[0];
                    RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess3[l];
                        IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];

                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout<<"Case 6: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 6: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }

                }


            //// Case 7: +-+-

                else if (DotProd1 > 0 && DotProd2 < 0 && DotProd3 > 0 && DotProd4 < 0)
                {
                    //// X on bottom edge

                    phi = DotProd1/(DotProd1 - DotProd2);
                    RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                    RidgeX[1] = point1[1];
                    RidgeX[2] = 0.;

                    //// X on top edge

                    ro = DotProd3/(DotProd3 - DotProd4);
                    RidgeY[0] = ( 1 - ro ) * point3[0] + ro * point4[0];
                    RidgeY[1] = point3[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                        IntHess2[l] = (1 - ro) * hess3[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {

                            std::cout<<"Case 7: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 7: "<< "X line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }


                }


            //// Case 8: +--+

                else if (DotProd1 > 0 && DotProd2 < 0 && DotProd3 < 0 && DotProd4 > 0)
                {

                    //// Interpolate the function bilinearly, we get the equation
                    //// f(x,y) = f_i,j * ( 1 - x ) * ( 1 - y ) + f_i+1,j * x * (1 - y) + f_i,j+1 ( 1 - x) * y + f_i+1,j+1 * x * y
                    //// f(x,y) = A*x*y + B*x + C*y + D

                    //// A = f_i,j - f_i+1,j - f_i,j+1 + f_i+1,j+1
                    double A = DotProd1 - DotProd2 - DotProd3 + DotProd4;

                    //// B = - f_i,j + f_i+1,j
                    double B = DotProd2 - DotProd1;

                    //// C = - f_i,j + f_i,j+1
                    double C = DotProd3 - DotProd1;

                    //// D = f_i,j
                    double D = DotProd1;

                    double AsympDecider = D - (B * C) / A;

                    double c_contour = - 10 * ( -1 * C / A - 0.3 ) * ( -1 * B / A - 0.5 ) + 4.5;

                    //// id D - BC/A > c, (c = 0), choose left case, else right

                    if ( AsympDecider > c_contour )
                    {

                        //// line 1: X component on bottom edge and Y component on left edge
                        phi = DotProd1 / (DotProd1 - DotProd2);
                        RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                        RidgeX[1] = point1[1];
                        RidgeX[2] = 0.;

                        ro = DotProd1 / (DotProd1 - DotProd3);
                        RidgeY[0] = point1[0];
                        RidgeY[1] = (1 - ro) * point1[1] + ro * point3[1];
                        RidgeY[2] = 0.;


                        ////line 2: X component on top edge and Y component on right edge
                        phi2 = DotProd3 / (DotProd3 - DotProd4);
                        RidgeX2[0] = (1 - phi2) * point3[0] + phi2 * point4[0];
                        RidgeX2[1] = point3[1];
                        RidgeX2[2] = 0.;

                        ro2 = DotProd2 / (DotProd2 - DotProd4);
                        RidgeY2[0] = point2[0];
                        RidgeY2[1] = (1 - ro2) * point2[1] + ro2 * point4[1];
                        RidgeY2[2] = 0.;


                        for (int l = 0; l < 4; l++)
                        {
                            IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                            IntHess2[l] = (1 - ro) * hess1[l] + ro * hess3[l];
                            IntHess3[l] = (1 - phi2) * hess3[l] + phi2 * hess4[l];
                            IntHess4[l] = (1 - ro2) * hess2[l] + ro2 * hess4[l];

                        }

                        if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                        {

                            RidgePoints->InsertPoint(X_Id, RidgeX);
                            RidgePoints->InsertPoint(Y_Id, RidgeY);

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {
                                std::cout<<"Case 8 first: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                                std::cout<<"Case 8 first: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;

                            }*/



                        }

                        if (RidgeEigenValue(IntHess3)  && RidgeEigenValue(IntHess4))
                        {
                            RidgePoints->InsertPoint(X_Id, RidgeX2 );
                            RidgePoints->InsertPoint(Y_Id, RidgeY2 );

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {
                                std::cout<<"Case 8 second: "<< "X line " << RidgeX2[0] << " " << RidgeX2[1] << " "<<RidgeX2[2] <<std::endl;
                                std::cout<<"Case 8 second: "<< "Y line " << RidgeY2[0] << " " << RidgeY2[1] << " "<<RidgeY2[2] <<std::endl;

                            }*/


                        }

                    }

                    else
                    {
                        //// line 1: X component on bottom edge and Y component on right edge
                        phi = DotProd1 / (DotProd1 - DotProd2);
                        RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                        RidgeX[1] = point1[1];
                        RidgeX[2] = 0.;

                        ro = DotProd2 / (DotProd2 - DotProd4);
                        RidgeY[0] = point2[0];
                        RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                        RidgeY[2] = 0.;

                        ////line 2: X component on top edge and Y component on left edge
                        phi2 = DotProd3 / (DotProd3 - DotProd4);
                        RidgeX2[0] = (1 - phi2) * point3[0] + phi2 * point4[0];
                        RidgeX2[1] = point3[1];
                        RidgeX2[2] = 0.;

                        ro2 = DotProd1 / (DotProd1 - DotProd3);
                        RidgeY2[0] = point1[0];
                        RidgeY2[1] = (1 - ro2) * point1[1] + ro2 * point3[1];
                        RidgeY2[2] = 0.;


                        for (int l = 0; l < 4; l++)
                        {
                            IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                            IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                            IntHess3[l] = (1 - phi2) * hess3[l] + phi2 * hess4[l];
                            IntHess4[l] = (1 - ro2) * hess1[l] + ro2 * hess3[l];

                        }

                        if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                        {

                            RidgePoints->InsertPoint(X_Id, RidgeX);
                            RidgePoints->InsertPoint(Y_Id, RidgeY);

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {

                                std::cout<<"Case 8 first: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                                std::cout<<"Case 8 first: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                            }*/
                        }

                        if (RidgeEigenValue(IntHess3)  && RidgeEigenValue(IntHess4))
                        {
                            RidgePoints->InsertPoint(X_Id, RidgeX2 );
                            RidgePoints->InsertPoint(Y_Id, RidgeY2 );

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {
                                std::cout<<"Case 8 second: "<< "X line " << RidgeX2[0] << " " << RidgeX2[1] << " "<<RidgeX2[2] <<std::endl;
                                std::cout<<"Case 8 second: "<< "Y line " << RidgeY2[0] << " " << RidgeY2[1] << " "<<RidgeY2[2] <<std::endl;

                            }*/

                        }


                    }


                    //std::cout<<"Case 8: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                    //std::cout<<"Case 8: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                    //std::cout<<"Case 8: "<< "X line " << RidgeX2[0] << " " << RidgeX2[1] << " "<<RidgeX2[2] <<std::endl;
                    //std::cout<<"Case 8: "<< "Y line " << RidgeY2[0] << " " << RidgeY2[1] << " "<<RidgeY2[2] <<std::endl;

                }


            //// Case 9: -++-

                else if (DotProd1 < 0 && DotProd2 > 0 && DotProd3 > 0 && DotProd4 < 0)
                {
                    //// Interpolate the function bilinearly, we get the equation
                    //// f(x,y) = f_i,j * ( 1 - x ) * ( 1 - y ) + f_i+1,j * x * (1 - y) + f_i,j+1 ( 1 - x) * y + f_i+1,j+1 * x * y
                    //// f(x,y) = A*x*y + B*x + C*y + D

                    //// A = f_i,j - f_i+1,j - f_i,j+1 + f_i+1,j+1
                    double A = DotProd1 - DotProd2 - DotProd3 + DotProd4;

                    //// B = - f_i,j + f_i+1,j
                    double B = DotProd2 - DotProd1;

                    //// C = - f_i,j + f_i,j+1
                    double C = DotProd3 - DotProd1;

                    //// D = f_i,j
                    double D = DotProd1;

                    //// Asymptotic Decider
                    double AsympDecider = D - (B * C) / A;

                    double c_contour = - 10 * ( -1 * C / A - 0.3 ) * ( -1 * B / A - 0.5 ) + 4.5;

                    //// id D - BC/A > c, (c = 0), choose left case, else right
                    if ( AsympDecider > c_contour )
                    {
                        //// line 1: X component on bottom edge and Y component on left edge
                        phi = DotProd1 / (DotProd1 - DotProd2);
                        RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                        RidgeX[1] = point1[1];
                        RidgeX[2] = 0.;

                        ro = DotProd1 / (DotProd1 - DotProd3);
                        RidgeY[0] = point1[0];
                        RidgeY[1] = (1 - ro) * point1[1] + ro * point3[1];
                        RidgeY[2] = 0.;

                        ////line 2: X component on top edge and Y component on right edge
                        phi2 = DotProd3 / (DotProd3 - DotProd4);
                        RidgeX2[0] = (1 - phi2) * point3[0] + phi2 * point4[0];
                        RidgeX2[1] = point3[1];
                        RidgeX2[2] = 0.;

                        ro2 = DotProd2 / (DotProd2 - DotProd4);
                        RidgeY2[0] = point2[0];
                        RidgeY2[1] = (1 - ro2) * point2[1] + ro2 * point4[1];
                        RidgeY2[2] = 0.;



                        for (int l = 0; l < 4; l++)
                        {
                            IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                            IntHess2[l] = (1 - ro) * hess1[l] + ro * hess3[l];
                            IntHess3[l] = (1 - phi2) * hess3[l] + phi2 * hess4[l];
                            IntHess4[l] = (1 - ro2) * hess2[l] + ro2 * hess4[l];

                        }

                        if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                        {

                            RidgePoints->InsertPoint(X_Id, RidgeX);
                            RidgePoints->InsertPoint(Y_Id, RidgeY);

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {

                                std::cout<<"Case 9 first: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                                std::cout<<"Case 9 first: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                            }*/
                        }

                        if (RidgeEigenValue(IntHess3)  && RidgeEigenValue(IntHess4))
                        {
                            RidgePoints->InsertPoint(X_Id, RidgeX2 );
                            RidgePoints->InsertPoint(Y_Id, RidgeY2 );

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {
                                std::cout<<"Case 9 second: "<< "X line " << RidgeX2[0] << " " << RidgeX2[1] << " "<<RidgeX2[2] <<std::endl;
                                std::cout<<"Case 9 second: "<< "Y line " << RidgeY2[0] << " " << RidgeY2[1] << " "<<RidgeY2[2] <<std::endl;
                            }*/

                        }

                    }

                    else
                    {
                        //// line 1: X component on bottom edge and Y component on right edge
                        phi = DotProd1 / (DotProd1 - DotProd2);
                        RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                        RidgeX[1] = point1[1];
                        RidgeX[2] = 0.;

                        ro = DotProd2 / (DotProd2 - DotProd4);
                        RidgeY[0] = point2[0];
                        RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                        RidgeY[2] = 0.;

                        ////line 2: X component on top edge and Y component on left edge
                        phi2 = DotProd3 / (DotProd3 - DotProd4);
                        RidgeX2[0] = (1 - phi2) * point3[0] + phi2 * point4[0];
                        RidgeX2[1] = point3[1];
                        RidgeX2[2] = 0.;

                        ro2 = DotProd1 / (DotProd1 - DotProd3);
                        RidgeY2[0] = point1[0];
                        RidgeY2[1] = (1 - ro2) * point1[1] + ro2 * point3[1];
                        RidgeY2[2] = 0.;


                        for (int l = 0; l < 4; l++)
                        {
                            IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                            IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                            IntHess3[l] = (1 - phi2) * hess3[l] + phi2 * hess4[l];
                            IntHess4[l] = (1 - ro2) * hess1[l] + ro2 * hess3[l];

                        }

                        if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                        {

                            RidgePoints->InsertPoint(X_Id, RidgeX);
                            RidgePoints->InsertPoint(Y_Id, RidgeY);



                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {
                                std::cout<<"Case 9 first: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                                std::cout<<"Case 9 first: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                            }*/

                        }

                        if (RidgeEigenValue(IntHess3)  && RidgeEigenValue(IntHess4))
                        {
                            RidgePoints->InsertPoint(X_Id, RidgeX2 );
                            RidgePoints->InsertPoint(Y_Id, RidgeY2 );

                            X_Id = X_Id + 2;
                            Y_Id = Y_Id + 2;

                            /*
                            if ( i == dims[0] - 1 || j == dims[1] - 1) {
                                std::cout<<"Case 9 second: "<< "X line " << RidgeX2[0] << " " << RidgeX2[1] << " "<<RidgeX2[2] <<std::endl;
                                std::cout<<"Case 9 second: "<< "Y line " << RidgeY2[0] << " " << RidgeY2[1] << " "<<RidgeY2[2] <<std::endl;
                            }*/
                        }

                    }

                    //std::cout<<"Case 9: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                    //std::cout<<"Case 9: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                    //std::cout<<"Case 9: "<< "X line " << RidgeX2[0] << " " << RidgeX2[1] << " "<<RidgeX2[2] <<std::endl;
                    //std::cout<<"Case 9: "<< "Y line " << RidgeY2[0] << " " << RidgeY2[1] << " "<<RidgeY2[2] <<std::endl;

                }


            //// Case 10: -+-+

                else if (DotProd1 < 0 && DotProd2 > 0 && DotProd3 < 0 && DotProd4 > 0)
                {
                    //// X on bottom edge

                    phi = DotProd1/(DotProd1 - DotProd2);
                    RidgeX[0] = (1 - phi) * point1[0] + phi * point2[0];
                    RidgeX[1] = point1[1];
                    RidgeX[2] = 0.;


                    //// X on top edge

                    ro = DotProd3/(DotProd3 - DotProd4);
                    RidgeY[0] = ( 1 - ro ) * point3[0] + ro * point4[0];
                    RidgeY[1] = point3[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                        IntHess2[l] = (1 - ro) * hess3[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout<<"Case 10: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 10: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }


                }


            //// Case 11: --++

                else if (DotProd1 < 0 && DotProd2 < 0 && DotProd3 > 0 && DotProd4 > 0)
                {
                    //// Y on left edge

                    phi = DotProd1 / ( DotProd1 - DotProd3 );
                    RidgeX[0] = point1[0];
                    RidgeX[1] = (1 - phi) * point1[1] + phi * point3[1];
                    RidgeX[2] = 0.;

                    //// Y on right edge

                    ro = DotProd2/( DotProd2 - DotProd4 );
                    RidgeY[0] = point2[0];
                    RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess3[l];
                        IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {

                            std::cout << "Case 11: " << "X line " << RidgeX[0] << " " << RidgeX[1] << " " << RidgeX[2]
                                      << std::endl;
                            std::cout << "Case 11: " << "Y line " << RidgeY[0] << " " << RidgeY[1] << " " << RidgeY[2]
                                      << std::endl;
                        }*/

                    }

                }


            //// Case 12: -+++

                else if (DotProd1 < 0 && DotProd2 > 0 && DotProd3 > 0 && DotProd4 > 0)
                {
                    //// Calculate x coordinate of ridge line on bottom edge

                    phi = DotProd1 / (DotProd1 - DotProd2);
                    RidgeX[0] = ( 1 - phi ) * point1[0] + phi * point2[0];
                    RidgeX[1] = point1[1];
                    RidgeX[2] = 0.;

                    //// Calculate y coordinate of ridge line on left edge

                    ro = DotProd1 / (DotProd1 - DotProd3);
                    RidgeY[0] = point1[0];
                    RidgeY[1] = ( 1 - ro ) * point1[1] + ro * point3[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                        IntHess2[l] = (1 - ro) * hess1[l] + ro * hess3[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout << "Case 12: " << "X line " << RidgeX[0] << " " << RidgeX[1] << " " << RidgeX[2]
                                      << std::endl;
                            std::cout << "Case 12: " << "Y line " << RidgeY[0] << " " << RidgeY[1] << " " << RidgeY[2]
                                      << std::endl;
                        }*/
                    }


                }


            //// Case 13: +-++

                else if (DotProd1 > 0 && DotProd2 < 0 && DotProd3 > 0 && DotProd4 > 0)
                {
                    //// X on bottom edge

                    phi = DotProd1 / (DotProd1 - DotProd2);
                    RidgeX[0] = ( 1 - phi ) * point1[0] + phi * point2[0];
                    RidgeX[1] = point2[1];
                    RidgeX[2] = 0.;


                    //// Y on Right edge

                    ro = DotProd2 / (DotProd2 - DotProd4);
                    RidgeY[0] = point2[0];
                    RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess1[l] + phi * hess2[l];
                        IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout << "Case 13: " << "X line " << RidgeX[0] << " " << RidgeX[1] << " " << RidgeX[2]
                                      << std::endl;
                            std::cout << "Case 13: " << "Y line " << RidgeY[0] << " " << RidgeY[1] << " " << RidgeY[2]
                                      << std::endl;
                        }*/
                    }

                }


                    //


            //// Case 14: ++-+

                else if (DotProd1 > 0 && DotProd2 > 0 && DotProd3 < 0 && DotProd4 > 0)
                {
                    //// X on top edge

                    phi = DotProd3 / (DotProd3 - DotProd4);
                    RidgeX[0] = ( 1 - phi ) * point3[0] + phi * point4[0];
                    RidgeX[1] = point3[1];
                    RidgeX[2] = 0.;


                    //// Y on left edge

                    ro = DotProd1 / (DotProd1 - DotProd3);
                    RidgeY[0] = point3[0];
                    RidgeY[1] = ( 1 - ro ) * point1[1] + ro * point3[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess3[l] + phi * hess4[l];
                        IntHess2[l] = (1 - ro) * hess1[l] + ro * hess3[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 1 || j == dims[1] - 1) {
                            std::cout<<"Case 14: "<< "X line " << RidgeX[0] << " " << RidgeX[1] << " "<<RidgeX[2] <<std::endl;
                            std::cout<<"Case 14: "<< "Y line " << RidgeY[0] << " " << RidgeY[1] << " "<<RidgeY[2] <<std::endl;
                        }*/
                    }


                }


            //// Case 15: +++-

                else if (DotProd1 > 0 && DotProd2 > 0 && DotProd3 > 0 && DotProd4 < 0)
                {

                    //// X on top edge

                    phi = DotProd3 / (DotProd3 - DotProd4);
                    RidgeX[0] = ( 1 - phi ) * point3[0] + phi * point4[0];
                    RidgeX[1] = point4[1];
                    RidgeX[2] = 0.;


                    //// Y on right edge
                    ro = DotProd2 / (DotProd2 - DotProd4);
                    RidgeY[0] = point4[0];
                    RidgeY[1] = (1 - ro) * point2[1] + ro * point4[1];
                    RidgeY[2] = 0.;

                    for (int l = 0; l < 4; l++)
                    {
                        IntHess1[l] = (1 - phi) * hess3[l] + phi * hess4[l];
                        IntHess2[l] = (1 - ro) * hess2[l] + ro * hess4[l];
                    }


                    if (RidgeEigenValue(IntHess1)  && RidgeEigenValue(IntHess2))
                    {

                        RidgePoints->InsertPoint(X_Id, RidgeX);
                        RidgePoints->InsertPoint(Y_Id, RidgeY);

                        X_Id = X_Id + 2;
                        Y_Id = Y_Id + 2;

                        /*
                        if ( i == dims[0] - 2 || j == dims[1] - 2) {
                            std::cout << "Case 15: " << "X line " << RidgeX[0] << " " << RidgeX[1] << " " << RidgeX[2]
                                      << std::endl;
                            std::cout << "Case 15: " << "Y line " << RidgeY[0] << " " << RidgeY[1] << " " << RidgeY[2]
                                      << std::endl;
                        }*/
                    }


                }


            //// Case 16: ++++, no ridges

                else if (DotProd1 > 0 && DotProd2 > 0 && DotProd3 > 0 && DotProd4 > 0)
                {
                    continue;
                }



            }
        }
    }

    std::cout<<"RidgePoints calculation finished" <<std::endl;
    //std::cout<<"number of ridge points in loop " <<number <<std::endl;

    return;
}

bool RidgeEigenValue (double IntHessian[4])
{

    Eigen::Matrix2d IntHessianmatrix(2,2);
    IntHessianmatrix(0, 0) = IntHessian[0];
    IntHessianmatrix(0, 1) = IntHessian[1];
    IntHessianmatrix(1, 0) = IntHessian[2];
    IntHessianmatrix(1, 1) = IntHessian[3];


    ////calculating eigenvectors and eigenvalues using Eigen
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(IntHessianmatrix, true);


    ////Get minimum eigenvalue and the corresponding eigenvector
    Eigen::MatrixXd IntValueReal(2,1);
    Eigen::MatrixXd IntVectorReal(2,2);

    IntValueReal = es.eigenvalues().real();     ///to get the real part of Evalues and Evectors
    IntVectorReal = es.eigenvectors().real();


    double lambda[2];
    double RealEV[2];

    lambda[0] = IntValueReal(0,0);
    lambda[1] = IntValueReal(1,0);


    if (lambda[0] < limit || lambda[1] < limit)
    {
        //std::cout << "eigenvalues " << lambda[0] << ", " << lambda[1] << ", " <<  limit  << std::endl;

        return true;
    }

    else
    {
        return false;
    }

}
