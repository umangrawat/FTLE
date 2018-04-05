#include "vtkColorMap.h"

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
#include <vtkPoints.h>
#include <vtkDataObject.h>

#include "vtkMath.h"
#include "vtkNew.h"

#include <vtkCellArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkFieldData.h>
#include <vtkLine.h>
#include <vtkUnsignedCharArray.h>

#include <vtkAssignAttribute.h>
#include <vtkCellData.h>

#include <random>
#include <vector>
#include <limits>
#include <iomanip>

using namespace std;


vtkStandardNewMacro(vtkColorMap);

//-----------------------------------------------------------------------------
vtkColorMap::vtkColorMap()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                 vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
vtkColorMap::~vtkColorMap()
{

}

//----------------------------------------------------------------------------
int vtkColorMap::FillInputPortInformation( int port, vtkInformation* info )
{

    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
        return 1;
    }
    return 0;
}


//----------------------------------------------------------------------------
int vtkColorMap::FillOutputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 )
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }
    if( port == 1){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );
        return 1;
    }

    return 0;
}


//----------------------------------------------------------------------------
int vtkColorMap::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    return 1;

}

//----------------------------------------------------------------------------
int vtkColorMap::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    //// Get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkImageData *output1 = vtkImageData::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

    int* dims = input->GetDimensions();
    double* spacing = input->GetSpacing();

    this->seedGrid = vtkSmartPointer<vtkImageData>::New();
    this->seedGrid->SetDimensions(input->GetDimensions());
    this->seedGrid->SetSpacing(input->GetSpacing());

    vtkIdType CellId = vtkIdType(this->cellId);

    //// GEt the required input arrays
    vtkSmartPointer<vtkDoubleArray> cellIndexArray = vtkSmartPointer<vtkDoubleArray>::New();
    cellIndexArray = vtkDoubleArray::SafeDownCast(input->GetCellData()->GetAbstractArray("Indices"));

    vtkSmartPointer<vtkDoubleArray> StartcellPointsArray = vtkSmartPointer<vtkDoubleArray>::New();
    StartcellPointsArray = vtkDoubleArray::SafeDownCast(input->GetCellData()->GetAbstractArray("Start Points"));

    vtkSmartPointer<vtkDoubleArray> EndcellPointsArray = vtkSmartPointer<vtkDoubleArray>::New();
    EndcellPointsArray = vtkDoubleArray::SafeDownCast(input->GetCellData()->GetAbstractArray("End Points"));

    ////Get the array of indices for the particular cell
    double cellIndex[static_cast<int>(cellIndexArray->GetNumberOfComponents())] = {0.};
    cellIndexArray->GetTuple(CellId, cellIndex);


    ////Get the starting points for the particular cell
    double StartcellPoints[static_cast<int>(StartcellPointsArray->GetNumberOfComponents())] = {0.};
    StartcellPointsArray->GetTuple(CellId, StartcellPoints);

    vtkSmartPointer<vtkDoubleArray> cellPoints = vtkSmartPointer<vtkDoubleArray>::New();
    cellPoints->SetNumberOfComponents(3);

    for (int c = 0; c < StartcellPointsArray->GetNumberOfComponents(); c++)
    {
        double cellpt[2] = {StartcellPoints[c*2], StartcellPoints[c*2+1]};
        vtkIdType idc = vtkIdType(c);
        cellPoints->InsertTuple3(idc, cellpt[0],cellpt[1],0.);
    }


    vtkSmartPointer<vtkPoints> cell = vtkSmartPointer<vtkPoints>::New();

    vtkIdType idPoint;

    //// Just taking the first point for locating cell
    int xIndex = static_cast<int>( StartcellPoints[0] / spacing[0] );
    int yIndex = static_cast<int>( StartcellPoints[1] / spacing[1] );

    ////Getting the bottom point for the cell
    double cellLeft[3] = {xIndex * spacing[0], yIndex * spacing[1], 0.};
    double cellRight[3] = {xIndex * spacing[0] + spacing[0], yIndex * spacing[1], 0.};
    double cellLeftTop[3] = {xIndex * spacing[0], yIndex * spacing[1] + spacing[1], 0.};
    double cellRightTop[3] = {xIndex * spacing[0] + spacing[0], yIndex * spacing[1] + spacing[1], 0.};

    //// Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(4);

    vtkIdType nextPoint1 = cell->InsertNextPoint(cellLeft);
    lines->InsertCellPoint(nextPoint1);
    vtkIdType nextPoint2 = cell->InsertNextPoint(cellRight);
    lines->InsertCellPoint(nextPoint2);
    vtkIdType nextPoint3 = cell->InsertNextPoint(cellRightTop);
    lines->InsertCellPoint(nextPoint3);
    vtkIdType nextPoint4 = cell->InsertNextPoint(cellLeftTop);
    lines->InsertCellPoint(nextPoint4);



    ////Pairing and sorting

    std::vector<pair<int, int>> vect;
    int n = sizeof(cellIndex)/sizeof(*cellIndex);

    int arr0[n] = {0};
    int arr1[n] = {0};

    for (int i = 0; i < n; i++)
    {
        arr0[i] = cellIndex[i];
        arr1[i] = i;

    }

    for (int m = 0; m < n; m++) {

        vect.push_back( make_pair (arr0[m], arr1[m]));
    }

    std::sort(vect.begin(), vect.end());

    for (int p=0; p<n; p++)
    {
        //std::cout << vect[p].first << " " << vect[p].second << std::endl;
    }


    vtkSmartPointer<vtkUnsignedCharArray> colorsPts = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colorsPts->SetNumberOfComponents(4);
    colorsPts->SetName ("Colors Points");

    vtkSmartPointer<vtkUnsignedCharArray> colorsCells = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colorsCells->SetNumberOfComponents(4);
    colorsCells->SetNumberOfTuples((dims[0] - 1)*(dims[1] - 1));
    colorsCells->SetName ("Colors Cells");

    ///setting intial values to 0

    for (int idzero = 0; idzero < colorsCells->GetNumberOfTuples(); idzero++)
    {
        unsigned char zeropt[4] = {255,255,255, 0};
        colorsCells->InsertTypedTuple(idzero,zeropt);
    }

    vtkSmartPointer<vtkPoints> outputPoints1 = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> outputCells = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();


    int idvec = 0;

    while ( vect.size() > 0 )
    {
        int index1 = static_cast<int>(vect[0].first);

        int x, y = 0;

        double point[3]={0.};
        cellPoints->GetTuple(idvec, point);
        //std::cout<<"index " << index1<< " " << point[0] << " "<< point[1]<<std::endl;
        outputPoints1->InsertNextPoint(point);

        idvec+= 1;
        int delIter = 0;

        unsigned char tempColor[4] = {(int)std::rand()%255,(int)std::rand()%255 ,(int)std::rand()%255, 255};

        cellArray->InsertNextCell(4);
        cellLocator(index1,dims,x,y);

        double cellpoint1[3] = {static_cast<double>(x) * spacing[0],static_cast<double>(y)*spacing[1], 0.};
        double cellpoint2[3] = {static_cast<double>(x) * spacing[0] + spacing[0], static_cast<double>(y) * spacing[1], 0};
        double cellpoint3[3] = {static_cast<double>(x) * spacing[0], static_cast<double>(y) * spacing[1] + spacing[1], 0};
        double cellpoint4[3] = {static_cast<double>(x) * spacing[0] + spacing[0], static_cast<double>(y) * spacing[1] + spacing[1], 0};

        vtkIdType pt1 = outputCells->InsertNextPoint(cellpoint1);
        cellArray->InsertCellPoint(pt1);
        vtkIdType pt2 = outputCells->InsertNextPoint(cellpoint2);
        cellArray->InsertCellPoint(pt2);
        vtkIdType pt3 = outputCells->InsertNextPoint(cellpoint3);
        cellArray->InsertCellPoint(pt3);
        vtkIdType pt4 = outputCells->InsertNextPoint(cellpoint4);
        cellArray->InsertCellPoint(pt4);

        colorsPts->InsertNextTypedTuple(tempColor);
        //colorsCells->InsertNextTypedTuple(tempColor);
        colorsCells->SetTypedTuple(index1,tempColor);

        for ( int k = 1; k < vect.size(); k++)
        {
            int index2 = static_cast<int>(vect[k].first);

            if ( index1 == index2 )
            {
                double point1[3] = { 0. };
                cellPoints->GetTuple(idvec, point1);
                //std::cout<<"index " << index2<< " "<< point1[0] << " "<< point1[1]<<std::endl;
                outputPoints1->InsertNextPoint(point1);
                colorsPts->InsertNextTypedTuple(tempColor);
                delIter += 1;
                idvec += 1;
            }
        }

        vect.erase( vect.begin(), vect.begin() + delIter + 1 );

        //std::cout <<"Size " << vect.size() << std::endl;
    }

    std::cout<<"output points size " << outputPoints1->GetNumberOfPoints() <<std::endl;

    output->SetPoints(outputPoints1);
    output->GetPointData()->SetScalars(colorsPts);

    //output1->SetPoints(outputCells);

    this->seedGrid->GetCellData()->AddArray(colorsCells);
    output1->ShallowCopy(this->seedGrid);

    //output->SetPoints(cell);
    //output->SetLines(lines);


    return 1;
}

void vtkColorMap::PrintSelf(ostream &os, vtkIndent indent)
{
}


void vtkColorMap:: cellLocator( int cellId, int* dims, int& x, int& y)
{
    y = cellId / (dims[0] - 1);
    x = cellId % (dims[0] - 1);

    return;
}
