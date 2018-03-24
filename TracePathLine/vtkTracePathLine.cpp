#include "vtkTracePathLine.h"

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

#include <vtkAssignAttribute.h>
#include <vtkCellData.h>

#include <random>
#include <vector>
#include <limits>
#include <iomanip>

using namespace std;

//int cellLocator( double Location[2], double spacing[3], int* resolution );
int getIndex(int x, int y, int* resolution)
{
    return ( y  * resolution[1]   + x );
}

vtkStandardNewMacro(vtkTracePathLine);

//-----------------------------------------------------------------------------
vtkTracePathLine::vtkTracePathLine()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                 vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
vtkTracePathLine::~vtkTracePathLine()
{

}

//----------------------------------------------------------------------------
int vtkTracePathLine::FillInputPortInformation( int port, vtkInformation* info )
{

    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
        return 1;
    }
    return 0;
}


//----------------------------------------------------------------------------
int vtkTracePathLine::FillOutputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 )
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );
        return 1;
    }

    return 0;
}

//----------------------------------------------------------------------------
int vtkTracePathLine::RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector)
{

    vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkInformation> inInfo = inputVector[0]->GetInformationObject(0);

    int resolution[6];
    resolution[0] = 0;
    resolution[1] = static_cast<int>(this->cellsNumber[0]);
    resolution[2] = 0;
    resolution[3] = static_cast<int>(this->cellsNumber[1]);
    resolution[4] = 0;
    resolution[5] = static_cast<int>(this->cellsNumber[2]);

    // Request numInTimes many time steps, starting from inTimes. inTimes is a pointer to the first time step (inTimes+3) would be a pointer to the fourth time step.

    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolution, 6);


    return 1;
}


//----------------------------------------------------------------------------
int vtkTracePathLine::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
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


    this->spacing[0] = this->bounds[0] / static_cast<double>(resolution[1]);
    this->spacing[1] = this->bounds[1] / static_cast<double>(resolution[3]);
    this->spacing[2] = this->bounds[2] / static_cast<double>(resolution[5]);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolution,6);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), resolution,6);
    outInfo->Set(vtkDataObject::SPACING(), this->spacing,3);
    outInfo->Set(vtkDataObject::ORIGIN(),this->origin,3);


    return 1;

}


//----------------------------------------------------------------------------
int vtkTracePathLine::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

    //// Get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkImageData* output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    this->seedGrid = vtkSmartPointer<vtkImageData>::New();

    int resolution[6];
    resolution[0] = 0;
    resolution[1] = static_cast<int>(this->cellsNumber[0]);
    resolution[2] = 0;
    resolution[3] = static_cast<int>(this->cellsNumber[1]);
    resolution[4] = 0;
    resolution[5] = static_cast<int>(this->cellsNumber[2]);

    int dimension[3] = { resolution[1] + 1, resolution[3] + 1, resolution[5] + 1  };

    ///Set Size and extent of the seedng grid
    this->seedGrid->SetDimensions(dimension);
    this->seedGrid->SetOrigin(this->origin);

    ///Calculate and set spacing of the seeding grid
    this->spacing[0] = this->bounds[0] / (static_cast<double>(resolution[1]));
    this->spacing[1] = this->bounds[1] / (static_cast<double>(resolution[3]));
    this->spacing[2] = this->bounds[2] / (static_cast<double>(resolution[5]));
    this->seedGrid->SetSpacing(this->spacing);

    ////Getting the point values from input data
    vtkSmartPointer<vtkDoubleArray> startPoints = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("StartPoints"));
    startPoints->SetNumberOfComponents(2);
    vtkSmartPointer<vtkDoubleArray> endPoints = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetAbstractArray("EndPoints"));
    endPoints->SetNumberOfComponents(2);

    double pointsPerCell = startPoints->GetNumberOfTuples() / ( resolution[1] * resolution[3] * 2 );

    std::cout<< " Number of Points per cell " << pointsPerCell << std::endl;

    ////Save the cell indices of both start and end locations
    vtkSmartPointer<vtkIntArray> LocationIndex = vtkSmartPointer<vtkIntArray>::New();
    LocationIndex->SetNumberOfComponents(2);


    //// allocate indices to all the startpoint and endpoint based on their local celllocation

    for (int i = 0; i < startPoints->GetNumberOfTuples() / 2; i++)              ////need to adjust for 2D. divided by 2
    {
        double startLocation[2], endLocation[2] = {0.};

        vtkIdType id = vtkIdType(i);
        startPoints->GetTuple( i, startLocation );
        endPoints->GetTuple( i, endLocation );

        int indexStart = cellLocator( startLocation, this->spacing, resolution );
        int indexEnd = cellLocator( endLocation, this->spacing, resolution );

        LocationIndex->InsertTuple2( id, indexStart, indexEnd );


    }

    std::cout<< "Location array size " << LocationIndex->GetNumberOfTuples() << std::endl;

    ////Need to sort the arrays, remove duplicates in order to calculate the number of cells pathlines ended at

    ////stores the respective cell number
    std::vector<vector<int>> TraceCells;
    //std::vector<int> pathlinevec;
    int numberOfPoints = 0;

    ///store the number of cells pathlines end from particular cell
    vtkSmartPointer<vtkDoubleArray> numberofCellsPathLinesEnd = vtkSmartPointer<vtkDoubleArray>::New();
    numberofCellsPathLinesEnd->SetNumberOfComponents(1);
    numberofCellsPathLinesEnd->SetNumberOfTuples( resolution[1] * resolution[3] );
    numberofCellsPathLinesEnd->SetName("Number Of Cells");


    for (int k = 0; k < resolution[1] * resolution[3]; k++)      ////number of cells
    {
        ////stores the cell's number pathlines traverse to
        std::vector<int> indexCells;
        vtkIdType idEnd = vtkIdType(k);

        for (int j = 0; j < LocationIndex->GetNumberOfTuples(); j++) {
            vtkIdType idLoc = vtkIdType(j);
            double indices[2] = {0.};
            LocationIndex->GetTuple(idLoc, indices);

            if (k == indices[0]) {
                indexCells.push_back(indices[1]);
            }
        }

        numberOfPoints += indexCells.size();

        std::sort(indexCells.begin(), indexCells.end());

        vector<int>::iterator ip;
        ip = std::unique(indexCells.begin(), indexCells.end());
        indexCells.resize(std::distance(indexCells.begin(), ip));

        /*
        if (indexCells.size() > pointsPerCell) {
            std::cout << "From cell " << k << std::endl;
            std::cout << "Pathlines end up in cells " << std::endl;
            for (int l = 0; l < indexCells.size(); l++) {
                std::cout << "  " << indexCells.at(l) << std::endl;
            }
        }
        std::cout << "________________________________________________________________________" << std::endl;
        */

        numberofCellsPathLinesEnd->InsertTuple1(idEnd, indexCells.size());
        //pathlinevec.push_back(indexCells.size());
        TraceCells.push_back( indexCells );
    }

    std::cout<< " Total number of points afterwards " << numberOfPoints << std::endl;

    /*
    auto old_count = pathlinevec.size();
    pathlinevec.resize(2 * old_count);
    std::copy_n(pathlinevec.begin(), old_count, pathlinevec.begin() + old_count);

    for (int pt = 0; pt < pathlinevec.size(); pt++)
    {
        vtkIdType idPt = vtkIdType(pt);
        int size = pathlinevec.at(pt);
        std::cout << pt <<" " << size << std::endl;
        numberofCellsPathLinesEnd->InsertTuple1( idPt, size );
    }
*/

    //std::cout <<"Number Of Tuples " << numberofCellsPathLinesEnd->GetNumberOfTuples()<<std::endl;

    this->seedGrid->GetCellData()->AddArray(numberofCellsPathLinesEnd);

    output->ShallowCopy(this->seedGrid);

    return 1;
}



/////check for the more points

void vtkTracePathLine::PrintSelf(ostream &os, vtkIndent indent)
{
}

int vtkTracePathLine:: cellLocator( double Location[2], double spacing[3], int* resolution )
{

    int xIndex, yIndex = 0;


    ///points which are on or more than x-boundary
    if ( Location[0] >= this->bounds[0] && Location[1] < this->bounds[1] )
    {
        xIndex = static_cast<int>( Location[0] / spacing[0] ) - 1;
        yIndex = static_cast<int>( Location[1] / spacing[1] );
    }

    ///points which are on or more than y-boundary
    else if ( Location[0] < this->bounds[0] && Location[1] >= this->bounds[1] )
    {
        xIndex = static_cast<int>( Location[0] / spacing[0] );
        yIndex = static_cast<int>( Location[1] / spacing[1] ) - 1;
    }

    ////both x and y boundaries
    else if ( Location[0] >= this->bounds[0] && Location[1] >= this->bounds[1] )
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