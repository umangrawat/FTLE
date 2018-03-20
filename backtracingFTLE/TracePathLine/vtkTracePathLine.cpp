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
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }

    return 0;
}

//----------------------------------------------------------------------------
int vtkTracePathLine::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    /*

    int resolutionSeedGrid[6];
    resolutionSeedGrid[0] = 0;
    resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0])-1;
    resolutionSeedGrid[2] = 0;
    resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1])-1;
    resolutionSeedGrid[4] = 0;
    resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2])-1;


    this->spacingSeedGrid[0] = this->boundsSeedGrid[0] / static_cast<double>(resolutionSeedGrid[1]);
    this->spacingSeedGrid[1] = this->boundsSeedGrid[1] / static_cast<double>(resolutionSeedGrid[3]);
    this->spacingSeedGrid[2] = this->boundsSeedGrid[2] / static_cast<double>(resolutionSeedGrid[5]);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid,6);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), resolutionSeedGrid,6);
    outInfo->Set(vtkDataObject::SPACING(), this->spacingSeedGrid,3);
    outInfo->Set(vtkDataObject::ORIGIN(),this->originSeedGrid,3);
*/

    return 1;

}



int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}