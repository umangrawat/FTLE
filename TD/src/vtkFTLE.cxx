#include "vtkFTLE.h"
#include <vtkCell.h>

//#include "stdfunc.h"
#include <iomanip>

#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPixel.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangle.h"
#include <vtkProbeFilter.h>
#include <vtkDoubleArray.h>
#include <vtkCubeSource.h>
#include <vtkPlaneSource.h>
#include <vtkStreamLine.h>
#include <vtkStreamTracer.h>
#include <vtkResampleToImage.h>
#include <vtkIndent.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCallbackCommand.h>
#include "Integrator.h"
#include "chrono"
#include "cudaIntegrator.h"
#include "vtkPVExtractVOI.h"
#include "vtkMultiBlockDataSet.h"
#include <vtkDataObject.h>
#include <vtkCompositeDataSet.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataPipeline.h>
#include <vtkAppendPolyData.h>
#include <vtkImageAppend.h>

using namespace std;
using namespace std::chrono;



/// Get the index of a 3 dim imagedata (uniform grid)

int getIndex(int z,int y, int x, int* dataDims)
{
    return (z*dataDims[1]*dataDims[0]+ y*dataDims[0]+x);
}

/// fy refers to the current function value
/// fx and fz to the i-1 and i+1 function value
double centralDiff(double fx, double fz, double dist)
{
    return (fz-fx)/(2*dist);
}

double forwardDiff(double fy, double fz, double dist)
{
    return (fz-fy)/dist;
}
/// Changed because of small bug in finiteDiff, was to lazy changing 12 function calls
double backwardDiff(double fx, double fy, double dist)
{
    return (fx-fy)/dist;

}


void progressFunction ( vtkObject* caller,
                        long unsigned int vtkNotUsed(eventId),
                        void* vtkNotUsed(clientData),
                        void* vtkNotUsed(callData) )
{

    vtkStreamTracer* streamTracer = static_cast<vtkStreamTracer*>(caller);

    std::cout << "Progress: " << streamTracer->GetProgress()  * 100.0 <<"%" << '\r';

}


vtkStandardNewMacro(vtkFTLE)

//-----------------------------------------------------------------------------

vtkFTLE::vtkFTLE()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(3);                                                                                     //// changed to 3
    ////changed back to 0,0,0 from 0,0
    // by default process active point vectors
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                 vtkDataSetAttributes::VECTORS);

}

//-----------------------------------------------------------------------------
vtkFTLE::~vtkFTLE()
{
}

//----------------------------------------------------------------------------
int vtkFTLE::FillInputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 )
    {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");   ////changed frim vtkdataset
        return 1;
    }
    return 0;
}

//----------------------------------------------------------------------------
int vtkFTLE::FillOutputPortInformation( int port, vtkInformation* info )
{
    /*
    if ( port == 0 )
    {
        //TODO remove if not needed
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );
        return 1;
    }
    if( port == 1){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );    ////changed it from vtkdataobject
        return 1;
    }*/

    if ( port == 0 )
    {
        //TODO remove if not needed
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );
        return 1;
    }
    if( port == 1){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );    ////changed it from vtkdataobject
        return 1;
    }

    if( port == 2){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );    ////changed it from vtkdataobject
        return 1;
    }

    return 0;
}

//----------------------------------------------------------------------------
int vtkFTLE::RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector)
{

    vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkInformation> outInfo1 = outputVector->GetInformationObject(1);           ////
    vtkSmartPointer<vtkInformation> outInfo2 = outputVector->GetInformationObject(2);           ////added this
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) ||
           outInfo1->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
        // get the update times
        //double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());           ///commented out
        double upTime = outInfo1->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

        double *inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        int numInTimes = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), inTimes, numInTimes);


        int resolutionSeedGrid[6];
        /*
        resolutionSeedGrid[0] = 0;
        resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0]) - 1;
        resolutionSeedGrid[2] = 0;
        resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1]) - 1;
        resolutionSeedGrid[4] = 0;
        resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2]) - 1;
        */
        // Request numInTimes many time steps, starting from inTimes. inTimes is a pointer to the first time step (inTimes+3) would be a pointer to the fourth time step.

        auto tmp = inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
        inInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), tmp, 6);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), tmp, 6);

        std::cout << "request update extent " << tmp[1] << " " << tmp[3] << " " << tmp[5] << std::endl;
        //outInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid, 6);

    }

    return 1;
}


//----------------------------------------------------------------------------
int vtkFTLE::RequestInformation(vtkInformation *vtkNotUsed(request),
                                vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkInformation *outInfo2 = outputVector->GetInformationObject(2);                       ///// new addition
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);



    int resolutionSeedGrid[6];
    resolutionSeedGrid[0] = 0;
    resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0])-1;
    resolutionSeedGrid[2] = 0;
    resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1])-1;
    resolutionSeedGrid[4] = 0;
    resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2])-1;



    this->spacingSeedGrid[0] = this->boundsSeedGrid[0] /static_cast<double>(resolutionSeedGrid[1]);
    this->spacingSeedGrid[1] = this->boundsSeedGrid[1] /static_cast<double>(resolutionSeedGrid[3]);
    this->spacingSeedGrid[2] = this->boundsSeedGrid[2] /static_cast<double>(resolutionSeedGrid[5]);

    outInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid,6);     ////changed to 1
    outInfo1->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), resolutionSeedGrid,6);
    outInfo1->Set(vtkDataObject::SPACING(), this->spacingSeedGrid,3);
    outInfo1->Set(vtkDataObject::ORIGIN(),this->originSeedGrid,3);

    std::cout << "request info " << resolutionSeedGrid[1] << " " <<resolutionSeedGrid[3] << " "<<resolutionSeedGrid[5]<<std::endl;



    /*
    //outInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid,6);
    //outInfo1->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), resolutionSeedGrid,6);
    //outInfo1->Set(vtkDataObject::SPACING(), this->spacingSeedGrid,3);
    //outInfo1->Set(vtkDataObject::ORIGIN(),this->originSeedGrid,3);
*/
    return 1;


}


//----------------------------------------------------------------------------
int vtkFTLE::RequestDataObject( vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    if (this->GetNumberOfInputPorts() == 0 || this->GetNumberOfOutputPorts() == 0)
    {
        return 1;
    }

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    if (!inInfo)
    {
        return 0;
    }
    vtkDataObject *input = inInfo->Get(vtkDataObject::DATA_OBJECT());

    if (input)
    {
        // for each output
        for(int i=0; i < this->GetNumberOfOutputPorts(); ++i)
        {
            vtkInformation* info = outputVector->GetInformationObject(i);
            vtkDataSet *output = vtkDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
            vtkDataSet *output1 = vtkDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
            vtkPolyData *output2 = vtkPolyData::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));    ////changed here

            if (!output || !output->IsA(input->GetClassName()))
            {
                vtkSmartPointer<vtkDataObject> newOutput = input->NewInstance();
                info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
                newOutput->Delete();
            }

        }
        return 1;
    }

    return 0;
}



//----------------------------------------------------------------------------
int vtkFTLE::RequestData(vtkInformation *vtkNotUsed(request),
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    vtkSmartPointer <vtkMultiBlockDataSet> inData = vtkMultiBlockDataSet::New();

    if (!inInfo->Has(vtkDataObject::DATA_OBJECT()))
    {
        cout << "Input has no data object. No calculation done." << endl;
        return 1;
    }
    else
    {
        inData = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    }

    int numTimeSteps = inData->GetNumberOfBlocks();

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkInformation *outInfo2 = outputVector->GetInformationObject(2);                                                       ////new

    vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet *output1 = vtkDataSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output2 = vtkPolyData::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));

    double upTime = outInfo1->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());                                    ////changed here
    //double* times = outputVector->GetInformationObject(0)->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());        //// different timesteps
    //int numInTimes = outputVector->GetInformationObject(0)->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());    //// number of timesteps
    double* times = outputVector->GetInformationObject(1)->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());        //// different timesteps
    int numInTimes = outputVector->GetInformationObject(1)->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

    this->seedGrid = vtkSmartPointer<vtkImageData>::New();

    int resolutionSeedGrid[6];
        resolutionSeedGrid[0] = 0;
        resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0]) - 1;
        resolutionSeedGrid[2] = 0;
        resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1]) - 1;
        resolutionSeedGrid[4] = 0;
        resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2]) - 1;


        ///Set Size and extent of the seedng grid
        this->seedGrid->SetExtent(resolutionSeedGrid);
        this->seedGrid->SetOrigin(this->originSeedGrid);


        ///Calculate and set spacing of the seeding grid
        this->spacingSeedGrid[0] = this->boundsSeedGrid[0] / (static_cast<double>(resolutionSeedGrid[1])) ;
        this->spacingSeedGrid[1] = this->boundsSeedGrid[1] / (static_cast<double>(resolutionSeedGrid[3]));
        this->spacingSeedGrid[2] = this->boundsSeedGrid[2] / (static_cast<double>(resolutionSeedGrid[5]));
        this->seedGrid->SetSpacing(this->spacingSeedGrid);


    #if VTK_MAJOR_VERSION <= 5
        this->seedGrid->SetNumberOfScalarComponents(1);
        this->seedGrid->SetScalarTypeToDouble();
        this->seedGrid->AllocateScalars();
    #else
        this->seedGrid->AllocateScalars(VTK_DOUBLE, 1);
        //this->nextseedGrid->AllocateScalars(VTK_DOUBLE, 1);
    #endif


        this->printSeedGridProperties();
        vtkSmartPointer<vtkDoubleArray> flowMap = vtkSmartPointer<vtkDoubleArray>::New();
        flowMap->SetNumberOfComponents(2);
        flowMap->SetName("FlowMap");
        flowMap->SetNumberOfTuples(seedGrid->GetNumberOfPoints());
        vtkSmartPointer<vtkDoubleArray> seeded = vtkSmartPointer<vtkDoubleArray>::New();
        seeded->SetNumberOfComponents(1);
        seeded->SetName("seeded");
        seeded->SetNumberOfTuples(seedGrid->GetNumberOfPoints());

        cout << "Fill seed grid with zeros" << endl;
        int *dims = this->seedGrid->GetDimensions();
        double *dummyPtr;
        for (int z = 0; z < 1; z++) {
            for (int y = 0; y < dims[1] - 1; y++) {
                for (int x = 0; x < dims[0] - 1; x++) {

                    auto *ptr = static_cast<double *>(this->seedGrid->GetScalarPointer(x, y, z));
                    ptr[0] = 0;
                    ptr[1] = 0;
                    int index = getIndex(z, y, x, dims);
                    vtkIdType idx = vtkIdType(index);
                    double temp = 0.0;
                    dummyPtr = &temp;
                    flowMap->SetTuple(idx, dummyPtr);
                }
            }
        }


    /// used for in between every 2 time steps
    vtkSmartPointer<vtkPolyData> streamLines = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCallbackCommand> callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    streamLines = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkDoubleArray> ftleField;

    vtkSmartPointer<vtkAppendPolyData> StreamLinesappendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    //// final pathlines
    vtkSmartPointer<vtkPolyData> pathLines = vtkSmartPointer<vtkPolyData>::New();


    ////iterators
    int InitialTime, FinalTime;
    std::vector<vector<double>> previousPts;
    std::vector<vector<double>> StartPointsofPathline;
    std::vector<vector<double>> EndPointsofPathline;


    ///Forward FTLE

    if (this->tauFinal > 0) {

        if (this->tauFinal + this->tauInitial > times[numInTimes - 1] || this->tauInitial < times[0]) {
            std::cerr << "Incorrect time values" << std::endl;
            return false;
        } else {
            for (int k = 0; k < numInTimes - 1; k++) {
                if (this->tauInitial >= times[k] && this->tauInitial < times[k + 1]) {
                    InitialTime = k;
                }
            }

            for (int l = 0; l < numInTimes - 1; l++) {
                if (this->tauFinal + this->tauInitial > times[l] && this->tauFinal + this->tauInitial <= times[l + 1]) {
                    FinalTime = l + 1;
                }
            }

        }

        std::cout << "time iterators" << InitialTime << " " <<FinalTime << std::endl;

        for ( int i = InitialTime; i < FinalTime;  i++ )
        {

            std::cout << "_________________________________________________________________" << std::endl;

            ////storing the consecutive time step data for computations
            vtkSmartPointer<vtkImageData> data = vtkImageData::SafeDownCast(inData->GetBlock(i));
            vtkSmartPointer<vtkImageData> nextdata = vtkImageData::SafeDownCast(inData->GetBlock(i + 1));


            vtkIndent indentMyAss(4);
            this->inputGrid = data;
            this->inputGridSpacing = data->GetSpacing();
            this->inputGridDimensions = data->GetDimensions();

            this->nextinputGrid = nextdata;
            this->nextinputGridSpacing = nextdata->GetSpacing();
            this->nextinputGridDimensions = nextdata->GetDimensions();


            callBack->SetCallback(progressFunction);

            cout << "Integrate stream lines" << endl;

            ///Needed to tell the stream tracer which vectors to use
            this->inputGrid->GetPointData()->SetActiveVectors(this->inputGrid->GetPointData()->GetArrayName(0));
            this->nextinputGrid->GetPointData()->SetActiveVectors(this->nextinputGrid->GetPointData()->GetArrayName(0));

            if (this->useGPU) {

                int temp = 1;

                cout << "Start Integration" << endl;


                Integrator *customIntegrator = new Integrator(temp, this->integrationStepSize, fabs(this->tauInitial),
                                                              fabs(this->tauFinal), LIMIT_DOUBLE,
                                                              this->inputGrid, this->seedGrid, this->nextinputGrid, i,
                                                              InitialTime, FinalTime, times, previousPts);
                cout << "init integrator" << endl;
                this->inputGrid->GetOrigin(customIntegrator->origin);
                customIntegrator->setOriginSource(this->originSeedGrid);
                customIntegrator->calcStreamlines = this->showStreamLines;
                customIntegrator->maxNumberOfSteps = this->maxNumberOfSteps;

                ////storing the values of enpoints of every iteration in a vector array, to be used for next iteration
                previousPts = customIntegrator->PathlineEndPoints;

                ///storing the starting point for each pathline
                if ( i == InitialTime)
                {
                    StartPointsofPathline = customIntegrator->StartPoint;
                }

                ////freeing up the memory
                customIntegrator->PathlineEndPoints.erase(customIntegrator->PathlineEndPoints.begin(),
                                                          customIntegrator->PathlineEndPoints.end());
                customIntegrator->NewEndPoints.erase(customIntegrator->NewEndPoints.begin(),
                                                     customIntegrator->NewEndPoints.end());


                customIntegrator->integrateRK4GPUWrapper();
                this->seedGrid = customIntegrator->source;
                delete customIntegrator;

            }
            else
            {

                int temp = 1;

                cout << "Start Integration" << endl;

                Integrator *customIntegrator = new Integrator(temp, this->integrationStepSize, fabs(this->tauInitial),
                                                              fabs(this->tauFinal), LIMIT_DOUBLE,
                                                              this->inputGrid, this->seedGrid, this->nextinputGrid, i,
                                                              InitialTime, FinalTime, times, previousPts);
                cout << "init integrator" << endl;
                this->inputGrid->GetOrigin(customIntegrator->origin);
                customIntegrator->setOriginSource(this->originSeedGrid);
                customIntegrator->calcStreamlines = this->showStreamLines;
                customIntegrator->maxNumberOfSteps = this->maxNumberOfSteps;

                streamLines = customIntegrator->integrateRK4();

                ////storing the values of enpoints of every iteration in a vector array, to be used for next iteration
                previousPts = customIntegrator->PathlineEndPoints;

                ///storing the starting point for each pathline
                if ( i == InitialTime)
                {
                    StartPointsofPathline = customIntegrator->StartPoint;
                }

                std::cout << "previous pts size " << previousPts.size() << std::endl;
                std::cout << "Start pts size " << StartPointsofPathline.size() << std::endl;

                ////freeing up the memory
                customIntegrator->PathlineEndPoints.erase(customIntegrator->PathlineEndPoints.begin(),
                                                          customIntegrator->PathlineEndPoints.end());
                customIntegrator->NewEndPoints.erase(customIntegrator->NewEndPoints.begin(),
                                                     customIntegrator->NewEndPoints.end());

                delete customIntegrator;

            }

            cout << "Number of stream line cells " << streamLines->GetNumberOfCells() << endl;
            cout << "Number of stream lines " << streamLines->GetNumberOfLines() << endl;

            StreamLinesappendFilter->AddInputData(streamLines);
            StreamLinesappendFilter->Update();

            std::cout<< "Current Progress " << int((double(i + 1 - InitialTime)/double(FinalTime - InitialTime))* 100.0)<< std::endl;
        }
    }

    ////Backward FTLE
    else {

        if ( this->tauInitial > times[numInTimes - 1] || this->tauInitial + this->tauFinal < times[0]) {
            std::cerr << "Incorrect time values" << std::endl;
            return false;
        } else {
            for (int k = 0; k < numInTimes - 1; k++) {
                if (this->tauInitial > times[k] && this->tauInitial <= times[k + 1]) {
                    InitialTime = k + 1;
                }

            }

            for (int l = 0; l < numInTimes - 1; l++) {
                if ( this->tauInitial + this->tauFinal >= times[l] && this->tauInitial + this->tauFinal < times[l + 1]) {
                    FinalTime = l;
                }
            }

        }


        std::cout << "time iterators" << InitialTime << " " <<FinalTime << std::endl;

        for ( int i = InitialTime; i > FinalTime;  i-- )
        {

            std::cout << "_________________________________________________________________" << std::endl;

            ////storing the consecutive time step data for computations
            vtkSmartPointer<vtkImageData> data = vtkImageData::SafeDownCast(inData->GetBlock(i));
            vtkSmartPointer<vtkImageData> nextdata = vtkImageData::SafeDownCast(inData->GetBlock(i - 1));


            vtkIndent indentMyAss(4);
            this->inputGrid = data;
            this->inputGridSpacing = data->GetSpacing();
            this->inputGridDimensions = data->GetDimensions();

            this->nextinputGrid = nextdata;
            this->nextinputGridSpacing = nextdata->GetSpacing();
            this->nextinputGridDimensions = nextdata->GetDimensions();


            callBack->SetCallback(progressFunction);

            cout << "Integrate stream lines" << endl;

            ///Needed to tell the stream tracer which vectors to use
            this->inputGrid->GetPointData()->SetActiveVectors(this->inputGrid->GetPointData()->GetArrayName(0));
            this->nextinputGrid->GetPointData()->SetActiveVectors(this->nextinputGrid->GetPointData()->GetArrayName(0));

            if (this->useGPU) {

                int temp = -1;                                      ////Backward FTLE

                cout << "Start Integration" << endl;


                Integrator *customIntegrator = new Integrator(temp, this->integrationStepSize, this->tauInitial,
                                                              this->tauFinal, LIMIT_DOUBLE,
                                                              this->inputGrid, this->seedGrid, this->nextinputGrid, i,
                                                              InitialTime, FinalTime, times, previousPts);
                cout << "init integrator" << endl;
                this->inputGrid->GetOrigin(customIntegrator->origin);
                customIntegrator->setOriginSource(this->originSeedGrid);
                customIntegrator->calcStreamlines = this->showStreamLines;
                customIntegrator->maxNumberOfSteps = this->maxNumberOfSteps;

                ////storing the values of enpoints of every iteration in a vector array, to be used for next iteration
                previousPts = customIntegrator->PathlineEndPoints;

                ///storing the starting point for each pathline
                if ( i == InitialTime)
                {
                    StartPointsofPathline = customIntegrator->StartPoint;
                }

                std::cout << "Start pts size " << StartPointsofPathline.size() << std::endl;

                ////freeing up the memory
                customIntegrator->PathlineEndPoints.erase(customIntegrator->PathlineEndPoints.begin(),
                                                          customIntegrator->PathlineEndPoints.end());
                customIntegrator->NewEndPoints.erase(customIntegrator->NewEndPoints.begin(),
                                                     customIntegrator->NewEndPoints.end());


                customIntegrator->integrateRK4GPUWrapper();
                this->seedGrid = customIntegrator->source;
                delete customIntegrator;

            }
            else
            {

                int temp = -1;

                cout << "Start Integration" << endl;

                Integrator *customIntegrator = new Integrator(temp, this->integrationStepSize, this->tauInitial, this->tauFinal, LIMIT_DOUBLE,
                                                              this->inputGrid, this->seedGrid, this->nextinputGrid, i,
                                                              InitialTime, FinalTime, times, previousPts);
                cout << "init integrator" << endl;
                this->inputGrid->GetOrigin(customIntegrator->origin);
                customIntegrator->setOriginSource(this->originSeedGrid);
                customIntegrator->calcStreamlines = this->showStreamLines;
                customIntegrator->maxNumberOfSteps = this->maxNumberOfSteps;

                streamLines = customIntegrator->integrateRK4();

                ////storing the values of enpoints of every iteration in a vector array, to be used for next iteration
                previousPts = customIntegrator->PathlineEndPoints;

                ///storing the starting point for each pathline
                if ( i == InitialTime)
                {
                    StartPointsofPathline = customIntegrator->StartPoint;
                }

                std::cout << "previous pts size " << previousPts.size() << std::endl;
                std::cout << "Start pts size " << StartPointsofPathline.size() << std::endl;

                ////freeing up the memory
                customIntegrator->PathlineEndPoints.erase(customIntegrator->PathlineEndPoints.begin(),
                                                          customIntegrator->PathlineEndPoints.end());
                customIntegrator->NewEndPoints.erase(customIntegrator->NewEndPoints.begin(),
                                                     customIntegrator->NewEndPoints.end());

                delete customIntegrator;

            }

            cout << "Number of stream line cells " << streamLines->GetNumberOfCells() << endl;
            cout << "Number of stream lines " << streamLines->GetNumberOfLines() << endl;

            StreamLinesappendFilter->AddInputData(streamLines);
            StreamLinesappendFilter->Update();

            std::cout<< "Current Progress " << int((double( InitialTime - i + 1)/double(InitialTime - FinalTime))* 100.0)<< std::endl;
        }
    }


    pathLines->ShallowCopy(StreamLinesappendFilter->GetOutput());
    std::cout << "pathlines " << pathLines->GetNumberOfCells() << std::endl;


    if (!this->checkBounds(this->boundsSeedGrid))
    {
        vtkErrorMacro("The user input bounds extent the data bounds.");
        //output->DeepCopy(seedGrid);
        //return 1;
    }

    ////set the output in the dataset, to be used later to append the previous results
        //newpoints->SetBlock(i, streamLines);





    if (0 == streamLines->GetNumberOfCells() && !useGPU) {
            //output->DeepCopy(this->seedGrid);
            return 1;
    }

    if (this->useGPU) {

    } else {


/*
        if (this->seedGrid->GetNumberOfPoints() == streamLines->GetNumberOfCells()) {      ////
            double point[2];
            double seed[1];

            //int numberStreamlines = streamLines->GetNumberOfCells();
            //int numberPathlines = pathLines->GetNumberOfCells();

            //std::cout << numberPathlines << " " << numberStreamlines << std::endl;

            //for (int indexLine =  numberPathlines - numberStreamlines; indexLine < numberPathlines; ++indexLine) {        ////
            for (int indexLine =  0; indexLine < streamLines->GetNumberOfCells(); ++indexLine) {        ////
                int lastPointIndex = streamLines->GetCell(indexLine)->GetNumberOfPoints() - 1;

                //std::cout <<"lastPointIndex " << lastPointIndex << std::endl;
                if (lastPointIndex >= 0) {
                    streamLines->GetCell(indexLine)->GetPoints()->GetPoint(lastPointIndex, point);         ////
                    flowMap->SetTuple(indexLine, point);
                    seed[0] = 1;
                    seeded->SetTuple(indexLine, seed);
                } else {
                    std::cerr << "stream Line not seeded" << std::endl;
                }
            }

        }*/


        if (this->seedGrid->GetNumberOfPoints() == previousPts.size()) {      ////
            double point[2];
            std::vector<double> pathlineEndpoint;
            double seed[1];

            for (int indexLine = 0; indexLine < previousPts.size(); ++indexLine) {        ////

                //int lastPointIndex = streamLines->GetCell(indexLine)->GetNumberOfPoints() - 1;
                //std::cout <<"lastPointIndex " << lastPointIndex << std::endl;
                //if (lastPointIndex >= 0) {
                pathlineEndpoint = previousPts.at(indexLine);
                    //streamLines->GetCell(indexLine)->GetPoints()->GetPoint(lastPointIndex, point);         ////
                point[0] = pathlineEndpoint[0];
                point[1] = pathlineEndpoint[1];
                flowMap->SetTuple(indexLine, point);
                seed[0] = 1;
                seeded->SetTuple(indexLine, seed);

                //} else {
                    //std::cerr << "stream Line not seeded" << std::endl;
                //}
            }

        }

    }

    vtkSmartPointer<vtkDoubleArray> startPoints = vtkSmartPointer<vtkDoubleArray>::New();
    startPoints->SetNumberOfComponents(2);
    startPoints->SetNumberOfTuples(StartPointsofPathline.size());
    startPoints->SetName("StartPoints");

    vtkSmartPointer<vtkDoubleArray> endPoints = vtkSmartPointer<vtkDoubleArray>::New();
    endPoints->SetNumberOfComponents(2);
    endPoints->SetNumberOfTuples(previousPts.size());
    endPoints->SetName("EndPoints");


    for (int pointIter = 0; pointIter < previousPts.size(); pointIter++)
    {
        std::vector<double> start;
        std::vector<double> end;
        start = StartPointsofPathline.at(pointIter);
        end = previousPts.at(pointIter);

        vtkIdType idPoint = vtkIdType(pointIter);

        startPoints->SetTuple2(idPoint, start[0], start[1]);
        endPoints->SetTuple2(idPoint, end[0], end[1]);
    }



    cout << "Get end positions of stream line integration " << endl;
    cout << "Number of seed grid points " << this->seedGrid->GetNumberOfPoints() << endl;
    cout << "Number of seed cells " << this->seedGrid->GetNumberOfCells() << endl;

    cout << "Copy to output" << endl;

    if (this->useGPU) {
        //ftleField->SetName("FTLE");
        //ftleField->Print(std::cout);
        //this->seedGrid->GetPointData()->AddArray(ftleField);
        vtkSmartPointer<vtkDoubleArray> ftle = this->computeDerivatives(seedGrid, flowMap->GetName());
        ftle->SetName("FTLECorrect");

        this->seedGrid->GetPointData()->AddArray(ftle);
        }
    else
    {
        this->seedGrid->GetPointData()->AddArray(flowMap);
        vtkSmartPointer<vtkDoubleArray> ftle = this->computeDerivatives(seedGrid, flowMap->GetName());
        ftle->SetName("FTLE");
        this->seedGrid->GetPointData()->AddArray(ftle);
        this->seedGrid->GetPointData()->AddArray(seeded);
        this->seedGrid->GetPointData()->AddArray(startPoints);
        this->seedGrid->GetPointData()->AddArray(endPoints);
    }



        ///Translate SeedGrid to correct origin
        double newOrigin[3];
        this->inputGrid->GetOrigin(newOrigin);
        newOrigin[0] += this->originSeedGrid[0];
        newOrigin[1] += this->originSeedGrid[1];
        newOrigin[2] += this->originSeedGrid[2];
        this->seedGrid->SetOrigin(newOrigin);

    vtkSmartPointer<vtkImageData> copy = vtkImageData::SafeDownCast(inData->GetBlock(0));

    output->ShallowCopy(copy);
    output1->ShallowCopy(this->seedGrid);
    output2->ShallowCopy(pathLines);


    return 1;
}

//----------------------------------------------------------------------------
void vtkFTLE::getData(vtkSmartPointer<vtkImageData> in_grid, double* &in_data, int numComponents, bool &allocated_in_data)
{
    switch (in_grid->GetScalarType())
    {
    case VTK_FLOAT:
    {
        ////changed 3D to 2D
        const size_t size = dim[0]*dim[1]*numComponents;
        std::vector<float> tmp_data(size);
        memcpy(tmp_data.data(), in_grid->GetScalarPointer(), size*sizeof(float));
        std::vector<double> tmp_double(tmp_data.begin(),tmp_data.end());
        in_data = new double[size];
        allocated_in_data = true;
        memcpy(in_data, tmp_double.data(), size*sizeof(double));
    }
    break;
    case VTK_DOUBLE:
    default:
        in_data = static_cast<double *>(in_grid->GetScalarPointer());
        break;
    }
}

void vtkFTLE::PrintSelf(ostream &os, vtkIndent indent)
{
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkDoubleArray> vtkFTLE::computeDerivatives(vtkSmartPointer<vtkImageData> data, char* arrayName)
{


    vtkSmartPointer<vtkDoubleArray> ftle = vtkSmartPointer<vtkDoubleArray>::New();

    ftle->SetNumberOfComponents(1);
    ftle->SetNumberOfTuples(data->GetPointData()->GetNumberOfTuples());
    cout<<"Try to get array by name"<<endl;
    vtkSmartPointer<vtkDataArray> field;
    if(true) {
        field = vtkDoubleArray::SafeDownCast(
                data->GetPointData()->GetAbstractArray(arrayName));
    }
    else {
        field = vtkFloatArray::SafeDownCast(
                data->GetPointData()->GetAbstractArray(arrayName));
    }
    field->GetNumberOfTuples();
    cout<<"Test returned array succesful"<<endl;
    int* dataDims = data->GetDimensions();
    double* spac = data->GetSpacing();

    ////changed 3D to 2D
    mat2 jacobi;
    vec2 du;
    vec2 dv;


    for (int z=0; z< dataDims[2]; z++)                                        ////changed to 1
    {
        for (int y=0; y<dataDims[1]; y++)
        {
            for (int x=0; x<dataDims[0]; x++)
            {
                ////changed (z,y,x) to (y,x)                                                                                               
                int index = getIndex(z,y,x,dataDims);
                vtkIdType id = vtkIdType(index);
                ////changed 3 to 2
                double t1[2], t2[2];


                /// y Component
                if(y==0)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y+1,x,dataDims);
                    field->GetTuple(id,t2);
                    du[1] = forwardDiff(t1[0],t2[0],spac[0]);
                    dv[1] = forwardDiff(t1[1],t2[1],spac[1]);
                    ////dw[1] = forwardDiff(t1[2],t2[2],spac[2]);
                }
                else if (y == dataDims[1]-1)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y-1,x,dataDims);
                    field->GetTuple(id,t2);
                    du[1] = backwardDiff(t1[0],t2[0],spac[0]);
                    dv[1] = backwardDiff(t1[1],t2[1],spac[1]);
                    ////dw[1] = backwardDiff(t1[2],t2[2],spac[2]);

                }
                else
                {
                    id = getIndex(z,y-1,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y+1,x,dataDims);
                    field->GetTuple(id,t2);
                    du[1] = centralDiff(t1[0],t2[0],spac[0]);
                    dv[1] = centralDiff(t1[1],t2[1],spac[1]);
                    ////dw[1] = centralDiff(t1[2],t2[2],spac[2]);
                }


                ///x Component
                if(x==0)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y,x+1,dataDims);
                    field->GetTuple(id,t2);
                    du[0] = forwardDiff(t1[0],t2[0],spac[0]);
                    dv[0] = forwardDiff(t1[1],t2[1],spac[1]);
                    ////dw[0] = forwardDiff(t1[2],t2[2],spac[2]);
                }
                else if (x == dataDims[0]-1)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y,x-1,dataDims);
                    field->GetTuple(id,t2);
                    du[0] = backwardDiff(t1[0],t2[0],spac[0]);
                    dv[0] = backwardDiff(t1[1],t2[1],spac[1]);
                    ////dw[0] = backwardDiff(t1[2],t2[2],spac[2]);

                }
                else
                {
                    id = getIndex(z,y,x-1,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y,x+1,dataDims);
                    field->GetTuple(id,t2);
                    du[0] = centralDiff(t1[0],t2[0],spac[0]);
                    dv[0] = centralDiff(t1[1],t2[1],spac[1]);
                    ////dw[0] = centralDiff(t1[2],t2[2],spac[2]);

                }

                ////changed 3D to 2D
                mat2setrows(jacobi, du,dv);
                mat2 jacobiT;
                mat2 cgTensor;
                mat2trp(jacobi,jacobiT);
                mat2mul(jacobiT,jacobi, cgTensor);
                double eMax = 0;
                vec2 eigenV;
                int realEigen = mat2eigenvalues(cgTensor, eigenV);
                ////changed max(max()) to max()
                eMax = max(eigenV[0],eigenV[1]);

                ////changed 3 to 2
                if(realEigen !=2 || fabs(eMax)< LIMIT_DOUBLE)
                {
                    //std::cerr<<"Eigenvalues of Cauchy Green Tensor are not real"<<endl;
                    //out<<realEigen<< " eigenvalues"<< eigenV[0]<< " " <<eigenV[1]<< " "<< eigenV[2]<<" "<<LIMIT_DOUBLE<<endl;
                    eMax = 0;
                    id = getIndex(z,y,x,dataDims);
                    /// Possible error source because eMax is only a scalar
                    ftle->SetTuple(id,&eMax);
                    continue;        
                }

                eMax = 1.0/ fabs((this->tauFinal)) * log(sqrt(eMax));
                id = getIndex(z,y,x,dataDims);
                /// Possible error source because eMax is only a scalar
                ftle->SetTuple(id,&eMax);

            }
        }
    }


    return ftle;

}


//----------------------------------------------------------------------------
void vtkFTLE::printSeedGridProperties()
{
    int* ext1 = this->seedGrid->GetExtent();
    double* space1 = this->seedGrid->GetSpacing();
    int* dimensions = this->seedGrid->GetDimensions();

    cout<<"_________________________________________________________________________"<<endl;
    //cout<<"Extent: "<<ext1[0]/100.<<" "<<ext1[1]/100.<<" "<<ext1[2]/100.<<" "<<ext1[3]/100.<<" "<<ext1[4]/100.<<" "<<ext1[5]/100.<<" "<<endl;
    cout<<"Extent: "<<ext1[0]<<" "<<ext1[1]*space1[0]<<" "<<ext1[2]<<" "<<ext1[3]*space1[1]<<" "<<ext1[4]<<" "<<ext1[5]*space1[2]<<" "<<endl;
    cout<<"Dimensions: "<<dimensions[0]<<" "<<dimensions[1]<<" "<<dimensions[2]<<" "<<endl;
    //cout<<this->dimSeedGrid[0]<<" "<<this->dimSeedGrid[1]<<" "<<this->dimSeedGrid[2]<<" "<<endl;
    cout<<"Spacing: "<<space1[0]<<" "<<space1[1]<<" "<<space1[2]<<" "<<endl;


}


///Check if the input dimensions exceed the data dimensions given to the filter
bool vtkFTLE::checkBounds(double* dimensions)
{
    int* inputDim = this->inputGrid->GetDimensions();
    double* spacing= this->inputGrid->GetSpacing();
    double bounds[3];
    bounds[0] = inputDim[0] / spacing[0];        ////changed from / to *
    bounds[1] = inputDim[1] / spacing[1];
    bounds[2] = inputDim[2] / spacing[2];

    if(this->boundsSeedGrid[0] > bounds[0] || this->boundsSeedGrid[1] > bounds[1] || this->boundsSeedGrid[2] > bounds[2])
    {
        return false;
    }

    return true;
}