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

using namespace std;
using namespace std::chrono;

//unsigned int i = 0;

vtkSmartPointer<vtkMultiBlockDataSet> newpoints = vtkSmartPointer<vtkMultiBlockDataSet>::New();

/// Get the index of a 3 dim imagedata (uniform grid)

int getIndex(int z,int y, int x, int* dataDims)
{
    return (z*dataDims[1]*dataDims[0]+ y*dataDims[0]+x);
}


/*
/// Get the index of a 2 dim imagedata (uniform grid)
int getIndex(int y, int x, int* dataDims)
{
    return (y*dataDims[0]+x);
}
*/


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
/*
vtkFTLE::vtkFTLE()
{
    this->SetNumberOfOutputPorts(2);

    // by default process active point vectors
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                 vtkDataSetAttributes::VECTORS);

}
*/

vtkFTLE::vtkFTLE()
{
    this->SetNumberOfInputPorts(1);      ////new addition
    this->SetNumberOfOutputPorts(2);
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
    if ( port == 0 )
    {
        //TODO remove if not needed
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
        return 1;
    }
    if( port == 1){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );    ////changed it from vtkdataobject
        return 1;
    }
    return 0;
}

int vtkFTLE::RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector)
{

    vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkInformation> outInfo1 = outputVector->GetInformationObject(1);           ////
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) ||
            outInfo1->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
        // get the update times
        double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

        double *inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        int numInTimes = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), inTimes, numInTimes);


        int resolutionSeedGrid[6];
        resolutionSeedGrid[0] = 0;
        resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0]) - 1;
        resolutionSeedGrid[2] = 0;
        resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1]) - 1;
        resolutionSeedGrid[4] = 0;
        resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2]) - 1;

        // Request numInTimes many time steps, starting from inTimes. inTimes is a pointer to the first time step (inTimes+3) would be a pointer to the fourth time step.

        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid, 6);

    }

    return 1;
}


/*
int vtkFTLE::RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector)
{
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

    ////changed 6 to 4

    int extentSeedGrid[4];
    extentSeedGrid[0] = 0;
    extentSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0])-1;
    extentSeedGrid[2] = 0;
    extentSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1])-1;
    ////extentSeedGrid[4] = 0;
    ////extentSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2])-1;

    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extentSeedGrid, 4);
    return 1;
}
*/



int vtkFTLE::RequestInformation(vtkInformation *vtkNotUsed(request),
                                vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);



    int resolutionSeedGrid[6];
    resolutionSeedGrid[0] = 0;
    resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0])-1;
    resolutionSeedGrid[2] = 0;
    resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1])-1;
    resolutionSeedGrid[4] = 0;
    resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2])-1;



    this->spacingSeedGrid[0] = this->boundsSeedGrid[0] *100 /static_cast<double>(resolutionSeedGrid[1]);
    this->spacingSeedGrid[1] = this->boundsSeedGrid[1] *100/static_cast<double>(resolutionSeedGrid[3]);
    this->spacingSeedGrid[2] = this->boundsSeedGrid[2] *100/static_cast<double>(resolutionSeedGrid[5]);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid,6);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), resolutionSeedGrid,6);
    outInfo->Set(vtkDataObject::SPACING(), this->spacingSeedGrid,3);
    outInfo->Set(vtkDataObject::ORIGIN(),this->originSeedGrid,3);

    outInfo1->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), resolutionSeedGrid,6);
    outInfo1->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), resolutionSeedGrid,6);
    outInfo1->Set(vtkDataObject::SPACING(), this->spacingSeedGrid,3);
    outInfo1->Set(vtkDataObject::ORIGIN(),this->originSeedGrid,3);
    //outInfo->Set( vtkCompositeDataPipeline::COMPOSITE_DATA_META_DATA(), inInfo);

    return 1;


}

/*
int vtkFTLE::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);


    //// changed 6 to 4
    int extentSeedGrid[4];
    extentSeedGrid[0] = 0;
    extentSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0])-1;
    extentSeedGrid[2] = 0;
    extentSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1])-1;
    ////extentSeedGrid[4] = 0;
    ////extentSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2])-1;



    this->spacingSeedGrid[0] = this->boundsSeedGrid[0]/static_cast<double>(extentSeedGrid[1]);
    this->spacingSeedGrid[1] = this->boundsSeedGrid[1]/static_cast<double>(extentSeedGrid[3]);
    ////this->spacingSeedGrid[2] = this->boundsSeedGrid[2]/static_cast<double>(extentSeedGrid[5]);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extentSeedGrid,4);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extentSeedGrid,4);

    ////changed 3 to 2
    outInfo->Set(vtkDataObject::SPACING(), this->spacingSeedGrid,2);
    //outInfo->Set(vtkDataObject::ORIGIN(),this->originSeedGrid,3);




    return 1;


}
*/

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
            vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
            vtkMultiBlockDataSet *output1 = vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));


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

    double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    double* times = outputVector->GetInformationObject(0)->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());        //// different timesteps
    int numInTimes = outputVector->GetInformationObject(0)->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());    //// number of timesteps


    vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkMultiBlockDataSet *output1 = vtkMultiBlockDataSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

    ////iterators
    int InitialTime, FinalTime;

    //if ( this->tauFinal >= numTimeSteps )
    if ( this->tauFinal > times[numInTimes - 1] || this->tauInitial < times[0])
    {
        std::cerr << "Incorrect time values" << std::endl;
        return false;
    }
        /*
    else
    {
        InitialTime = static_cast<unsigned int>(this->tauInitial);

        auto timePos = static_cast<unsigned int>(this->tauFinal);

        if ( (this->tauFinal) > timePos )
        {
            FinalTime = timePos + 1;
        }
        else if ((this->tauFinal) == timePos)
        {
            FinalTime = timePos;
        }
    }
         */

    else
    {
        for ( int k = 0; k < numInTimes - 1; k++)
        {
            if ( this->tauInitial >= times[k] && this->tauInitial < times[k + 1])
            {
                InitialTime = k;
            }

        }

        for (  int l = 0; l < numInTimes - 1; l++)
        {
            if ( this->tauFinal > times[l] && this->tauFinal <= times[l + 1] )
            {
                FinalTime = l + 1;
            }
        }

    }

    for (int a = 0; a < numInTimes; a ++)
    {
        std::cout <<"time steps" << times[a] << std::endl;
    }

    std::cout << "time iterators" << InitialTime << " " <<FinalTime << std::endl;

    std::vector<vector<double>> previousPts;

    //for (unsigned int i = 0; i < numTimeSteps - 1; i ++)
    for ( int i = InitialTime; i < FinalTime;  i ++)
    {
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

        this->seedGrid = vtkSmartPointer<vtkImageData>::New();



    int resolutionSeedGrid[6];
        resolutionSeedGrid[0] = 0;
        resolutionSeedGrid[1] = static_cast<int>(this->dimensionSeedGrid[0]) - 1;
        resolutionSeedGrid[2] = 0;
        resolutionSeedGrid[3] = static_cast<int>(this->dimensionSeedGrid[1]) - 1;
        resolutionSeedGrid[4] = 0;
        resolutionSeedGrid[5] = static_cast<int>(this->dimensionSeedGrid[2]) - 1;

        ///Set Size and extent of the seedng grid
        seedGrid->SetExtent(resolutionSeedGrid);
        seedGrid->SetOrigin(this->originSeedGrid);

        //double seedgridExtent[3];
        //seedgridExtent[0] = boundsSeedGrid[0] / this->inputGridSpacing[0];
        //seedgridExtent[1] = boundsSeedGrid[1] / this->inputGridSpacing[1];
        //seedgridExtent[2] = boundsSeedGrid[2] / this->inputGridSpacing[2];

        ///Calculate and set spacing of the seeding grid
        this->spacingSeedGrid[0] = this->boundsSeedGrid[0] *100 / (static_cast<double>(resolutionSeedGrid[1])) ;
        this->spacingSeedGrid[1] = this->boundsSeedGrid[1] *100/ (static_cast<double>(resolutionSeedGrid[3]));
        this->spacingSeedGrid[2] = this->boundsSeedGrid[2] *100/ (static_cast<double>(resolutionSeedGrid[5]));
       // this->spacingSeedGrid[0] = seedgridExtent[0] / (static_cast<double>(resolutionSeedGrid[1])) ;
       // this->spacingSeedGrid[1] = seedgridExtent[1]/ (static_cast<double>(resolutionSeedGrid[3]));
       // this->spacingSeedGrid[2] = seedgridExtent[2] / (static_cast<double>(resolutionSeedGrid[5]));
        this->seedGrid->SetSpacing(this->spacingSeedGrid);



    #if VTK_MAJOR_VERSION <= 5
        this->seedGrid->SetNumberOfScalarComponents(1);
        this->seedGrid->SetScalarTypeToDouble();
        this->seedGrid->AllocateScalars();
    #else
        this->seedGrid->AllocateScalars(VTK_DOUBLE, 1);
        //this->nextseedGrid->AllocateScalars(VTK_DOUBLE, 1);
    #endif

        if (!this->checkBounds(this->boundsSeedGrid))
        {
            vtkErrorMacro("The user input bounds extent the data bounds.");
            //output->DeepCopy(seedGrid);
            //return 1;
        }


        this->printSeedGridProperties();
        vtkSmartPointer<vtkDoubleArray> flowMap = vtkSmartPointer<vtkDoubleArray>::New();
        flowMap->SetNumberOfComponents(3);
        flowMap->SetName("FlowMap");
        flowMap->SetNumberOfTuples(seedGrid->GetNumberOfPoints());
        vtkSmartPointer<vtkDoubleArray> seeded = vtkSmartPointer<vtkDoubleArray>::New();
        seeded->SetNumberOfComponents(1);
        seeded->SetName("seeded");
        seeded->SetNumberOfTuples(seedGrid->GetNumberOfPoints());

        cout << "Fill seed grid with zeros" << endl;
        int *dims = this->seedGrid->GetDimensions();
        double *dummyPtr;
        for (int z = 0; z < dims[2] - 1; z++) {
            for (int y = 0; y < dims[1] - 1; y++) {
                for (int x = 0; x < dims[0] - 1; x++) {

                    auto *ptr = static_cast<double *>(this->seedGrid->GetScalarPointer(x, y, z));
                    ptr[0] = 0;
                    int index = getIndex(z, y, x, dims);                        ///changed y,x
                    vtkIdType idx = vtkIdType(index);
                    double temp = 0.0;
                    dummyPtr = &temp;
                    flowMap->SetTuple(idx, dummyPtr);
                }
            }
        }


        vtkSmartPointer<vtkPolyData> streamLines = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkStreamTracer> streamTracer = vtkSmartPointer<vtkStreamTracer>::New();
        vtkSmartPointer<vtkCallbackCommand> callBack = vtkSmartPointer<vtkCallbackCommand>::New();
        streamLines = vtkSmartPointer<vtkPolyData>::New();

        callBack->SetCallback(progressFunction);


        cout << "Integrate stream lines" << endl;

        ///Needed to tell the stream tracer which vectors to use
        this->inputGrid->GetPointData()->SetActiveVectors(this->inputGrid->GetPointData()->GetArrayName(0));
        this->nextinputGrid->GetPointData()->SetActiveVectors(this->nextinputGrid->GetPointData()->GetArrayName(0));
        vtkSmartPointer<vtkDoubleArray> ftleField;

        if (this->useGPU) {

            int temp = 1;
            if (this->tauFinal < 0)
                temp = -1;
            cout << "Start Integration" << endl;


            Integrator *customIntegrator = new Integrator(temp, this->integrationStepSize, fabs(this->tauInitial), fabs(this->tauFinal), LIMIT_DOUBLE,
                                                          this->inputGrid, this->seedGrid, this->nextinputGrid, i, InitialTime, FinalTime, times, previousPts);
            cout << "init integrator" << endl;
            this->inputGrid->GetOrigin(customIntegrator->origin);
            customIntegrator->setOriginSource(this->originSeedGrid);
            customIntegrator->calcStreamlines = this->showStreamLines;
            customIntegrator->maxNumberOfSteps = this->maxNumberOfSteps;

            ////storing the values of enpoints of every iteration in a vector array, to be used for next iteration
            previousPts = customIntegrator->PathlineEndPoints;

            ////freeing up the memory
            customIntegrator->PathlineEndPoints.erase (customIntegrator->PathlineEndPoints.begin(),customIntegrator->PathlineEndPoints.end());
            customIntegrator->NewEndPoints.erase (customIntegrator->NewEndPoints.begin(),customIntegrator->NewEndPoints.end());


            //streamLines->PrintSelf(std::cout, indentMyAss);
            clock_t begin = std::clock();
            high_resolution_clock::time_point t1 = high_resolution_clock::now();

            customIntegrator->integrateRK4GPUWrapper();
            this->seedGrid = customIntegrator->source;
            //ftleField->Print(std::cout);
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            clock_t end = clock();
            //double duration = duration_cast<microseconds>( t2 - t1 ).count() ;
            double duration = double(end - begin) / CLOCKS_PER_SEC;
            std::cout << "Time for integration " << duration << " seconds" << std::endl;
            delete customIntegrator;

        }
        else
        {
            int temp = 1;
            if (this->tauFinal < 0)
                temp = -1;
            cout << "Start Integration" << endl;


            Integrator *customIntegrator = new Integrator(temp, this->integrationStepSize, fabs(this->tauInitial), fabs(this->tauFinal), LIMIT_DOUBLE,
                                                          this->inputGrid, this->seedGrid, this->nextinputGrid, i, InitialTime, FinalTime, times, previousPts);
            cout << "init integrator" << endl;
            this->inputGrid->GetOrigin(customIntegrator->origin);
            customIntegrator->setOriginSource(this->originSeedGrid);
            customIntegrator->calcStreamlines = this->showStreamLines;
            customIntegrator->maxNumberOfSteps = this->maxNumberOfSteps;


            //streamLines->PrintSelf(std::cout, indentMyAss);
            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            clock_t begin = std::clock();

            streamLines = customIntegrator->integrateRK4();
            clock_t end = clock();
            high_resolution_clock::time_point t2 = high_resolution_clock::now();

            ////storing the values of enpoints of every iteration in a vector array, to be used for next iteration
            previousPts = customIntegrator->PathlineEndPoints;

            std::cout<< "previous pts size " <<previousPts.size()<<std::endl;

            ////freeing up the memory
            customIntegrator->PathlineEndPoints.erase (customIntegrator->PathlineEndPoints.begin(),customIntegrator->PathlineEndPoints.end());
            customIntegrator->NewEndPoints.erase (customIntegrator->NewEndPoints.begin(),customIntegrator->NewEndPoints.end());

            double duration = double(end - begin) / CLOCKS_PER_SEC;
            std::cout << "Time for integration " << duration / 1000.0 << " seconds" << std::endl;
            delete customIntegrator;
        }

        ////set the output in the dataset, to be used later to append the previous results
        newpoints->SetBlock(i, streamLines);

        cout << "Get end positions of stream line integration " << endl;
        cout << "Number of seed grid points " << this->seedGrid->GetNumberOfPoints() << endl;
        cout << "Number of seed cells " << this->seedGrid->GetNumberOfCells() << endl;
        cout << "Number of stream line cells " << streamLines->GetNumberOfCells() << endl;
        cout << "Number of stream lines " << streamLines->GetNumberOfLines() << endl;


        if (0 == streamLines->GetNumberOfCells() && !useGPU) {
            //output->DeepCopy(this->seedGrid);
            return 1;
        }

        if (this->useGPU) {

        } else {

            if (this->seedGrid->GetNumberOfPoints() == streamLines->GetNumberOfCells()) {      ////
                double point[3];
                double seed[1];

                for (int indexLine = 0; indexLine < streamLines->GetNumberOfCells(); ++indexLine) {        ////
                    int lastPointIndex = streamLines->GetCell(indexLine)->GetNumberOfPoints() - 1;

                    //std::cout <<"lastPointIndex " << lastPointIndex << std::endl;
                    if (lastPointIndex >= 0)
                    {
                        streamLines->GetCell(indexLine)->GetPoints()->GetPoint(lastPointIndex, point);         ////
                        flowMap->SetTuple(indexLine, point);
                        seed[0] = 1;
                        seeded->SetTuple(indexLine, seed);
                    } else
                    {
                        std::cerr << "stream Line not seeded" << std::endl;
                    }

                }

            }


        }

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
        }


        vtkSmartPointer<vtkPolyData> newlines = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    //// pass all the arrays here as points, then draw the line through points

        for (int j = InitialTime; j <= i; j++) {

            vtkSmartPointer<vtkDataObject> addData = newpoints->GetBlock(j);

            vtkSmartPointer<vtkPolyData> streamlineData = vtkSmartPointer<vtkPolyData>::New();
            streamlineData->ShallowCopy(addData);

            appendFilter->AddInputData(streamlineData);
            appendFilter->Update();
        }

        newlines->ShallowCopy(appendFilter->GetOutput());

        /*

        vtkIdType numCells = streamlineData->GetNumberOfCells();

        std::cout<<"number of cells "<< numCells <<std::endl;

        vtkSmartPointer<vtkPoints> pointstoadd = vtkSmartPointer<vtkPoints>::New();

        //vtkIdType numberpoints = streamlineData->GetNumberOfPoints();
        //std::cout<< "number of points in the polydata "<< numberpoints <<std::endl;


        //vtkSmartPointer<vtkPoints> pointstoadd = streamlineData->GetPoints();
        //vtkIdType numPoints = pointstoadd->GetNumberOfPoints();               ////or
        //std::cout<< "number of points to add" <<numPoints<<std::endl;         ////or



        //outputlines->InsertNextCell(numPoints);                               ////or
        //outputlines= streamlineData->GetLines();



        for (int k = 0; k < numPoints; k++)
        {
            vtkIdType id = vtkIdType(k);
            pointstoadd->GetPoint(id,newpoint);
            vtkIdType nextPoint1 = outputpoints->InsertNextPoint(newpoint);
            outputlines->InsertCellPoint(nextPoint1);
        }


        for (int k = 0; k < numCells; k++)
        {
            vtkIdType id = vtkIdType(k);
            vtkSmartPointer<vtkCell> celldata;
            celldata=streamlineData->GetCell(k);

            //pointstoadd->GetPoint(id,newpoint);                           ////
            pointstoadd = celldata->GetPoints();

            vtkIdType numpointstoadd = pointstoadd->GetNumberOfPoints();
            //check += numpointstoadd;

            //std::cout<<"num of points to add " << numpointstoadd <<std::endl;
            //vtkIdType nextPoint1 = outputpoints->InsertNextPoint(newpoint);
            //outputlines->InsertCellPoint(pointstoadd->GetNumberOfPoints());
            outputlines->InsertNextCell(celldata);

            for ( int l = 0; l < pointstoadd->GetNumberOfPoints(); l++ )
            {
                double newpoint[2];

                vtkIdType id1 = vtkIdType(l);
                pointstoadd->GetPoint(id1, newpoint);
                outputpoints->InsertNextPoint(newpoint);
            }

        }
        */


        vtkIdType numoutputpoints = newlines->GetNumberOfPoints();
        std::cout << "num output points" << numoutputpoints << std::endl;

        vtkIdType numOfCells = newlines->GetNumberOfCells();
        std::cout << "num Of cells" << numOfCells << std::endl;


        ///Translate SeedGrid to correct origin
        double newOrigin[3];
        this->inputGrid->GetOrigin(newOrigin);
        newOrigin[0] += this->originSeedGrid[0];
        newOrigin[1] += this->originSeedGrid[1];
        newOrigin[2] += this->originSeedGrid[2];
        this->seedGrid->SetOrigin(newOrigin);

        output->SetBlock(i, this->seedGrid);
        output1->SetBlock(i, newlines);

        std::cout<< "Current Progress " << int((double(i + 1 - InitialTime)/double(FinalTime - InitialTime))* 100.0)<< std::endl;

        }
    /*
        i += 1;

        if (i == numTimeSteps - 1)
        {
            i = 0;
            //break;

        }
*/

    return 1;
}


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

/*

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
    mat3 jacobi;
    vec3 du;
    vec3 dv;
    vec3 dw;


    for (int z=0; z<dataDims[2]; z++)
    {
        for (int y=0; y<dataDims[1]; y++)
        {
            for (int x=0; x<dataDims[0]; x++)
            {
                int index = getIndex(z,y,x,dataDims);
                //double grad[3];
                vtkIdType id = vtkIdType(index);
                double t1[3], t2[3];

                if(z==0)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z+1,y,x,dataDims);
                    field->GetTuple(id,t2);
                    du[2] = forwardDiff(t1[0],t2[0],spac[0]);
                    dv[2] = forwardDiff(t1[1],t2[1],spac[1]);
                    dw[2] = forwardDiff(t1[2],t2[2],spac[2]);

                }
                else if (z == dataDims[2]-1)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z-1,y,x,dataDims);
                    field->GetTuple(id,t2);
                    du[2] = backwardDiff(t1[0],t2[0],spac[0]);
                    dv[2] = backwardDiff(t1[1],t2[1],spac[1]);
                    dw[2] = backwardDiff(t1[2],t2[2],spac[2]);

                }
                else
                {
                    id = getIndex(z-1,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z+1,y,x,dataDims);
                    field->GetTuple(id,t2);
                    du[2] = centralDiff(t1[0],t2[0],spac[0]);
                    dv[2] = centralDiff(t1[1],t2[1],spac[1]);
                    dw[2] = centralDiff(t1[2],t2[2],spac[2]);
                }


                /// y Component
                if(y==0)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y+1,x,dataDims);
                    field->GetTuple(id,t2);
                    du[1] = forwardDiff(t1[0],t2[0],spac[0]);
                    dv[1] = forwardDiff(t1[1],t2[1],spac[1]);
                    dw[1] = forwardDiff(t1[2],t2[2],spac[2]);
                }
                else if (y == dataDims[1]-1)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y-1,x,dataDims);
                    field->GetTuple(id,t2);
                    du[1] = backwardDiff(t1[0],t2[0],spac[0]);
                    dv[1] = backwardDiff(t1[1],t2[1],spac[1]);
                    dw[1] = backwardDiff(t1[2],t2[2],spac[2]);

                }
                else
                {
                    id = getIndex(z,y-1,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y+1,x,dataDims);
                    field->GetTuple(id,t2);
                    du[1] = centralDiff(t1[0],t2[0],spac[0]);
                    dv[1] = centralDiff(t1[1],t2[1],spac[1]);
                    dw[1] = centralDiff(t1[2],t2[2],spac[2]);
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
                    dw[0] = forwardDiff(t1[2],t2[2],spac[2]);
                }
                else if (x == dataDims[0]-1)
                {
                    id = getIndex(z,y,x,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y,x-1,dataDims);
                    field->GetTuple(id,t2);
                    du[0] = backwardDiff(t1[0],t2[0],spac[0]);
                    dv[0] = backwardDiff(t1[1],t2[1],spac[1]);
                    dw[0] = backwardDiff(t1[2],t2[2],spac[2]);

                }
                else
                {
                    id = getIndex(z,y,x-1,dataDims);
                    field->GetTuple(id,t1);
                    id = getIndex(z,y,x+1,dataDims);
                    field->GetTuple(id,t2);
                    du[0] = centralDiff(t1[0],t2[0],spac[0]);
                    dv[0] = centralDiff(t1[1],t2[1],spac[1]);
                    dw[0] = centralDiff(t1[2],t2[2],spac[2]);

                }

                mat3setrows(jacobi, du,dv,dw);
                mat3 jacobiT;
                mat3 cgTensor;
                mat3trp(jacobi,jacobiT);
                mat3mul(jacobiT,jacobi, cgTensor);
                double eMax = 0;
                vec3 eigenV;
                int realEigen = mat3eigenvalues(cgTensor, eigenV);
                eMax = max(max(eigenV[0],eigenV[1]),eigenV[2]);


                if(realEigen !=3 || fabs(eMax)< LIMIT_DOUBLE)
                {
                    //std::cerr<<"Eigenvalues of Cauchy Green Tensor are not real"<<endl;
                    //out<<realEigen<< " eigenvalues"<< eigenV[0]<< " " <<eigenV[1]<< " "<< eigenV[2]<<" "<<LIMIT_DOUBLE<<endl;
                    eMax = 0;
                    id = getIndex(z,y,x,dataDims);
                    /// Possible error source because eMax is only an scalar
                    ftle->SetTuple(id,&eMax);
                    continue;
     ///Translate SeedGrid to correct origin           }



                eMax = 1.0/ fabs(this->tau) * log(sqrt(eMax));
                id = getIndex(z,y,x,dataDims);
                /// Possible error source because eMax is only an scalar
                ftle->SetTuple(id,&eMax);

            }
        }
    }


    return ftle;

}

*/

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


    for (int z=0; z<dataDims[2]; z++)
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
                    ////changed (z,y,x) to (y,x)										
                    id = getIndex(z,y,x,dataDims);
                    /// Possible error source because eMax is only an scalar
                    ftle->SetTuple(id,&eMax);
                    continue;        
                }

                eMax = 1.0/ fabs((this->tauFinal) - (this->tauInitial)) * log(sqrt(eMax));
                ////changed (z,y,x) to (y,x)											
                id = getIndex(z,y,x,dataDims);
                /// Possible error source because eMax is only an scalar
                ftle->SetTuple(id,&eMax);

            }
        }
    }


    return ftle;

}



void vtkFTLE::printSeedGridProperties()
{
    int* ext1 = this->seedGrid->GetExtent();
    double* space1 = this->seedGrid->GetSpacing();
    int* dimensions = this->seedGrid->GetDimensions();

    cout<<"_________________________________________________________________________"<<endl;
    //cout<<"Extent: "<<ext1[0]/100.<<" "<<ext1[1]/100.<<" "<<ext1[2]/100.<<" "<<ext1[3]/100.<<" "<<ext1[4]/100.<<" "<<ext1[5]/100.<<" "<<endl;
    cout<<"Extent: "<<ext1[0]<<" "<<ext1[1]*space1[0]/100.<<" "<<ext1[2]<<" "<<ext1[3]*space1[1]/100.<<" "<<ext1[4]<<" "<<ext1[5]*space1[2]/100.<<" "<<endl;
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
    bounds[0] = inputDim[0] / spacing[0];
    bounds[1] = inputDim[1] / spacing[1];
    bounds[2] = inputDim[2] / spacing[2];
    if(this->boundsSeedGrid[0] > bounds[0] || this->boundsSeedGrid[1] > bounds[1] || this->boundsSeedGrid[2] > bounds[2])
    {
        return false;
    }

    return true;
}