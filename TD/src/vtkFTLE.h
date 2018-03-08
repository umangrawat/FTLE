#ifndef __vtkFTLE_h
#define __vtkFTLE_h


#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkImageAlgorithm.h"
#include "vtkSmartPointer.h" // compiler errors if this is forward declared
#include "vtkPointSet.h"
#include <vtkMultiTimeStepAlgorithm.h>
#include <vtkMultiBlockDataSet.h>

using namespace std;
typedef unsigned int uint;


class vtkPolyData;
class vtkImageData;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;

////changed 3D to 2D
////#define DUMPVECTOR(x) (cout << "(" << x[0] << "," << x[1] << "," << x[2] << ")" << endl)
#define DUMPVECTOR(x) (cout << "(" << x[0] << "," << x[1] << ")" << endl)



class vtkFTLE : public vtkMultiTimeStepAlgorithm
{
public:
    static vtkFTLE *New();
    vtkTypeMacro(vtkFTLE, vtkMultiTimeStepAlgorithm)
    void PrintSelf(ostream &os, vtkIndent indent);

    vtkSetMacro(tauInitial, double);
    vtkGetMacro(tauInitial, double);

    vtkSetMacro(tauFinal, double);
    vtkGetMacro(tauFinal, double);

    vtkSetMacro(maxNumberOfSteps, int);
    vtkGetMacro(maxNumberOfSteps, int);

    vtkSetMacro(integrationStepSize, double);
    vtkGetMacro(integrationStepSize, double);
    
    vtkSetVector3Macro(originSeedGrid, double);
    vtkGetVector3Macro(originSeedGrid, double);

    vtkSetVector3Macro(boundsSeedGrid, double);
    vtkGetVector3Macro(boundsSeedGrid, double);

    vtkSetVector3Macro(dimensionSeedGrid, double);
    vtkGetVector3Macro(dimensionSeedGrid, double);

    vtkSetMacro(useGPU, bool);
    vtkGetMacro(useGPU, bool);

    vtkSetMacro(showStreamLines, bool);
    vtkGetMacro(showStreamLines, bool);


    vtkSmartPointer<vtkImageData> inputGrid;
    vtkSmartPointer<vtkImageData> seedGrid;
    double* inputGridSpacing;
    int* inputGridDimensions;

    vtkSmartPointer<vtkImageData> nextinputGrid;
    double* nextinputGridSpacing;
    int* nextinputGridDimensions;




protected:
    vtkFTLE();
    ~vtkFTLE();
    uint dim[3];
    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
    // Generate output
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
    int RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector);
    int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    void getData(vtkSmartPointer<vtkImageData> in_grid, double* &in_data, int numComponents, bool &allocated_in_data);
    vtkSmartPointer<vtkDoubleArray> computeDerivatives(vtkSmartPointer<vtkImageData> data, char* arrayName);
    bool checkBounds(double* );

    double tauInitial;
    double tauFinal;

    double originSeedGrid[3];
    double boundsSeedGrid[3];
    double dimensionSeedGrid[3];
    double spacingSeedGrid[3];
    int maxNumberOfSteps;

    double integrationStepSize;
    bool useGPU;
    bool showStreamLines;

    vtkSmartPointer<vtkImageData> m_evalImageData;
    public:
    void printSeedGridProperties();
};
#endif



/*
////changed 3D to 2D
class vtkFTLE : public vtkImageAlgorithm
{
public:
    static vtkFTLE *New();
    vtkTypeMacro(vtkFTLE, vtkImageAlgorithm)
    void PrintSelf(ostream &os, vtkIndent indent);

    vtkSetMacro(tau, double);
    vtkGetMacro(tau, double);
    
    vtkSetMacro(maxNumberOfSteps, int);
    vtkGetMacro(maxNumberOfSteps, int);

    vtkSetMacro(integrationStepSize, double);
    vtkGetMacro(integrationStepSize, double);
    
    vtkSetVector2Macro(originSeedGrid, double);
    vtkGetVector2Macro(originSeedGrid, double);

    vtkSetVector2Macro(boundsSeedGrid, double);
    vtkGetVector2Macro(boundsSeedGrid, double);

    vtkSetVector2Macro(dimensionSeedGrid, double);
    vtkGetVector2Macro(dimensionSeedGrid, double);

    vtkSetMacro(useGPU, bool);
    vtkGetMacro(useGPU, bool);

    vtkSetMacro(showStreamLines, bool);
    vtkGetMacro(showStreamLines, bool);

 	




    vtkSmartPointer<vtkImageData> inputGrid;
    vtkSmartPointer<vtkImageData> seedGrid;
    double* inputGridSpacing;
    int* inputGridDimensions;


protected:
    vtkFTLE();
    ~vtkFTLE();
    uint dim[2];
    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
    // Generate output
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
    int RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector);

    void getData(vtkSmartPointer<vtkImageData> in_grid, double* &in_data, int numComponents, bool &allocated_in_data);
    vtkSmartPointer<vtkDoubleArray> computeDerivatives(vtkSmartPointer<vtkImageData> data, char* arrayName);
    bool checkBounds(double* );

    double tau;


    double originSeedGrid[2];
    double boundsSeedGrid[2];
    double dimensionSeedGrid[2];
    double spacingSeedGrid[2];
    int maxNumberOfSteps;

    double integrationStepSize;
    bool useGPU;
    bool showStreamLines;






    vtkSmartPointer<vtkImageData> m_evalImageData;
    public:
    void printSeedGridProperties();
};
#endif
*/
