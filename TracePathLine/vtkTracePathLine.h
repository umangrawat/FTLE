#ifndef __vtkTracePathLine_h
#define __vtkTracePathLine_h

#include <vtkDataSetAlgorithm.h>
#include "vtkSmartPointer.h"
#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkImageAlgorithm.h"
#include "vtkPointSet.h"

using namespace std;
typedef unsigned int uint;


class vtkInformation;
class vtkInformationVector;
class vtkDataSet;
class vtkPolyData;
class vtkImageData;

class vtkTracePathLine : public vtkDataSetAlgorithm
{
public:
    static vtkTracePathLine *New();

    vtkTypeMacro(vtkTracePathLine, vtkDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    vtkSetVector3Macro(origin, double);
    vtkGetVector3Macro(origin, double);

    vtkSetVector3Macro(bounds, double);
    vtkGetVector3Macro(bounds, double);

    vtkSetVector3Macro(dimension, int);
    vtkGetVector3Macro(dimension, int);

    vtkSmartPointer<vtkImageData> inputGrid;
    vtkSmartPointer<vtkImageData> seedGrid;
    double* inputGridSpacing;
    int* inputGridDimensions;




protected:
    vtkTracePathLine();
    ~vtkTracePathLine();

    // Usual data generation method
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int FillOutputPortInformation( int port, vtkInformation* info );
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
    virtual int RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector);
    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    // Generate output

    double origin[3];
    double bounds[3];
    int dimension[3];
    double spacing[3];

private:

};

#endif