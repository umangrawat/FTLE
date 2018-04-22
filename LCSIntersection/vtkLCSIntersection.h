#ifndef __vtkLCSIntersection_h
#define __vtkLCSIntersection_h

#include <vtkDataSetAlgorithm.h>
#include "vtkSmartPointer.h"
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>

using namespace std;
typedef unsigned int uint;


class vtkInformation;
class vtkInformationVector;
class vtkDataSet;
class vtkPolyData;
class vtkImageData;

class vtkLCSIntersection : public vtkDataSetAlgorithm
{
public:
    static vtkLCSIntersection *New();
    vtkTypeMacro(vtkLCSIntersection, vtkDataSetAlgorithm);

    void PrintSelf(ostream& os, vtkIndent indent);
    vtkSetVector3Macro(cellsNumber, int);
    vtkGetVector3Macro(cellsNumber, int);


protected:
    vtkLCSIntersection();
  ~vtkLCSIntersection();

  // Usual data generation method
    virtual int FillInputPortInformation(int port, vtkInformation *info);

    virtual int FillOutputPortInformation( int port, vtkInformation* info );
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual int cellLocator( double Location[2], double bounds[6], double spacing[3], int* resolution );
    //virtual int RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector);
    //virtual int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


    //bool RidgeEigenValue (double IntHessian[4]);
    // Generate output
    int cellsNumber[3];

private:

};

#endif
