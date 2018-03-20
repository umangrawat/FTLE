#ifndef __vtkTracePathLine_h
#define __vtkTracePathLine_h

#include <vtkDataSetAlgorithm.h>
#include "vtkSmartPointer.h"

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

    //vtkSetMacro(LimitEigenValue, double);
    //vtkGetMacro(LimitEigenValue, double);
    //double LimitEigenValue;

protected:
    vtkTracePathLine();
    ~vtkTracePathLine();

    // Usual data generation method
    virtual int FillInputPortInformation(int port, vtkInformation *info);

    virtual int FillOutputPortInformation( int port, vtkInformation* info );
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    //bool RidgeEigenValue (double IntHessian[4]);
    // Generate output


private:

};

#endif