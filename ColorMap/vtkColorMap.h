#ifndef __vtkColorMap_h
#define __vtkColorMap_h

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

class vtkColorMap : public vtkDataSetAlgorithm
{
public:
    static vtkColorMap *New();


    vtkTypeMacro(vtkColorMap, vtkDataSetAlgorithm);

    void PrintSelf(ostream& os, vtkIndent indent);

    vtkSetMacro(cellId, double);
    vtkGetMacro(cellId, double);
    double cellId;

    vtkSmartPointer<vtkImageData> seedGrid;

protected:
    vtkColorMap();
    ~vtkColorMap();

    // Usual data generation method
    virtual int FillInputPortInformation(int port, vtkInformation *info);

    virtual int FillOutputPortInformation( int port, vtkInformation* info );
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    void cellLocator( int cellId, int* dims, int& x, int& y);
    // Generate output


private:

};

#endif
