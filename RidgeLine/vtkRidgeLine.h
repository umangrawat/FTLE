#ifndef __vtkRidgeLine_h
#define __vtkRidgeLine_h

#include <vtkDataSetAlgorithm.h>
#include "vtkSmartPointer.h"

using namespace std;
typedef unsigned int uint;


class vtkInformation;
class vtkInformationVector;
class vtkDataSet;
class vtkPolyData;
class vtkImageData;

class vtkRidgeLine : public vtkDataSetAlgorithm
{
public:
  static vtkRidgeLine *New();


  vtkTypeMacro(vtkRidgeLine, vtkDataSetAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

    vtkSetMacro(LimitEigenValue, double);
    vtkGetMacro(LimitEigenValue, double);
    double LimitEigenValue;



protected:
  vtkRidgeLine();
  ~vtkRidgeLine();

  // Usual data generation method
   virtual int FillInputPortInformation(int port, vtkInformation *info);

   virtual int FillOutputPortInformation( int port, vtkInformation* info );
   virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

   virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    // Generate output


private:

};

#endif
