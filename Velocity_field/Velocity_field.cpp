#include <vtkVersion.h>
#include <vtkArrowSource.h>
#include <vtkCellArray.h>
#include <vtkGlyph2D.h>
#include <vtkImageData.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkInteractorStyleImage.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkFieldData.h>
#include <stdio.h>

void double_gyre (float x, float y, float &u, float &v);
float f(float x);

int getIndex(int z, int y, int x, int* dataDims)
{
    return (z*dataDims[1]*dataDims[0]+ y*dataDims[0]+x);
}

int main()
{
    float x, y, u, v;

    ///initialising image data and dimensions
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(200,100,2);
    //imageData->SetDimensions(8,6,2);



    #if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(3);
    imageData->SetScalarTypeToFloat();
    imageData->AllocateScalars();
    #else
    imageData->AllocateScalars(VTK_FLOAT,3);
    #endif


    int* dims = imageData->GetDimensions();

    //vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
    //scalars->SetNumberOfComponents(3);
    //scalars->SetName("Scalars");

    //vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for ( float k = 0.; k < dims[2]; k++)
    {
    	for ( float j = 0.; j < dims[1]; j++)
    	{

             for ( float i = 0.; i < dims[0]; i++)
             {

                 //y = (- 1 * dims[1]/2 + j) / 100;
                 x = 1.0;

                 ///storing the values in pixel

        	     //std::cout<< i << ", " << j << " " <<  y <<  std::endl;
        	     float* pixel = static_cast<float*>(imageData->GetScalarPointer(i,j,k));
	             //std::cout<<"Ptr: " << pixel << endl;
                    pixel[0] = x;
                    pixel[1] = 0.0;
                    pixel[2] = 0.0;

             }
    	}
    }


    //imageData->GetPointData()->SetScalars(scalars);
    imageData->GetPointData()->SetActiveVectors("ImageScalars");


    ///wrting imagedata to file

    vtkSmartPointer<vtkXMLImageDataWriter> writer =vtkSmartPointer<vtkXMLImageDataWriter>::New();
   #if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
   #else
    writer->SetInputData(imageData);
   #endif
    writer->SetFileName("Xvelocityfield_20000.vti");
    writer->Write();


    return EXIT_SUCCESS;


}

/*
///Double Gyre calculations

void double_gyre (float x, float y, float&u, float&v)
{
    //float A = 0.1;   ////changed from 0.1 to 1000000
    //float h = 0.01;

    //float df = ( f( x + h ) - f( x ) )/ h;

    //velocity field
    //u = -1. * vtkMath::Pi() * A * sin ( vtkMath::Pi() * f(x) ) * cos ( vtkMath::Pi() * y );
    //v = vtkMath::Pi() * A * cos ( vtkMath::Pi() * f(x) ) * sin ( vtkMath::Pi() * y ) * df;

    u = - 1 * vtkMath::Pi() * sin ( vtkMath::Pi() * x ) * cos ( vtkMath::Pi() * y );
    v = vtkMath::Pi() * cos ( vtkMath::Pi() * x ) * sin ( vtkMath::Pi() * y );

    return;

}

float f (float x)
{
    float a, b, f;
    float t = 10.;
    float eps = 0.25;                ////changed from 0.25 to 2500000
    float omega = 2. * vtkMath::Pi() / 10.;

    a = eps * sin (omega * t);
    b = 1. - 2. * eps * sin (omega * t);

    f = a * x * x + b * x;

    return f;

}
*/