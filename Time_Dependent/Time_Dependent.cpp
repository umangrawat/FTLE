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
#include <string.h>

void double_gyre (float x, float y, float t, float &u, float &v);
float f(float x, float t);

int main()
{
    float t;
    int ind = 0;
    int limit;
    std::cout<< "Enter value of limit" << std::endl;
    std::cin>>limit;
    while ( t < limit )
    {
        float x, y, u, v;


        ///initialising image data and dimensions
        vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
        imageData->SetDimensions(201,101,2);
        //imageData->SetSpacing(0.01,0.01,0.01);
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

        for (float k = 0.; k < dims[2]; k++)
        {
            for (float j = 0.; j < dims[1]; j++)
            {
                for (float i = 0.; i < dims[0]; i++)
                {

                    x = i / 100;
                    y = j / 100;

                    double_gyre(x, y, t, u, v);
                    //points->InsertNextPoint( x, y, 0 );
                    //scalars->InsertTuple3(idx1, u, v, 0 );

                    ///storing the values in pixel

                    //std::cout<< u << ", " << v << " " << x <<" " << y <<  std::endl;
                    float *pixel = static_cast<float *>(imageData->GetScalarPointer(i, j, k));
                    //std::cout<<"Ptr: " << pixel << endl;
                    pixel[0] = u;
                    pixel[1] = v;
                    pixel[2] = 0.0;

                }
            }
        }

        //imageData->GetPointData()->SetScalars(scalars);
        imageData->GetPointData()->SetActiveVectors("ImageScalars");

        std::stringstream index;
        index << ind;
        std::string filename = "doubleGyre_";
        std::string filename2 = ".vti";
        std::string path = filename + std::string(index.str()) + filename2;
        //filename = filename + string (index.str());

        ///wrting imagedata to file

        vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageData->GetProducerPort());
        #else
        writer->SetInputData(imageData);
        #endif
        writer->SetFileName(path.c_str());
        //writer->SetFileName("vel_0009.vti");
        writer->Write();

        ind++;
        t = t + 0.05;
    }

    return EXIT_SUCCESS;


}

///Double Gyre calculations

void double_gyre (float x, float y, float t, float&u, float&v)
{
    float A = 0.1;                  ////changed from 0.1 to 1000000
    float h = 0.01;

    float df = ( f( x + h, t ) - f( x, t ) )/ h;

    //velocity field
    u = -1. * vtkMath::Pi() * A * sin ( vtkMath::Pi() * f(x, t) ) * cos ( vtkMath::Pi() * y );
    v = vtkMath::Pi() * A * cos ( vtkMath::Pi() * f(x, t) ) * sin ( vtkMath::Pi() * y ) * df;

    //u = - 1 * vtkMath::Pi() * sin ( vtkMath::Pi() * x ) * cos ( vtkMath::Pi() * y );
    //v = vtkMath::Pi() * cos ( vtkMath::Pi() * x ) * sin ( vtkMath::Pi() * y );

    return;

}

float f (float x, float t)
{
    float a, b, f;
    //float t = 2.;
    float eps = 0.25;                ////changed from 0.25 to 2500000
    float omega = 2. * vtkMath::Pi() / 10.;

    a = eps * sin (omega * t);
    b = 1. - 2. * eps * sin (omega * t);

    f = a * x * x + b * x;

    return f;

}