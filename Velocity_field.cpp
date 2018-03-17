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


int getIndex(int z, int y, int x, int* dataDims)
{
    return (z*dataDims[1]*dataDims[0]+ y*dataDims[0]+x);
}

int main()
{

    float t, limit, accelLeft, accelRight, boundaryVelLeft, boundaryVelRight = 0.;

    int ind = 0;

    std::cout<< "Enter value of time limit" << std::endl;
    std::cin>>limit;

    std::cout<< "Enter value of acceleration from Left " << std::endl;
    std::cin>>accelLeft;

    std::cout<< "Enter value of acceleration from Right " << std::endl;
    std::cin>>accelRight;

    //// multiplying with -1 since the acceleration is in opposite direction
    accelRight = accelRight * ( - 1.);

    int dims[3] = {201, 101, 2};
    double spacing[3] = {0.01, 0.01, 0.01};

    double extent[3] = {0.};
    extent[0] = (dims[0] - 1)* spacing[0];
    extent[1] = (dims[1] - 1)* spacing[1];
    extent[2] = (dims[2] - 1)* spacing[2];

    float x, y = 0.;

    while ( t < limit )
    {

        ///initialising image data and dimensions
        vtkSmartPointer <vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
        imageData->SetDimensions(dims);
        imageData->SetSpacing(spacing);

/*
    #if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(3);
    imageData->SetScalarTypeToDouble();
    imageData->AllocateScalars();
    #else
    imageData->AllocateScalars(VTK_FLOAT,3);
    #endif
*/

        vtkSmartPointer <vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
        int numPoints = dims[0] * dims[1] * dims[2];
        scalars->SetNumberOfComponents(3);
        scalars->SetNumberOfTuples(numPoints);
        scalars->SetName("Scalars");

        //float jacobian[4] = {0, 0, 0, 2.};

        for (float k = 0.; k < dims[2]; k++) {
            for (float j = 0.; j < dims[1]; j++) {
                for (float i = 0.; i < dims[0]; i++) {

                    int index1;
                    index1 = getIndex(k, j, i, dims);

                    double interpolator = i / ( dims[0] - 1 );
                    x = ( accelLeft * ( 1 - interpolator ) + accelRight *  interpolator ) * t;                            ////have to interpolate using acceleration from both sides

                    y = ( -1 * ( dims[1] - 1 ) / 2 + j ) * 2. / 100;
                    //std::cout << "velocity " << x <<" " << y << std::endl;

                    ///storing the values in pixel

                    float tuple[3];
                    tuple[0] = x;
                    tuple[1] = y;
                    tuple[2] = 0.0;
                    assert(index1 < dims[0] * dims[1] * dims[2]);

                    scalars->SetTuple(index1, tuple);

                    /*
                    //std::cout<< i << ", " << j << " " <<  y <<  std::endl;
                    float* pixel = static_cast<float*>(imageData->GetScalarPointer(i,j,k));
                    //std::cout<<"Ptr: " << pixel << endl;
                       pixel[0] = x;
                       pixel[1] = 0.0;
                       pixel[2] = 0.0;
                    */
                }
            }
        }

        imageData->GetPointData()->AddArray(scalars);

        std::stringstream index;
        index << ind;
        std::string filename = "VelocityField_";
        std::string filename2 = ".vti";
        std::string path = filename + std::string(index.str()) + filename2;

        ///wrting imagedata to file

        vtkSmartPointer <vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageData->GetProducerPort());
        #else
        writer->SetInputData(imageData);
        #endif
        writer->SetFileName(path.c_str());
        writer->Write();

        boundaryVelLeft = accelLeft * t;
        boundaryVelRight = accelRight * t;

        extent[0] = extent[0] - ( boundaryVelLeft - boundaryVelRight ) * t;
        spacing[0] = extent[0] / (dims[0] - 1);

        //std::cout << extent[0] << " " << spacing[0] <<std::endl;

        if ( extent[0] <= 0 )
        {
            std::cerr<<" Wrong acceleration values " <<std::endl;
            break;
        }

        ind++;
        t = t + 0.05;

    }

    return EXIT_SUCCESS;

}