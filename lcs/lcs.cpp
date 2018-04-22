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
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include "vtkIntersectionPolyDataFilter.h"
#include <vtkAppendPolyData.h>
#include <vtkXMLPolyDataWriter.h>

using namespace std;

int main(int argc, char *argv[])
{
    std::string filename = "15to10Ridge";
    std::string filename2 = ".vtp";
    std::string path = filename + filename2;

    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(path.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> data =  vtkSmartPointer<vtkPolyData>::New();
    data->ShallowCopy(reader->GetOutput());
    
    std::string filename1 = "10to15Ridge";
    std::string path1 = filename1 + filename2;
    vtkSmartPointer<vtkXMLPolyDataReader> reader1 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader1->SetFileName(path1.c_str());
    reader1->Update();

    vtkSmartPointer<vtkPolyData> nextdata =  vtkSmartPointer<vtkPolyData>::New();
    nextdata->ShallowCopy(reader1->GetOutput());

    //vtkSmartPointer<vtkBooleanOperationPolyDataFilter> boolFilter = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
    //boolFilter->SetOperationToUnion();
    //boolFilter->SetInputData(0, data);
    //boolFilter->SetInputData(1, nextdata);
    //boolFilter->Update();

    //vtkSmartPointer<vtkIntersectionPolyDataFilter> intersectionPolyDataFilter = vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
    //intersectionPolyDataFilter->SetInputConnection( 0, reader->GetOutputPort() );
    //intersectionPolyDataFilter->SetInputConnection( 1, reader1->GetOutputPort() );
    //intersectionPolyDataFilter->Update();




    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    appendFilter->AddInputData( data );
    appendFilter->Update();
    appendFilter->AddInputData( nextdata );
    appendFilter->Update();

    vtkSmartPointer<vtkPolyData> result =  vtkSmartPointer<vtkPolyData>::New();
    result->ShallowCopy(appendFilter->GetOutput());

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("10to15.vtp");
    #if VTK_MAJOR_VERSION <= 5
    writer->SetInput(result);
    #else
    writer->SetInputData(result);
    #endif
    writer->Write();

    
    return EXIT_SUCCESS;
    
}

