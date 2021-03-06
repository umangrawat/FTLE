// Loadable modules
//
// Generated by /export/home/urawat/Downloads/paraview_build/bin/vtkkwProcessXML-pv5.3
//
#ifndef __vtkSMXML_RidgeLine_h
#define __vtkSMXML_RidgeLine_h

#include <string.h>


// From file /export/home/urawat/Desktop/RidgeLine/RidgeLine.xml
static const char* const RidgeLineRidgeLineInterfaces0 =
"<ServerManagerConfiguration>\n"
"    <ProxyGroup name=\"filters\">\n"
"        <SourceProxy\n"
"                name=\"RidgeLine\"\n"
"                class=\"vtkRidgeLine\"\n"
"                label=\"RidgeLine\">	<!-- name as it will appear in ParaView under \"Filters\" menu -->\n"
"\n"
"            <Documentation\n"
"                    long_help=\"Extracts ridge lines on the given patch of the dataset\"\n"
"                    short_help=\"Extracts ridge lines\">\n"
"                Extracts ridge lines on uniform grids.\n"
"            </Documentation>\n"
"\n"
"            <InputProperty\n"
"                    name=\"Input\"\n"
"                    command=\"AddInputConnection\"\n"
"                    clean_command=\"RemoveAllInputs\">\n"
"                <ProxyGroupDomain name=\"groups\">\n"
"                    <Group name=\"sources\"/> <!-- the input might come either from a reader-->\n"
"                    <Group name=\"filters\"/> <!-- ... or from another filter -->\n"
"                </ProxyGroupDomain>\n"
"                <DataTypeDomain name=\"input_type\">\n"
"                    <DataType value=\"vtkImageData\"/> <!-- We require that the data has type vtkDataSet; might as well be more specific (e.g., vtkRectilinearGrid) -->\n"
"                </DataTypeDomain>\n"
"            </InputProperty>\n"
"\n"
"\n"
"            <DoubleVectorProperty\n"
"                    name=\"SetEigenvalueLimit\"\n"
"                    command=\"SetLimitEigenValue\"\n"
"                    number_of_elements=\"1\"\n"
"                    default_values=\"0.0\">\n"
"                <DoubleRangeDomain name=\"range\" min=\"-10000.0\" max=\"10000.0\" />\n"
"                <Documentation>Set the minimum limit of eigenvalue of the data.</Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"        </SourceProxy>\n"
"    </ProxyGroup>\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* RidgeLineRidgeLineInterfaces()
{
  size_t len = ( 0
    + strlen(RidgeLineRidgeLineInterfaces0) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, RidgeLineRidgeLineInterfaces0);
  return res;
}



#endif
