// Loadable modules
//
// Generated by /export/home/urawat/Downloads/paraview_build/bin/vtkkwProcessXML-pv5.3
//
#ifndef __vtkSMXML_ColorMap_h
#define __vtkSMXML_ColorMap_h

#include <string.h>


// From file /export/home/urawat/Desktop/ColorMap/ColorMap.xml
static const char* const ColorMapColorMapInterfaces0 =
"<ServerManagerConfiguration>\n"
"    <ProxyGroup name=\"filters\">\n"
"        <SourceProxy\n"
"                name=\"ColorMap\"\n"
"                class=\"vtkColorMap\"\n"
"                label=\"ColorMap\">	<!-- name as it will appear in ParaView under \"Filters\" menu -->\n"
"\n"
"            <Documentation\n"
"                    long_help=\"Maps the pathlines to respective cells\"\n"
"                    short_help=\"Maps the points\">\n"
"                Maps the points back to cells of origin.\n"
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
"                    name=\"Set Cell Id\"\n"
"                    command=\"SetcellId\"\n"
"                    number_of_elements=\"1\"\n"
"                    default_values=\"0.0\">\n"
"                <Documentation>Set the CellId for evaluation.</Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <OutputPort name=\"Points\" index=\"0\" id=\"port0\"/>\n"
"            <OutputPort name=\"Cells\" index=\"1\" id=\"port1\"/>\n"
"\n"
"        </SourceProxy>\n"
"    </ProxyGroup>\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* ColorMapColorMapInterfaces()
{
  size_t len = ( 0
    + strlen(ColorMapColorMapInterfaces0) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, ColorMapColorMapInterfaces0);
  return res;
}



#endif
