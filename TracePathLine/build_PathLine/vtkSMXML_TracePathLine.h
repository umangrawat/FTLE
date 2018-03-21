// Loadable modules
//
// Generated by /export/home/urawat/Downloads/paraview_build/bin/vtkkwProcessXML-pv5.3
//
#ifndef __vtkSMXML_TracePathLine_h
#define __vtkSMXML_TracePathLine_h

#include <string.h>


// From file /export/home/urawat/Desktop/TracePathLine/TracePathLine.xml
static const char* const TracePathLineTracePathLineInterfaces0 =
"<ServerManagerConfiguration>\n"
"    <ProxyGroup name=\"filters\">\n"
"        <SourceProxy\n"
"            name=\"TracePathLine\"\n"
"            class=\"vtkTracePathLine\"\n"
"            label=\"Trace PathLine\">	<!-- name as it will appear in ParaView under \"Filters\" menu -->\n"
"\n"
"            <Documentation\n"
"                long_help=\"Traces the pathline end points to cells\"\n"
"                short_help=\"Traces the pathline end points to cells\">\n"
"                Locates pathline end points to cells.\n"
"            </Documentation>\n"
"\n"
"            <InputProperty\n"
"                    name=\"Input\"\n"
"                    command=\"AddInputConnection\"\n"
"                    clean_command=\"RemoveAllInputs\">\n"
"                <ProxyGroupDomain name=\"groups\">\n"
"                        <Group name=\"sources\"/> <!-- the input might come either from a reader-->\n"
"                        <Group name=\"filters\"/> <!-- ... or from another filter -->\n"
"                </ProxyGroupDomain>\n"
"                <DataTypeDomain name=\"input_type\">\n"
"                     <DataType value=\"vtkImageData\"/> <!-- We require that the data has type vtkDataSet; might as well be more specific (e.g., vtkRectilinearGrid) -->\n"
"                </DataTypeDomain>\n"
"            </InputProperty>\n"
"\n"
"            <DoubleVectorProperty\n"
"                    name=\"Seed grid origin, relative to input\"\n"
"                    command=\"Setorigin\"\n"
"                    number_of_elements=\"3\"\n"
"                    default_values=\"0 0 0\">\n"
"                <Documentation>\n"
"                    Origin of seeding grid.\n"
"                </Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <DoubleVectorProperty\n"
"                    name=\"Seed grid physical extent\"\n"
"                    command=\"Setbounds\"\n"
"                    number_of_elements=\"3\"\n"
"                    default_values=\"2.0 1.0 0.0\">\n"
"                <Documentation>\n"
"                    Size of seeding grid.\n"
"                </Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <IntVectorProperty\n"
"                    name=\"Seed grid resolution\"\n"
"                    command=\"Setdimension\"\n"
"                    number_of_elements=\"3\"\n"
"                    default_values=\"20 10 2\">\n"
"                <Documentation>\n"
"                    Resolution of seeding grid.\n"
"                </Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"        </SourceProxy>\n"
"    </ProxyGroup>\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* TracePathLineTracePathLineInterfaces()
{
  size_t len = ( 0
    + strlen(TracePathLineTracePathLineInterfaces0) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, TracePathLineTracePathLineInterfaces0);
  return res;
}



#endif