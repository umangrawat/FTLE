// Loadable modules
//
// Generated by /export/home/urawat/Downloads/paraview_build/bin/vtkkwProcessXML-pv5.3
//
#ifndef __vtkSMXML_FTLECUDACPUv5_h
#define __vtkSMXML_FTLECUDACPUv5_h

#include <string.h>


// From file /export/home/urawat/Desktop/TD/src/FTLE_Server.xml
static const char* const FTLECUDACPUv5FTLE_ServerInterfaces0 =
"<ServerManagerConfiguration>\n"
"    <ProxyGroup name=\"filters\">\n"
"        <!-- ================================================================== -->\n"
"        <SourceProxy name=\"FTLECUDA\" class=\"vtkFTLE\" label=\"FTLE CUDA\">\n"
"            <Documentation\n"
"                    long_help=\"Computes FTLE on the given patch of the dataset\"\n"
"                    short_help=\"Computes FTLE\">\n"
"                Computes Finite-Time Lyapunov Exponents on uniform grids.\n"
"            </Documentation>\n"
"\n"
"            <InputProperty\n"
"                    name=\"Input\"\n"
"                    command=\"SetInputConnection\">\n"
"                <!-- clean_command=\"RemoveAllInputs\"> -->\n"
"                <ProxyGroupDomain name=\"groups\">\n"
"                    <Group name=\"sources\"/>\n"
"                    <Group name=\"filters\"/>\n"
"                </ProxyGroupDomain>\n"
"                <DataTypeDomain name=\"input_type\">\n"
"                    <DataType value=\"vtkImageData\"/>\n"
"                </DataTypeDomain>\n"
"\n"
"                <InputArrayDomain attribute_type=\"point\"\n"
"                                  name=\"input_vectors\"\n"
"                                  number_of_components=\"3\"\n"
"                                  optional=\"1\"/>\n"
"\n"
"                <Documentation>\n"
"                    Compute Finite-Time Lyapunov exponent, which characterizes the rate of seperatation of massless\n"
"                    tracer particles within a flow.\n"
"                </Documentation>\n"
"            </InputProperty>\n"
"\n"
"            <DoubleVectorProperty name=\"Start time\"\n"
"                                  command=\"SettauInitial\"\n"
"                                  number_of_elements=\"1\"\n"
"                                  default_values=\"0.0\">\n"
"                <Documentation>Set initial timestep for the flow map computation</Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <DoubleVectorProperty name=\"Advection time\"\n"
"                                  command=\"SettauFinal\"\n"
"                                  number_of_elements=\"1\"\n"
"                                  default_values=\"15.0\">\n"
"                <Documentation>Set final timestep for the flow map computation</Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"\n"
"            <IntVectorProperty\n"
"                name=\"Maximum Number of Steps\"\n"
"                command=\"SetmaxNumberOfSteps\"\n"
"                number_of_elements=\"1\"\n"
"                default_values=\"10000\">\n"
"                <Documentation>The maximum number of steps for the integrator.</Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <DoubleVectorProperty\n"
"                    name=\"Seed grid origin, relative to input\"\n"
"                    command=\"SetoriginSeedGrid\"\n"
"                    number_of_elements=\"3\"\n"
"                    default_values=\"0 0 0\">\n"
"                <Documentation>\n"
"                    Origin of seeding grid.\n"
"                </Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <DoubleVectorProperty\n"
"                    name=\"Seed grid physical extent\"\n"
"                    command=\"SetboundsSeedGrid\"\n"
"                    number_of_elements=\"3\"\n"
"                    default_values=\"2.0 1.0 0.0\">\n"
"                <Documentation>\n"
"                    Size of seeding grid.\n"
"                </Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <DoubleVectorProperty\n"
"                    name=\"Seed grid resolution\"\n"
"                    command=\"SetdimensionSeedGrid\"\n"
"                    number_of_elements=\"3\"\n"
"                    default_values=\"20.0 10.0 2.0\">\n"
"                <Documentation>\n"
"                    Resolution of seeding grid.\n"
"                </Documentation>\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <IntVectorProperty name=\"UseGPU\"\n"
"                               command=\"SetuseGPU\"\n"
"                               number_of_elements=\"1\"\n"
"                               default_values=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>Use GPU for stream/path line integration, not yet implemented.</Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <IntVectorProperty name=\"ShowStreamLines(CPU only)\"\n"
"                               command=\"SetshowStreamLines\"\n"
"                               number_of_elements=\"1\"\n"
"                               default_values=\"0\">\n"
"                <BooleanDomain name=\"bool\"/>\n"
"                <Documentation>Return streamlines within the multiBlock</Documentation>\n"
"            </IntVectorProperty>\n"
"\n"
"            <DoubleVectorProperty name=\"Integration step size\"\n"
"                                  command=\"SetintegrationStepSize\"\n"
"                                  number_of_elements=\"1\"\n"
"                                  default_values=\"0.1\">\n"
"                <Documentation>Set step size for the fixed step size RK4 integration scheme.</Documentation>\n"
"\n"
"            </DoubleVectorProperty>\n"
"\n"
"            <OutputPort name=\"Ignore\" index=\"0\" id=\"port0\"/>\n"
"            <OutputPort name=\"FTLE Field\" index=\"1\" id=\"port1\"/>\n"
"            <OutputPort name=\"Trajectories(CPU ONLY)\" index=\"2\" id=\"port2\"/>\n"
"\n"
"        </SourceProxy>\n"
"    </ProxyGroup>\n"
"    <!-- End Filters Group -->\n"
"</ServerManagerConfiguration>\n"
"\n";
// Get single string
char* FTLECUDACPUv5FTLE_ServerInterfaces()
{
  size_t len = ( 0
    + strlen(FTLECUDACPUv5FTLE_ServerInterfaces0) );
  char* res = new char[ len + 1];
  res[0] = 0;
  strcat(res, FTLECUDACPUv5FTLE_ServerInterfaces0);
  return res;
}



#endif
