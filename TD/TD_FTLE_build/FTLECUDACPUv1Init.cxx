/* #undef PARAVIEW_USE_UNIFIED_BINDINGS */
#define NO_PYTHON_BINDINGS_AVAILABLE
#ifdef NO_PYTHON_BINDINGS_AVAILABLE
#undef PARAVIEW_USE_UNIFIED_BINDINGS
#endif
#ifdef PARAVIEW_USE_UNIFIED_BINDINGS
#include "vtkPython.h"
#include "vtkPythonInterpreter.h"
#endif

#include "vtkClientServerInterpreter.h"

#ifndef PARAVIEW_BUILD_SHARED_LIBS
#define PARAVIEW_BUILD_SHARED_LIBS
#endif
#if defined(PARAVIEW_BUILD_SHARED_LIBS) && defined(_WIN32)
# define VTK_WRAP_CS_EXPORT __declspec(dllexport)
#else
# define VTK_WRAP_CS_EXPORT
#endif

#ifdef PARAVIEW_USE_UNIFIED_BINDINGS
extern "C" void real_initFTLECUDACPUv1Python(const char *modulename);

void initFTLECUDACPUv1Python()
{
  static const char modulename[] = "FTLECUDACPUv1Python";
  real_initFTLECUDACPUv1Python(modulename);
}
#endif

extern void vtkFTLE_Init(vtkClientServerInterpreter* csi);


extern "C" void VTK_WRAP_CS_EXPORT FTLECUDACPUv1_Initialize(
  vtkClientServerInterpreter *csi)
{
#ifdef PARAVIEW_USE_UNIFIED_BINDINGS
  if (!vtkPythonInterpreter::IsInitialized())
    {
    vtkPythonInterpreter::Initialize();
    }

  static bool initialized = false;

  if (!initialized)
    {
    initialized = true;
    PyImport_AppendInittab("FTLECUDACPUv1Python", initFTLECUDACPUv1Python);
    }

  csi->Load("FTLECUDACPUv1");
#endif

  vtkFTLE_Init(csi);

}
