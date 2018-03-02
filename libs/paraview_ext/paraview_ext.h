// helpers for Paraview
// filip sadlo
// cgl eth 2007

#ifndef _PARAVIEW_EXT_H_
#define _PARAVIEW_EXT_H_

#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"

using namespace std;

vtkFloatArray* castToVtkFloatArray(vtkDataArray* inArray);

// Basic
void passInterpolatePointData(vtkUnstructuredGrid *input,
                              vtkPointSet *output);

void passLineFieldVertexData(vtkStructuredGrid *input,
                             std::vector<int> *usedNodes,
                             std::vector<int> *usedNodesVertCnt,
                             vtkPointSet *output);

vtkUnstructuredGrid *generateUniformUSG(const char *name,
                                        float originX, float originY, float originZ,
                                        int cellsX, int cellsY, int cellsZ,
                                        float cellSize);

char *getLabelName(const char *widgetName);

bool arrayExists(vtkFloatArray *optionalScalar, vtkFloatArray *requiredVector);

#endif // _PARAVIEW_EXT_H_
