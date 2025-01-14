#include <ttkHelloWorld.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <ImplicitTriangulation.h>

const std::array<ttk::SimplexId, 64 * 14 * 3> offsetsLUT{
  0,  -1, -1, 1,  -1, -1, 0,  0,  -1, 1,  0,  -1, 0,  -1, 0,  1,  -1, 0,  1,
  0,  0,  -1, 0,  1,  0,  0,  1,  -1, 0,  0,  -1, 1,  0,  0,  1,  0,  -1, 1,
  1,  0,  1,  1,  0,  -1, -1, 1,  -1, -1, 0,  0,  -1, 1,  0,  -1, 0,  -1, 0,
  1,  -1, 0,  1,  0,  0,  0,  1,  0,  0,  1,  1,  0,  0,  1,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0,  0,  -1, 1,  0,  0,  1,  0,  -1, 0,
  1,  0,  0,  1,  -1, 1,  1,  0,  1,  1,  0,  0,  -1, 0,  -1, -1, 0,  -1, 0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0,  0,
  -1, 1,  0,  0,  1,  0,  -1, 0,  1,  0,  0,  1,  -1, 1,  1,  0,  1,  1,  1,
  0,  0,  1,  0,  -1, 0,  0,  -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  1,  0,  0,  -1, 1,  0,  -1,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  -1, 0,  0,  -1, 1,  0,  0,  1,  0,  -1, 0,  1,  0,  0,
  1,  -1, 1,  1,  0,  1,  1,  0,  0,  -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, -1, 1,  -1, -1,
  0,  0,  -1, 1,  0,  -1, 0,  -1, 0,  1,  -1, 0,  1,  0,  0,  -1, 0,  0,  -1,
  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1,
  -1, 1,  -1, -1, 0,  0,  -1, 1,  0,  -1, 0,  -1, 0,  1,  -1, 0,  1,  0,  0,
  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  -1, 0,  -1, 0,  0,  -1, 0,  1,  0,  0,  1,  0,  -1, -1, 0,  0,
  -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  -1, 0,  0,  -1, 1,  0,  0,  1,  0,  -1, 0,  1,
  0,  0,  1,  -1, 1,  1,  0,  1,  1,  0,  -1, 0,  1,  -1, 0,  1,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0,  1,  -1, 0,  1,  0,
  0,  0,  0,  1,  0,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0,  0,  -1,
  1,  0,  0,  1,  0,  -1, 0,  1,  0,  0,  1,  -1, 1,  1,  0,  1,  1,  0,  -1,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  -1, 0,  0,  -1, 1,  0,  0,  1,  0,  -1, 0,  1,  0,  0,  1,
  -1, 1,  1,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,
  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0,  0,  -1, 1,  0,  0,
  1,  0,  -1, 0,  1,  0,  0,  1,  -1, 1,  1,  0,  1,  1,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  -1, 0,  -1, 0,  0,  -1, 0,  1,  0,  0,  1,  1,  -1, 0,  1,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  -1, 0,  1,  -1, 0,  1,  0,  0,  0,  0,  1,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 0,  -1, 0,  0,  -1, 0,  1,  0,
  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, -1, 1,  -1,
  -1, 0,  0,  -1, 1,  0,  -1, 0,  -1, 0,  1,  -1, 0,  1,  0,  0,  -1, 0,  0,
  -1, 1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  -1, -1, 1,  -1, -1, 0,  0,  -1, 1,  0,  -1, 0,  -1, 0,  1,  -1, 0,  1,  0,
  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  -1, -1, 0,  0,  -1, 1,  0,  0,  1,  0,  0,  -1, -1, 0,
  -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, -1, 0,  0,  -1, 1,
  0,  0,  1,  0,  1,  0,  -1, 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, 1,
  0,  -1, 1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  -1, -1, 0,  0,  -1, 1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, -1, 1,  -1, -1, 0,  0,  -1, 1,  0,
  -1, 0,  -1, 0,  1,  -1, 0,  1,  0,  0,  -1, 0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, -1, 1,  -1, -1, 0,
  0,  -1, 1,  0,  -1, 0,  -1, 0,  1,  -1, 0,  1,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  -1, -1,
  0,  0,  -1, 0,  -1, 0,  -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0};
std::array<ttk::SimplexId, 64 * 14> offsetsLUT2;

const std::array<ttk::SimplexId, 64> nNeighborsLUT{
  14, 10, 10, 0, 10, 6, 8, 0, 10, 8, 6,  0, 0, 0, 0, 0, 10, 6, 8, 0, 8, 4,
  7,  0,  6,  4, 4,  0, 0, 0, 0,  0, 10, 8, 6, 0, 6, 4, 4,  0, 8, 7, 4, 0,
  0,  0,  0,  0, 0,  0, 0, 0, 0,  0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0};

struct MyImplicitTriangulationBase {

  virtual ttk::SimplexId getNumberOfVertices() const {
    return 0;
  };

  virtual ttk::SimplexId
    getVertexNeighborNumber(const ttk::SimplexId idx) const {
    return 0;
  }

  virtual void getVertexNeighbor(const ttk::SimplexId idx,
                                 const ttk::SimplexId n,
                                 ttk::SimplexId &nIdx) const {
  }

  virtual void preconditionVertexNeighbors() {
  }
};

struct MyImplicitTriangulation : MyImplicitTriangulationBase {

  ttk::SimplexId dim[3];
  ttk::SimplexId dimM1[3];
  float dimM1F[3];
  ttk::SimplexId dimXY;

  void setDimension(int *dim_) {
    this->dim[0] = dim_[0];
    this->dim[1] = dim_[1];
    this->dim[2] = dim_[2];
    this->dimM1[0] = dim_[0] - 1;
    this->dimM1[1] = dim_[1] - 1;
    this->dimM1[2] = dim_[2] - 1;

    this->dimM1F[0] = (float)this->dimM1[0];
    this->dimM1F[1] = (float)this->dimM1[1];
    this->dimM1F[2] = (float)this->dimM1[2];

    this->dimXY = this->dim[0] * this->dim[1];
  }

  ttk::SimplexId getNumberOfVertices() const final {
    return this->dim[0] * this->dim[1] * this->dim[2];
  }

  ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId idx) const final {
    ttk::SimplexId xyz[3];
    ttk::SimplexId &x = xyz[0];
    ttk::SimplexId &y = xyz[1];
    ttk::SimplexId &z = xyz[2];
    xyz[2] = idx / this->dimXY;
    const auto idx2 = idx - (xyz[2] * this->dimXY);
    xyz[1] = idx2 / this->dim[0];
    xyz[0] = idx2 % this->dim[0];

    int key = (x == 0                ? 1
               : x == this->dimM1[0] ? 2
                                     : 0)
              + (y == 0                ? 4
                 : y == this->dimM1[1] ? 8
                                       : 0)
              + (z == 0                ? 16
                 : z == this->dimM1[2] ? 32
                                       : 0);

    return nNeighborsLUT[key];
  }

  void getVertexNeighbor(const ttk::SimplexId idx,
                         const ttk::SimplexId n,
                         ttk::SimplexId &nIdx) const final {
    ttk::SimplexId xyz[3];
    ttk::SimplexId &x = xyz[0];
    ttk::SimplexId &y = xyz[1];
    ttk::SimplexId &z = xyz[2];
    xyz[2] = idx / this->dimXY;
    const auto idx2 = idx - (xyz[2] * this->dimXY);
    xyz[1] = idx2 / this->dim[0];
    xyz[0] = idx2 % this->dim[0];

    int key = (x == 0                ? 1
               : x == this->dimM1[0] ? 2
                                     : 0)
              + (y == 0                ? 4
                 : y == this->dimM1[1] ? 8
                                       : 0)
              + (z == 0                ? 16
                 : z == this->dimM1[2] ? 32
                                       : 0);

    nIdx = idx + offsetsLUT2[key * 14 + n];
  }

  void preconditionVertexNeighbors() final {
    for(int i = 0, j = 0; i < 64 * 14; i++, j += 3) {
      const auto &dx = offsetsLUT[j + 0];
      const auto &dy = offsetsLUT[j + 1];
      const auto &dz = offsetsLUT[j + 2];

      offsetsLUT2[i] = dx + dy * this->dim[0] + dz * this->dimXY;
    }
  }
};

// struct MyImplicitTriangulation : MyImplicitTriangulationBase {

//   ttk::SimplexId dim[3];
//   ttk::SimplexId dimM1[3];
//   float dimM1F[3];
//   ttk::SimplexId dimXY;

//   void setDimension(int* dim_){
//     this->dim[0] = dim_[0];
//     this->dim[1] = dim_[1];
//     this->dim[2] = dim_[2];
//     this->dimM1[0] = dim_[0]-1;
//     this->dimM1[1] = dim_[1]-1;
//     this->dimM1[2] = dim_[2]-1;

//     this->dimM1F[0] = (float)this->dimM1[0];
//     this->dimM1F[1] = (float)this->dimM1[1];
//     this->dimM1F[2] = (float)this->dimM1[2];

//     this->dimXY = this->dim[0]*this->dim[1];
//   }

//   ttk::SimplexId getNumberOfVertices() const override {
//     return this->dim[0]*this->dim[1]*this->dim[2];
//   }

//   ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId idx) const
//   override {
//     ttk::SimplexId xyz[3];
//     ttk::SimplexId& x = xyz[0];
//     ttk::SimplexId& y = xyz[1];
//     ttk::SimplexId& z = xyz[2];
//     xyz[2] = idx / this->dimXY;
//     const auto idx2 = idx - (xyz[2] * this->dimXY);
//     xyz[1] = idx2 / this->dim[0];
//     xyz[0] = idx2 % this->dim[0];

//     int key =
//       (x==0?1:x==this->dimM1[0]?2:0)+
//       (y==0?4:y==this->dimM1[1]?8:0)+
//       (z==0?16:z==this->dimM1[2]?32:0)
//     ;

//     return nNeighborsLUT[key];
//   }

//   void getVertexNeighbor(const ttk::SimplexId idx, const ttk::SimplexId n,
//   ttk::SimplexId& nIdx) const override {
//     ttk::SimplexId xyz[3];
//     ttk::SimplexId& x = xyz[0];
//     ttk::SimplexId& y = xyz[1];
//     ttk::SimplexId& z = xyz[2];
//     xyz[2] = idx / this->dimXY;
//     const auto idx2 = idx - (xyz[2] * this->dimXY);
//     xyz[1] = idx2 / this->dim[0];
//     xyz[0] = idx2 % this->dim[0];

//     int key =
//       (x==0?1:x==this->dimM1[0]?2:0)+
//       (y==0?4:y==this->dimM1[1]?8:0)+
//       (z==0?16:z==this->dimM1[2]?32:0)
//     ;

//     nIdx = idx + offsetsLUT2[key*14+n];
//   }

//   void preconditionVertexNeighbors() override {
//     for(int i=0,j=0; i<64*14; i++,j+=3){
//       const auto& dx =  offsetsLUT[j+0];
//       const auto& dy =  offsetsLUT[j+1];
//       const auto& dz =  offsetsLUT[j+2];

//       offsetsLUT2[i] = dx + dy*this->dim[0] + dz*this->dimXY;
//     }
//   }
// };

// struct MyImplicitTriangulation : MyImplicitTriangulationBase {

//   ttk::SimplexId dim[3];
//   ttk::SimplexId dimM1[3];
//   float dimM1F[3];
//   ttk::SimplexId dimXY;

//   void setDimension(int* dim_){
//     this->dim[0] = dim_[0];
//     this->dim[1] = dim_[1];
//     this->dim[2] = dim_[2];
//     this->dimM1[0] = dim_[0]-1;
//     this->dimM1[1] = dim_[1]-1;
//     this->dimM1[2] = dim_[2]-1;

//     this->dimM1F[0] = (float)this->dimM1[0];
//     this->dimM1F[1] = (float)this->dimM1[1];
//     this->dimM1F[2] = (float)this->dimM1[2];

//     this->dimXY = this->dim[0]*this->dim[1];
//   }

//   ttk::SimplexId getNumberOfVertices() const override {
//     return this->dim[0]*this->dim[1]*this->dim[2];
//   }

//   ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId idx) const
//   override {
//     ttk::SimplexId xyz[3];
//     ttk::SimplexId& x = xyz[0];
//     ttk::SimplexId& y = xyz[1];
//     ttk::SimplexId& z = xyz[2];
//     xyz[2] = idx / this->dimXY;
//     const auto idx2 = idx - (xyz[2] * this->dimXY);
//     xyz[1] = idx2 / this->dim[0];
//     xyz[0] = idx2 % this->dim[0];

//     // build case key
//     // float xd = float(x)/this->dimM1F[0];
//     // float yd = float(y)/this->dimM1F[1];
//     // float zd = float(z)/this->dimM1F[2];
//     // int key = static_cast<int>(
//     //       (1.-ceil(xd))+
//     //   2.0*floor(xd)+
//     //   4.0*(1.-ceil(yd))+
//     //   8.0*floor(yd)+
//     //   16.0*(1.-ceil(zd))+
//     //   32.0*floor(zd)
//     // );
//     // float xd = float(x)/this->dimM1F[0];
//     // float yd = float(y)/this->dimM1F[1];
//     // float zd = float(z)/this->dimM1F[2];
//     // int key =
//     //     static_cast<int>(1.0-xd)+
//     //   2*static_cast<int>(xd)+
//     //   4*static_cast<int>(1.0-yd)+
//     //   8*static_cast<int>(yd)+
//     //   16*static_cast<int>(1.0-zd)+
//     //   32*static_cast<int>(zd)
//     // ;

//     int key =
//       (x==0?1:x==this->dimM1[0]?2:0)+
//       (y==0?4:y==this->dimM1[1]?8:0)+
//       (z==0?16:z==this->dimM1[2]?32:0)
//     ;

//     return nNeighborsLUT[key];
//   }

//   void getVertexNeighbor(const ttk::SimplexId idx, const ttk::SimplexId n,
//   ttk::SimplexId& nIdx) const override {
//     ttk::SimplexId xyz[3];
//     ttk::SimplexId& x = xyz[0];
//     ttk::SimplexId& y = xyz[1];
//     ttk::SimplexId& z = xyz[2];
//     xyz[2] = idx / this->dimXY;
//     const auto idx2 = idx - (xyz[2] * this->dimXY);
//     xyz[1] = idx2 / this->dim[0];
//     xyz[0] = idx2 % this->dim[0];

//     int key =
//       (x==0?1:x==this->dimM1[0]?2:0)+
//       (y==0?4:y==this->dimM1[1]?8:0)+
//       (z==0?16:z==this->dimM1[2]?32:0)
//     ;
//     // float xd = float(x)/this->dimM1F[0];
//     // float yd = float(y)/this->dimM1F[1];
//     // float zd = float(z)/this->dimM1F[2];
//     // int key = static_cast<int>(
//     //       (1.-ceil(xd))+
//     //   2.0*floor(xd)+
//     //   4.0*(1.-ceil(yd))+
//     //   8.0*floor(yd)+
//     //   16.0*(1.-ceil(zd))+
//     //   32.0*floor(zd)
//     // );
//     // float xd = float(x)/this->dimM1F[0];
//     // float yd = float(y)/this->dimM1F[1];
//     // float zd = float(z)/this->dimM1F[2];
//     // int key =
//     //     static_cast<int>(1.0-xd)+
//     //   2*static_cast<int>(xd)+
//     //   4*static_cast<int>(1.0-yd)+
//     //   8*static_cast<int>(yd)+
//     //   16*static_cast<int>(1.0-zd)+
//     //   32*static_cast<int>(zd)
//     ;

//     nIdx = idx + offsetsLUT2[key*14+n];
//   }

//   void preconditionVertexNeighbors() override {
//     for(int i=0,j=0; i<64*14; i++,j+=3){
//       const auto& dx =  offsetsLUT[j+0];
//       const auto& dy =  offsetsLUT[j+1];
//       const auto& dz =  offsetsLUT[j+2];

//       offsetsLUT2[i] = dx + dy*this->dim[0] + dz*this->dimXY;
//     }
//   }
// };

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkHelloWorld);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkHelloWorld::ttkHelloWorld() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkHelloWorld::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkHelloWorld::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkHelloWorld::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  // Get input array that will be processed
  //
  // Note: VTK provides abstract functionality to handle array selections, but
  //       this essential functionality is unfortunately not well documented.
  //       Before you read further, please keep in mind the the TTK developer
  //       team is not responsible for the existing VTK Api ;-)
  //
  //       In a nutshell, prior to the RequestData execution one has to call
  //
  //           SetInputArrayToProcess (
  //               int idx,
  //               int port,
  //               int connection,
  //               int fieldAssociation,
  //               const char *name
  //            )
  //
  //       The parameter 'idx' is often misunderstood: lets say the filter
  //       requires n arrays, then idx enumerates them from 0 to n-1.
  //
  //       The 'port' is the input port index at which the object is connected
  //       from which we want to get the array.
  //
  //       The 'connection' is the connection index at that port (we have to
  //       specify this because VTK allows multiple connections at the same
  //       input port).
  //
  //       The 'fieldAssociation' integer specifies if the array should be taken
  //       from 0: point data, 1: cell data, or 2: field data.
  //
  //       The final parameter is the 'name' of the array.
  //
  //       Example: SetInputArrayToProcess(3,1,0,1,"EdgeLength") will store that
  //                for the 3rd array the filter needs the cell data array named
  //                "EdgeLength" that it will retrieve from the vtkDataObject
  //                at input port 1 (first connection). During the RequestData
  //                method one can then actually retrieve the 3rd array it
  //                requires for its computation by calling
  //                GetInputArrayToProcess(3, inputVector)
  //
  //       If this filter is run within ParaView, then the UI will automatically
  //       call SetInputArrayToProcess (see HelloWorld.xml file).
  //
  //       During the RequestData execution one can then retrieve an actual
  //       array with the method "GetInputArrayToProcess".
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

  // Create an output array that has the same data type as the input array
  // Note: vtkSmartPointers are well documented
  //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  vtkSmartPointer<vtkDataArray> const outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->SetName(this->OutputArrayName.data()); // set array name
  outputArray->SetNumberOfComponents(1); // only one component per tuple
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class

  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  // MyImplicitTriangulation<1025,1025,1025,1024,1024,1024,1025*1025> trian;
  MyImplicitTriangulation trian;
  int dim[3];
  ((vtkImageData *)inputDataSet)->GetDimensions(dim);
  trian.setDimension(dim);
  trian.preconditionVertexNeighbors();
  ttkVtkTemplateMacro(
    inputArray->GetDataType(), triangulation->getType(),
    (status = this->computeAverages<VTK_TT, MyImplicitTriangulation>(
       (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),
       (VTK_TT *)ttkUtils::GetVoidPointer(inputArray), &trian)));
  // int status = 0; // this integer checks if the base code returns an error
  ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
                      (status = this->computeAverages<VTK_TT, TTK_TT>(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),
                         (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         (TTK_TT *)triangulation->getData())));
  ttkVtkTemplateMacro(
    inputArray->GetDataType(), triangulation->getType(),
    (status = this->computeAverages<VTK_TT, ttk::ImplicitTriangulation>(
       (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),
       (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
       (ttk::ImplicitTriangulation *)triangulation->getData())));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}
