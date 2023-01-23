#include <ttkScalarFieldCriticalPoints2.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkScalarFieldCriticalPoints2);

ttkScalarFieldCriticalPoints2::ttkScalarFieldCriticalPoints2() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkScalarFieldCriticalPoints2::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkScalarFieldCriticalPoints2::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    return 1;
  }
  return 0;
}

int ttkScalarFieldCriticalPoints2::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  auto inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  auto inputOrder = this->GetOrderArray(inputDataSet, 0);
  if(!inputOrder)
    return 0;

  auto ascManifold = inputDataSet->GetPointData()->GetArray("AscendingManifold");
  auto desManifold = inputDataSet->GetPointData()->GetArray("DescendingManifold");
  if(!ascManifold || !desManifold) {
    this->printErr("Unable to retrieve ascending or descending manifold.");
    return 0;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  std::vector<ttk::ScalarFieldCriticalPoints2::CriticalPoint> criticalPoints;
  int status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (
      status = this->computeCritialPoints<T0>(
        criticalPoints,
        ttkUtils::GetPointer<ttk::SimplexId>(inputOrder),
        ttkUtils::GetPointer<ttk::SimplexId>(desManifold),
        ttkUtils::GetPointer<ttk::SimplexId>(ascManifold),
        static_cast<const T0 *>(triangulation->getData())
      )
    )
  );

  if(status != 1)
    return 0;

  // Components Output
  {
    const size_t nCriticalPoints = criticalPoints.size();
    auto output = vtkPolyData::GetData(outputVector, 0);

    auto typeArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
    typeArray->SetName("Type");
    typeArray->SetNumberOfTuples(nCriticalPoints);

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToFloat();
    points->SetNumberOfPoints(nCriticalPoints);

    auto connectivityArray = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivityArray->SetNumberOfTuples(nCriticalPoints);

    auto offsetArray = vtkSmartPointer<vtkIdTypeArray>::New();
    offsetArray->SetNumberOfTuples(nCriticalPoints + 1);

    ttkTypeMacroT(
      triangulation->getType(),
      (
        status = this->computeCellAndPointArray<vtkIdType, T0>(
          criticalPoints,
          ttkUtils::GetPointer<float>(points->GetData()),
          ttkUtils::GetPointer<unsigned char>(typeArray),
          ttkUtils::GetPointer<vtkIdType>(offsetArray),
          ttkUtils::GetPointer<vtkIdType>(connectivityArray),
          static_cast<const T0 *>(triangulation->getData())
        )
      )
    );

    if(status != 1)
      return 0;

    output->SetPoints(points);
    auto pd = output->GetPointData();
    // pd->AddArray(sizeArray);
    pd->AddArray(typeArray);
    {
      size_t nPdArrays = inputDataSet->GetPointData()->GetNumberOfArrays();
      for(size_t a = 0; a < nPdArrays; a++) {
        auto array = inputDataSet->GetPointData()->GetArray(a);
        auto oArray = vtkSmartPointer<vtkDataArray>::Take(array->NewInstance());
        oArray->SetName(array->GetName());
        oArray->SetNumberOfValues(nCriticalPoints);
        // TODO: handle non 1-component arrays
        for(size_t p = 0; p < nCriticalPoints; p++) {
          auto cp = criticalPoints[p];
          oArray->SetTuple(p, cp.idx, array);
        }
        output->GetPointData()->AddArray(oArray);
      }
    }
    // Copy Field Data
    output->GetFieldData()->ShallowCopy(inputDataSet->GetFieldData());
  }

  // return success
  return 1;
}
