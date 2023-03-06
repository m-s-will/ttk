#include <ttkScalarFieldCriticalPoints2.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkMultiBlockDataSet.h>
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

  auto inputOrder = this->GetOrderArray(inputDataSet, 0);
  auto ascManifold
    = inputDataSet->GetPointData()->GetArray(ttk::MorseSmaleAscendingName);
  auto desManifold
    = inputDataSet->GetPointData()->GetArray(ttk::MorseSmaleDescendingName);

  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(!inputOrder)
    return 0;
  if(!ascManifold || !desManifold) {
    this->printErr("Unable to retrieve ascending or descending manifold.");
    return 0;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  // ttkTypeMacroT(
  //   triangulation->getType(),
  //   (
  //     this->computeLookUpTable(
  //       static_cast<const T0 *>(triangulation->getData())
  //     )
  //   )
  // );
  // return 1;

   // type -> thread -> cp
  std::array<std::vector<std::vector<ttk::SimplexId>>,4> criticalPoints;

  int status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (
      status = this->computeCritialPoints<T0>(
        criticalPoints,
        ttkUtils::GetPointer<ttk::SimplexId>(inputOrder),
        ttkUtils::GetPointer<ttk::SimplexId>(ascManifold),
        ttkUtils::GetPointer<ttk::SimplexId>(desManifold),
        static_cast<const T0 *>(triangulation->getData())
      )
    )
  );

  if(status != 1)
    return 0;

  // generate VTK output
  {
    ttk::Timer timer;

    const std::string msg{"Generating VTK Output"};
    this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

    // auto output = vtkMultiBlockDataSet::GetData(outputVector, 0);
    auto output = vtkPolyData::GetData(outputVector, 0);
    auto of = output->GetFieldData();
    for(int b=0; b<4; b++){
      auto ttkVertexScalarFieldArray = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
      ttkVertexScalarFieldArray->SetName(("cp"+std::to_string(b)+"id").data());
      of->AddArray(ttkVertexScalarFieldArray);
    }

    #pragma omp parallel for schedule(static,1) num_threads(this->threadNumber_)
    for(int b=0; b<4; b++){
      const auto& cp = criticalPoints[b];
      const size_t nVectors = cp.size();
      size_t nCriticalPoints = 0;
      for(size_t i=0; i<nVectors; i++)
        nCriticalPoints += cp[i].size();

      auto array = of->GetArray(("cp"+std::to_string(b)+"id").data());
      array->SetNumberOfTuples(nCriticalPoints);

      this->computeIdArray(
        ttkUtils::GetPointer<ttk::SimplexId>(array),
        cp
      );
    }

    this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
  }

  return 1;
}
