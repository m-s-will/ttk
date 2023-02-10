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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
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

  // type -> thread -> cp
  std::vector<std::vector<std::vector<CriticalPoint>>> criticalPoints(4);

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

  // generate VTK output
  {
    ttk::Timer timer;

    const std::string msg{"Generating VTK Output"};
    this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

    auto output = vtkMultiBlockDataSet::GetData(outputVector, 0);

    #pragma omp parallel for schedule(static,1) num_threads(this->threadNumber_)
    for(int b=0; b<4; b++){
      auto block = vtkSmartPointer<vtkPolyData>::New();
      const auto& cp = criticalPoints[b];
      const size_t nVectors = cp.size();
      size_t nCriticalPoints = 0;
      for(size_t i=0; i<nVectors; i++)
        nCriticalPoints += cp[i].size();

      auto ttkVertexScalarFieldArray = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
      ttkVertexScalarFieldArray->SetName("ttkVertexScalarField");
      ttkVertexScalarFieldArray->SetNumberOfTuples(nCriticalPoints);

      // auto outputOrder = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
      // outputOrder->SetName(inputOrder->GetName());
      // outputOrder->SetNumberOfTuples(nCriticalPoints);

      // points and cells
      {
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
            this->computeCellAndPointArray<vtkIdType, T0>(
              ttkUtils::GetPointer<float>(points->GetData()),
              ttkUtils::GetPointer<ttk::SimplexId>(ttkVertexScalarFieldArray),
              // ttkUtils::GetPointer<ttk::SimplexId>(outputOrder),
              ttkUtils::GetPointer<vtkIdType>(offsetArray),
              ttkUtils::GetPointer<vtkIdType>(connectivityArray),
              cp,
              // ttkUtils::GetPointer<ttk::SimplexId>(inputOrder),
              static_cast<const T0 *>(triangulation->getData())
            )
          )
        );

        block->SetPoints(points);
        auto cellArray = vtkSmartPointer<vtkCellArray>::New();
        cellArray->SetData(offsetArray, connectivityArray);
        block->SetVerts(cellArray);
      }

      // Point and Field Data
      {
        auto pd = block->GetPointData();
        pd->AddArray(ttkVertexScalarFieldArray);
        // pd->AddArray(outputOrder);

        // const ttk::SimplexId* ids = ttkUtils::GetPointer<ttk::SimplexId>(ttkVertexScalarFieldArray);
        // size_t nPdArrays = inputDataSet->GetPointData()->GetNumberOfArrays();
        // for(size_t a = 0; a < nPdArrays; a++) {
        //   auto array = inputDataSet->GetPointData()->GetArray(a);
        //   auto oArray = vtkSmartPointer<vtkDataArray>::Take(array->NewInstance());
        //   oArray->SetName(array->GetName());
        //   oArray->SetNumberOfValues(nCriticalPoints);
        //   for(size_t p = 0; p < nCriticalPoints; p++) {
        //     oArray->SetTuple(p, ids[p], array);
        //   }
        //   pd->AddArray(oArray);
        // }
        // Copy Field Data
        block->GetFieldData()->ShallowCopy(inputDataSet->GetFieldData());
      }
      output->SetBlock(b,block);
    }

    this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
  }

  return 1;
}
