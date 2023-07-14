#include <ttkMacros.h>
#include <ttkSimplifiedMT.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkUnstructuredGrid.h>

namespace {
  template <typename vtkArrayType, typename vectorType>
  void setArray(vtkArrayType &vtkArray, vectorType &vector) {
    ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
  }
} // namespace

vtkStandardNewMacro(ttkSimplifiedMT);

ttkSimplifiedMT::ttkSimplifiedMT() {
  this->setDebugMsgPrefix("SimplifiedMT");
  this->SetNumberOfInputPorts(2);
  //this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkSimplifiedMT::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  // domain
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  // persistence diagram
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkSimplifiedMT::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkSimplifiedMT::dispatch(vtkDataArray *const inputScalars,
                                    vtkPolyData *const outputSeparators,
                                    vtkDataArray *const persistenceScalars,
                                    const triangulationType &triangulation) {

  const auto scalars = ttkUtils::GetPointer<scalarType>(inputScalars);
  this->printMsg("Getting persistence scalars");
  this->printMsg("inputscalars type:" + std::to_string(inputScalars->GetDataType()));
  this->printMsg("persistenceScalars type:" + std::to_string(persistenceScalars->GetDataType()));
  const auto persScalars = ttkUtils::GetPointer<scalarType>(persistenceScalars);
  ttk::SimplexId numberOfPersistent = persistenceScalars->GetNumberOfTuples();
  this->printMsg("Got persistence scalars");

  const int dim = triangulation.getDimensionality();

  output_points_.clear();
  output_cells_labels_.clear();
  output_cells_connectivity_.clear();

  this->printMsg("Starting execeute");
  const int status
    = this->execute<scalarType, triangulationType>(scalars, persScalars, numberOfPersistent, triangulation);

  if(status != 0)
    return !this->printErr("SimplifiedMT.execute() error");

  vtkNew<vtkFloatArray> pointsCoords{};
  pointsCoords->SetNumberOfComponents(3);
  setArray(pointsCoords, output_points_);

  vtkNew<ttkSimplexIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(output_numberOfCells_ + 1);
  connectivity->SetNumberOfComponents(1);
  setArray(connectivity, output_cells_connectivity_);

  vtkNew<vtkUnsignedLongLongArray> hashArr{};
  hashArr->SetNumberOfComponents(1);
  hashArr->SetName("Hash");
  setArray(hashArr, output_cells_labels_);

  if(dim == 2 || dim == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < output_numberOfCells_ + 1; ++i) {
      offsets->SetTuple1(i, dim * i);
    }
  }

  vtkNew<vtkPoints> points{};
  points->SetData(pointsCoords);
  outputSeparators->SetPoints(points);

  vtkNew<vtkCellArray> cells{};
#ifndef TTK_ENABLE_64BIT_IDS
  cells->Use32BitStorage();
#endif // TTK_ENABLE_64BIT_IDS
  cells->SetData(offsets, connectivity);
  if(dim == 3) {
    outputSeparators->SetPolys(cells);
  } else {
    outputSeparators->SetLines(cells);
  }

  auto cellData = outputSeparators->GetCellData();
  cellData->AddArray(hashArr);

  return 1;
}

int ttkSimplifiedMT::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  const auto persistence = vtkDataSet::GetData(inputVector[1]);
  auto outputSeparators = vtkPolyData::GetData(outputVector, 0);

  if(!input)
    return !this->printErr("Input pointer is NULL.");

  if(input->GetNumberOfPoints() == 0)
    return !this->printErr("Input has no point.");

  if(!outputSeparators)
    return !this->printErr("Output pointers are NULL.");

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(triangulation == nullptr)
    return !this->printErr("Triangulation is null");


  const auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

  if(inputScalars == nullptr)
    return !this->printErr("wrong scalars.");

  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "'...");
  for (size_t i = 0; i < persistence->GetPointData()->GetNumberOfArrays(); i++) {
    this->printMsg("Array " + std::to_string(i) + " name: " + persistence->GetPointData()->GetArray(i)->GetName());
  }
  const auto persistenceScalars = persistence->GetPointData()->GetArray('ttkVertexScalarField');
  if(persistenceScalars == nullptr)
    return !this->printErr("no persistence scalars.");

  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();

  if(!numberOfVertices)
    return !this->printErr("Input has no vertices.");

  int status{};

  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (status = dispatch<VTK_TT, TTK_TT>(
                         inputScalars, outputSeparators, persistenceScalars,
                         *static_cast<TTK_TT *>(triangulation->getData()))));

  return status;
}
