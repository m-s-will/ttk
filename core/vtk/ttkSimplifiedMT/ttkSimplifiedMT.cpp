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
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
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
                                    vtkDataArray *const inputManifold,
                                    vtkPolyData *const outputSeparators,
                                    vtkDataArray *const persistenceScalars,
                                    const triangulationType &triangulation) {

  const auto scalars = ttkUtils::GetPointer<scalarType>(inputScalars);
  const auto manifold = ttkUtils::GetPointer<SimplexId>(inputManifold);
  const auto persScalars = ttkUtils::GetPointer<double>(persistenceScalars);
  ttk::SimplexId numberOfPersistent = persistenceScalars->GetNumberOfTuples();
  this->printMsg("Number of persistent extrema: " + std::to_string(numberOfPersistent));
  const int dim = triangulation.getDimensionality();

  output_points_.clear();
  output_cells_labels_.clear();
  output_cells_connectivity_.clear();

  this->printMsg("Starting execute");
  const int status
    = this->execute<scalarType, triangulationType>(scalars, manifold, persScalars, numberOfPersistent, triangulation);

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

  const auto input = vtkDataSet::GetData(inputVector[0], 0);
  const auto persistence = vtkDataSet::GetData(inputVector[1], 0);
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
  const auto inputManifold = this->GetInputArrayToProcess(1, inputVector);

  if(inputScalars == nullptr || inputManifold == nullptr)
    return !this->printErr("wrong scalars.");

  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "' with manifold `" + std::string(inputManifold->GetName()) + "'");
  for (int i = 0; i < persistence->GetPointData()->GetNumberOfArrays(); i++) {
    this->printMsg("Array " + std::to_string(i) + " name: " + persistence->GetPointData()->GetArray(i)->GetName());
  }

  const auto persistenceScalars = persistence->GetPointData()->GetArray(0); // scalar values of start vertices of persistence pairs
  const auto persistenceScalars2 = persistence->GetPointData()->GetArray(1); // scalar valuse of end vertices of persistence pairs
  vtkDoubleArray* persistenceScalarsCombined = vtkDoubleArray::New();
  persistenceScalarsCombined->SetName("persistenceScalarsCombined");
  persistenceScalarsCombined->SetNumberOfComponents(1);
  for (int i = 0; i < persistenceScalars->GetNumberOfTuples(); ++i)
  {
    double value = persistenceScalars->GetTuple1(i);
    persistenceScalarsCombined->InsertNextTypedTuple(&value);
  }
  for (int i = 0; i < persistenceScalars2->GetNumberOfTuples(); ++i)
  {
    double value = persistenceScalars2->GetTuple1(i);
    persistenceScalarsCombined->InsertNextTypedTuple(&value);
  }

  //if(persistenceScalars == nullptr)
  //  return !this‚->printErr("no persistence scalars.");
  if(persistenceScalarsCombined == nullptr)
    return !this->printErr("no persistence scalars.");

  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();

  if(!numberOfVertices)
    return !this->printErr("Input has no vertices.");

  int status{};

  ttkVtkTemplateMacro(inputScalars->GetDataType(), triangulation->getType(),
                      (status = dispatch<VTK_TT, TTK_TT>(
                         inputScalars, inputManifold, outputSeparators, persistenceScalarsCombined,
                         *static_cast<TTK_TT *>(triangulation->getData()))));


  return status;
}
