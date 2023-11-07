#include <ttkMacros.h>
#include <ttkConnectedComponentsPC.h>
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

vtkStandardNewMacro(ttkConnectedComponentsPC);

ttkConnectedComponentsPC::ttkConnectedComponentsPC() {
  this->setDebugMsgPrefix("ConnectedComponentsPC");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkConnectedComponentsPC::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkConnectedComponentsPC::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}



int ttkConnectedComponentsPC::RequestData(vtkInformation *ttkNotUsed(request),
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 0);

  if(!input)
    return !this->printErr("Input pointer is NULL.");

  if(input->GetNumberOfPoints() == 0)
    return !this->printErr("Input has no point.");

  if(!outputSegmentation)
    return !this->printErr("Output pointers are NULL.");

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);

  if(triangulation == nullptr)
    return !this->printErr("Triangulation is null");

  this->preconditionTriangulation(triangulation);

  const auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

  if(inputScalars == nullptr)
    return !this->printErr("No input scalars");


  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "'...");

  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();

  if(!numberOfVertices)
    return !this->printErr("Input has no vertices.");

  vtkNew<ttkSimplexIdTypeArray> segmentation{};

  if(!segmentation)
    return !this->printErr("Segmentation vtkDataArray allocation problem.");

  segmentation->SetNumberOfComponents(1);
  segmentation->SetNumberOfTuples(numberOfVertices);
  //segmentation->SetName(inputScalars->GetName() + "_Segmentation");
  segmentation->SetName("Segmentation");



  int ret{};
#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    this->printMsg("Running with MPI in vtk layer");
    auto pointData = input->GetPointData();
    auto globalIds = pointData->GetArray("GlobalPointIds");

ttkTemplateMacro(
    triangulation->getType(),
    (ret = this->execute<TTK_TT>(ttkUtils::GetPointer<SimplexId>(segmentation),this->IsoValue, ttkUtils::GetPointer<SimplexId>(inputScalars),
                            *static_cast<TTK_TT *>(triangulation->getData()), ttkUtils::GetPointer<ttk::SimplexId>(globalIds))));
  } else {
ttkTemplateMacro(
    triangulation->getType(),
    (ret = this->execute<TTK_TT>(ttkUtils::GetPointer<SimplexId>(segmentation),this->IsoValue, ttkUtils::GetPointer<SimplexId>(inputScalars),
                            *static_cast<TTK_TT *>(triangulation->getData()))));  }
#else
ttkTemplateMacro(
    triangulation->getType(),
    (ret = this->execute<TTK_TT>(ttkUtils::GetPointer<SimplexId>(segmentation),this->IsoValue, ttkUtils::GetPointer<SimplexId>(inputScalars),
                            *static_cast<TTK_TT *>(triangulation->getData()))));#endif // TTK_ENABLE_MPI


#endif
  if(ret != 0)
    return -1;

  outputSegmentation->ShallowCopy(input);

  vtkPointData *pointData = outputSegmentation->GetPointData();

  if(!pointData)
    return !this->printErr("outputSegmentation has no point data.");

  pointData->AddArray(segmentation);


  return 1;
}