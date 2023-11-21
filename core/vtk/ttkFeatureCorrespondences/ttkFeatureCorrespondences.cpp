#include <ttkFeatureCorrespondences.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkFeatureCorrespondences);

ttkFeatureCorrespondences::ttkFeatureCorrespondences() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkFeatureCorrespondences::~ttkFeatureCorrespondences() {
}

int ttkFeatureCorrespondences::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkFeatureCorrespondences::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkFeatureCorrespondences::RequestData(vtkInformation *ttkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {

  auto iImage = vtkImageData::GetData(inputVector[0]);
  if(!iImage)
    return !this->printErr("Unable to retrieve input data object.");

  int dim[3];
  iImage->GetDimensions(dim);

  auto iMatrix = this->GetInputArrayToProcess(0, inputVector);
  if(!iMatrix)
    return !this->printErr("Unable to retrieve input matrix.");

  if(this->GetInputArrayAssociation(0, inputVector) != 0)
    return !this->printErr("Input matrix needs to be a point data array.");

  if(iMatrix->GetNumberOfComponents() != 1)
    return !this->printErr("Input matrix needs to be a scalar array.");

  auto oMatrix = vtkSmartPointer<vtkDataArray>::Take(iMatrix->NewInstance());
  oMatrix->DeepCopy(iMatrix);

  int status = 0;

  // compute output matrix
  switch(this->OptimizationMethod) {
    case OPTIMIZATION_METHOD::N_SMALLEST_CORRESPONDENCES_PER_FEATURE:
    case OPTIMIZATION_METHOD::N_LARGEST_CORRESPONDENCES_PER_FEATURE: {
      ttkTypeMacroA(
        iMatrix->GetDataType(),
        (status = this->sortAndReduceCorrespondencesPerFeature<T0>(
           ttkUtils::GetPointer<T0>(oMatrix), ttkUtils::GetPointer<T0>(iMatrix),
           dim[0], dim[1], this->NumberOfLargestCorrespondencesPerFeature,
           this->OptimizationMethod
             == OPTIMIZATION_METHOD::N_SMALLEST_CORRESPONDENCES_PER_FEATURE)));
      break;
    }
    case OPTIMIZATION_METHOD::THRESHOLD_ABOVE: {
      ttkTypeMacroA(
        iMatrix->GetDataType(),
        (status = this->mapEachElement<T0>(
           ttkUtils::GetPointer<T0>(oMatrix), ttkUtils::GetPointer<T0>(iMatrix),
           dim[0], dim[1],
           [=](const T0 &v) { return v >= this->Threshold ? v : 0; })));
      break;
    }
    case OPTIMIZATION_METHOD::THRESHOLD_BELOW: {
      ttkTypeMacroA(
        iMatrix->GetDataType(),
        (status = this->mapEachElement<T0>(
           ttkUtils::GetPointer<T0>(oMatrix), ttkUtils::GetPointer<T0>(iMatrix),
           dim[0], dim[1],
           [=](const T0 &v) { return v <= this->Threshold ? v : 0; })));
      break;
    }
    default: {
      return !this->printErr("Unsupported Optimization Method");
    }
  }
  if(!status)
    return 0;

  auto oImage = vtkImageData::GetData(outputVector, 0);
  oImage->ShallowCopy(iImage);
  oImage->GetPointData()->AddArray(oMatrix);

  return 1;
}