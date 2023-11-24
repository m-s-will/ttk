#include <ttkScalarFieldFromPoints.h>

#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkScalarFieldFromPoints);

ttkScalarFieldFromPoints::ttkScalarFieldFromPoints() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkScalarFieldFromPoints::~ttkScalarFieldFromPoints() {
}

int ttkScalarFieldFromPoints::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;

  return 1;
}

int ttkScalarFieldFromPoints::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  else
    return 0;

  return 1;
}

int ttkScalarFieldFromPoints::RequestInformation(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector) {

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Extent
  int wholeExtent[6]
    = {0, (int)this->Resolution[0] - 1, 0, (int)this->Resolution[1] - 1,
       0, (int)this->Resolution[2] - 1};
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent, 6);

  return 1;
}

int ttkScalarFieldFromPoints::RequestData(vtkInformation *,
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {
  // Get input
  auto input = vtkPointSet::GetData(inputVector[0], 0);
  if(!input) {
    this->printErr("There is no input provided.");
    return 0;
  }

  const size_t nPoints = input->GetNumberOfPoints();

  // Get point attributes arrays from input
  auto ampArray = GetInputArrayToProcess(2, input);
  if(!ampArray) {
    this->printErr("No amplitude array was provided.");
    return 0;
  }

  auto varArray = GetInputArrayToProcess(3, input);
  if(!varArray) {
    this->printErr("No variance array was provided.");
    return 0;
  }

  auto idArray = GetInputArrayToProcess(4, input);
  if(!idArray) {
    this->printErr("No point id array was provided.");
    // return 0;
  }

  auto output = vtkImageData::GetData(outputVector);
  output->SetDimensions(
    this->Resolution[0], this->Resolution[1], this->Resolution[2]);
  output->SetOrigin(
    this->ImageBounds[0], this->ImageBounds[2], this->ImageBounds[4]);

  const double spacing[3]{
    this->Resolution[0] > 1 ? (this->ImageBounds[1] - this->ImageBounds[0])
                                / (this->Resolution[0] - 1)
                            : 0,
    this->Resolution[1] > 1 ? (this->ImageBounds[3] - this->ImageBounds[2])
                                / (this->Resolution[1] - 1)
                            : 0,
    this->Resolution[2] > 1 ? (this->ImageBounds[5] - this->ImageBounds[4])
                                / (this->Resolution[2] - 1)
                            : 0};
  output->SetSpacing(spacing);

  output->AllocateScalars(VTK_DOUBLE, 1);

  if(this->Resolution[0] < 1 || this->Resolution[1] < 1
     || this->Resolution[2] < 1)
    return !this->printErr("Resolution contains zeros.");

  // Get data array to put the maximum mixture results in
  auto maxScalarArray = output->GetPointData()->GetArray(0);
  maxScalarArray->SetName("MaxScalars");

  auto nPixels = maxScalarArray->GetNumberOfTuples();

  // Create array to store ids of maximum in
  auto maxIdArray = vtkSmartPointer<vtkIntArray>::New();
  maxIdArray->SetName("MaxId");
  maxIdArray->SetNumberOfTuples(nPixels);
  maxIdArray->SetNumberOfComponents(1);

  // Create array to store additive mixture in
  auto addScalarArray = vtkSmartPointer<vtkDoubleArray>::New();
  addScalarArray->SetName("AddScalars");
  addScalarArray->SetNumberOfTuples(nPixels);
  addScalarArray->SetNumberOfComponents(1);

  // Used to check base layer execution status
  int status = 0;

  // Execute either 2D or 3D case for the chosen kernel
  switch(this->Kernel) {
    case 0: {
      ttkTypeMacroR(
        input->GetPoints()->GetDataType(),
        (status
         = this->computeScalarField3D<ScalarFieldFromPoints::Gaussian, T0>(
           ttkUtils::GetPointer<double>(maxScalarArray),
           ttkUtils::GetPointer<double>(addScalarArray),
           ttkUtils::GetPointer<int>(maxIdArray),
           ttkUtils::GetPointer<T0>(input->GetPoints()->GetData()),
           ttkUtils::GetPointer<int>(idArray),
           ttkUtils::GetPointer<double>(ampArray),
           ttkUtils::GetPointer<double>(varArray), this->ImageBounds, spacing,
           this->Resolution, nPoints, nPixels)));
      break;
    }
  }

  // On error cancel filter execution
  if(status == 0)
    return 0;

  output->GetPointData()->AddArray(addScalarArray);
  output->GetPointData()->AddArray(maxIdArray);

  return 1;
}
