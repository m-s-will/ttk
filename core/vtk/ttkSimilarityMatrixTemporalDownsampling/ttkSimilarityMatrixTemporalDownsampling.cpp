#include <ttkSimilarityMatrixTemporalDownsampling.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityMatrixTemporalDownsampling);

ttkSimilarityMatrixTemporalDownsampling::
  ttkSimilarityMatrixTemporalDownsampling() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkSimilarityMatrixTemporalDownsampling::
  ~ttkSimilarityMatrixTemporalDownsampling() {
}

int ttkSimilarityMatrixTemporalDownsampling::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkSimilarityMatrixTemporalDownsampling::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkSimilarityMatrixTemporalDownsampling::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // get input as multiblock
  auto inputMB = vtkMultiBlockDataSet::GetData(inputVector[0]);
  int nTimesteps = inputMB->GetNumberOfBlocks();

  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  if(this->SamplingInterval < 1) {
    return !this->printErr("Sampling interval must be at least 1");
  } else if(this->SamplingInterval == 1) {
    outputMB->ShallowCopy(inputMB);
    return 1;
  }

  int interval = this->SamplingInterval;

  for(int i = 1; i < nTimesteps; i += interval) {
    if(i + (interval - 1) > nTimesteps)
      break;

    ttk::Timer timer;

    const std::string msg = "Multiplying " + std::to_string(interval)
                            + " matrices, ids: " + std::to_string(i - 1)
                            + " to " + std::to_string(i + (interval - 2));
    this->printMsg(
      msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

    auto prevMatrix = vtkImageData::SafeDownCast(inputMB->GetBlock(i - 1));
    int prevDims[3];
    int curDims[3];
    auto outMatrix = vtkSmartPointer<vtkImageData>::New();

    for(int t = i; t < i + (interval - 1); ++t) {
      auto curMatrix = vtkImageData::SafeDownCast(inputMB->GetBlock(t));

      auto prevData = (t == i) ? GetInputArrayToProcess(0, prevMatrix)
                               : prevMatrix->GetPointData()->GetArray(0);
      auto curData = GetInputArrayToProcess(0, curMatrix);
      if(!prevData || !curData)
        return !this->printErr(
          "Input data is missing array containing matrix data.");

      // Get dims
      prevMatrix->GetDimensions(prevDims);
      curMatrix->GetDimensions(curDims);

      outMatrix->SetDimensions(prevDims[0], curDims[1], 1);
      outMatrix->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
      auto matrixData = outMatrix->GetPointData()->GetArray(0);
      matrixData->SetName("Downscaled");

      int status = 0;
      ttkTypeMacroI(
        prevData->GetDataType(),
        status = this->multiplyMatrices<T0>(
          ttkUtils::GetPointer<unsigned char>(matrixData),
          ttkUtils::GetPointer<const T0>(prevData),
          ttkUtils::GetPointer<const T0>(curData), prevDims, curDims));
      if(!status)
        return 0;

      // Get the field data arrays with indices
      auto prevIds = GetInputArrayToProcess(2, prevMatrix);
      auto curIds = GetInputArrayToProcess(3, curMatrix);
      if(!prevIds || !curIds)
        return !this->printErr(
          "Input field data is missing arrays containing id data.");

      auto outFD = outMatrix->GetFieldData();
      outFD->AddArray(prevIds);
      outFD->AddArray(curIds);

      prevMatrix->DeepCopy(outMatrix);
    }

    auto timeArray = vtkSmartPointer<vtkIntArray>::New();
    timeArray->SetName("Timesteps");
    timeArray->SetNumberOfComponents(1);
    timeArray->SetNumberOfTuples(2);
    timeArray->SetTuple1(0, i - 1);
    timeArray->SetTuple1(1, i + (interval - 2));
    outMatrix->GetFieldData()->AddArray(timeArray);

    outputMB->SetBlock(outputMB->GetNumberOfBlocks(), outMatrix);

    this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
  }

  // return success
  return 1;
}
