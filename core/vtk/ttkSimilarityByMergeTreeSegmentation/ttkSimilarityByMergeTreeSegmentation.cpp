#include <ttkSimilarityByMergeTreeSegmentation.h>

#include <vtkObjectFactory.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkIntArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByMergeTreeSegmentation);

ttkSimilarityByMergeTreeSegmentation::ttkSimilarityByMergeTreeSegmentation() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkSimilarityByMergeTreeSegmentation::~ttkSimilarityByMergeTreeSegmentation() {
}

int ttkSimilarityByMergeTreeSegmentation::ComputeSimilarityMatrix(
  vtkImageData *similarityMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {

  // unpack inputs
  auto inputsAsMB0 = static_cast<vtkMultiBlockDataSet *>(inputDataObjects0);
  auto inputsAsMB1 = static_cast<vtkMultiBlockDataSet *>(inputDataObjects1);

  auto d0 = vtkDataSet::SafeDownCast(inputsAsMB0->GetBlock(0));
  auto d1 = vtkDataSet::SafeDownCast(inputsAsMB1->GetBlock(0));

  auto m0 = vtkDataSet::SafeDownCast(inputsAsMB0->GetBlock(1));
  auto m1 = vtkDataSet::SafeDownCast(inputsAsMB1->GetBlock(1));

  const int nNodes0 = m0->GetNumberOfPoints();
  const int nNodes1 = m1->GetNumberOfPoints();

  // extract arrays
  auto seg0 = vtkIntArray::SafeDownCast(d0->GetPointData()->GetArray("NodeId"));
  auto seg1 = vtkIntArray::SafeDownCast(d1->GetPointData()->GetArray("NodeId"));
  if(!seg0 || !seg1)
    return !this->printErr(
      "Unable to retrieve `NodeId` arrays from segmentations.");

  auto next0
    = vtkIntArray::SafeDownCast(m0->GetPointData()->GetArray("NextId"));
  auto next1
    = vtkIntArray::SafeDownCast(m1->GetPointData()->GetArray("NextId"));
  if(!next0 || !next1)
    return !this->printErr(
      "Unable to retrieve `NextId` arrays from merge trees.");

  auto scalars0 = this->GetInputArrayToProcess(0, m0);
  auto scalars1 = this->GetInputArrayToProcess(0, m1);
  if(!scalars0 || !scalars1)
    return !this->printErr("Unable to retrieve merge tree scalar arrays.");

  // initialize similarity matrix
  similarityMatrix->SetDimensions(nNodes0, nNodes1, 1);
  similarityMatrix->AllocateScalars(VTK_INT, 1);
  auto matrixData = similarityMatrix->GetPointData()->GetArray(0);
  matrixData->SetName("Overlap");

  // compute overlap of segments
  int status = 0;
  ttkTypeMacroA(
    scalars0->GetDataType(),
    (status = this->computeSegmentationOverlap<int, T0>(
       ttkUtils::GetPointer<int>(matrixData),

       ttkUtils::GetPointer<const int>(seg0),
       ttkUtils::GetPointer<const int>(seg1), seg0->GetNumberOfTuples(),
       ttkUtils::GetPointer<const int>(next0),
       ttkUtils::GetPointer<const int>(next1),
       ttkUtils::GetPointer<const T0>(scalars0),
       ttkUtils::GetPointer<const T0>(scalars1), nNodes0, nNodes1)));
  if(!status)
    return 0;

  status = ttkSimilarityAlgorithm::AddIndexIdMaps(
    similarityMatrix, this->GetInputArrayToProcess(1, m0),
    this->GetInputArrayToProcess(1, m1));
  if(!status)
    return 0;

  return 1;
}
