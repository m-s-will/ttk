#include <ttkSimilarityById.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityById);

ttkSimilarityById::ttkSimilarityById() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkSimilarityById::~ttkSimilarityById() {
}

int ttkSimilarityById::ComputeSimilarityMatrix(
  vtkImageData *similarityMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {
  // unpack input
  auto p0 = vtkDataSet::SafeDownCast(inputDataObjects0);
  auto p1 = vtkDataSet::SafeDownCast(inputDataObjects1);
  if(!p0 || !p1)
    return !this->printErr("Input data objects need to be vtkPointSets.");

  const int nPoints0 = p0->GetNumberOfPoints();
  const int nPoints1 = p1->GetNumberOfPoints();

  // Get ids
  auto ids0 = GetInputArrayToProcess(0, p0);
  auto ids1 = GetInputArrayToProcess(0, p1);
  if(!ids0 || !ids1)
    return !this->printErr("Input data is missing array containing ids.");

  auto uniqueIds0 = vtkSmartPointer<vtkDataArray>::Take(ids0->NewInstance());
  auto uniqueIds1 = vtkSmartPointer<vtkDataArray>::Take(ids1->NewInstance());
  uniqueIds0->SetName(ids0->GetName());
  uniqueIds1->SetName(ids1->GetName());
  uniqueIds0->SetNumberOfComponents(1);
  uniqueIds1->SetNumberOfComponents(1);
  uniqueIds0->SetNumberOfTuples(nPoints0);
  uniqueIds1->SetNumberOfTuples(nPoints1);
  int numUniqueIds0 = 0;
  int numUniqueIds1 = 0;

  // Traverse arrays for unique ids
  int status = 0;

  ttkTypeMacroI(ids0->GetDataType(),
                (status = this->computeUniqueIds<T0>(
                   ttkUtils::GetPointer<T0>(uniqueIds0),
                   ttkUtils::GetPointer<T0>(uniqueIds1), numUniqueIds0,
                   numUniqueIds1, ttkUtils::GetPointer<const T0>(ids0),
                   ttkUtils::GetPointer<const T0>(ids1), nPoints0, nPoints1)));
  if(!status)
    return 0;

  uniqueIds0->Resize(numUniqueIds0);
  uniqueIds1->Resize(numUniqueIds1);

  // initialize similarity matrix i.e., identity matrix
  similarityMatrix->SetDimensions(numUniqueIds0, numUniqueIds1, 1);
  similarityMatrix->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  auto matrixData = similarityMatrix->GetPointData()->GetArray(0);
  matrixData->SetName("Identity");

  // // compute indentity matrix
  ttkTypeMacroI(
    ids0->GetDataType(), (status = this->computeIdentityMatrix<T0>(
                            ttkUtils::GetPointer<unsigned char>(matrixData),
                            ttkUtils::GetPointer<const T0>(uniqueIds0),
                            ttkUtils::GetPointer<const T0>(uniqueIds1),
                            numUniqueIds0, numUniqueIds1)));
  if(!status)
    return 0;

  status = ttkSimilarityAlgorithm::AddIndexIdMaps(
    similarityMatrix, uniqueIds0, uniqueIds1);
  if(!status)
    return 0;

  return 1;
}
