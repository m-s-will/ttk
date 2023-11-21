#include <ttkMergeSimilarityMatrices.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMergeSimilarityMatrices);

ttkMergeSimilarityMatrices::ttkMergeSimilarityMatrices() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkMergeSimilarityMatrices::FillInputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkMergeSimilarityMatrices::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkMergeSimilarityMatrices::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Get input object from input vector
  vtkMultiBlockDataSet *input = vtkMultiBlockDataSet::GetData(inputVector[0]);
  if(!input)
    return 0;

  // Get first block as MultiPieceDataSet
  vtkMultiPieceDataSet *inputMatrices
    = vtkMultiPieceDataSet::SafeDownCast(input->GetBlock(0));

  // Number of matrices and arrays
  int nMatrices = inputMatrices->GetNumberOfPieces();
  this->printMsg("The matrix consists of " + std::to_string(nMatrices)
                 + " submatrices.");

  // Go through each piece and store the unique ids
  std::set<ttk::SimplexId> uniqIds_t_1{};
  std::set<ttk::SimplexId> uniqIds_t{};
  auto allIds_t_1 = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  auto allIds_t = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  allIds_t_1->SetNumberOfComponents(1);
  allIds_t->SetNumberOfComponents(1);

  for(int i = 0; i < nMatrices; i++) {
    auto matrix = vtkImageData::SafeDownCast(inputMatrices->GetPiece(i));
    auto ids_t_1 = ttkAlgorithm::GetInputArrayToProcess(0, matrix);
    auto ids_t = ttkAlgorithm::GetInputArrayToProcess(1, matrix);

    int nSubFeatures0 = ids_t_1->GetNumberOfTuples();
    int nSubFeatures1 = ids_t->GetNumberOfTuples();

    // get id names from arrays
    if(i == 0) {
      allIds_t_1->SetName(ids_t_1->GetName());
      allIds_t->SetName(ids_t->GetName());
    }

    for(int j = 0; j < nSubFeatures0; j++) {
      auto res = uniqIds_t_1.emplace(ids_t_1->GetTuple1(j));
      // new id
      if(res.second)
        allIds_t_1->InsertNextTuple1(ids_t_1->GetTuple1(j));
    }
    for(int j = 0; j < nSubFeatures1; j++) {
      auto res = uniqIds_t.emplace(ids_t->GetTuple1(j));
      // new id
      if(res.second)
        allIds_t->InsertNextTuple1(ids_t->GetTuple1(j));
    }
  }

  // Create mappings between ids and row/column indices
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> columnMap_t_1{};
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> rowMap_t{};

  int status
    = ttkSimilarityAlgorithm::BuildIdIndexMap(columnMap_t_1, allIds_t_1);
  status = ttkSimilarityAlgorithm::BuildIdIndexMap(rowMap_t, allIds_t);

  // Initialize matrix data
  ttk::SimplexId nFeatures0 = allIds_t_1->GetNumberOfTuples();
  ttk::SimplexId nFeatures1 = allIds_t->GetNumberOfTuples();
  auto similarityMatrix = vtkSmartPointer<vtkImageData>::New();
  similarityMatrix->SetDimensions(nFeatures0, nFeatures1, 1);
  auto fieldData = similarityMatrix->GetFieldData();
  fieldData->AddArray(allIds_t_1);
  fieldData->AddArray(allIds_t);

  // Add first matrix array as Scalars array and init + get data from
  // submatrices
  auto pdPiece0 = inputMatrices->GetPiece(0)->GetPointData();
  int nArrays = pdPiece0->GetNumberOfArrays();
  auto fArray = pdPiece0->GetArray(0);
  similarityMatrix->AllocateScalars(fArray->GetDataType(), 1);
  auto firstArray = similarityMatrix->GetPointData()->GetArray(0);
  firstArray->SetName(fArray->GetName());

  status = 0;
  ttkTypeMacroA(
    firstArray->GetDataType(),
    (status = this->initializeSimilarityMatrix(
       ttkUtils::GetPointer<T0>(firstArray), nFeatures0, nFeatures1)));
  if(!status)
    return 0;

  for(int j = 0; j < nMatrices; j++) {
    auto matrix = vtkImageData::SafeDownCast(inputMatrices->GetPiece(j));
    auto subMatrixData = matrix->GetPointData()->GetArray(0);
    auto ids_t_1 = ttkAlgorithm::GetInputArrayToProcess(0, matrix);
    auto ids_t = ttkAlgorithm::GetInputArrayToProcess(1, matrix);
    int nSubFeatures0 = ids_t_1->GetNumberOfTuples();
    int nSubFeatures1 = ids_t->GetNumberOfTuples();

    status = 0;
    ttkTypeMacroA(subMatrixData->GetDataType(),
                  (status = this->addSubMatrixToMatrix<T0, ttk::SimplexId>(
                     ttkUtils::GetPointer<T0>(firstArray),
                     ttkUtils::GetPointer<T0>(subMatrixData),
                     ttkUtils::GetPointer<ttk::SimplexId>(ids_t_1),
                     ttkUtils::GetPointer<ttk::SimplexId>(ids_t), nSubFeatures0,
                     nSubFeatures1, nFeatures0, columnMap_t_1, rowMap_t)));
    if(!status)
      return 0;
  }

  // Add all similarities if more than one is defined on the matrix
  for(int i = 1; i < nArrays; i++) {
    auto curArray = pdPiece0->GetArray(i);
    auto matrixData
      = vtkSmartPointer<vtkDataArray>::Take(curArray->NewInstance());
    matrixData->SetNumberOfTuples(firstArray->GetNumberOfTuples());
    matrixData->SetName(curArray->GetName());
    similarityMatrix->GetPointData()->AddArray(matrixData);

    ttkTypeMacroA(
      matrixData->GetDataType(),
      (status = this->initializeSimilarityMatrix(
         ttkUtils::GetPointer<T0>(matrixData), nFeatures0, nFeatures1)));
    if(!status)
      return 0;

    for(int j = 0; j < nMatrices; j++) {
      auto matrix = vtkImageData::SafeDownCast(inputMatrices->GetPiece(j));
      auto subMatrixData = matrix->GetPointData()->GetArray(i);
      auto ids_t_1 = ttkAlgorithm::GetInputArrayToProcess(0, matrix);
      auto ids_t = ttkAlgorithm::GetInputArrayToProcess(1, matrix);
      int nSubFeatures0 = ids_t_1->GetNumberOfTuples();
      int nSubFeatures1 = ids_t->GetNumberOfTuples();

      status = 0;
      ttkTypeMacroA(
        subMatrixData->GetDataType(),
        (status = this->addSubMatrixToMatrix<T0, ttk::SimplexId>(
           ttkUtils::GetPointer<T0>(matrixData),
           ttkUtils::GetPointer<T0>(subMatrixData),
           ttkUtils::GetPointer<ttk::SimplexId>(ids_t_1),
           ttkUtils::GetPointer<ttk::SimplexId>(ids_t), nSubFeatures0,
           nSubFeatures1, nFeatures0, columnMap_t_1, rowMap_t)));
      if(!status)
        return 0;
    }
  }

  // Add matrix to output
  vtkMultiBlockDataSet *outputDataSet
    = vtkMultiBlockDataSet::GetData(outputVector, 0);

  outputDataSet->SetBlock(0, similarityMatrix);

  // return success
  return 1;
}
