#include <ttkSimilarityByGradient.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>

#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByGradient);

ttkSimilarityByGradient::ttkSimilarityByGradient() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}
ttkSimilarityByGradient::~ttkSimilarityByGradient() {
}

ttkSimplexIdTypeArray *GetVertexIdArray(vtkDataSet *input) {
  return ttkSimplexIdTypeArray::SafeDownCast(
    input->GetPointData()->GetArray("ttkVertexScalarField"));
}

int ttkSimilarityByGradient::ComputeSimilarityMatrix(
  vtkImageData *similarityMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {
  auto inputs0AsMB = static_cast<vtkMultiBlockDataSet *>(inputDataObjects0);
  auto inputs1AsMB = static_cast<vtkMultiBlockDataSet *>(inputDataObjects1);

  auto domain0 = vtkDataSet::SafeDownCast(inputs0AsMB->GetBlock(0));
  auto domain1 = vtkDataSet::SafeDownCast(inputs1AsMB->GetBlock(0));
  auto seeds0 = vtkDataSet::SafeDownCast(inputs0AsMB->GetBlock(1));
  auto seeds1 = vtkDataSet::SafeDownCast(inputs1AsMB->GetBlock(1));

  int totalFeatures0 = seeds0->GetNumberOfElements(0);
  int totalFeatures1 = seeds1->GetNumberOfElements(0);

  auto triangulation = this->GetTriangulation(domain0);
  this->preconditionTriangulation(triangulation);
  int totalVertices = triangulation->getNumberOfVertices();

  auto manifoldSeg0 = ttkAlgorithm::GetInputArrayToProcess(2, domain0);
  auto manifoldSeg1 = ttkAlgorithm::GetInputArrayToProcess(2, domain1);

  auto vertexIds0 = GetVertexIdArray(seeds0);
  auto vertexIds1 = GetVertexIdArray(seeds1);

  if(!manifoldSeg0 || !manifoldSeg1)
    return !this->printErr("No manifold arrays.");

  // allocate similarity matrices
  similarityMatrix->SetDimensions(totalFeatures0, totalFeatures1, 1);
  similarityMatrix->AllocateScalars(VTK_INT, 1);

  auto forManifold = vtkSmartPointer<vtkDoubleArray>::New();
  forManifold->SetNumberOfTuples(totalFeatures0 * totalFeatures1);
  forManifold->SetNumberOfComponents(1);
  forManifold->SetName("ForwardManifold");
  similarityMatrix->GetPointData()->AddArray(forManifold);

  auto backManifold = vtkSmartPointer<vtkDoubleArray>::New();
  backManifold->SetNumberOfTuples(totalFeatures0 * totalFeatures1);
  backManifold->SetNumberOfComponents(1);
  backManifold->SetName("BackwardManifold");
  similarityMatrix->GetPointData()->AddArray(backManifold);

  // Create mappings between ids and row/column indices for matrix creation
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> map0{};
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> map1{};

  int status = ttkSimilarityAlgorithm::BuildIdIndexMap(map0, vertexIds0);
  status = ttkSimilarityAlgorithm::BuildIdIndexMap(map1, vertexIds1);

  // MANIFOLDS //
  std::map<ttk::SimplexId, std::map<ttk::SimplexId, int>> forwardMapManifold;
  std::map<ttk::SimplexId, std::map<ttk::SimplexId, int>> backwardMapManifold;
  status = 0;
  status = this->computeMapManifoldOverlap<ttk::SimplexId>(
    forwardMapManifold, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
    ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1), totalVertices);
  if(!status)
    return 0;
  status = 0;
  status = this->computeMapManifoldOverlap<ttk::SimplexId>(
    backwardMapManifold, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
    ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0), totalVertices);
  if(!status)
    return 0;

  // MANIFOLD SIZES
  std::map<ttk::SimplexId, int> manifoldSize0;
  std::map<ttk::SimplexId, int> manifoldSize1;
  status = 0;
  status = this->computeManifoldSize<ttk::SimplexId>(
    manifoldSize0, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
    totalVertices);
  status = 0;
  status = this->computeManifoldSize<ttk::SimplexId>(
    manifoldSize1, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
    totalVertices);

  status = 0;
  status = this->computeSimilarityManifoldMatrix<ttk::SimplexId>(
    ttkUtils::GetPointer<double>(forManifold), forwardMapManifold,
    manifoldSize0, ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0),
    totalFeatures0, totalFeatures1, map1,
    [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n, ttk::SimplexId) {
      return j * n + i;
    },
    "Computing Manifold Forward Similarity Matrix");
  if(!status)
    return 0;
  status = 0;
  status = this->computeSimilarityManifoldMatrix<ttk::SimplexId>(
    ttkUtils::GetPointer<double>(backManifold), backwardMapManifold,
    manifoldSize1, ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1),
    totalFeatures1, totalFeatures0, map0,
    [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId, ttk::SimplexId m) {
      return i * m + j;
    },
    "Computing Manifold Backward Similarity Matrix");
  if(!status)
    return 0;

  status = 0;
  status = this->AddIndexIdMaps(similarityMatrix,
                                this->GetInputArrayToProcess(1, seeds0),
                                this->GetInputArrayToProcess(1, seeds1));
  if(!status)
    return 0;

  return 1;
}