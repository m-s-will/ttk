#include <ttkSimilarityByGradientV3.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByGradientV3);

ttkSimilarityByGradientV3::ttkSimilarityByGradientV3() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}
ttkSimilarityByGradientV3::~ttkSimilarityByGradientV3() {
}

ttkSimplexIdTypeArray *GetVertexIdArray(vtkDataSet *input) {
  return ttkSimplexIdTypeArray::SafeDownCast(
    input->GetPointData()->GetArray("ttkVertexScalarField"));
}

int ttkSimilarityByGradientV3::ComputeSimilarityMatrix(
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

  auto orderArray0 = ttkAlgorithm::GetOrderArray(domain0, 0);
  auto orderArray1 = ttkAlgorithm::GetOrderArray(domain1, 0);

  auto vertexIds0 = GetVertexIdArray(seeds0);
  auto vertexIds1 = GetVertexIdArray(seeds1);

  if(!manifoldSeg0 || !manifoldSeg1)
    return !this->printErr("No manifold arrays.");

#ifdef TTK_ENABLE_MPI
  // The matches and the vertex id arrays that may contain features outside the
  // current rank
  std::vector<ttk::SimplexId> forwardMapValsLocal{};
  std::vector<ttk::SimplexId> backwardMapValsLocal{};
  auto newVertexIds0
    = vtkSmartPointer<ttkSimplexIdTypeArray>::Take(vertexIds0->NewInstance());
  newVertexIds0->DeepCopy(vertexIds0);
  auto newVertexIds1
    = vtkSmartPointer<ttkSimplexIdTypeArray>::Take(vertexIds1->NewInstance());
  newVertexIds1->DeepCopy(vertexIds1);

  if(ttk::isRunningWithMPI()) {
    // Sets that store matches outside the current Rank
    std::set<ttk::SimplexId> outsideMatches0{};
    std::set<ttk::SimplexId> outsideMatches1{};

    // Forward
    int status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeMap<ttk::SimplexId>(
         forwardMapValsLocal, outsideMatches1, (T0 *)triangulation->getData(),
         ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
         ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0), totalFeatures0)));
    if(!status)
      return 0;

    // Backward
    status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeMap<ttk::SimplexId, T0>(
         backwardMapValsLocal, outsideMatches0, (T0 *)triangulation->getData(),
         ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
         ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1), totalFeatures1)));
    if(!status)
      return 0;

    // Calculate number of total features needed to represent the matches in the
    // matrix and add these features to the newVertexId arrays
    totalFeatures0 += outsideMatches0.size();
    totalFeatures1 += outsideMatches1.size();

    for(auto vertexId : outsideMatches0) {
      newVertexIds0->InsertNextTuple1(vertexId);
    }
    for(auto vertexId : outsideMatches1) {
      newVertexIds1->InsertNextTuple1(vertexId);
    }
  }

#endif

  // allocate similarity matrices
  similarityMatrix->SetDimensions(totalFeatures0, totalFeatures1, 1);
  similarityMatrix->AllocateScalars(VTK_INT, 1);

  auto forwardCombi = similarityMatrix->GetPointData()->GetArray(0);
  forwardCombi->SetName("Forward");

  auto backwardCombi= vtkSmartPointer<vtkIntArray>::New();
  backwardCombi->DeepCopy(forwardCombi);
  backwardCombi->SetName("Backward");
  similarityMatrix->GetPointData()->AddArray(backwardCombi);

  auto forwardDist = vtkSmartPointer<vtkIntArray>::New();
  forwardDist->DeepCopy(forwardCombi);
  forwardDist->SetName("ForwardDist");
  similarityMatrix->GetPointData()->AddArray(forwardDist);

  auto backwardDist = vtkSmartPointer<vtkIntArray>::New();
  backwardDist->DeepCopy(forwardCombi);
  backwardDist->SetName("BackwardDist");
  similarityMatrix->GetPointData()->AddArray(backwardDist);

  auto forManifold = vtkSmartPointer<vtkDoubleArray>::New();
  forManifold->SetNumberOfTuples(totalFeatures0*totalFeatures1);
  forManifold->SetNumberOfComponents(1);
  forManifold->SetName("ForwardManifold");
  similarityMatrix->GetPointData()->AddArray(forManifold);

  auto backManifold = vtkSmartPointer<vtkDoubleArray>::New();
  backManifold->SetNumberOfTuples(totalFeatures0*totalFeatures1);
  backManifold->SetNumberOfComponents(1);
  backManifold->SetName("BackwardManifold");
  similarityMatrix->GetPointData()->AddArray(backManifold);

#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    int status = 0;
    status = this->computeSimilarityMatrixMPI<ttk::SimplexId>(
      ttkUtils::GetPointer<int>(forwardCombi), forwardMapValsLocal,
      ttkUtils::GetPointer<ttk::SimplexId>(newVertexIds1), totalFeatures0,
      totalFeatures1,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n, ttk::SimplexId) {
        return j * n + i;
      },
      "Computing Forward Similarity Matrix");
    if(!status)
      return 0;

    status = this->computeSimilarityMatrixMPI<ttk::SimplexId>(
      ttkUtils::GetPointer<int>(backwardCombi), backwardMapValsLocal,
      ttkUtils::GetPointer<ttk::SimplexId>(newVertexIds0), totalFeatures1,
      totalFeatures0,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId, ttk::SimplexId m) {
        return i * m + j;
      },
      "Computing Backward Similarity Matrix");
    if(!status)
      return 0;

    status = 0;
    status
      = this->AddIndexIdMaps(similarityMatrix, newVertexIds0, newVertexIds1);
    if(!status)
      return 0;
  }
  else {
    // combinatorial
    std::vector<std::map<ttk::SimplexId, int>> forwardCombiMap(totalFeatures0);
    std::vector<std::map<ttk::SimplexId, int>> backwardCombiMap(totalFeatures1);

    // Forward
    int status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeMapCombinatorial<ttk::SimplexId>(
          forwardCombiMap, (T0 *)triangulation->getData(),
          ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
          ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0), totalFeatures0, this->NeighbourhoodSize)));
    if(!status)
      return 0;

    // Backward
    status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeMapCombinatorial<ttk::SimplexId>(
          backwardCombiMap, (T0 *)triangulation->getData(),
          ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
          ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1), totalFeatures1, this->NeighbourhoodSize)));
    if(!status)
      return 0;

    // Create mappings between ids and row/column indices for matrix creation
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> map0{};
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> map1{};

    status = ttkSimilarityAlgorithm::BuildIdIndexMap(map0, vertexIds0);
    status = ttkSimilarityAlgorithm::BuildIdIndexMap(map1, vertexIds1);

    status = this->computeSimilarityMatrix<ttk::SimplexId>(
        ttkUtils::GetPointer<int>(forward),
        forwardCombiMap,
        totalFeatures0, totalFeatures1, map1,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
            ttk::SimplexId) { return j * n + i; },
        "Computing Combi Forward Similarity Matrix");
    if(!status)
      return 0;

    status = this->computeSimilarityMatrix<ttk::SimplexId>(
        ttkUtils::GetPointer<int>(backwardCombi),
        backwardCombiMap,
        totalFeatures1, totalFeatures0, map0,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
            ttk::SimplexId m) { return i * m + j; },
        "Computing Combi Backward Similarity Matrix");
    if(!status)
      return 0;

    // MANIFOLDS //
    std::map<ttk::SimplexId, std::map<ttk::SimplexId, int>> forwardMapManifold;
    std::map<ttk::SimplexId, std::map<ttk::SimplexId, int>> backwardMapManifold;
    status = 0;
    status = this->computeMapManifoldOverlap<ttk::SimplexId>(
        forwardMapManifold, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0), ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
        totalVertices);
    if(!status)
      return 0;
    status = 0;
    status = this->computeMapManifoldOverlap<ttk::SimplexId>(
        backwardMapManifold, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1), ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
        totalVertices);
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
        ttkUtils::GetPointer<double>(forManifold),
        forwardMapManifold,
        manifoldSize0,
        manifoldSize1,
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0),
        totalFeatures0, totalFeatures1, map1,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
            ttk::SimplexId) { return j * n + i; },
        "Computing Manifold Forward Similarity Matrix");
    if(!status)
      return 0;
    status = 0;
    status = this->computeSimilarityManifoldMatrix<ttk::SimplexId>(
        ttkUtils::GetPointer<double>(backManifold),
        backwardMapManifold,
        manifoldSize1,
        manifoldSize0,
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1),
        totalFeatures1, totalFeatures0, map0,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
            ttk::SimplexId m) { return i * m + j; },
        "Computing Manifold Backward Similarity Matrix");
    if(!status)
      return 0;
    

    status = 0;
    status = this->AddIndexIdMaps(similarityMatrix,
                                  this->GetInputArrayToProcess(1, seeds0),
                                  this->GetInputArrayToProcess(1, seeds1));
    if(!status)
      return 0;
  }
#else

  // Create mappings between ids and row/column indices for matrix creation
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> map0{};
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> map1{};

  int status = ttkSimilarityAlgorithm::BuildIdIndexMap(map0, vertexIds0);
  status = ttkSimilarityAlgorithm::BuildIdIndexMap(map1, vertexIds1);

  // combinatorial
  std::vector<std::map<ttk::SimplexId, int>> forwardCombiMap(totalFeatures0);
  std::vector<std::map<ttk::SimplexId, int>> backwardCombiMap(totalFeatures1);

  // combi - distribution
  std::vector<std::map<ttk::SimplexId, int>> forwardCombiDistrMap(totalFeatures0);
  std::vector<std::map<ttk::SimplexId, int>> backwardCombiDistrMap(totalFeatures1);

  // euclidean grid distance
  std::vector<std::map<ttk::SimplexId, int>> forwardDistMap(totalFeatures0);
  std::vector<std::map<ttk::SimplexId, int>> backwardDistMap(totalFeatures1);

  // Forward combi
  status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (status = this->computeMapCombinatorial<ttk::SimplexId>(
        forwardCombiMap, (T0 *)triangulation->getData(),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0), totalFeatures0, this->NeighbourhoodSize)));
  if(!status)
    return 0;

  status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (status = this->computeMapCombinatorialDistribution<ttk::SimplexId>(
        forwardCombiDistrMap, (T0 *)triangulation->getData(),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0), ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1), ttkUtils::GetPointer<ttk::SimplexId>(orderArray0), ttkUtils::GetPointer<ttk::SimplexId>(orderArray1), totalFeatures0, this->NeighbourhoodSize)));
  if(!status)
    return 0;

  // Backward combi
  status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (status = this->computeMapCombinatorial<ttk::SimplexId>(
        backwardCombiMap, (T0 *)triangulation->getData(),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1), totalFeatures1, this->NeighbourhoodSize)));
  if(!status)
    return 0;

  // Forward dist
  status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (status = this->computeMapDistance<ttk::SimplexId>(
        forwardDistMap, (T0 *)triangulation->getData(),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0), totalFeatures0, this->NeighbourhoodSize, this->NeighbourhoodDistance)));
  if(!status)
    return 0;

  // Backward dist
  status = 0;
  ttkTypeMacroT(
    triangulation->getType(),
    (status = this->computeMapDistance<ttk::SimplexId>(
        backwardDistMap, (T0 *)triangulation->getData(),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
        ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1), totalFeatures1, this->NeighbourhoodSize, this->NeighbourhoodDistance)));
  if(!status)
    return 0;

  status = this->computeSimilarityMatrix<ttk::SimplexId>(
      ttkUtils::GetPointer<int>(forwardCombi),
      forwardCombiMap,
      totalFeatures0, totalFeatures1, map1,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
          ttk::SimplexId) { return j * n + i; },
      "Computing Combi Forward Similarity Matrix");
  if(!status)
    return 0;

  status = this->computeSimilarityMatrix<ttk::SimplexId>(
      ttkUtils::GetPointer<int>(backwardCombi),
      backwardCombiMap,
      totalFeatures1, totalFeatures0, map0,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
          ttk::SimplexId m) { return i * m + j; },
      "Computing Combi Backward Similarity Matrix");
  if(!status)
    return 0;

  status = this->computeSimilarityMatrix<ttk::SimplexId>(
      ttkUtils::GetPointer<int>(forwardDist),
      forwardDistMap,
      totalFeatures0, totalFeatures1, map1,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
          ttk::SimplexId) { return j * n + i; },
      "Computing Distance Forward Similarity Matrix");
  if(!status)
    return 0;

  status = this->computeSimilarityMatrix<ttk::SimplexId>(
      ttkUtils::GetPointer<int>(backwardDist),
      backwardDistMap,
      totalFeatures1, totalFeatures0, map0,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
          ttk::SimplexId m) { return i * m + j; },
      "Computing Distance Backward Similarity Matrix");
  if(!status)
    return 0;

  // MANIFOLDS //
  std::map<ttk::SimplexId, std::map<ttk::SimplexId, int>> forwardMapManifold;
  std::map<ttk::SimplexId, std::map<ttk::SimplexId, int>> backwardMapManifold;
  status = 0;
  status = this->computeMapManifoldOverlap<ttk::SimplexId>(
      forwardMapManifold, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0), ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
      totalVertices);
  if(!status)
    return 0;
  status = 0;
  status = this->computeMapManifoldOverlap<ttk::SimplexId>(
      backwardMapManifold, ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1), ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
      totalVertices);
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
      ttkUtils::GetPointer<double>(forManifold),
      forwardMapManifold,
      manifoldSize0,
      manifoldSize1,
      ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0),
      totalFeatures0, totalFeatures1, map1,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
          ttk::SimplexId) { return j * n + i; },
      "Computing Manifold Forward Similarity Matrix");
  if(!status)
    return 0;
  status = 0;
  status = this->computeSimilarityManifoldMatrix<ttk::SimplexId>(
      ttkUtils::GetPointer<double>(backManifold),
      backwardMapManifold,
      manifoldSize1,
      manifoldSize0,
      ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1),
      totalFeatures1, totalFeatures0, map0,
      [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
          ttk::SimplexId m) { return i * m + j; },
      "Computing Manifold Backward Similarity Matrix");
  if(!status)
    return 0;
  

  status = 0;
  status = this->AddIndexIdMaps(similarityMatrix,
                                this->GetInputArrayToProcess(1, seeds0),
                                this->GetInputArrayToProcess(1, seeds1));
  if(!status)
    return 0;
#endif

  return 1;
}