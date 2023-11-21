#include <ttkSimilarityByGradientV2.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>

#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>



vtkStandardNewMacro(ttkSimilarityByGradientV2);

ttkSimilarityByGradientV2::ttkSimilarityByGradientV2() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}
ttkSimilarityByGradientV2::~ttkSimilarityByGradientV2() {
}

ttkSimplexIdTypeArray *GetVertexIdArray(vtkDataSet *input) {
  return ttkSimplexIdTypeArray::SafeDownCast(
    input->GetPointData()->GetArray("ttkVertexScalarField"));
}

int ttkSimilarityByGradientV2::ComputeSimilarityMatrix(
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

  auto manifoldSeg0 =  ttkAlgorithm::GetInputArrayToProcess(2, domain0);
  auto manifoldSeg1 =  ttkAlgorithm::GetInputArrayToProcess(2, domain1);

  if (!manifoldSeg0 || !manifoldSeg1)
    return !this->printErr("No manifold arrays.");

  #ifdef TTK_ENABLE_MPI
    // The matches and the vertex id arrays that may contain features outside the current rank 
    std::vector<ttk::SimplexId> forwardMapValsLocal{};
    std::vector<ttk::SimplexId> backwardMapValsLocal{};
    auto vertexIds0 = GetVertexIdArray(seeds0);
    auto vertexIds1 = GetVertexIdArray(seeds1);
    auto newVertexIds0 = vtkSmartPointer<ttkSimplexIdTypeArray>::Take(vertexIds0->NewInstance());
    newVertexIds0->DeepCopy(vertexIds0);
    auto newVertexIds1 = vtkSmartPointer<ttkSimplexIdTypeArray>::Take(vertexIds1->NewInstance());
    newVertexIds1->DeepCopy(vertexIds1);

    if (ttk::isRunningWithMPI()) {
      // Sets that store matches outside the current Rank
      std::set<ttk::SimplexId> outsideMatches0{};
      std::set<ttk::SimplexId> outsideMatches1{};

      // Forward
      int status = 0;
      ttkTypeMacroT(
        triangulation->getType(),
        (status = this->computeMap<ttk::SimplexId>(
          forwardMapValsLocal,
          outsideMatches1,
          (T0 *)triangulation->getData(),
          ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
          ttkUtils::GetPointer<ttk::SimplexId>(vertexIds0),
          totalFeatures0)));
      if(!status)
        return 0;

      // Backward
      status = 0;
      ttkTypeMacroT(
        triangulation->getType(),
        (status = this->computeMap<ttk::SimplexId, T0>(
          backwardMapValsLocal,
          outsideMatches0,
          (T0 *)triangulation->getData(),
          ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
          ttkUtils::GetPointer<ttk::SimplexId>(vertexIds1),
          totalFeatures1)));
      if(!status)
        return 0;

      // Calculate number of total features needed to represent the matches in the matrix and add these features to the newVertexId arrays
      totalFeatures0 += outsideMatches0.size();
      totalFeatures1 += outsideMatches1.size();

      for (auto vertexId : outsideMatches0) {
        newVertexIds0->InsertNextTuple1(vertexId);
      }
      for (auto vertexId : outsideMatches1) {
        newVertexIds1->InsertNextTuple1(vertexId);
      }

    }
  
  #endif

  // allocate similarity matrices
  similarityMatrix->SetDimensions(totalFeatures0, totalFeatures1, 1);
  similarityMatrix->AllocateScalars(VTK_INT, 1);

  auto forward = similarityMatrix->GetPointData()->GetArray(0);
  forward->SetName("Forward");

  auto backward = vtkSmartPointer<vtkIntArray>::New();
  backward->DeepCopy(forward);
  backward->SetName("Backward");
  similarityMatrix->GetPointData()->AddArray(backward);

  #ifdef TTK_ENABLE_MPI
    if (ttk::isRunningWithMPI()) {
      int status = 0;
      status = this->computeSimilarityMatrixMPI<ttk::SimplexId>(
        ttkUtils::GetPointer<int>(forward),
        forwardMapValsLocal,
        ttkUtils::GetPointer<ttk::SimplexId>(newVertexIds1),
        totalFeatures0, totalFeatures1,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
          ttk::SimplexId) { return j * n + i; },
          "Computing Forward Similarity Matrix");
      if(!status)
        return 0;

      status = this->computeSimilarityMatrixMPI<ttk::SimplexId>(
        ttkUtils::GetPointer<int>(backward),
        backwardMapValsLocal,
        ttkUtils::GetPointer<ttk::SimplexId>(newVertexIds0),
        totalFeatures1, totalFeatures0,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
          ttk::SimplexId m) { return i * m + j; },
          "Computing Backward Similarity Matrix");
      if(!status)
        return 0;
    }
  #else
    int status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeSimilarityMatrix<ttk::SimplexId>(
        ttkUtils::GetPointer<int>(forward),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg1),
        ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds0)),
        ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds1)),
        totalFeatures0, totalFeatures1,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
            ttk::SimplexId) { return j * n + i; },
            "Computing Forward Similarity Matrix")));
    if(!status)
      return 0;

    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeSimilarityMatrix<ttk::SimplexId>(
        ttkUtils::GetPointer<int>(backward),
        ttkUtils::GetPointer<ttk::SimplexId>(manifoldSeg0),
        ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds1)),
        ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds0)),
        totalFeatures0, totalFeatures1,
        [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId,
            ttk::SimplexId m) { return i * m + j; },
            "Computing Forward Similarity Matrix")));
    if(!status)
      return 0;
  #endif

  
  int status1 = 0;
  status1 = this->AddIndexIdMaps(similarityMatrix,
                                this->GetInputArrayToProcess(1, seeds0),
                                this->GetInputArrayToProcess(1, seeds1));
  if(!status)
    return 0;

  return 1;
}