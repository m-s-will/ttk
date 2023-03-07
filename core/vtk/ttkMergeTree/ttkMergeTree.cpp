#include <ttkMergeTree.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMergeTree);

ttkMergeTree::ttkMergeTree() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

ttkMergeTree::~ttkMergeTree() {
}

int ttkMergeTree::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMergeTree::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      return 1;
    case 1:
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
      return 1;
    default:
      return 0;
  }
}

int ttkMergeTree::RequestData(vtkInformation *,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  // Get the input
  auto input = vtkDataSet::GetData(inputVector[0]);
  if(!input) {
    this->printErr("Unable to retrieve input data object.");
    return 0;
  }
  const size_t nVertices = input->GetNumberOfPoints();

  // Get triangulation of the input object
  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  // Precondition triangulation
  this->PreconditionTriangulation(triangulation);

  // Get input array
  auto scalarArray = this->GetInputArrayToProcess(0, inputVector);
  if(!scalarArray) {
    this->printErr("Unable to retrieve scalar array.");
    return 0;
  }

  // Order Array
  auto orderArray = this->GetOrderArray(input, 0);
  const auto orderArrayData
    = ttkUtils::GetPointer<const ttk::SimplexId>(orderArray);

  // Init segmentation
  auto segmentationIds
    = vtkSmartPointer<vtkDataArray>::Take(orderArray->NewInstance());
  segmentationIds->SetName("BranchId");
  segmentationIds->SetNumberOfComponents(1);
  segmentationIds->SetNumberOfTuples(nVertices);

  // Compute merge tree segmentation
  std::vector<ttk::mt::Propagation<ttk::SimplexId>> propagations;
  std::vector<const ttk::mt::Propagation<ttk::SimplexId> *> sortedPropagations;
  {
    int status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = this->computeMergeTreeSegmentation<ttk::SimplexId, T0>(
         ttkUtils::GetPointer<ttk::SimplexId>(segmentationIds), propagations,

         static_cast<T0 *>(triangulation->getData()), orderArrayData,
         this->GetType())));
    if(!status)
      return 0;

    ttkTypeMacroAT(
      scalarArray->GetDataType(), triangulation->getType(),
      (status = this->finalizePropagations<T0, ttk::SimplexId, T1>(
         sortedPropagations,
         ttkUtils::GetPointer<ttk::SimplexId>(segmentationIds), propagations,

         static_cast<T1 *>(triangulation->getData()),
         ttkUtils::GetPointer<T0>(scalarArray))));
    if(!status)
      return 0;
  }

  // Finalize Output
  {
    ttk::Timer timer;
    this->printMsg(
      "Generating Output Data Objects", 0, 0, ttk::debug::LineMode::REPLACE);

    // Create segmentation output
    {
      auto segmentation = vtkDataSet::GetData(outputVector, 1);
      segmentation->ShallowCopy(input);

      auto segmentationPD = segmentation->GetPointData();
      segmentationPD->AddArray(segmentationIds);
      segmentationPD->AddArray(orderArray);
    }

    // Compute merge tree output
    {
      auto mergeTree = vtkUnstructuredGrid::GetData(outputVector, 0);

      const int nPropagations = propagations.size();

      std::unordered_map<int, std::pair<int, int>>
        vertexToNodeMap; // vertexId -> NodeId,BranchId

      // build map
      {
        for(int p = 0; p < nPropagations; p++) {
          vertexToNodeMap.insert(
            {sortedPropagations[p]->criticalPoints[0], {p, p}});
        }
        for(int p = 0, i = nPropagations; p < nPropagations; p++) {
          const auto &criticalPoints = sortedPropagations[p]->criticalPoints;
          const int nCriticalPoints = criticalPoints.size();
          for(int c = 1; c < nCriticalPoints; c++) {
            auto it = vertexToNodeMap.find(criticalPoints[c]);
            if(it == vertexToNodeMap.end())
              vertexToNodeMap.insert({criticalPoints[c], {i++, p}});
          }
        }
      }

      const int nNodes = vertexToNodeMap.size();
      const int nEdges = nNodes - 1;

      // nodes
      {
        auto pointCoords = vtkSmartPointer<vtkFloatArray>::New();
        pointCoords->SetNumberOfComponents(3);
        pointCoords->SetNumberOfTuples(nNodes);
        auto pointCoordData
          = static_cast<float *>(ttkUtils::GetVoidPointer(pointCoords));

        auto vertexId = vtkSmartPointer<vtkIntArray>::New();
        vertexId->SetName("VertexId");
        vertexId->SetNumberOfTuples(nNodes);
        auto vertexIdData
          = static_cast<int *>(ttkUtils::GetVoidPointer(vertexId));

        auto branchId = vtkSmartPointer<vtkIntArray>::New();
        branchId->SetName("BranchId");
        branchId->SetNumberOfTuples(nNodes);
        auto branchIdData
          = static_cast<int *>(ttkUtils::GetVoidPointer(branchId));

        auto outputScalars
          = vtkSmartPointer<vtkDataArray>::Take(scalarArray->NewInstance());
        outputScalars->SetName(scalarArray->GetName());
        outputScalars->SetNumberOfTuples(nNodes);

        for(const auto &it : vertexToNodeMap) {
          const auto &propId = it.second.second;
          const auto &nodeId = it.second.first;
          const auto &offset = nodeId * 3;

          vertexIdData[nodeId] = it.first;
          branchIdData[nodeId] = propId;

          outputScalars->SetTuple(nodeId, it.first, scalarArray);

          triangulation->getVertexPoint(it.first, pointCoordData[offset],
                                        pointCoordData[offset + 1],
                                        pointCoordData[offset + 2]);
        }

        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetData(pointCoords);
        mergeTree->SetPoints(points);

        auto mergeTreePD = mergeTree->GetPointData();
        mergeTreePD->AddArray(vertexId);
        mergeTreePD->AddArray(branchId);
        mergeTreePD->AddArray(outputScalars);
      }

      // edges
      {
        auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
        offsets->SetNumberOfTuples(nEdges + 1);
        auto offsetsData
          = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(offsets));
        for(int i = 0; i <= nEdges; i++)
          offsetsData[i] = i * 2;

        auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
        connectivity->SetNumberOfTuples(nEdges * 2);
        auto connectivityData
          = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(connectivity));

        auto edgebranchId = vtkSmartPointer<vtkIntArray>::New();
        edgebranchId->SetName("BranchId");
        edgebranchId->SetNumberOfTuples(nEdges);
        auto edgebranchIdData
          = static_cast<int *>(ttkUtils::GetVoidPointer(edgebranchId));

        auto edgePersistence = vtkSmartPointer<vtkDoubleArray>::New();
        edgePersistence->SetName("Persistence");
        edgePersistence->SetNumberOfTuples(nEdges);
        auto edgePersistenceData
          = static_cast<double *>(ttkUtils::GetVoidPointer(edgePersistence));

        auto nextId = vtkSmartPointer<vtkIntArray>::New();
        nextId->SetName("NextId");
        nextId->SetNumberOfTuples(nNodes);
        auto nextIdData = static_cast<int *>(ttkUtils::GetVoidPointer(nextId));

        auto type = vtkSmartPointer<vtkSignedCharArray>::New();
        type->SetName("Type");
        type->SetNumberOfTuples(nNodes);
        auto typeData
          = static_cast<signed char *>(ttkUtils::GetVoidPointer(type));

        auto pairId = vtkSmartPointer<vtkIntArray>::New();
        pairId->SetName("PairId");
        pairId->SetNumberOfTuples(nNodes);
        auto pairIdData = static_cast<int *>(ttkUtils::GetVoidPointer(pairId));

        for(int i = 0; i < nNodes; i++) {
          typeData[i] = 0;
          nextIdData[i] = -1;
        }

        for(int p = 0, i = 0, j = 0; p < nPropagations; p++) {
          const auto &criticalPoints = sortedPropagations[p]->criticalPoints;
          const int nCriticalPoints = criticalPoints.size();

          double persistence = std::abs(
            (*scalarArray->GetTuple(criticalPoints[0]))
            - (*scalarArray->GetTuple(criticalPoints[nCriticalPoints - 1])));

          typeData[vertexToNodeMap[criticalPoints[0]].first] = 1;

          for(int c = 1; c < nCriticalPoints; c++) {
            auto &v0it = vertexToNodeMap[criticalPoints[c - 1]];
            auto &v1it = vertexToNodeMap[criticalPoints[c]];
            edgebranchIdData[j] = v0it.second;
            edgePersistenceData[j++] = persistence;

            nextIdData[v0it.first] = v1it.first;

            connectivityData[i++] = v0it.first;
            connectivityData[i++] = v1it.first;
          }

          // persistence pair
          {
            auto &v0it = vertexToNodeMap[criticalPoints[0]];
            auto &v1it = vertexToNodeMap[criticalPoints[nCriticalPoints - 1]];
            pairIdData[v0it.first] = v1it.first;
            pairIdData[v1it.first] = v0it.first;
          }
        }

        auto cells = vtkSmartPointer<vtkCellArray>::New();
        cells->SetData(offsets, connectivity);
        mergeTree->SetCells(VTK_LINE, cells);
        auto mergeTreeCD = mergeTree->GetCellData();
        mergeTreeCD->AddArray(edgebranchId);
        mergeTreeCD->AddArray(edgePersistence);

        auto mergeTreePD = mergeTree->GetPointData();
        mergeTreePD->AddArray(nextId);
        mergeTreePD->AddArray(pairId);
        mergeTreePD->AddArray(type);
      }
    }

    this->printMsg("Generating Output Data Objects", 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1);
  }

  return 1;
}
