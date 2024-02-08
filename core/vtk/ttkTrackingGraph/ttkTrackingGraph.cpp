#include <ttkTrackingGraph.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkSimilarityAlgorithm.h>
#include <ttkUtils.h>

typedef std::vector<std::tuple<vtkFieldData *, vtkFieldData *, int>>
  ArrayMapSet;

vtkStandardNewMacro(ttkTrackingGraph);

ttkTrackingGraph::ttkTrackingGraph() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkTrackingGraph::~ttkTrackingGraph() {
}

int ttkTrackingGraph::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    if(port == 1)
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkTrackingGraph::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename CT, typename IT>
int countEdges(int &nEdges,
               const int *dim,
               const CT *similarities,
               const int nIds0 = 0,
               const int nIds1 = 0,
               const IT *ids0 = nullptr,
               const IT *ids1 = nullptr,
               const vtkDataArray *indexIdMap0 = nullptr,
               const vtkDataArray *indexIdMap1 = nullptr) {

  if(nIds0 > 0 && nIds1 > 0) {
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> idIndexMap0;
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> idIndexMap1;
    ttkSimilarityAlgorithm::BuildIdIndexMap(idIndexMap0, indexIdMap0);
    ttkSimilarityAlgorithm::BuildIdIndexMap(idIndexMap1, indexIdMap1);

    for(int i = 0; i < nIds0; i++) {
      for(int j = 0; j < nIds1; j++) {
        auto iId = static_cast<ttk::SimplexId>(ids0[i]);
        auto jId = static_cast<ttk::SimplexId>(ids1[j]);
        const auto &iIt = idIndexMap0.find(iId);
        const auto &jIt = idIndexMap1.find(jId);
        if(iIt == idIndexMap0.end() || jIt == idIndexMap1.end())
          continue;

        const auto &iIdx = iIt->second;
        const auto &jIdx = jIt->second;
        const int index = jIdx * dim[0] + iIdx;
        if(similarities[index] > 0)
          nEdges++;
      }
    }
  } else {
    for(int i = 0; i < dim[0]; i++)
      for(int j = 0; j < dim[1]; j++) {
        const auto &iIdx = i;
        const auto &jIdx = j;
        const int index = jIdx * dim[0] + iIdx;
        if(similarities[index] > 0)
          nEdges++;
      }
  }

  return 1;
}

template <typename CT, typename IT>
int generateEdges(vtkPolyData *output,
                  int &iEdge,
                  vtkFieldData *trackingGraphCD,
                  vtkFieldData *similaritiesPD,
                  const CT *similarities,
                  const int &offset0,
                  const int &offset1,
                  const int *dim,
                  const int nIds0 = 0,
                  const int nIds1 = 0,
                  const IT *ids0 = nullptr,
                  const IT *ids1 = nullptr,
                  const vtkDataArray *indexIdMap0 = nullptr,
                  const vtkDataArray *indexIdMap1 = nullptr) {

  std::vector<std::pair<vtkAbstractArray *, vtkAbstractArray *>> arrayMap;
  for(int a = 0; a < trackingGraphCD->GetNumberOfArrays(); a++) {
    auto oArray = trackingGraphCD->GetAbstractArray(a);
    arrayMap.push_back({oArray, similaritiesPD->GetArray(oArray->GetName())});
  }

  if(nIds0 > 0 && nIds1 > 0) {
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> idIndexMap0;
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> idIndexMap1;
    ttkSimilarityAlgorithm::BuildIdIndexMap(idIndexMap0, indexIdMap0);
    ttkSimilarityAlgorithm::BuildIdIndexMap(idIndexMap1, indexIdMap1);

    for(int i = 0; i < nIds0; i++) {
      for(int j = 0; j < nIds1; j++) {
        auto iId = static_cast<ttk::SimplexId>(ids0[i]);
        auto jId = static_cast<ttk::SimplexId>(ids1[j]);
        const auto &iIt = idIndexMap0.find(iId);
        const auto &jIt = idIndexMap1.find(jId);
        if(iIt == idIndexMap0.end() || jIt == idIndexMap1.end())
          continue;

        const auto &iIdx = iIt->second;
        const auto &jIdx = jIt->second;
        const int index = jIdx * dim[0] + iIdx;
        if(similarities[index] <= 0)
          continue;

        vtkIdType points[2]{offset0 + i, offset1 + j};
        output->InsertNextCell(VTK_LINE, 2, points);

        for(auto &it : arrayMap)
          it.first->SetTuple(iEdge, index, it.second);

        iEdge++;
      }
    }
  } else {
    for(int i = 0; i < dim[0]; i++) {
      for(int j = 0; j < dim[1]; j++) {
        const auto &iIdx = i;
        const auto &jIdx = j;
        const int index = jIdx * dim[0] + iIdx;
        if(similarities[index] <= 0)
          continue;

        const vtkIdType points[2]{offset0 + i, offset1 + j};
        output->InsertNextCell(VTK_LINE, 2, points);

        for(auto &it : arrayMap)
          it.first->SetTuple(iEdge, index, it.second);

        iEdge++;
      }
    }
  }

  return 1;
}

int ttkTrackingGraph::CountNodesAndEdges(int &nNodes,
                                         int &nEdges,
                                         std::vector<int> &nodeIdxOffsets,
                                         vtkMultiBlockDataSet *similarities,
                                         vtkMultiBlockDataSet *features
                                         = nullptr) {

  const int nSteps = similarities->GetNumberOfBlocks() + 1;

  nNodes = 0;
  nEdges = 0;
  nodeIdxOffsets.resize(nSteps);
  nodeIdxOffsets[0] = 0;

  if(features) {
    // nodes
    for(int t = 1; t < nSteps; t++) {
      auto f = static_cast<vtkPointSet *>(features->GetBlock(t - 1));
      const auto n = f->GetNumberOfPoints();
      nNodes += n;
      nodeIdxOffsets[t] = n + nodeIdxOffsets[t - 1];
    }

    // edges
    for(int t = 1; t < nSteps; t++) {
      auto c = static_cast<vtkImageData *>(similarities->GetBlock(t - 1));
      int dim[3];
      c->GetDimensions(dim);
      auto matrix = this->GetInputArrayToProcess(0, c);
      auto cFD = c->GetFieldData();

      auto f0 = static_cast<vtkPointSet *>(features->GetBlock(t - 1));
      auto f1 = static_cast<vtkPointSet *>(features->GetBlock(t));
      if(f0->GetNumberOfPoints() < 1 || f1->GetNumberOfPoints() < 1)
        continue;

      const auto idArrayName = ttkSimilarityAlgorithm::GetIdArrayName(cFD);
      auto l0 = f0->GetPointData()->GetArray(idArrayName.data());
      auto l1 = f1->GetPointData()->GetArray(idArrayName.data());

      vtkDataArray *indexIdMapP{nullptr};
      vtkDataArray *indexIdMapC{nullptr};
      if(!ttkSimilarityAlgorithm::GetIndexIdMaps(indexIdMapP, indexIdMapC, cFD))
        return !this->printErr("Unable to retrieve Index-Id-Maps.");

      int status = 0;
      ttkTypeMacroAI(matrix->GetDataType(), l0->GetDataType(),
                     status = countEdges(
                       nEdges, dim, ttkUtils::GetPointer<T0>(matrix),
                       l0->GetNumberOfValues(), l1->GetNumberOfValues(),
                       ttkUtils::GetPointer<T1>(l0),
                       ttkUtils::GetPointer<T1>(l1), indexIdMapP, indexIdMapC));
      if(!status)
        return 0;
    }
  } else {
    for(int t = 1; t < nSteps; t++) {
      auto c = static_cast<vtkImageData *>(similarities->GetBlock(t - 1));
      int dim[3];
      c->GetDimensions(dim);

      // nodes
      if(t == 1)
        nNodes += dim[0];
      nNodes += dim[1];

      nodeIdxOffsets[t] = dim[0] + nodeIdxOffsets[t - 1];

      if(dim[0] == 0 || dim[1] == 0)
        continue;

      // edges
      auto matrix = this->GetInputArrayToProcess(0, c);
      switch(matrix->GetDataType()) {
        vtkTemplateMacro((countEdges<VTK_TT, int>(
          nEdges, dim, ttkUtils::GetPointer<VTK_TT>(matrix))));
      }
    }
  }

  return 1;
};

int collectArrays(vtkFieldData *arrayMap,
                  vtkMultiBlockDataSet *objects,
                  int attribute) {
  const int nSteps = objects->GetNumberOfBlocks();
  for(int t = 0; t < nSteps; t++) {
    auto block = objects->GetBlock(t);
    auto data = block->GetAttributesAsFieldData(attribute);
    for(int a = 0; a < data->GetNumberOfArrays(); a++)
      arrayMap->AddArray(data->GetAbstractArray(a));
  }

  return 1;
}

int ttkTrackingGraph::Validate(vtkMultiBlockDataSet *similarities,
                               vtkMultiBlockDataSet *features) {
  ttk::Timer timer;
  const std::string msg = "Validating Input";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const std::string errmsg0 = "Similarity Matrix input must be a flat "
                              "vtkMultiBlockDataSet that contains "
                              "only vtkImageData objects.";
  const std::string errmsg1 = "Features must be vtkPointSet where each feature "
                              "is represented by a single point.";

  if(!similarities)
    return !this->printErr(errmsg0);

  const int nSteps = similarities->GetNumberOfBlocks();
  for(int t = 0; t < nSteps; t++) {
    auto block = vtkImageData::SafeDownCast(similarities->GetBlock(t));
    if(!block)
      return !this->printErr(errmsg0);

    auto matrix = this->GetInputArrayToProcess(0, block);
    if(!matrix && block->GetNumberOfPoints() > 0)
      return !this->printErr("Unable to retrieve similarity matrix.");
  }

  if(features) {
    if(features->GetNumberOfBlocks() - 1 != similarities->GetNumberOfBlocks())
      return !this->printMsg(
        std::vector<std::string>({"Number of feature sets (F) and "
                                  "similarity matrices (C) inconsistent.",
                                  "They must satisfy F = C + 1."}),
        ttk::debug::Priority::ERROR);

    // automatically determine indexIdMapName

    std::string idArrayName;
    for(int t = 0; t < nSteps; t++) {
      idArrayName = ttkSimilarityAlgorithm::GetIdArrayName(
        similarities->GetBlock(t)->GetFieldData());
      if(idArrayName.size() < 1)
        return !this->printMsg(
          std::vector<std::string>(
            {"Similarity Matrices not augmented with IndexIdMaps.",
             "Unable to perform FeatureId lookup."}),
          ttk::debug::Priority::ERROR);
    }

    for(int t = 0; t <= nSteps; t++) {
      auto f = vtkPointSet::SafeDownCast(features->GetBlock(t));
      if(!f)
        return !this->printErr(errmsg1);
      auto ids = f->GetPointData()->GetArray(idArrayName.data());
      if(f->GetNumberOfPoints() > 0 && !ids)
        return !this->printErr("Unable to retrieve feature ids.");
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::GenerateTrackingGraphFromFeatures(
  vtkPolyData *trackingGraph,
  vtkMultiBlockDataSet *similarities,
  vtkMultiBlockDataSet *features) {

  // initialize trackingGraph
  ttk::Timer timer;
  std::string msg = "Initializing Output";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nSteps = similarities->GetNumberOfBlocks() + 1;

  int nNodes, nEdges;
  std::vector<int> nodeIdxOffsets;
  if(!this->CountNodesAndEdges(
       nNodes, nEdges, nodeIdxOffsets, similarities, features))
    return 0;

  if(nNodes < 1)
    return this->printWrn("Empty Input");

  // collect arrays to augment tracking graph
  auto featuresPD = vtkSmartPointer<vtkFieldData>::New();
  collectArrays(featuresPD, features, 0);
  auto featuresFD = vtkSmartPointer<vtkFieldData>::New();
  collectArrays(featuresFD, features, 2);
  auto similaritiesPD = vtkSmartPointer<vtkFieldData>::New();
  collectArrays(similaritiesPD, similarities, 0);

  // allocating memory
  auto trackingGraphPD = trackingGraph->GetPointData();
  auto trackingGraphCD = trackingGraph->GetCellData();

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);

  trackingGraph->SetPoints(points);
  trackingGraph->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

  {
    ArrayMapSet oSet{{trackingGraphPD, featuresPD, nNodes},
                     {trackingGraphPD, featuresFD, nNodes},
                     {trackingGraphCD, similaritiesPD, nEdges}};
    for(auto it : oSet) {
      auto arrayTemplates = std::get<1>(it);
      for(int a = 0; a < arrayTemplates->GetNumberOfArrays(); a++) {
        auto arrayTemplate = arrayTemplates->GetAbstractArray(a);
        auto arrayInstance = vtkSmartPointer<vtkAbstractArray>::Take(
          arrayTemplate->NewInstance());
        arrayInstance->SetName(arrayTemplate->GetName());
        arrayInstance->SetNumberOfComponents(
          arrayTemplate->GetNumberOfComponents());
        arrayInstance->SetNumberOfTuples(std::get<2>(it));
        std::get<0>(it)->AddArray(arrayInstance);
      }
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());
  timer.reStart();

  // ---------------------------------------------------------------------------
  msg = "Generating Spatial Tracking Graph";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  // Nodes
  {
    auto pointCoords = points->GetData();
    for(int t = 0, q = 0; t < nSteps; t++) {
      auto f = vtkPointSet::SafeDownCast(features->GetBlock(t));
      const int n = f->GetNumberOfPoints();
      if(n < 1)
        continue;
      pointCoords->InsertTuples(q, n, 0, f->GetPoints()->GetData());

      auto fPD = f->GetPointData();
      auto fFD = f->GetFieldData();
      for(int a = 0; a < trackingGraphPD->GetNumberOfArrays(); a++) {
        auto oArray = trackingGraphPD->GetAbstractArray(a);
        const int nComponents = oArray->GetNumberOfComponents();
        auto iArray = fPD->GetAbstractArray(oArray->GetName());
        if(iArray) {
          oArray->InsertTuples(q, n, 0, iArray);
        } else {
          iArray = fFD->GetAbstractArray(oArray->GetName());
          for(int i = 0; i < n; i++)
            oArray->InsertTuples(q + i, nComponents, 0, iArray);
        }
      }

      q += n;
    }
  }

  // Edges
  {
    for(int t = 1, q = 0; t < nSteps; t++) {
      auto c = static_cast<vtkImageData *>(similarities->GetBlock(t - 1));
      int dim[3];
      c->GetDimensions(dim);
      auto matrix = this->GetInputArrayToProcess(0, c);
      auto cFD = c->GetFieldData();

      const auto idArrayName = ttkSimilarityAlgorithm::GetIdArrayName(cFD);

      auto f0 = static_cast<vtkPointSet *>(features->GetBlock(t - 1));
      auto f1 = static_cast<vtkPointSet *>(features->GetBlock(t));
      if(f0->GetNumberOfPoints() < 1 || f1->GetNumberOfPoints() < 1)
        continue;

      auto l0 = f0->GetPointData()->GetArray(idArrayName.data());
      auto l1 = f1->GetPointData()->GetArray(idArrayName.data());

      vtkDataArray *indexIdMapP{nullptr};
      vtkDataArray *indexIdMapC{nullptr};
      if(!ttkSimilarityAlgorithm::GetIndexIdMaps(indexIdMapP, indexIdMapC, cFD))
        return !this->printErr("Unable to retrieve Index-Id-Maps.");

      ttkTypeMacroAI(
        matrix->GetDataType(), l0->GetDataType(),
        (generateEdges<T0, T1>(
          trackingGraph, q, trackingGraphCD, c->GetPointData(),
          ttkUtils::GetPointer<T0>(matrix), nodeIdxOffsets[t - 1],
          nodeIdxOffsets[t], dim, f0->GetNumberOfPoints(),
          f1->GetNumberOfPoints(), ttkUtils::GetPointer<T1>(l0),
          ttkUtils::GetPointer<T1>(l1), indexIdMapP, indexIdMapC)));
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::GenerateTrackingGraphFromMatrix(
  vtkPolyData *trackingGraph, vtkMultiBlockDataSet *similarities) {

  // initialize output
  ttk::Timer timer;
  std::string msg = "Initializing Output";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nSteps = similarities->GetNumberOfBlocks();

  int nNodes, nEdges;
  std::vector<int> nodeIdxOffsets;
  if(!this->CountNodesAndEdges(nNodes, nEdges, nodeIdxOffsets, similarities))
    return 0;

  auto similaritiesPD = vtkSmartPointer<vtkFieldData>::New();
  collectArrays(similaritiesPD, similarities, 0);

  auto similaritiesFD = vtkSmartPointer<vtkFieldData>::New();
  collectArrays(similaritiesFD, similarities, 2);

  // allocate memory
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);

  auto trackingGraphPD = trackingGraph->GetPointData();
  auto trackingGraphCD = trackingGraph->GetCellData();
  for(int a = 0; a < similaritiesPD->GetNumberOfArrays(); a++) {
    auto arrayTemplate = similaritiesPD->GetAbstractArray(a);
    auto arrayInstance
      = vtkSmartPointer<vtkAbstractArray>::Take(arrayTemplate->NewInstance());
    arrayInstance->SetName(arrayTemplate->GetName());
    arrayInstance->SetNumberOfComponents(
      arrayTemplate->GetNumberOfComponents());
    arrayInstance->SetNumberOfTuples(nEdges);
    trackingGraphCD->AddArray(arrayInstance);
  }

  vtkSmartPointer<vtkDataArray> ids;
  {
    vtkDataArray *indexIdMapP{nullptr};
    vtkDataArray *indexIdMapC{nullptr};
    if(!ttkSimilarityAlgorithm::GetIndexIdMaps(
         indexIdMapP, indexIdMapC, similaritiesFD))
      return !this->printErr("Unable to retrieve Index-Id-Maps.");

    std::string idArrayName
      = ttkSimilarityAlgorithm::GetIdArrayName(similaritiesFD);

    ids = vtkSmartPointer<vtkDataArray>::Take(indexIdMapP->NewInstance());
    ids->SetName(idArrayName.data());
    ids->SetNumberOfTuples(nNodes);
    trackingGraphPD->AddArray(ids);
  }

  auto timeIdx = vtkSmartPointer<vtkIntArray>::New();
  {
    timeIdx->SetName("TimeIdx");
    timeIdx->SetNumberOfTuples(nNodes);
    trackingGraphPD->AddArray(timeIdx);
  }

  trackingGraph->SetPoints(points);
  trackingGraph->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

  this->printMsg(msg, 1, timer.getElapsedTime());

  // Generate Tracking Graph
  timer.reStart();
  msg = "Generating Tracking Graph";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  auto pointCoordsData = ttkUtils::GetPointer<float>(points->GetData());
  auto timeIdxData = ttkUtils::GetPointer<int>(timeIdx);

  for(int t = 0, nodeIdx = 0, nodeIdx3 = 0, edgeIdx = 0; t < nSteps; t++) {
    auto c = static_cast<vtkImageData *>(similarities->GetBlock(t));
    int dim[3];
    c->GetDimensions(dim);

    auto cFD = c->GetFieldData();
    vtkDataArray *indexIdMapP{nullptr};
    vtkDataArray *indexIdMapC{nullptr};
    if(!ttkSimilarityAlgorithm::GetIndexIdMaps(indexIdMapP, indexIdMapC, cFD))
      return !this->printErr("Unable to retrieve Index-Id-Maps.");

    // nodes
    if(t == 0) {
      ids->InsertTuples(0, dim[0], 0, indexIdMapP);
      for(int i = 0; i < dim[0]; i++) {
        pointCoordsData[nodeIdx3++] = 0;
        pointCoordsData[nodeIdx3++] = i;
        pointCoordsData[nodeIdx3++] = 0;

        timeIdxData[nodeIdx++] = 0;
      }
    }

    ids->InsertTuples(nodeIdx, dim[1], 0, indexIdMapC);
    for(int i = 0; i < dim[1]; i++) {
      pointCoordsData[nodeIdx3++] = t + 1;
      pointCoordsData[nodeIdx3++] = i;
      pointCoordsData[nodeIdx3++] = 0;

      timeIdxData[nodeIdx++] = t + 1;
    }

    // edges
    if(dim[0] < 1 || dim[1] < 1)
      continue;

    auto cMatrix = this->GetInputArrayToProcess(0, c);
    if(!cMatrix)
      return !this->printErr("Unable to retrieve similarity matrix.");

    switch(cMatrix->GetDataType()) {
      vtkTemplateMacro((generateEdges<VTK_TT, int>(
        trackingGraph, edgeIdx, trackingGraphCD, c->GetPointData(),
        ttkUtils::GetPointer<VTK_TT>(cMatrix), nodeIdxOffsets[t + 0],
        nodeIdxOffsets[t + 1], dim)));
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::RequestData(vtkInformation *,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  // get input / output
  auto similarities = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto features = vtkMultiBlockDataSet::GetData(inputVector[1]);
  auto trackingGraph = vtkPolyData::GetData(outputVector);

  if(!this->Validate(similarities, features))
    return 0;

  if(features) {
    if(!this->GenerateTrackingGraphFromFeatures(
         trackingGraph, similarities, features))
      return 0;
  } else {
    if(!this->GenerateTrackingGraphFromMatrix(trackingGraph, similarities))
      return 0;
  }

  return 1;
}