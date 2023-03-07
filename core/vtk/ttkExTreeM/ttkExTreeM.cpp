#include <ttkExTreeM.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkPolyData.h>

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

#include <PathCompression.h>
#include <ScalarFieldCriticalPoints2.h>

vtkStandardNewMacro(ttkExTreeM);

ttkExTreeM::ttkExTreeM() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

ttkExTreeM::~ttkExTreeM() {
}

int ttkExTreeM::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkExTreeM::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      return 1;
    case 1:
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
      return 1;
    default:
      return 0;
  }
}

int ttkExTreeM::RequestData(vtkInformation *,
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
  this->preconditionTriangulation(triangulation);

  // Get input array
  auto scalarArray = this->GetInputArrayToProcess(0, inputVector);
  if(!scalarArray) {
    this->printErr("Unable to retrieve scalar array.");
    return 0;
  }

  // Order Array
  auto orderArray = this->GetOrderArray(input, 0);
  auto orderArrayData = ttkUtils::GetPointer<ttk::SimplexId>(orderArray);

  vtkNew<ttkSimplexIdTypeArray> ascendingManifold{};
  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(nVertices);
  ascendingManifold->SetName(ttk::MorseSmaleAscendingName);

  vtkNew<ttkSimplexIdTypeArray> descendingManifold{};
  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(nVertices);
  descendingManifold->SetName(ttk::MorseSmaleDescendingName);

  vtkNew<ttkSimplexIdTypeArray> segmentationId{};
  segmentationId->SetNumberOfComponents(1);
  segmentationId->SetNumberOfTuples(nVertices);
  segmentationId->SetName("SegmentationId");
  auto segmentationIdData = ttkUtils::GetPointer<ttk::SimplexId>(segmentationId);

  // compute path compression
  {
    ttk::PathCompression subModule;
    subModule.setThreadNumber(this->threadNumber_);
    subModule.setDebugLevel(this->debugLevel_);
    subModule.setComputeSegmentation(true,true,false);

    ttk::PathCompression::OutputManifold om{
      ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
      ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
      nullptr
    };

    int status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (status = subModule.execute<T0>(
        om,
        ttkUtils::GetPointer<const ttk::SimplexId>(orderArray),
        *static_cast<T0*>(triangulation->getData())
      ))
    );
    if(status!=0)
      return 0;
  }

  // compute critical points
  // type -> thread -> cp
  std::array<std::vector<std::vector<ttk::ScalarFieldCriticalPoints2::CP>>,4> criticalPoints;
  {
    ttk::ScalarFieldCriticalPoints2 subModule;
    subModule.setThreadNumber(this->threadNumber_);
    subModule.setDebugLevel(this->debugLevel_);

    int status = 0;
    ttkTypeMacroT(
      triangulation->getType(),
      (
        status = subModule.computeCritialPoints<T0>(
          criticalPoints,
          ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
          ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
          ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
          static_cast<const T0 *>(triangulation->getData())
        )
      )
    );
    if(status != 1)
      return 0;
  }

  // debug: build extremum graph
  {
    auto tree = vtkPolyData::GetData(outputVector, 0);

    // points
    ttk::SimplexId maxO = nVertices-1;
    ttk::SimplexId maxId = 0;
    ttk::SimplexId minId = 0;
    {
      int nPoints = 0;
      for(int i=2; i<4; i++){
        const auto& cp = criticalPoints[i];
        const int nThreads = cp.size();
        for(int j=0; j<nThreads; j++){
          nPoints += cp[j].size();
        }
      }
      nPoints++; // make room for global minimum

      // init point data
      vtkNew<ttkSimplexIdTypeArray> treeOrder;
      treeOrder->SetNumberOfTuples(nPoints);
      treeOrder->SetName(orderArray->GetName());
      auto treeOrderData = ttkUtils::GetPointer<ttk::SimplexId>(treeOrder);
      tree->GetPointData()->AddArray(treeOrder);

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetDataTypeToFloat();
      points->SetNumberOfPoints(nPoints);
      tree->SetPoints(points);
      auto pointCoords = static_cast<float*>(points->GetData()->GetVoidPointer(0));
      int pointCursor = 0;

      auto addPoint = [=](int& pc, const ttk::SimplexId& v){
        treeOrderData[pc] = orderArrayData[v];
        auto pc3 = pc*3;
        segmentationIdData[v] = pc;
        triangulation->getVertexPoint(
          v,
          pointCoords[pc3],
          pointCoords[pc3+1],
          pointCoords[pc3+2]
        );
        pc++;
      };

      // maxima
      {
        const auto& cp = criticalPoints[3];
        const int nThreads = cp.size();
        for(int j=0; j<nThreads; j++){
          const auto& cp_ = cp[j];
          const int n = cp_.size();
          for(int k=0; k<n; k++){
            const auto& v = cp_[k].id;
            if(orderArrayData[v]==maxO)
              maxId = v;
            addPoint(pointCursor, v);
          }
        }
      }

      // saddles
      {
        const auto& cp = criticalPoints[2];
        const int nThreads = cp.size();
        for(int j=0; j<nThreads; j++){
          const auto& cp_ = cp[j];
          const int n = cp_.size();
          for(int k=0; k<n; k++){
            addPoint(pointCursor, cp_[k].id);
          }
        }
      }

      // find and add global minimum
      {
        const auto& cp = criticalPoints[0];
        const int nThreads = cp.size();
        [=](int& pc, ttk::SimplexId& id){
          for(int j=0; j<nThreads; j++){
            const auto& cp_ = cp[j];
            const int n = cp_.size();
            for(int k=0; k<n; k++){
              const auto& v = cp_[k].id;
              if(orderArrayData[v]!=0)
                continue;
              id = v;
              addPoint(pc, v);
              return;
            }
          }
        }(pointCursor,minId);
      }
    }

    // edges
    {
      int nEdges = 0;
      {
        const auto& cp = criticalPoints[2];
        const int nThreads = cp.size();
        for(int j=0; j<nThreads; j++){
          const auto& cp_ = cp[j];
          const int n = cp_.size();
          for(int k=0; k<n; k++){
            nEdges+=cp_[k].neighbors.size();
          }
        }
      }
      nEdges++; // make room for main branch
      tree->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

      // add edges
      {
        // int edgeCursor = 0;
        const auto& cp = criticalPoints[2];
        const int nThreads = cp.size();
        for(int j=0; j<nThreads; j++){
          const auto& cp_ = cp[j];
          const int n = cp_.size();
          for(int k=0; k<n; k++){
            const auto& v = segmentationIdData[cp_[k].id];
            const auto& neighbors = cp_[k].neighbors;
            const int nNeighbors = neighbors.size();
            for(int l=0; l<nNeighbors; l++){
              vtkIdType points[2]{
                v,
                segmentationIdData[neighbors[l]]
              };
              tree->InsertNextCell(VTK_LINE, 2, points);
            }
          }
        }
      }

      // add main branch
      {
        vtkIdType points[2]{
          segmentationIdData[minId],
          segmentationIdData[maxId]
        };
        tree->InsertNextCell(VTK_LINE, 2, points);
      }
    }


  }




  // // Init segmentation
  // auto segmentationIds
  //   = vtkSmartPointer<vtkDataArray>::Take(orderArray->NewInstance());
  // segmentationIds->SetName("BranchId");
  // segmentationIds->SetNumberOfComponents(1);
  // segmentationIds->SetNumberOfTuples(nVertices);

  // // Compute merge tree segmentation
  // std::vector<ttk::mt::Propagation<ttk::SimplexId>> propagations;
  // std::vector<const ttk::mt::Propagation<ttk::SimplexId> *> sortedPropagations;
  // {
  //   int status = 0;
  //   ttkTypeMacroT(
  //     triangulation->getType(),
  //     (status = this->computeMergeTreeSegmentation<ttk::SimplexId, T0>(
  //       ttkUtils::GetPointer<ttk::SimplexId>(segmentationIds), propagations,

  //       static_cast<T0 *>(triangulation->getData()), orderArrayData,
  //       this->GetType())));
  //   if(!status)
  //     return 0;

  //   ttkTypeMacroAT(
  //     scalarArray->GetDataType(), triangulation->getType(),
  //     (status = this->finalizePropagations<T0, ttk::SimplexId, T1>(
  //       sortedPropagations,
  //       ttkUtils::GetPointer<ttk::SimplexId>(segmentationIds), propagations,

  //       static_cast<T1 *>(triangulation->getData()),
  //       ttkUtils::GetPointer<T0>(scalarArray))));
  //   if(!status)
  //     return 0;
  // }

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
      segmentationPD->AddArray(ascendingManifold);
      segmentationPD->AddArray(descendingManifold);
      segmentationPD->AddArray(segmentationId);
    }

    this->printMsg("Generating Output Data Objects", 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1);
  }

  return 1;
}
