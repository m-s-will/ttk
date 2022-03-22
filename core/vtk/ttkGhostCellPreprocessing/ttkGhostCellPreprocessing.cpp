#include <ttkGhostCellPreprocessing.h>
#include <ttkMacros.h>
#include <ttkUtils.h>
#include <mpi.h>

#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkMPIController.h>
#include <unordered_map>

vtkStandardNewMacro(ttkGhostCellPreprocessing);

ttkGhostCellPreprocessing::ttkGhostCellPreprocessing() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->setDebugMsgPrefix("GhostCellPreprocessing");

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkGhostCellPreprocessing::Modified);
}

int ttkGhostCellPreprocessing::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkGhostCellPreprocessing::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}



int ttkGhostCellPreprocessing::RequestData(vtkInformation *ttkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkDataSet::GetData(outputVector);
  ttk::Timer tm{};

  if(input == nullptr || output == nullptr) {
    return 0;
  }

  output->ShallowCopy(input);

  auto pointData = input->GetPointData();
  int nVertices = input->GetNumberOfPoints();

  

  auto vtkGlobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  if (vtkGlobalPointIds != nullptr && vtkGhostCells != nullptr){
    vtkMPIController *controller = vtkMPIController::SafeDownCast(vtkMPIController::GetGlobalController());
    int numProcs = controller->GetNumberOfProcesses();
    int rank = controller->GetLocalProcessId();
    if (rank == 0) this->printMsg("Global Point Ids and Ghost Cells exist, therefore we can continue!");
    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));
    this->printMsg("#Points: " + std::to_string(nVertices));

    vtkNew<vtkIntArray> rankArray{};
    rankArray->SetName("RankArray");
    rankArray->SetNumberOfComponents(1);
    rankArray->SetNumberOfTuples(nVertices);
    std::vector<int> currentRankUnknownIds;
    std::vector<std::vector<int>> allUnknownIds(numProcs);
    std::unordered_map<int, int> idToRankMap;
    std::unordered_map<int, int> gIdToLocalMap;
    for (int i = 0; i < nVertices; i++){
      int ghostCellVal = vtkGhostCells->GetComponent(i, 0);
      int globalId = vtkGlobalPointIds->GetComponent(i, 0);
      if (ghostCellVal == 0){
        // if the ghost cell value is 0, then this vertex mainly belongs to this rank
        rankArray->SetComponent(i, 0, rank);
        idToRankMap[globalId] = rank;
      } else{
        // otherwise the vertex belongs to another rank and we need to find out to which one
        // this needs to be done by broadcasting the global id and hoping for some other rank to answer
        currentRankUnknownIds.push_back(globalId);
        gIdToLocalMap[globalId] = i;
      }
    }
    allUnknownIds[rank] = currentRankUnknownIds;
    int sizeOfCurrentRank;
    // first each rank gets the information which rank needs which globalid
    for (int r = 0; r < numProcs; r++){
      if (r == rank) sizeOfCurrentRank = currentRankUnknownIds.size();
      MPI_Bcast(&sizeOfCurrentRank, 1, MPI_INT, r, MPI_COMM_WORLD);
      allUnknownIds[r].resize(sizeOfCurrentRank);
      //this->printMsg("This is rank " + std::to_string(rank) + ", rank " + std::to_string(r) + " needs values for " + std::to_string(sizeOfCurrentRank) + " vertices.");
      MPI_Bcast(allUnknownIds[r].data(), sizeOfCurrentRank, MPI_INT, r, MPI_COMM_WORLD);
    }
    
    // then we check if the needed globalid values are present in the local globalid map
    // if so, we send the rank value to the requesting rank
    for (int r = 0; r < numProcs; r++){
      for (int gId : allUnknownIds[r]){
        if (idToRankMap.count(gId)){
          //send the value back
          MPI_Send(&gId, 1, MPI_INT, r, 101, MPI_COMM_WORLD);
        }
      }
    }

    // when we receive values, fill the map further, to get values for other ranks more quickly?
    // only exactly 1 send / recv operation per value?
    // MPI_Recv (value blala)
    // save where the value should have been? Struct?
    // rankArray->SetComponent(which one?, 0, value);
    // idToRankMap[globalId] = value;

    for (size_t i = 0; i < allUnknownIds[rank].size(); i++){
      int receivedGlobal;
      MPI_Status status;
      MPI_Recv(&receivedGlobal, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      int sourceRank = status.MPI_SOURCE;
      int localVal = gIdToLocalMap[receivedGlobal];
      rankArray->SetComponent(localVal, 0, sourceRank);
    }


    output->GetPointData()->AddArray(rankArray);



    return 1;
  } else{
    this->printMsg("Either Global Point Ids or Ghost Cells don't exist, returning.");
    return 1;
  }

}

 