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
#include <unordered_set>

vtkStandardNewMacro(ttkGhostCellPreprocessing);
using IT = long long int;

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
  IT nVertices = input->GetNumberOfPoints();

  

  auto vtkGlobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  if (vtkGlobalPointIds != nullptr && vtkGhostCells != nullptr){
    vtkMPIController *controller = vtkMPIController::SafeDownCast(vtkMPIController::GetGlobalController());
    int numProcs = controller->GetNumberOfProcesses();
    int rank = controller->GetLocalProcessId();
    if (rank == 0) this->printMsg("Global Point Ids and Ghost Cells exist, therefore we can continue!");
    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));
    this->printMsg("#Points: " + std::to_string(nVertices));
    MPI_Datatype MIT = MPI_LONG_LONG_INT;

    vtkNew<vtkIntArray> rankArray{};
    rankArray->SetName("RankArray");
    rankArray->SetNumberOfComponents(1);
    rankArray->SetNumberOfTuples(nVertices);
    std::vector<IT> currentRankUnknownIds;
    std::vector<std::vector<IT>> allUnknownIds(numProcs);
    std::unordered_set<IT> gIdSet;
    std::unordered_map<IT, IT> gIdToLocalMap;
    for (int i = 0; i < nVertices; i++){
      IT ghostCellVal = vtkGhostCells->GetComponent(i, 0);
      IT globalId = vtkGlobalPointIds->GetComponent(i, 0);
      if (ghostCellVal == 0){
        // if the ghost cell value is 0, then this vertex mainly belongs to this rank
        rankArray->SetComponent(i, 0, rank);
        gIdSet.insert(globalId);
      } else{
        // otherwise the vertex belongs to another rank and we need to find out to which one
        // this needs to be done by broadcasting the global id and hoping for some other rank to answer
        currentRankUnknownIds.push_back(globalId);
        gIdToLocalMap[globalId] = i;
      }
    }
    allUnknownIds[rank] = currentRankUnknownIds;
    IT sizeOfCurrentRank;
    // first each rank gets the information which rank needs which globalid
    for (int r = 0; r < numProcs; r++){
      if (r == rank) sizeOfCurrentRank = currentRankUnknownIds.size();
      MPI_Bcast(&sizeOfCurrentRank, 1, MIT, r, MPI_COMM_WORLD);
      allUnknownIds[r].resize(sizeOfCurrentRank);
      //this->printMsg("This is rank " + std::to_string(rank) + ", rank " + std::to_string(r) + " needs values for " + std::to_string(sizeOfCurrentRank) + " vertices.");
      MPI_Bcast(allUnknownIds[r].data(), sizeOfCurrentRank, MIT, r, MPI_COMM_WORLD);
    }
    
    // then we check if the needed globalid values are present in the local globalid map
    // if so, we send the rank value to the requesting rank
    for (int r = 0; r < numProcs; r++){
      std::vector<IT> gIdsToSend; 
      for (IT gId : allUnknownIds[r]){
        if (gIdSet.count(gId)){
          //add the value to the vector which will be sent
          gIdsToSend.push_back(gId);
        }
      }
      //send whole vector of data
      MPI_Send(gIdsToSend.data(), gIdsToSend.size(), MIT, r, 101, MPI_COMM_WORLD);

    }

    // receive a variable amount of values from different ranks
    size_t i = 0;
    while (i < allUnknownIds[rank].size()){
      std::vector<IT> receivedGlobals;
      receivedGlobals.resize(allUnknownIds[rank].size());
      MPI_Status status;
      int amount;
      MPI_Recv(receivedGlobals.data(), allUnknownIds[rank].size(), MIT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      int sourceRank = status.MPI_SOURCE;
      MPI_Get_count(&status, MIT, &amount);
      receivedGlobals.resize(amount);
      for (IT receivedGlobal : receivedGlobals){
        IT localVal = gIdToLocalMap[receivedGlobal];
        rankArray->SetComponent(localVal, 0, sourceRank);
        i++;
      }
    }


    output->GetPointData()->AddArray(rankArray);

    this->printMsg("Preprocessed RankArray", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);


    return 1;
  } else{
    this->printMsg("Either Global Point Ids or Ghost Cells don't exist, returning.");
    return 1;
  }

}

 