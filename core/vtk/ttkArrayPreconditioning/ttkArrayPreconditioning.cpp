#include <OrderDisambiguation.h>
#include <ttkArrayPreconditioning.h>
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
#include <regex>
#include <unordered_map>

vtkStandardNewMacro(ttkArrayPreconditioning);

ttkArrayPreconditioning::ttkArrayPreconditioning() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->setDebugMsgPrefix("ArrayPreconditioning");

  // ensure that modifying the selection re-triggers the filter
  // (c.f. vtkPassSelectedArrays.cxx)
  this->ArraySelection->AddObserver(
    vtkCommand::ModifiedEvent, this, &ttkArrayPreconditioning::Modified);
}

int ttkArrayPreconditioning::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkArrayPreconditioning::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}


void ttkArrayPreconditioning::ReceiveAndAddToVector(int burstSize, MPI_Datatype mpi_values, int rankFrom, int tag, std::vector<std::vector<ttk::value>> &unsortedReceivedValues){
  std::vector<ttk::value> receivedValues;
  // be prepared to receive burstsize of elements, resize after receiving to the correct size
  receivedValues.resize(burstSize, {0,0,0});
  MPI_Status status;
  int amount;
  MPI_Recv(receivedValues.data(), burstSize, mpi_values, rankFrom, tag, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, mpi_values, &amount);
  //this->printMsg("Received " + std::to_string(amount) + " values from rank " + std::to_string(i));
  receivedValues.resize(amount, {0,0,0});
  unsortedReceivedValues[rankFrom] = receivedValues;
}

int ttkArrayPreconditioning::RequestData(vtkInformation *ttkNotUsed(request),
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

  std::vector<vtkDataArray *> scalarArrays{};

  if(SelectFieldsWithRegexp) {
    // select all input point data arrays whose name is matching the regexp
    const auto n = pointData->GetNumberOfArrays();
    for(int i = 0; i < n; ++i) {
      auto array = pointData->GetArray(i);
      if(array != nullptr && array->GetName() != nullptr
         && std::regex_match(array->GetName(), std::regex(RegexpString))) {
        scalarArrays.emplace_back(array);
      }
    }
  } else {
    // get all selected input point data arrays
    for(int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
      auto array = pointData->GetArray(i);
      if(array != nullptr && array->GetName() != nullptr
         && ArraySelection->ArrayIsEnabled(array->GetName())) {
        scalarArrays.emplace_back(array);
      }
    }
  }

  auto vtkGlobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  if (vtkGlobalPointIds != nullptr && vtkGhostCells != nullptr){
    vtkMPIController *controller = vtkMPIController::SafeDownCast(vtkMPIController::GetGlobalController());
    MPI_Datatype mpi_values;
    const int nitems = 4;
    int blocklengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_FLOAT, MPI_LONG, MPI_INT, MPI_INT};
    MPI_Aint offsets[4];
    offsets[0] = offsetof(ttk::value, scalar);
    offsets[1] = offsetof(ttk::value, globalId);
    offsets[2] = offsetof(ttk::value, localId);
    offsets[3] = offsetof(ttk::value, ordering);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_values);
    MPI_Type_commit(&mpi_values);

    int numProcs = controller->GetNumberOfProcesses();
    int rank = controller->GetLocalProcessId();
    int intTag = 100;
    int structTag = 101;
    int boolTag = 102;
    if (rank == 0) this->printMsg("Global Point Ids and Ghost Cells exist, therefore we are in distributed mode!");
    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));

    // add the order array for every scalar array, except the ghostcells and the global ids
    for(auto scalarArray : scalarArrays) {
      std::string arrayName = std::string(scalarArray->GetName());
      if (arrayName != "GlobalPointIds" && arrayName != "vtkGhostType"){
                
        if (rank == 0) this->printMsg("Arrayname: " + arrayName);
        if (rank == 0) this->printMsg("Arraytype: " + std::to_string(scalarArray->GetDataType()));
        if (rank == 0) this->printMsg("#Points: " + std::to_string(nVertices));
        ttk::Timer fillAndSortTimer;
        std::vector<ttk::value> sortingValues;
        sortingValues = ttk::populateVector(nVertices,
                        ttkUtils::GetPointer<float>(scalarArray),
                        ttkUtils::GetPointer<long int>(vtkGlobalPointIds),
                        ttkUtils::GetPointer<char>(vtkGhostCells));
          
             
        // sort the scalar array distributed first by the scalar value itself, then by the global id
        ttk::sortVerticesDistributed(sortingValues);

        // when all are done sorting, rank 0 requests the highest values and merges them
        controller->Barrier();
        if (rank == 0) {
          this->printMsg("Filling vector and sorting for each rank done, starting merge.", 
                          1,
                          fillAndSortTimer.getElapsedTime());
        }

        std::vector<ttk::value> orderedValuesForRank;  	  	    
        std::vector<ttk::value> finalValues;
        ttk::Timer mergeTimer;
        size_t totalSize = 0;
        if (rank == 0){
            this->printMsg("Rank 0 starts merging");
            // get the nVertices from each rank, add them to get the complete size of the dataset 
            for (int i = 0; i < numProcs; i++){
              if (i == 0){
                totalSize +=sortingValues.size();
              } else{
                int receivedSize;
                vtkIdType values = 1;
                controller->Receive(&receivedSize, values, i, intTag);
                totalSize += receivedSize;
              }
            }
            MPI_Bcast(&totalSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
            this->printMsg("Total amount of distributed points: " + std::to_string(totalSize)); 
            int currentOrder = totalSize;
            std::vector<std::vector<ttk::value>> unsortedReceivedValues;
            unsortedReceivedValues.resize(numProcs);
            std::vector<std::vector<ttk::value>> orderResendValues;
            orderResendValues.resize(numProcs);


            // receive the first batch of values
            for (int i = 0; i < numProcs; i++){
              if (i == 0){
                std::vector<ttk::value> ownValues = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
                unsortedReceivedValues[i] = ownValues;
              } else{
                this->ReceiveAndAddToVector(BurstSize, mpi_values, i, structTag, unsortedReceivedValues);
              }
            }
            
            while (finalValues.size() < totalSize){
              //take the current maximum scalar over all ranks 
              int rankIdOfMaxScalar = -1;
              float maxScalar = -FLT_MAX;
              for (int i = 0; i < numProcs; i++){
                if (unsortedReceivedValues[i].size() > 0){
                  int thisId = i;
                  float thisScalar = unsortedReceivedValues[i].back().scalar;
                  if (thisScalar > maxScalar){
                    maxScalar = thisScalar;
                    rankIdOfMaxScalar = thisId;
                  }
                }
              }
              if (rankIdOfMaxScalar == -1){
                this->printMsg("All vectors are empty, but out final vector is not complete yet. Either something went wrong or some rank didn't send their values yet.");
                return 0;
              }

              // move the struct from the unsortedReceivedValues subvector to the finalValues vector to get an ordering
              ttk::value currentValue = unsortedReceivedValues[rankIdOfMaxScalar].back();
              currentValue.ordering = currentOrder;
              currentOrder--;
              orderResendValues[rankIdOfMaxScalar].push_back(currentValue);
              finalValues.push_back(currentValue);
              unsortedReceivedValues[rankIdOfMaxScalar].pop_back();
              if (unsortedReceivedValues[rankIdOfMaxScalar].size() == 0){
                //this->printMsg("Vector for Rank " + std::to_string(rankIdOfMaxScalar) + " is empty, we either need to receive more or are done");
                
                if (rankIdOfMaxScalar == 0){
                  // append the ordered values to the correct vector
                  orderedValuesForRank.insert(orderedValuesForRank.end(), orderResendValues[rankIdOfMaxScalar].begin(), orderResendValues[rankIdOfMaxScalar].end());
                  orderResendValues[rankIdOfMaxScalar].clear();
                  if (sortingValues.size() > 0){
                    std::vector<ttk::value> ownValues = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
                    unsortedReceivedValues[rankIdOfMaxScalar] = ownValues;
                  } else {
                    this->printMsg("We are done with rank 0!");
                  }
                } else {
                  // receive more values from rank, send ordering to the rank
                  // send to the finished rank that we want more
                  MPI_Send(orderResendValues[rankIdOfMaxScalar].data(), orderResendValues[rankIdOfMaxScalar].size(), mpi_values, rankIdOfMaxScalar, structTag, MPI_COMM_WORLD);
                  orderResendValues[rankIdOfMaxScalar].clear();
                  //check if there are more values to be received. If so, receive them
                  bool moreVals = false;
                  MPI_Recv(&moreVals, 1, MPI_CXX_BOOL, rankIdOfMaxScalar, boolTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  if (moreVals){
                    //this->printMsg("ANOTHER ONE " + std::to_string(rankIdOfMaxScalar));
                    this->ReceiveAndAddToVector(BurstSize, mpi_values, rankIdOfMaxScalar, structTag, unsortedReceivedValues);
                  } else {
                    this->printMsg("We are done with rank " + std::to_string(rankIdOfMaxScalar));
                  }
                }
              }
            }

            this->printMsg("Finished with sorting, max value is " + std::to_string(finalValues[0].scalar) + ", min value is " + std::to_string(finalValues.back().scalar));
        } else {
          int nValues = sortingValues.size();
          controller->Send(&nValues, 1, 0, intTag);
          MPI_Bcast(&totalSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
          finalValues.resize(totalSize, {0,0,0});

          // send the next burstsize values and then wait for an answer from the root rank
          while (sortingValues.size() > 0){
            std::vector<ttk::value> sendValues = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
            MPI_Send(sendValues.data(), sendValues.size(), mpi_values, 0, structTag, MPI_COMM_WORLD);
            std::vector<ttk::value> receivedValues;

            // be prepared to receive burstsize of elements, resize after receiving to the correct size
            receivedValues.resize(BurstSize, {0,0,0});
            MPI_Status status;
            int amount;
            MPI_Recv(receivedValues.data(), BurstSize, mpi_values, 0, structTag, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, mpi_values, &amount);
            receivedValues.resize(amount, {0,0,0});
            orderedValuesForRank.insert(orderedValuesForRank.end(), receivedValues.begin(), receivedValues.end());

            // afterwards send to root if there are still values to be sent
            bool moreVals = sortingValues.size() > 0;
            MPI_Send(&moreVals, 1, MPI_CXX_BOOL, 0, boolTag, MPI_COMM_WORLD);
          }


        }

        // all ranks do the following
        controller->Barrier();
        if (rank == 0) {
          this->printMsg("Merging done, sending results to ranks and constructing order array.", 
                          1,
                          mergeTimer.getElapsedTime());
        }
        ttk::Timer sendTimer;        
        MPI_Bcast(finalValues.data(), totalSize, mpi_values, 0, MPI_COMM_WORLD);
        this->printMsg("Sent results, generating order array for rank " + std::to_string(rank), 
                   1,
                   sendTimer.getElapsedTime());
        ttk::Timer orderTimer;        
        std::unordered_map<float, int> orderMap = ttk::buildOrderMap(finalValues, totalSize);
        // every rank now has an orderedValuesForRank array with the points sorted in descending order and their correct order
        // now we need to transform this to a correct vtk orderarray and append it

        vtkNew<ttkSimplexIdTypeArray> orderArray{};
        orderArray->SetName(this->GetOrderArrayName(scalarArray).data());
        orderArray->SetNumberOfComponents(1);
        orderArray->SetNumberOfTuples(nVertices);
        for (int i = 0; i < nVertices; i++){
          float scalarVal = scalarArray->GetComponent(i, 0);
          int orderVal = orderMap[scalarVal];
          orderArray->SetComponent(i, 0, orderVal);
        }
        output->GetPointData()->AddArray(orderArray);
        this->printMsg("Generated order array for scalar array `"
                   + std::string{scalarArray->GetName()} + "', rank " + std::to_string(rank), 
                   1,
                   orderTimer.getElapsedTime());
      }

    }
    this->printMsg("Preconditioned selected scalar arrays", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);
    return 1;
  }


  

  for(auto scalarArray : scalarArrays) {
    vtkNew<ttkSimplexIdTypeArray> orderArray{};
    orderArray->SetName(this->GetOrderArrayName(scalarArray).data());
    orderArray->SetNumberOfComponents(1);
    orderArray->SetNumberOfTuples(nVertices);

    switch(scalarArray->GetDataType()) {
      vtkTemplateMacro(ttk::sortVertices(
        nVertices, static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
        static_cast<int *>(nullptr),
        static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(orderArray)),
        this->threadNumber_));
    }

    output->GetPointData()->AddArray(orderArray);
    this->printMsg("Generated order array for scalar array `"
                   + std::string{scalarArray->GetName()} + "'");
  }

  this->printMsg("Preconditioned selected scalar arrays", 1.0,
                 tm.getElapsedTime(), this->threadNumber_);

  return 1;
}
