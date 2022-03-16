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

  auto vtkglobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  if (vtkglobalPointIds != nullptr && vtkGhostCells != nullptr){
    vtkMPIController *controller = vtkMPIController::SafeDownCast(vtkMPIController::GetGlobalController());
    MPI_Datatype mpi_values;
    const int nitems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_INT, MPI_INT};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(ttk::value, scalar);
    offsets[1] = offsetof(ttk::value, globalId);
    offsets[2] = offsetof(ttk::value, localId);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_values);
    MPI_Type_commit(&mpi_values);

    int numProcs = controller->GetNumberOfProcesses();
    int rank = controller->GetLocalProcessId();
    int intTag = 100;
    int structTag = 101;
    if (rank == 0) this->printMsg("Global Point Ids and Ghost Cells exist, therefore we are in distributed mode!");

    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));

    // add the order array for every scalar array, except the ghostcells and the global ids
    for(auto scalarArray : scalarArrays) {
      std::string arrayName = std::string(scalarArray->GetName());
      if (arrayName != "GlobalPointIds" && arrayName != "vtkGhostType"){
                
        if (rank == 0) this->printMsg("Arrayname: " + arrayName);
        if (rank == 0) this->printMsg("Arraytype: " + std::to_string(scalarArray->GetDataType()));
        this->printMsg("#Points: " + std::to_string(nVertices));
        // sort the scalar array distributed first by the scalar value itself, then by the global id
        //std::vector<std::tuple<scalarArray->GetDataType(), int, int>> sortingValues;
        
        std::vector<ttk::value> sortingValues;
        /*
        ttkTypeMacroA(scalarArray->GetDataType(), 
                      (sortingValues = ttk::populateVector<T0>(nVertices,
                          ttkUtils::GetPointer<float>(scalarArray),
                          ttkUtils::GetPointer<ttk::SimplexId>(vtkglobalPointIds),
                          ttkUtils::GetPointer<char>(vtkGhostCells))));*/

        sortingValues = ttk::populateVector(nVertices,
                          ttkUtils::GetPointer<float>(scalarArray),
                          ttkUtils::GetPointer<ttk::SimplexId>(vtkglobalPointIds),
                          ttkUtils::GetPointer<char>(vtkGhostCells));
                      
        /*
        switch(scalarArray->GetDataType()) {
          vtkTemplateMacro(ttk::populateVector(
            nVertices, 
            static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(scalarArray)),
            static_cast<int *>(ttkUtils::GetVoidPointer(vtkglobalPointIds)),
            static_cast<int *>(ttkUtils::GetVoidPointer(vtkGhostCells)),
            sortingValues));
        }*/


        /*
        for (int i = 0; i < nVertices; i++){
          if (vtkGhostCells->GetComponent(i, 0) == 0){
            double scalarValue = scalarArray->GetComponent(i, 0);
            int globalId = vtkglobalPointIds->GetComponent(i, 0);
            int localId = i;
            sortingValues.emplace_back(scalarValue, globalId, localId);
          }
        }*/
        auto element0 = sortingValues[0];
        this->printMsg("Rank " +  std::to_string(rank)
                      + ", #Elements in Vector " + std::to_string(sortingValues.size()) 
                      + ", first values before sort are " + std::to_string(element0.scalar) + " " 
                      + std::to_string(element0.globalId) + " " + std::to_string(element0.localId));

        // sort the vector of this rank, first by their scalarvalue and if they are the same, by their globalId
        /*
        std::sort(sortingValues.begin(), sortingValues.end(),
        [&](const std::tuple<double, int, int> a, const std::tuple<double, int, int> b) {
          return (std::get<0>(a) < std::get<0>(b))
                 || (std::get<0>(a) == std::get<0>(a) && std::get<1>(a) < std::get<1>(b));
        });
        */

        //ttkTypeMacroA(scalarArray->GetDataType(), ttk::sortVerticesDistributed<T0>(sortingValues));
        ttk::sortVerticesDistributed(sortingValues);

        element0 = sortingValues[0];
        this->printMsg("Rank " +  std::to_string(rank)
                      + ", #Elements in Vector " + std::to_string(sortingValues.size()) 
                      + ", first values after sort are " + std::to_string(element0.scalar) + " " 
                      + std::to_string(element0.globalId) + " " + std::to_string(element0.localId));

        // when all are done sorting, rank 0 requests the highest values and merges them
        controller->Barrier();
        //int BurstSize = SetBurstSize();
        if (rank == 0){
            this->printMsg("Rank 0 starts merging");
            int totalSize = 0;
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
            this->printMsg("Total amount of distributed points: " + std::to_string(totalSize)); 
  	  	    std::vector<std::vector<ttk::value>> finalValues;
            //finalValues.resize(nVertices, {0,0,0});


            for (int i = 0; i < numProcs; i++){
              if (i == 0){
                std::vector<ttk::value> sendValues = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
                finalValues.push_back(sendValues);
              } else{
                std::vector<ttk::value> receivedValues;
                receivedValues.resize(BurstSize, {0,0,0});
                MPI_Recv(receivedValues.data(), BurstSize, mpi_values, i, structTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                finalValues.push_back(receivedValues);
              }
            }
            this->printMsg("Size of finalValues: " + std::to_string(finalValues.size()));
            this->printMsg("Size of first vector in finalValues: " + std::to_string(finalValues[0].size()));

            // initally gather #BurstSize of elements from each rank
            //controller->Gather(sendValues.data(), finalValues.data(), BurstSize, 0);
            //MPI_Gather(sendValues.data(), BurstSize, mpi_values, finalValues.data(), BurstSize, mpi_values, 0, MPI_COMM_WORLD);
            //finalValues.resize(BurstSize * numProcs);

            /*
            while (finalValues.size() < totalSize){
              //request points from the different ranks 

            }*/
          


        } else {
          vtkIdType values = 1;
          int nValues = sortingValues.size();
          controller->Send(&nValues, values, 0, intTag);
          std::vector<ttk::value> sendValues = ttk::returnVectorForBurstsize(sortingValues, BurstSize);
          MPI_Send(sendValues.data(), BurstSize, mpi_values, 0, structTag, MPI_COMM_WORLD);
          //ttkTypeMacroA(scalarArray->GetDataType(), sendValues = ttk::returnVectorForBurstsize<T0>(sortingValues, BurstSize));

        }

      }
        
      //output->GetPointData()->AddArray(orderArray);
      //this->printMsg("Generated order array for scalar array `" + std::string{scalarArray->GetName()} + "'");

    }

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
