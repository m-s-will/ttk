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

  

  auto vtkglobalPointIds = pointData->GetGlobalIds();
  auto vtkGhostCells = pointData->GetArray("vtkGhostType");
  if (vtkglobalPointIds != nullptr && vtkGhostCells != nullptr){
    vtkMPIController *controller = vtkMPIController::SafeDownCast(vtkMPIController::GetGlobalController());
    int numProcs = controller->GetNumberOfProcesses();
    int rank = controller->GetLocalProcessId();
    if (rank == 0) this->printMsg("Global Point Ids and Ghost Cells exist, therefore we can continue!");
    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));

    vtkNew<vtkIntArray> rankArray{};
    rankArray->SetName("RankArray");
    rankArray->SetNumberOfComponents(1);
    rankArray->SetNumberOfTuples(nVertices);
    for (int i = 0; i < nVertices; i++){
      int ghostCellVal = vtkGhostCells->GetComponent(i, 0);
      if (ghostCellVal == 0){
        // if the ghost cell value is 0, then this vertex mainly belongs to this rank
        rankArray->SetComponent(i, 0, rank);
      } else{
        // otherwise the vertex belongs to another rank and we need to find out to which one
        // this needs to be done by broadcasting the global id and hoping for some other rank to answer
      }
    }
    output->GetPointData()->AddArray(orderArray);



    return 1;
  } else{
    if (rank == 0) this->printMsg("Either Global Point Ids or Ghost Cells don't exist, returning.");
    return 1;
  }

}

 