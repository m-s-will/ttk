#include <OrderDisambiguation.h>
#include <ttkArrayPreconditioning.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

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
  auto nVertices = input->GetNumberOfPoints();

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

  bool globalPointIdsExist = false;
  bool ghostCellsExist = false;
  for (auto scalarArray : scalarArrays) {
    std::string arrayName = std::string(scalarArray->GetName());
    if (arrayName == "GlobalPointIds") globalPointIdsExist = true;
    if (arrayName == "vtkGhostType") ghostCellsExist = true;
  }

  if (globalPointIdsExist && ghostCellsExist){
    this->printMsg("Global Point Ids and Ghost Cells exist, therefore we are in distributed mode!");
    for(auto scalarArray : scalarArrays) {
      std::string arrayName = std::string(scalarArray->GetName());
      if (arrayName != "GlobalPointIds" && arrayName != "vtkGhostType"){
        auto vtkglobalPointIds = pointData->GetGlobalIds();
        vtkMPIController *controller = vtkMPIController::SafeDownCast(vtkMPIController::GetGlobalController());
        int numProcs = controller->GetNumberOfProcesses();
        int rank = controller->GetLocalProcessId();
        this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));

        // sort the scalar array distributed first by the scalar value itself, then by the global id

        std::vector<std::tuple<double, int, int>> sortingValues;

        for (int i = 0; i < nVertices; i++){
          double scalarValue = scalarArray->GetComponent(i, 0);
          int globalId = vtkglobalPointIds->GetComponent(i, 0);
          int localId = i;
          sortingValues.emplace_back(scalarValue, globalId, localId);
        }
        auto element0 = sortingValues[0];
        this->printMsg("#Elements in Vector " + std::to_string(sortingValues.size()) + ", first values are " + std::to_string(std::get<0>(element0)) + " " + std::to_string(std::get<1>(element0)) + " " + std::to_string(std::get<2>(element0)));

        // sort the vector of this rank, first by their scalarvalue and if they are the same, by their globalId
        std::sort(sortingValues.begin(), sortingValues.end(),
        [&](const std::tuple<double, int, int> a, const std::tuple<double, int, int> b) {
          return (std::get<0>(a) < std::get<0>(b))
                 || (std::get<0>(a) == std::get<0>(a) && std::get<1>(a) < std::get<1>(b));
        });

        this->printMsg("#Elements in Vector " + std::to_string(sortingValues.size()) + ", first values are " + std::to_string(std::get<0>(element0)) + " " + std::to_string(std::get<1>(element0)) + " " + std::to_string(std::get<2>(element0)));

        // when all are done sorting, rank 0 requests the highest values and merges them

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
