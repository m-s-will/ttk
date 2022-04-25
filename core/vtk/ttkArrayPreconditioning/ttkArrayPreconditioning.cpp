#include <OrderDisambiguation.h>
#include <ttkArrayPreconditioning.h>
#include <ttkMacros.h>
#include <ttkUtils.h>
#include <mpi.h>

#include <regex>
#include <unordered_map>
#include <vtkCommand.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

vtkStandardNewMacro(ttkArrayPreconditioning);

ttkArrayPreconditioning::ttkArrayPreconditioning() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

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
  size_t nVertices = input->GetNumberOfPoints();

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
  auto rankArray = pointData->GetArray("RankArray");
  if(vtkGlobalPointIds != nullptr && vtkGhostCells != nullptr
     && rankArray != nullptr) {
#ifdef TTK_ENABLE_MPI
    // add the order array for every scalar array, except the ghostcells, the
    // rankarray and the global ids
    for(auto scalarArray : scalarArrays) {
      int status;
      std::string arrayName = std::string(scalarArray->GetName());
      if(arrayName != "GlobalPointIds" && arrayName != "vtkGhostType"
         && arrayName != "RankArray") {
        this->printMsg("Arrayname: " + arrayName);
        vtkNew<ttkSimplexIdTypeArray> orderArray{};
        orderArray->SetName(this->GetOrderArrayName(scalarArray).data());
        orderArray->SetNumberOfComponents(1);
        orderArray->SetNumberOfTuples(nVertices);

        ttkTypeMacroAI(
          scalarArray->GetDataType(), vtkGlobalPointIds->GetDataType(),
          (status = processScalarArray<T0, T1>(
             ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
             ttkUtils::GetPointer<T0>(scalarArray),
             ttkUtils::GetPointer<T1>(vtkGlobalPointIds),
             ttkUtils::GetPointer<int>(rankArray),
             ttkUtils::GetPointer<char>(vtkGhostCells), nVertices, BurstSize)));

        // On error cancel filter execution
        if(status != 1)
          return 0;
        output->GetPointData()->AddArray(orderArray);
      }
    }
    this->printMsg("Preconditioned selected scalar arrays", 1.0,
                   tm.getElapsedTime(), this->threadNumber_);
    return 1;
#else
    this->printMsg(
      "Necessary arrays are present, but TTK is not build with MPI support");
    return 0;
#endif
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
