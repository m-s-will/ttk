#include <ttkBranchDecomposition.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <TrackingGraph.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkBranchDecomposition);

ttkBranchDecomposition::ttkBranchDecomposition() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkBranchDecomposition::~ttkBranchDecomposition() {
}

int ttkBranchDecomposition::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkBranchDecomposition::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkBranchDecomposition::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto input = vtkPointSet::GetData(inputVector[0]);
  if(!input)
    return 0;

  auto vtkTrackingGraph = vtkPointSet::GetData(outputVector);
  vtkTrackingGraph->ShallowCopy(input);

  if(this->GetInputArrayAssociation(0, inputVector) != 0)
    return !this->printErr("Time needs to be a point data array.");
  auto timeArray = this->GetInputArrayToProcess(0, inputVector);
  if(!timeArray)
    return !this->printErr("Unable to retrieve time array.");

  auto attributeArrayAssociation
    = this->GetInputArrayAssociation(1, inputVector);
  if(attributeArrayAssociation != 0 && attributeArrayAssociation != 1)
    return !this->printErr("Attribute needs to be a point or cell data array.");
  auto attributeArray = this->GetInputArrayToProcess(1, inputVector);
  if(!attributeArray)
    return !this->printErr("Unable to retrieve attribute array.");

  if(timeArray->GetNumberOfComponents() != 1
     || attributeArray->GetNumberOfComponents() != 1)
    return !this->printErr("Input arrays need to be scalar arrays.");

  const int nNodes = vtkTrackingGraph->GetNumberOfPoints();
  auto connectivityList
    = vtkTrackingGraph->IsA("vtkPolyData")
        ? static_cast<vtkPolyData *>(vtkTrackingGraph)
            ->GetLines()
            ->GetConnectivityArray()
      : vtkTrackingGraph->IsA("vtkUnstructuredGrid")
        ? static_cast<vtkUnstructuredGrid *>(vtkTrackingGraph)
            ->GetCells()
            ->GetConnectivityArray()
        : nullptr;

  if(!connectivityList)
    return !this->printErr("Unable to retrieve connectivity list.");

  const int nEdges = connectivityList->GetNumberOfValues() / 2;

  ttk::TrackingGraph ttkTrackingGraph;
  ttkTrackingGraph.setDebugLevel(this->debugLevel_);
  ttkTypeMacroI(
    connectivityList->GetDataType(),
    ttkTrackingGraph.preconditionInOutEdges<T0>(
      nNodes, nEdges, ttkUtils::GetConstPointer<const T0>(connectivityList)));

  auto branchIdP = vtkSmartPointer<vtkIntArray>::New();
  branchIdP->SetName("BranchId");
  branchIdP->SetNumberOfComponents(1);
  branchIdP->SetNumberOfTuples(nNodes);
  vtkTrackingGraph->GetPointData()->AddArray(branchIdP);

  auto branchIdC = vtkSmartPointer<vtkIntArray>::New();
  branchIdC->SetName("BranchId");
  branchIdC->SetNumberOfComponents(1);
  branchIdC->SetNumberOfTuples(nEdges);
  vtkTrackingGraph->GetCellData()->AddArray(branchIdC);

  int status = 0;
  ttkTypeMacroAA(timeArray->GetDataType(), attributeArray->GetDataType(),
                 (status = this->computeBranchDecompositionByAttribute<T0, T1>(
                    ttkUtils::GetPointer<int>(branchIdP),
                    ttkUtils::GetPointer<int>(branchIdC), ttkTrackingGraph,
                    ttkUtils::GetConstPointer<const T0>(timeArray),
                    ttkUtils::GetConstPointer<const T1>(attributeArray),
                    attributeArrayAssociation)));
  return status;
}