#include <ttkPairExtrema.h>

#include <vtkInformation.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPairExtrema);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkPairExtrema::ttkPairExtrema() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkPairExtrema::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkPairExtrema::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <class triangulationType>
int ttkPairExtrema::getSkeletonArcs(
  vtkUnstructuredGrid *outputSkeletonArcs,
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
  const ttk::SimplexId *order,
  const triangulationType *triangulation) {
  vtkNew<vtkUnstructuredGrid> skeletonArcs{};
  ttk::SimplexId pointIds[2];
  ttk::SimplexId pointOrders[2];
  vtkNew<vtkPoints> points{};
  vtkNew<vtkLongLongArray> data{};
  data->SetNumberOfComponents(1);
  data->SetName("Order");
  vtkNew<vtkLongLongArray> gIdArray{};
  gIdArray->SetNumberOfComponents(1);
  gIdArray->SetName("GlobalPointIds");
  float point[3];
  std::map<ttk::SimplexId, ttk::SimplexId> addedPoints;
  ttk::SimplexId currentId = 0;
  for(auto const &p : persistencePairs) {
    pointIds[0] = p.first;
    pointIds[1] = p.second;
    pointOrders[0] = order[p.first];
    pointOrders[1] = order[p.second];
    // add each point only once to the vtkPoints
    // addedPoints.insert(x).second inserts x and is true if x was not in
    // addedPoints beforehand
    if(addedPoints.insert({pointIds[0], currentId}).second) {
      // this->printMsg("point " + std::to_string(pointIds[0]));
      triangulation->getVertexPoint(pointIds[0], point[0], point[1], point[2]);
      points->InsertNextPoint(point);
      data->InsertNextTuple1(pointOrders[0]);
      gIdArray->InsertNextTuple1(pointIds[0]);
      currentId++;
    }
    if(addedPoints.insert({pointIds[1], currentId}).second) {
      // this->printMsg("point " + std::to_string(pointIds[1]));
      triangulation->getVertexPoint(pointIds[1], point[0], point[1], point[2]);
      points->InsertNextPoint(point);
      data->InsertNextTuple1(pointOrders[1]);
      gIdArray->InsertNextTuple1(pointIds[1]);
      currentId++;
    }
    // this->printMsg("Join Tree Arc: " + std::to_string(pointIds[0]) + " "
    //                + std::to_string(pointIds[1]));
    pointIds[0] = addedPoints.at(pointIds[0]);
    pointIds[1] = addedPoints.at(pointIds[1]);
    skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
  }
  skeletonArcs->SetPoints(points);
  outputSkeletonArcs->ShallowCopy(skeletonArcs);
  outputSkeletonArcs->GetPointData()->AddArray(data);
  outputSkeletonArcs->GetPointData()->AddArray(gIdArray);

  return 1;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkPairExtrema::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  vtkPointSet *inputCriticalPoints = vtkPointSet::GetData(inputVector[1]);
  if(!inputDataSet || !inputCriticalPoints)
    return 0;

  auto manifoldPointData = inputDataSet->GetPointData();
  auto ascendingManifold = manifoldPointData->GetArray("AscendingManifold");
  auto tempArray = manifoldPointData->GetArray("DescendingManifold");
  auto criticalPointData = inputCriticalPoints->GetPointData();
  auto criticalGlobalIds = criticalPointData->GetArray("GlobalPointIds");
  auto criticalType = criticalPointData->GetArray("CriticalType");
  auto order = ttkAlgorithm::GetOrderArray(inputDataSet, 0);
  if(!ascendingManifold | !criticalType | !order | !tempArray
     | !criticalGlobalIds) {
    this->printErr("Unable to retrieve input arrays.");
    return 0;
  }

  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.

  if((ascendingManifold->GetNumberOfComponents() != 1) |
     (criticalType->GetNumberOfComponents() != 1) |
     (order->GetNumberOfComponents() != 1) ) {
    this->printErr("Input arrays needs to be a scalar arrays.");
    return 0;
  }

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");

  ttk::SimplexId nCriticalPoints = criticalType->GetNumberOfTuples();
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> persistencePairs{};

  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class

  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  // construct the tree
  ttkTypeMacroT(triangulation->getType(),
                (status = this->computePairs<T0>(
                   persistencePairs, ttkUtils::GetPointer<char>(criticalType),
                   ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
                   ttkUtils::GetPointer<ttk::SimplexId>(tempArray),
                   ttkUtils::GetPointer<ttk::SimplexId>(order),
                   (T0 *)triangulation->getData(),
                   ttkUtils::GetPointer<ttk::SimplexId>(criticalGlobalIds),
                   nCriticalPoints)));

  // On error cancel filter execution
  if(status != 1)
    return 0;

  //  Construct output, see ttkFTMTree.cpp
  auto outputSkeletonArcs = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 1);

  ttkTypeMacroT(
    triangulation->getType(),
    status = getSkeletonArcs<T0>(outputSkeletonArcs, persistencePairs,
                                 ttkUtils::GetPointer<ttk::SimplexId>(order),
                                 (T0 *)triangulation->getData()));

  outputSegmentation->ShallowCopy(inputDataSet);


  // return success
  return 1;
}
