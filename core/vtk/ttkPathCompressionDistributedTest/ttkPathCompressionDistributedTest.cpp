#include <ttkPathCompressionDistributedTest.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkIntArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPathCompressionDistributedTest);

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
ttkPathCompressionDistributedTest::ttkPathCompressionDistributedTest() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPathCompressionDistributedTest::~ttkPathCompressionDistributedTest() {
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkPathCompressionDistributedTest::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
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
int ttkPathCompressionDistributedTest::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }

  return 0;
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
int ttkPathCompressionDistributedTest::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  this->printMsg("Request data called!");
  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;



  /*

  dataSetWithGlobalIds->GetPointData()->AddArray(ttkBoundaryVertices);

  // create a ghost cell layer from VTK ghostcells (level 1 and level 2 ghost cells), this extends the dataset
  auto generator = vtkSmartPointer<vtkGhostCellsGenerator>::New();
  generator->SetNumberOfGhostLayers(1);
  generator->BuildIfRequiredOff();
  generator->SetInputData(dataSetWithGlobalIds);
  generator->Update();

  vtkDataSet *vtkGhostLayer = vtkDataSet::SafeDownCast(generator->GetOutput());
  vtkPointData *vtkGhostPoints = vtkGhostLayer->GetPointData();
  vtkDataArray *ttkboundaryPointValues = vtkGhostPoints->GetScalars("boundaryVertices");


  

  */

  
  // Get input array that will be processed
  //
  // Note: VTK provides abstract functionality to handle array selections, but
  //       this essential functionality is unfortunately not well documented.
  //       Before you read further, please keep in mind the the TTK developer
  //       team is not responsible for the existing VTK Api ;-)
  //
  //       In a nutshell, prior to the RequestData execution one has to call
  //
  //           SetInputArrayToProcess (
  //               int idx,
  //               int port,
  //               int connection,
  //               int fieldAssociation,
  //               const char *name
  //            )
  //
  //       The parameter 'idx' is often missunderstood: lets say the filter
  //       requires n arrays, then idx enumerates them from 0 to n-1.
  //
  //       The 'port' is the input port index at which the object is connected
  //       from which we want to get the array.
  //
  //       The 'connection' is the connection index at that port (we have to
  //       specify this because VTK allows multiple connections at the same
  //       input port).
  //
  //       The 'fieldAssociation' integer specifies if the array should be taken
  //       from 0: point data, 1: cell data, or 2: field data.
  //
  //       The final parameter is the 'name' of the array.
  //
  //       Example: SetInputArrayToProcess(3,1,0,1,"EdgeLength") will store that
  //                for the 3rd array the filter needs the cell data array named
  //                "EdgeLength" that it will retrieve from the vtkDataObject
  //                at input port 1 (first connection). During the RequestData
  //                method one can then actually retrieve the 3rd array it
  //                requires for its computation by calling
  //                GetInputArrayToProcess(3, inputVector)
  //
  //       If this filter is run within ParaView, then the UI will automatically
  //       call SetInputArrayToProcess (see PathCompressionDistributedTest.xml file).
  //
  //       During the RequestData execution one can then retrieve an actual
  //       array with the method "GetInputArrayToProcess".
  //vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);

  auto order = ttkAlgorithm::GetOrderArray(
    inputDataSet, 0);
  if(!order) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  
  if(order->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation!");
  this->printMsg("  Scalar Array: " + std::string(order->GetName()));


  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation){   
    return 0;
  } else {
    this->printMsg("Triangulation type:" + std::to_string(static_cast<int>(triangulation->getType())));
  }

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class
  this->printMsg("Triangulation type:" + std::to_string(static_cast<int>(triangulation->getType())));


  vtkSmartPointer<vtkIntArray> descendingManifold
    = vtkSmartPointer<vtkIntArray>::New();
  descendingManifold->SetName("DescendingManifold"); // set array name
  descendingManifold->SetNumberOfComponents(1); // only one component per tuple
  descendingManifold->SetNumberOfTuples(order->GetNumberOfTuples());

  vtkSmartPointer<vtkIntArray> ascendingManifold
    = vtkSmartPointer<vtkIntArray>::New();
  ascendingManifold->SetName("AscendingManifold"); // set array name
  ascendingManifold->SetNumberOfComponents(1); // only one component per tuple
  ascendingManifold->SetNumberOfTuples(order->GetNumberOfTuples());

  this->printMsg("  Output Array 1: " + std::string(descendingManifold->GetName()));
  this->printMsg("  Output Array 2: " + std::string(ascendingManifold->GetName()));

  auto pointData = inputDataSet->GetPointData();
  auto rankArray = pointData->GetArray("rankArray");
  auto globalIds = pointData->GetGlobalIds();

  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  ttkTypeMacroIT(order->GetDataType(), triangulation->getType(),
                      (status = this->computeCompression<T0, T1>(
                         ttkUtils::GetPointer<int>(descendingManifold),
                         ttkUtils::GetPointer<int>(ascendingManifold),
                         ttkUtils::GetPointer<T0>(order),
                         ttkUtils::GetPointer<int>(rankArray),
                         ttkUtils::GetPointer<T0>(globalIds),
                         (T1 *)triangulation->getData())));



  // On error cancel filter execution
  if(status != 1)
    return 0;
  

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(descendingManifold);
  outputDataSet->GetPointData()->AddArray(ascendingManifold);


  
  // return success
  return 1;
}
