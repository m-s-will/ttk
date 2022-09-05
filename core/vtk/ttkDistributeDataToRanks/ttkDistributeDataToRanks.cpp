#include <ttkDistributeDataToRanks.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkIntArray.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <limits>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkDistributeDataToRanks);

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
ttkDistributeDataToRanks::ttkDistributeDataToRanks() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkDistributeDataToRanks::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
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
int ttkDistributeDataToRanks::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
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
int ttkDistributeDataToRanks::RequestData(vtkInformation *ttkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting distribution...");

#ifdef TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI()) {
    vtkMPICommunicator *comm = vtkMPICommunicator::SafeDownCast(
      vtkMultiProcessController::GetGlobalController()->GetCommunicator());
    comm->Broadcast(inputDataSet, 0);

  } else {
    this->printMsg("MPI is not initialized, please run with mpirun!");
  }
#else
  this->printMsg("TTK is not build with TTK_ENABLE_MPI, this filter does "
                 "nothing in sequential mode.");
#endif // TTK_ENABLE_MPI

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  // outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}
