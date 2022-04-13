#include <ttkBoundingBoxNeighborDetector.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>

#include <mpi.h>
#include <limits>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkBoundingBoxNeighborDetector);

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
ttkBoundingBoxNeighborDetector::ttkBoundingBoxNeighborDetector() {
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
int ttkBoundingBoxNeighborDetector::FillInputPortInformation(int port, vtkInformation *info) {
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
int ttkBoundingBoxNeighborDetector::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

// returns true if they intersect, false if not
bool checkForIntersection(double* myBB, double* theirBB){
  return !(myBB[0] > theirBB[1]  // my left side is right of their right side
        || myBB[1] < theirBB[0]  // my right side is left of their left side
        || myBB[2] > theirBB[3]  // my bottom side is above their top side
        || myBB[3] < theirBB[2]  // my top side is under their bottom side
        || myBB[4] > theirBB[5]  // my front side is behind their back side
        || myBB[5] < theirBB[4]  // my back side is in front of their front side
      );
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
int ttkBoundingBoxNeighborDetector::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  //auto pointData = inputDataSet->GetPointData();
  //int nVertices = inputDataSet->GetNumberOfPoints();

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");



  int flag_i;
  MPI_Initialized(&flag_i);
  if (flag_i){
    int numProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) this->printMsg("MPI is initialized, therefore we are in distributed mode!");
    this->printMsg("#Ranks " + std::to_string(numProcs) + ", this is rank " + std::to_string(rank));
    double* boundingBox = inputDataSet->GetBounds();
    std::vector<double *> rankBoundingBoxes(numProcs);
    rankBoundingBoxes[rank] = boundingBox;
    for (int r = 0; r < numProcs; r++){
      if (r != rank) rankBoundingBoxes[r] = (double *)malloc(6 * sizeof(double));
      MPI_Bcast(rankBoundingBoxes[r], 6, MPI_DOUBLE, r, MPI_COMM_WORLD);
    } 

    
    //double epsilon = std::numeric_limits<double>::epsilon();
    double epsilon = 0.0001;
    // inflate our own bounding box by epsilon
    for (int i = 0; i < 6; i++){
      if (i % 2 == 0) boundingBox[i]-=epsilon;
      if (i % 2 == 1) boundingBox[i]+=epsilon;
    }
    std::vector<int> neighbors;
    // the first value of the neighbors-array is the rank from which is originating
    neighbors.push_back(rank);
    for (int i = 0; i < numProcs; i++){
      if (i != rank){
        double* theirBoundingBox = rankBoundingBoxes[i];
        if (checkForIntersection(boundingBox, theirBoundingBox)){
          this->printMsg("Rank " + std::to_string(rank) + " is neighbors with Rank " + std::to_string(i));
          neighbors.push_back(i);
        }
        
      }

      // all ranks send data to rank 0 and that rank adds to fielddata?
      vtkNew<vtkIntArray> neighborsArray{};
      neighborsArray->SetName("Neighbors");
      neighborsArray->SetNumberOfComponents(neighbors.size());
      neighborsArray->InsertNextTuple((float*)neighbors.data());
      inputDataSet->GetFieldData()->AddArray(neighborsArray);
    }
    
  } else {
    this->printMsg("MPI is not initialized, please run with mpirun!");
  }












  // Create an output array that has the same data type as the input array
  // Note: vtkSmartPointers are well documented
  //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  /*vtkSmartPointer<vtkDataArray> outputArray
    = vtkSmartPointer<vtkDataArray>::New(inputArray->NewInstance());
  outputArray->SetName(this->OutputArrayName.data()); // set array name
  outputArray->SetNumberOfComponents(1); // only one component per tuple
  outputArray->SetNumberOfTuples(nVertices);
  */
  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  /*
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
 // this->preconditionTriangulation(triangulation); // implemented in base class
  
  
  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
                      (status = this->computeAverages<VTK_TT, TTK_TT>(
                         (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),
                         (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
                         (TTK_TT *)triangulation->getData())));

  // On error cancel filter execution
  if(status != 1)
    return 0;
  */
  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  //outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}
