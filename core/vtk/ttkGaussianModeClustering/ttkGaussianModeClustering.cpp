#include <ttkGaussianModeClustering.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>

#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkGaussianModeClustering);

ttkGaussianModeClustering::ttkGaussianModeClustering() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

ttkGaussianModeClustering::~ttkGaussianModeClustering() {
}

int ttkGaussianModeClustering::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkGaussianModeClustering::ComputeSimilarityMatrix(
  vtkImageData *similarityMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {
  // unpack input
  auto p0 = vtkPointSet::SafeDownCast(inputDataObjects0);
  auto p1 = vtkPointSet::SafeDownCast(inputDataObjects1);
  if(!p0 || !p1)
    return !this->printErr("Input data objects need to be vtkPointSets.");

  const int nPoints0 = p0->GetNumberOfPoints();
  const int nPoints1 = p1->GetNumberOfPoints();

  // Get ids
  auto ids0 = GetInputArrayToProcess(0, p0);
  auto ids1 = GetInputArrayToProcess(0, p1);
  if(!ids0 || !ids1)
    return !this->printErr("Input data is missing array containing ids.");

  // Get point amplitudes
  auto amps0 = GetInputArrayToProcess(2, p0);
  auto amps1 = GetInputArrayToProcess(2, p1);
  if(!amps0 || !amps1)
    return !this->printErr(
      "Input data is missing array containing point amplitudes.");

  // Get point variances
  auto vars0 = GetInputArrayToProcess(3, p0);
  auto vars1 = GetInputArrayToProcess(3, p1);
  if(!vars0 || !vars1)
    return !this->printErr(
      "Input data is missing array containing point variances.");

  // Get point coords
  auto coords0 = p0->GetPoints()->GetData();
  auto coords1 = p1->GetPoints()->GetData();

  if(coords0->GetDataType() != coords1->GetDataType())
    return !this->printErr("Input vtkPointSet need to have same precision.");

  // Check if both data objects umbrella clusters should be calculated
  if(clustersPerTimestep_.size() == 0) {
    std::map<int, std::vector<int>> clusters0;
    std::map<int, std::vector<int>> clusters1;

    int status = 0;

    ttkTypeMacroR(
      coords0->GetDataType(),
      (status = this->computeClusters<T0>(
         clusters0, nClusters0_, ttkUtils::GetPointer<const T0>(coords0),
         ttkUtils::GetPointer<const double>(amps0),
         ttkUtils::GetPointer<const double>(vars0), nPoints0)));
    if(!status)
      return 0;

    clustersPerTimestep_.push_back(clusters0);
  }
  {
    std::map<int, std::vector<int>> clusters1;

    int status = 0;

    ttkTypeMacroR(
      coords1->GetDataType(),
      (status = this->computeClusters<T0>(
         clusters1, nClusters1_, ttkUtils::GetPointer<const T0>(coords1),
         ttkUtils::GetPointer<const double>(amps1),
         ttkUtils::GetPointer<const double>(vars1), nPoints1)));

    if(!status)
      return 0;

    clustersPerTimestep_.push_back(clusters1);
  }

  auto &clusters0 = clustersPerTimestep_[clustersPerTimestep_.size() - 2];
  auto &clusters1 = clustersPerTimestep_[clustersPerTimestep_.size() - 1];

  // initialize similarity matrix
  similarityMatrix->SetDimensions(nClusters0_, nClusters1_, 1);
  similarityMatrix->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  auto matrixData = similarityMatrix->GetPointData()->GetArray(0);
  matrixData->SetName("ClusterCorrespondence");

  // compute matrix
  int status = 0;
  ttkTypeMacroI(
    ids0->GetDataType(),
    status = this->computeMatrix<T0>(
      ttkUtils::GetPointer<unsigned char>(matrixData), clusters0, clusters1,
      ttkUtils::GetPointer<const T0>(ids0),
      ttkUtils::GetPointer<const T0>(ids1), nClusters0_, nClusters1_));

  if(!status)
    return 0;

  // create id vector for umbrellas connected to the values in the matrix
  auto clusterIds0 = vtkSmartPointer<vtkDataArray>::Take(ids0->NewInstance());
  auto clusterIds1 = vtkSmartPointer<vtkDataArray>::Take(ids1->NewInstance());
  clusterIds0->SetName(ids0->GetName());
  clusterIds1->SetName(ids1->GetName());
  clusterIds0->SetNumberOfComponents(1);
  clusterIds1->SetNumberOfComponents(1);
  clusterIds0->SetNumberOfTuples(nClusters0_);
  clusterIds1->SetNumberOfTuples(nClusters1_);

  int index = 0;
  for(const auto &i : clusters0) {
    if(i.second[0] != -1) {
      clusterIds0->SetTuple1(index, ids0->GetTuple1(i.first));
      index++;
    }
  }

  index = 0;
  for(const auto &i : clusters1) {
    if(i.second[0] != -1) {
      clusterIds1->SetTuple1(index, ids1->GetTuple1(i.first));
      index++;
    }
  }

  // Create arrays for time data
  auto timeArray = vtkSmartPointer<vtkIntArray>::New();
  timeArray->SetName("Timesteps");
  timeArray->SetNumberOfComponents(1);
  timeArray->SetNumberOfTuples(2);
  timeArray->SetTuple1(0, clustersPerTimestep_.size() - 2);
  timeArray->SetTuple1(1, clustersPerTimestep_.size() - 1);
  similarityMatrix->GetFieldData()->AddArray(timeArray);

  status = ttkSimilarityAlgorithm::AddIndexIdMaps(
    similarityMatrix, clusterIds0, clusterIds1);
  if(!status)
    return 0;

  // Propagate
  nClusters0_ = nClusters1_;
  nClusters1_ = 0;

  return 1;
}

int ttkGaussianModeClustering::AddClusterIds(vtkDataObject *inputDataObjects,
                                             const size_t t) {
  auto p0 = vtkPointSet::SafeDownCast(inputDataObjects);
  if(!p0)
    return !this->printErr("No points");

  // Get number of points
  const int nPoints = p0->GetNumberOfPoints();

  // Do nothing if there are no points
  if(nPoints == 0)
    return 1;

  // Get ids
  auto ids = GetInputArrayToProcess(0, p0);
  if(!ids)
    return !this->printErr("Input data is missing array containing ids.");

  // Add umbrella and tresholded umbrella ids to point output
  auto clusterIds = vtkSmartPointer<vtkIntArray>::New();
  clusterIds->SetName("ClusterId");
  clusterIds->SetNumberOfComponents(1);
  clusterIds->SetNumberOfTuples(nPoints);

  for(const auto &u : clustersPerTimestep_[t]) {
    for(long unsigned int j = 0; j < u.second.size(); j++) {
      if(u.second[j] == -1)
        clusterIds->SetTuple1(u.first, -1);
      else
        clusterIds->SetTuple1(u.second[j], ids->GetTuple1(u.first));
    }
  }

  p0->GetPointData()->AddArray(clusterIds);

  return 1;
}

int ttkGaussianModeClustering::FormatClusters(vtkDataObject *inputDataObjects,
                                              vtkPolyData *outputPoints,
                                              const size_t t) {
  auto inPointSet = vtkPointSet::SafeDownCast(inputDataObjects);
  if(!inPointSet)
    return !this->printErr("No points");

  auto inPD = inPointSet->GetPointData();
  auto inCD = inPointSet->GetCellData();

  int nPoints = inPointSet->GetNumberOfPoints();
  // Do nothing if there are no points
  if(nPoints == 0)
    return 1;

  // Get ids
  auto ids = GetInputArrayToProcess(0, inPointSet);
  if(!ids)
    return !this->printErr("Input data is missing array containing ids.");

  // Create new points and allocate point data
  auto newPoints = vtkSmartPointer<vtkPoints>::New();
  auto outPD = outputPoints->GetPointData();
  auto outCD = outputPoints->GetCellData();
  outPD->CopyAllocate(inPD);
  outCD->CopyAllocate(inCD);
  outputPoints->Allocate(inPointSet->GetNumberOfCells());
  newPoints->SetDataType(inPointSet->GetPoints()->GetDataType());

  // Copy all data from points to the new clusters
  int nNewPoints = 0;
  for(const auto &cluster : clustersPerTimestep_[t]) {
    if(cluster.second[0] != -1) {
      nNewPoints++;
    }
  }

  newPoints->SetNumberOfPoints(nNewPoints);
  auto newPointsData = newPoints->GetData();
  auto inPointsData = inPointSet->GetPoints()->GetData();
  int p = 0;
  auto vertId = vtkSmartPointer<vtkIdList>::New();
  for(const auto &cluster : clustersPerTimestep_[t]) {
    if(cluster.second[0] != -1) {
      newPointsData->SetTuple(p, cluster.first, inPointsData);
      vertId->InsertId(0, p);
      int c = outputPoints->InsertNextCell(VTK_VERTEX, vertId);
      outPD->CopyData(inPD, cluster.first, p);
      outCD->CopyData(inCD, cluster.first, c);
      p++;
    }
  }

  outputPoints->SetPoints(newPoints);
  outputPoints->Squeeze();
  outputPoints->GetFieldData()->DeepCopy(inPointSet->GetFieldData());

  return 1;
}

int ttkGaussianModeClustering::RequestData(vtkInformation *request,
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector) {

  // Clear previous data
  clustersPerTimestep_.clear();
  nClusters0_ = 0;
  nClusters1_ = 0;

  // Set class members in base class
  this->setClusteringType(this->ClusteringType);
  this->setThreshold(this->ScalarThreshold);

  // Calculate matrices
  this->ttkSimilarityAlgorithm::RequestData(request, inputVector, outputVector);

  // Copy input
  auto input = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto pointsOut = vtkMultiBlockDataSet::GetData(outputVector, 1);
  pointsOut->DeepCopy(input);

  /* FORMATTING OUTPUT */
  ttk::Timer timer;
  const std::string msg = "Add Feature Arrays to Points";
  this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);
  for(size_t t = 0; t < pointsOut->GetNumberOfBlocks(); t++) {
    if(!AddClusterIds(pointsOut->GetBlock(t), t))
      return 0;
  }
  this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

  const std::string msg2 = "Create Cluster Output";
  this->printMsg(
    msg2, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // Create and format all clusters via function
  auto clustersOut = vtkMultiBlockDataSet::GetData(outputVector, 2);
  for(size_t t = 0; t < pointsOut->GetNumberOfBlocks(); t++) {
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    if(!FormatClusters(pointsOut->GetBlock(t), polyData, t))
      return 0;
    clustersOut->SetBlock(t, polyData);
  }
  this->printMsg(msg2, 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}
