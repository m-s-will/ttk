#include <ttkInputPointAdvection.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// std includes
#include <random>

vtkStandardNewMacro(ttkInputPointAdvection);

ttkInputPointAdvection::ttkInputPointAdvection() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkInputPointAdvection::~ttkInputPointAdvection() {
}

int ttkInputPointAdvection::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;

  return 1;
}

int ttkInputPointAdvection::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkInputPointAdvection::rampFunction(const int t,
                                         const int lifetime,
                                         double &y) {
  // variables for ramp function, deciding slope, cutoff for constant value
  double cut = std::min(1.0 / 3, 10.0 / lifetime);
  double offset = 0.001;
  double invCut = (1 - offset) / cut;

  double x = double(t) / lifetime;

  if(x < cut || x == cut) {
    y = offset + invCut * x;
  } else if(x > cut && x < (1 - cut)) {
    y = 1.0;
  } else {
    y = offset + invCut - invCut * x;
  }

  return 1;
}

int ttkInputPointAdvection::formatOutput(vtkPolyData *pd,
                                         std::vector<int> &aliveIds,
                                         const int timestep) {
  // Function for formatting data arrays
  auto prepArray
    = [](vtkDataArray *array, std::string name, int nTuples, int nComponents) {
        array->SetName(name.data());
        array->SetNumberOfComponents(nComponents);
        array->SetNumberOfTuples(nTuples);
        return ttkUtils::GetVoidPointer(array);
      };

  auto dataPoints = vtkSmartPointer<vtkPoints>::New();

  // Set point initial data
  int nPs = aliveIds.size();
  dataPoints->SetDataType(VTK_DOUBLE);
  dataPoints->SetNumberOfPoints(nPs);

  // Create cell arrays off offset array and connectivity array.
  // Our cells are just vertices
  auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
  offsetArray->SetNumberOfTuples(nPs + 1);
  auto offsetArrayData
    = static_cast<int *>(ttkUtils::GetVoidPointer(offsetArray));

  auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
  connectivityArray->SetNumberOfTuples(nPs);
  auto connectivityArrayData
    = static_cast<int *>(ttkUtils::GetVoidPointer(connectivityArray));

  auto cellArray = vtkSmartPointer<vtkCellArray>::New();
  cellArray->SetData(offsetArray, connectivityArray);

  // Create arrays for point data
  auto idArray = vtkSmartPointer<vtkIntArray>::New();
  auto idArrayData = static_cast<int *>(prepArray(idArray, "PointId", nPs, 1));

  auto birthArray = vtkSmartPointer<vtkIntArray>::New();
  auto birthArrayData
    = static_cast<int *>(prepArray(birthArray, "Birth", nPs, 1));
  auto deathArray = vtkSmartPointer<vtkIntArray>::New();
  auto deathArrayData
    = static_cast<int *>(prepArray(deathArray, "Death", nPs, 1));

  auto amplitudeArray = vtkSmartPointer<vtkDoubleArray>::New();
  auto amplitudeArrayData
    = static_cast<double *>(prepArray(amplitudeArray, "Amplitude", nPs, 1));
  auto varianceArray = vtkSmartPointer<vtkDoubleArray>::New();
  auto varianceArrayData
    = static_cast<double *>(prepArray(varianceArray, "Variance", nPs, 1));

  // Create arrays for field data
  auto timeArray = vtkSmartPointer<vtkDoubleArray>::New();
  auto timeArrayData
    = static_cast<double *>(prepArray(timeArray, "Time", 1, 1));
  timeArrayData[0] = timestep * this->TimeInterval;

  // Add data to points and arrays
  int idx = 0;
  for(int j = 0; j < nPs; j++) {
    auto p = this->allPoints[aliveIds[j]];

    // Set position
    double pos[3] = {p.x, p.y, p.z};
    dataPoints->SetPoint(idx, pos);

    // Calculate the point amplitude in current timestep
    double pw = 0.0;
    rampFunction(timestep - p.birth, p.death - p.birth, pw);

    // Set data arrays
    idArrayData[idx] = p.pointId;
    birthArrayData[idx] = p.birth;
    deathArrayData[idx] = p.death;
    amplitudeArrayData[idx] = pw * p.amplitude;
    varianceArrayData[idx] = p.variance;

    // Set connectivity and offset array
    connectivityArrayData[idx] = idx;
    offsetArrayData[idx] = idx;
    idx++;
  }

  // Set offset array last index to the number of elements in the connectivity
  // array
  offsetArrayData[nPs] = nPs;

  // Format the output data structure into a dataset of vtkPolyData
  pd->SetPoints(dataPoints);
  pd->SetVerts(cellArray);

  auto pointData = pd->GetPointData();
  pointData->AddArray(idArray);
  pointData->AddArray(birthArray);
  pointData->AddArray(deathArray);

  pointData->AddArray(amplitudeArray);
  pointData->AddArray(varianceArray);

  pd->GetFieldData()->AddArray(timeArray);

  return 1;
}

int ttkInputPointAdvection::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  // Get input
  auto input = vtkPointSet::GetData(inputVector[0], 0);
  if(!input) {
    this->printErr("There is no input provided.");
    return 0;
  }

  const size_t nPoints = input->GetNumberOfPoints();
  auto inputPoints = input->GetPoints();

  // Create output
  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  /* Create sampling distributions */

  // Get random number engine and seed it with user provided seed
  std::mt19937 randGen(RandomSeed);

  // Position distributions
  std::uniform_real_distribution<> disX(0.0, 1.0);
  std::uniform_real_distribution<> disY(0.0, 1.0);
  std::uniform_real_distribution<> disZ(0.0, 1.0);

  // Time data distributions
  std::uniform_int_distribution<> disLifetime(
    this->Lifetime[0], this->Lifetime[1]);
  std::uniform_int_distribution<> disRespawnTime(
    this->RespawnTime[0], this->RespawnTime[1]);
  std::uniform_real_distribution<> disRate(1, 1);

  // Attributes distributions
  std::uniform_real_distribution<> disAmplitude(
    this->Amplitude[0], this->Amplitude[1]);
  std::uniform_real_distribution<> disVariance(
    this->Variance[0], this->Variance[1]);

  /* Create initial points */

  // Delete old data and create space for new
  this->allPoints.clear();
  this->allPoints.resize(nPoints);

  // Create all new points
  for(size_t i = 0; i < nPoints; i++) {
    auto &p = this->allPoints[i];
    double pos[3] = {0, 0, 0};
    inputPoints->GetPoint(i, pos);

    // Positions
    p.x = pos[0];
    p.y = pos[1];
    p.z = pos[2];

    // Point data
    p.pointId = i;

    p.timestep = 0;
    p.birth = 0 + disRespawnTime(randGen);
    p.death = p.birth + disLifetime(randGen);
    p.rate = disRate(randGen);

    p.amplitude = disAmplitude(randGen);
    p.variance = disVariance(randGen);
  }

  // maximum point id
  int maxPointId = nPoints - 1;

  /* set base layer variables */

  // Determine what vector field to use
  InputPointAdvection::VectorField vf
    = InputPointAdvection::VectorField::PerlinPerturbed;

  switch(this->VecField) {
    case 0: {
      vf = InputPointAdvection::VectorField::PerlinPerturbed;
      break;
    }
    case 1: {
      vf = InputPointAdvection::VectorField::PerlinGradient;
      break;
    }
    case 2: {
      vf = InputPointAdvection::VectorField::PosDiagonal;
      break;
    }
    case 3: {
      vf = InputPointAdvection::VectorField::PosX;
      break;
    }
    case 4: {
      vf = InputPointAdvection::VectorField::Sink;
      break;
    }
    case 5: {
      vf = InputPointAdvection::VectorField::Saddle;
      break;
    }
  }

  this->setVariables(this->StepLength, this->PerlinScaleFactor, vf);

  /* Loop through all timesteps and advect points */
  std::vector<int> alivePointsIds;
  ttk::Timer timer;
  const std::string msg
    = "Advecting " + std::to_string(nPoints) + " points in vector field for "
      + std::to_string(this->NumberOfTimesteps) + " timesteps";
  this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);
  for(int t = 0; t < this->NumberOfTimesteps; t++) {
    alivePointsIds.clear();

    // Go through all points
    for(size_t i = 0; i < nPoints; i++) {
      auto &p = this->allPoints[i];

      // Check what points are alive
      if(t >= p.birth && t <= p.death) {
        alivePointsIds.push_back(i);
      } else if(t > p.death) {
        // re-spawn
        ++maxPointId;

        // Positions
        p.x = disX(randGen);
        p.y = disY(randGen);
        p.z = disZ(randGen);

        // Point data
        p.pointId = maxPointId;

        p.timestep = t;
        p.birth = t + disRespawnTime(randGen);
        p.death = p.birth + disLifetime(randGen);
        p.rate = disRate(randGen);

        p.amplitude = disAmplitude(randGen);
        p.variance = disVariance(randGen);

        // check case re-spawn time is zero, so the point is immediately added
        if(p.birth == t) {
          alivePointsIds.push_back(i);
        }
      } else if(t < p.birth) {
        // do nothing
        p.timestep = t;
      }
    }

    // Format output
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    this->formatOutput(polyData, alivePointsIds, t);

    // Advect points
    this->advect(this->allPoints, alivePointsIds, t, this->TimeInterval);

    // Set data to a block in the output dataset
    size_t nBlocks = outputMB->GetNumberOfBlocks();
    outputMB->SetBlock(nBlocks, polyData);
  }

  this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

  // return success
  return 1;
}
