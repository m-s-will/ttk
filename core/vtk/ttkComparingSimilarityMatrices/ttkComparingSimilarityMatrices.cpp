#include <ttkComparingSimilarityMatrices.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTable.h>

#include <ttkMacros.h>
#include <ttkSimilarityAlgorithm.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkComparingSimilarityMatrices);

ttkComparingSimilarityMatrices::ttkComparingSimilarityMatrices() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkComparingSimilarityMatrices::~ttkComparingSimilarityMatrices() {
}

int ttkComparingSimilarityMatrices::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkComparingSimilarityMatrices::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkComparingSimilarityMatrices::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Get input objects from input vector
  auto alg = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto truth = vtkMultiBlockDataSet::GetData(inputVector[1]);

  // Check algorithm and truth has the same number of blocks
  if(alg->GetNumberOfBlocks() != truth->GetNumberOfBlocks()) {
    return !this->printErr("The number of matrices is different.");
  }

  // Clear old data and set number of input
  this->matrixEvents.clear();
  int nMatrices = alg->GetNumberOfBlocks();
  this->matrixEvents.resize(nMatrices);

  for(int i = 0; i < nMatrices; i++) {
    auto aMatrix = alg->GetBlock(i);
    auto algVals = GetInputArrayToProcess(0, aMatrix);
    if(!algVals)
      return !this->printErr(
        "Algorithm matrices is missing array containing matrix data.");

    auto gtMatrix = truth->GetBlock(i);
    auto gtVals = GetInputArrayToProcess(1, gtMatrix);
    if(!gtVals)
      return !this->printErr(
        "Ground truth matrices is missing array containing matrix data.");

    vtkDataArray *indexIdMapAlgP{nullptr};
    vtkDataArray *indexIdMapAlgC{nullptr};
    if(!ttkSimilarityAlgorithm::GetIndexIdMaps(
         indexIdMapAlgP, indexIdMapAlgC, aMatrix->GetFieldData()))
      return !this->printErr(
        "Unable to retrieve Index-Id-Maps for algorithm data.");
    vtkDataArray *indexIdMapGTP{nullptr};
    vtkDataArray *indexIdMapGTC{nullptr};
    if(!ttkSimilarityAlgorithm::GetIndexIdMaps(
         indexIdMapGTP, indexIdMapGTC, gtMatrix->GetFieldData()))
      return !this->printErr(
        "Unable to retrieve Index-Id-Maps for algorithm data.");

    long long int dimsAlg[2] = {
      indexIdMapAlgP->GetNumberOfTuples(), indexIdMapAlgC->GetNumberOfTuples()};
    long long int dimsGT[2] = {
      indexIdMapGTP->GetNumberOfTuples(), indexIdMapGTC->GetNumberOfTuples()};

    if(dimsAlg[0] != dimsGT[0] || dimsAlg[1] != dimsGT[1])
      return !this->printErr("Dimensions mismatch, filter assums both matrices "
                             "has the same features.");

    int status = 0;
    ttkTypeMacroAA(algVals->GetDataType(), indexIdMapAlgP->GetDataType(),
                   (status = this->compareEvents<T0, T1>(
                      matrixEvents[i], ttkUtils::GetPointer<T0>(algVals),
                      ttkUtils::GetPointer<T0>(gtVals),
                      ttkUtils::GetPointer<T1>(indexIdMapAlgP),
                      ttkUtils::GetPointer<T1>(indexIdMapAlgC),
                      ttkUtils::GetPointer<T1>(indexIdMapGTP),
                      ttkUtils::GetPointer<T1>(indexIdMapGTC),
                      dimsGT)));
    if(!status)
      return 0;
  }

  auto output = vtkTable::GetData(outputVector, 0);

  // Function for formatting data arrays
  auto prepArray
    = [](vtkDataArray *array, std::string name, int nTuples, int nComponents) {
        array->SetName(name.data());
        array->SetNumberOfComponents(nComponents);
        array->SetNumberOfTuples(nTuples);
        return ttkUtils::GetVoidPointer(array);
      };

  auto contsAlg = vtkSmartPointer<vtkIntArray>::New();
  auto contsAlgData
    = static_cast<int *>(prepArray(contsAlg, "ContsAlg", nMatrices + 1, 1));
  output->AddColumn(contsAlg);
  auto contsGt = vtkSmartPointer<vtkIntArray>::New();
  auto contsGtData
    = static_cast<int *>(prepArray(contsGt, "ContsGT", nMatrices + 1, 1));
  output->AddColumn(contsGt);
  auto contsOk = vtkSmartPointer<vtkIntArray>::New();
  auto contsOkData
    = static_cast<int *>(prepArray(contsOk, "ContsCorrect", nMatrices + 1, 1));
  output->AddColumn(contsOk);

  auto birthsAlg = vtkSmartPointer<vtkIntArray>::New();
  auto birthsAlgData
    = static_cast<int *>(prepArray(birthsAlg, "BirthsAlg", nMatrices + 1, 1));
  output->AddColumn(birthsAlg);
  auto birthsGt = vtkSmartPointer<vtkIntArray>::New();
  auto birthsGtData
    = static_cast<int *>(prepArray(birthsGt, "BirthsGT", nMatrices + 1, 1));
  output->AddColumn(birthsGt);
  auto birthsOk = vtkSmartPointer<vtkIntArray>::New();
  auto birthsOkData = static_cast<int *>(
    prepArray(birthsOk, "BirthsCorrect", nMatrices + 1, 1));
  output->AddColumn(birthsOk);

  auto deathsAlg = vtkSmartPointer<vtkIntArray>::New();
  auto deathsAlgData
    = static_cast<int *>(prepArray(deathsAlg, "DeathsAlg", nMatrices + 1, 1));
  output->AddColumn(deathsAlg);
  auto deathsGt = vtkSmartPointer<vtkIntArray>::New();
  auto deathsGtData
    = static_cast<int *>(prepArray(deathsGt, "DeathsGT", nMatrices + 1, 1));
  output->AddColumn(deathsGt);
  auto deathsOk = vtkSmartPointer<vtkIntArray>::New();
  auto deathsOkData = static_cast<int *>(
    prepArray(deathsOk, "DeathsCorrect", nMatrices + 1, 1));
  output->AddColumn(deathsOk);

  auto splitsAlg = vtkSmartPointer<vtkIntArray>::New();
  auto splitsAlgData
    = static_cast<int *>(prepArray(splitsAlg, "SplitsAlg", nMatrices + 1, 1));
  output->AddColumn(splitsAlg);
  auto splitsGt = vtkSmartPointer<vtkIntArray>::New();
  auto splitsGtData
    = static_cast<int *>(prepArray(splitsGt, "SplitsGT", nMatrices + 1, 1));
  output->AddColumn(splitsGt);
  auto splitsOk = vtkSmartPointer<vtkIntArray>::New();
  auto splitsOkData = static_cast<int *>(
    prepArray(splitsOk, "SplitsCorrect", nMatrices + 1, 1));
  output->AddColumn(splitsOk);

  auto mergesAlg = vtkSmartPointer<vtkIntArray>::New();
  auto mergesAlgData
    = static_cast<int *>(prepArray(mergesAlg, "MergesAlg", nMatrices + 1, 1));
  output->AddColumn(mergesAlg);
  auto mergesGt = vtkSmartPointer<vtkIntArray>::New();
  auto mergesGtData
    = static_cast<int *>(prepArray(mergesGt, "MergesGT", nMatrices + 1, 1));
  output->AddColumn(mergesGt);
  auto mergesOk = vtkSmartPointer<vtkIntArray>::New();
  auto mergesOkData = static_cast<int *>(
    prepArray(mergesOk, "MergesCorrect", nMatrices + 1, 1));
  output->AddColumn(mergesOk);

  //int nColumns = 15;
  //int totals[nColumns] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int totals[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for(int i = 0; i < nMatrices; i++) {
    auto me = matrixEvents[i];

    contsAlgData[i] = me.continuations[0];
    totals[0] += me.continuations[0];
    contsGtData[i] = me.continuations[1];
    totals[1] += me.continuations[1];
    contsOkData[i] = me.continuations[2];
    totals[2] += me.continuations[2];

    birthsAlgData[i] = me.births[0];
    totals[3] += me.births[0];
    birthsGtData[i] = me.births[1];
    totals[4] += me.births[1];
    birthsOkData[i] = me.births[2];
    totals[5] += me.births[2];

    deathsAlgData[i] = me.deaths[0];
    totals[6] += me.deaths[0];
    deathsGtData[i] = me.deaths[1];
    totals[7] += me.deaths[1];
    deathsOkData[i] = me.deaths[2];
    totals[8] += me.deaths[2];

    splitsAlgData[i] = me.splits[0];
    totals[9] += me.splits[0];
    splitsGtData[i] = me.splits[1];
    totals[10] += me.splits[1];
    splitsOkData[i] = me.splits[2];
    totals[11] += me.splits[2];

    mergesAlgData[i] = me.merges[0];
    totals[12] += me.merges[0];
    mergesGtData[i] = me.merges[1];
    totals[13] += me.merges[1];
    mergesOkData[i] = me.merges[2];
    totals[14] += me.merges[2];
  }

  contsAlgData[nMatrices] = totals[0];
  contsGtData[nMatrices] = totals[1];
  contsOkData[nMatrices] = totals[2];

  birthsAlgData[nMatrices] = totals[3];
  birthsGtData[nMatrices] = totals[4];
  birthsOkData[nMatrices] = totals[5];

  deathsAlgData[nMatrices] = totals[6];
  deathsGtData[nMatrices] = totals[7];
  deathsOkData[nMatrices] = totals[8];

  splitsAlgData[nMatrices] = totals[9];
  splitsGtData[nMatrices] = totals[10];
  splitsOkData[nMatrices] = totals[11];

  mergesAlgData[nMatrices] = totals[12];
  mergesGtData[nMatrices] = totals[13];
  mergesOkData[nMatrices] = totals[14];

  return 1;
}
