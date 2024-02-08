#include <ttkSimilarityByPersistencePairs.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkSimilarityByPersistencePairs);

ttkSimilarityByPersistencePairs::ttkSimilarityByPersistencePairs() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkSimilarityByPersistencePairs::~ttkSimilarityByPersistencePairs() {
}

int ttkSimilarityByPersistencePairs::ComputeSimilarityMatrix(
  vtkImageData *similarityMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {
  int status = 0;

  // unpack input
  vtkUnstructuredGrid *p0 = vtkUnstructuredGrid::SafeDownCast(
    inputDataObjects0); // set of persist. pairs from previous timestep
  vtkUnstructuredGrid *p1 = vtkUnstructuredGrid::SafeDownCast(
    inputDataObjects1); // PPs from current timestep

  if(!p0 || !p1)
    return !this->printErr(
      "Input data objects need to be vtkUnstructuredGrids.");

  // get point coordinates
  auto coords0 = p0->GetPoints()->GetData();
  auto coords1 = p1->GetPoints()->GetData();

  // get and format persistence diagrams
  using dataType = double;
  using dT
    = std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType, int,
                 dataType, float, float, float, dataType, float, float, float>;
  using mT = std::tuple<int, int, double>;

  std::vector<dT> CTDiagram0;
  std::vector<dT> CTDiagram1;

  const double spacing
    = 0; // 2d tracking -> height spacing parameter for 3d display
  status = getDiagram<double>(CTDiagram0, p0, spacing, 0);
  if(status < 0)
    return !this->printErr(
      "Could not extract diagram from first input data-set");

  status = getDiagram<double>(CTDiagram1, p1, spacing, 1);
  if(status < 0)
    return !this->printErr(
      "Could not extract diagram from second input data-set");

  if(coords0->GetDataType() != coords1->GetDataType())
    return !this->printErr("Input diagrams need to have the same data type.");

  // compute similarity matrix dimensions
  const int nFeatures0 = CTDiagram0.size();
  const int nFeatures1 = CTDiagram1.size();

  // initialize similarity matrix i.e., distance matrix
  similarityMatrix->SetDimensions(nFeatures0, nFeatures1, 1);
  similarityMatrix->AllocateScalars(VTK_FLOAT, 1); // matching output = float

  auto similaritiesArray = similarityMatrix->GetPointData()->GetArray(0);
  similaritiesArray->SetName("LiftedWassersteinDistance");
  auto similarityMatrixData = ttkUtils::GetPointer<float>(similaritiesArray);
  for(int i = 0; i < nFeatures0; ++i)
    for(int j = 0; j < nFeatures1; ++j) {
      similarityMatrixData[j * nFeatures0 + i] = 0;
    }

  // get metric parameters
  const std::string algorithm = DistanceAlgorithm;
  const std::string wasserstein = WassersteinMetric;
  const int pvAlgorithm = PVAlgorithm;
  double px;
  double py;
  double pz;
  double ps;
  double pe;
  const double maxJump = MaxJump;

  // Using advanced parameters if the user has tweaked them.
  if(PX != 0. || PY != 0. || PZ != 0. || PE != 1. || PS != 1.) {
    px = PX;
    py = PY;
    pz = PZ;
    pe = PE;
    ps = PS;
  } else // Using lifting parameter instead.
  {
    // 0 -> persistence; 1 -> geometry
    double geometricalLift = Lifting / 100.0;
    double persistenceLift = 1.0 - geometricalLift;
    px = geometricalLift;
    py = geometricalLift;
    pz = geometricalLift;
    pe = persistenceLift;
    ps = persistenceLift;
  }

  // compute similarities in basecode
  std::vector<mT> matchings;
  switch(coords0->GetDataType()) {
    vtkTemplateMacro(status = this->computeDistanceMatrix<double>(
                       CTDiagram0, CTDiagram1, matchings, px, py, pz, ps, pe,
                       algorithm, wasserstein, pvAlgorithm, maxJump));
  }
  if(status < 0) {
    this->printErr("Error computing distance matrix.");
    return -1;
  }

  // build matrix
  auto matchingsSize = matchings.size();
  if(matchingsSize > 0) {
    for(int i = 0; i < matchingsSize; ++i) {
      vtkIdType ids[2];
      mT t = matchings.at((unsigned long)i);
      auto n1 = (int)std::get<0>(t); // diagram 0
      auto n2 = (int)std::get<1>(t); // diagram 1
      if(n1 >= nFeatures0 || n2 >= nFeatures1)
        return !this->printErr(
          "Invalid indexing: feature index > feature number.");

      similarityMatrixData[n2 * nFeatures0 + n1] = (float)1; // std::get<2>(t);
    }
  }

  status = ttkSimilarityAlgorithm::AddIndexIdMaps(
    similarityMatrix, this->GetInputArrayToProcess(0, p0),
    this->GetInputArrayToProcess(0, p1));
  if(!status)
    return 0;

  return 1;
}

template <typename dataType>
int ttkSimilarityByPersistencePairs::getDiagram(
  std::vector<std::tuple<int,
                         ttk::CriticalType,
                         int,
                         ttk::CriticalType,
                         dataType,
                         int,
                         dataType,
                         float,
                         float,
                         float,
                         dataType,
                         float,
                         float,
                         float>> &diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber) {
  auto pointData = CTPersistenceDiagram_->GetPointData();
  auto cellData = CTPersistenceDiagram_->GetCellData();

  if(pointData == nullptr || cellData == nullptr) {
    return -1;
  }

  auto vertexIdentifierScalars = ttkSimplexIdTypeArray::SafeDownCast(
    pointData->GetArray(ttk::VertexScalarFieldName));

  auto nodeTypeScalars
    = vtkIntArray::SafeDownCast(pointData->GetArray("CriticalType"));
  auto pairIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray("PairIdentifier"));
  auto extremumIndexScalars
    = vtkIntArray::SafeDownCast(cellData->GetArray("PairType"));
  auto persistenceScalars
    = vtkDoubleArray::SafeDownCast(cellData->GetArray("Persistence"));
  auto birthScalars
    = vtkDoubleArray::SafeDownCast(pointData->GetArray("Birth"));
  auto deathScalars
    = vtkDoubleArray::SafeDownCast(pointData->GetArray("Death"));

  vtkPoints *points = CTPersistenceDiagram_->GetPoints();
  if(!pairIdentifierScalars)
    return -2;

  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  // Continuous indexing (no gap in indices)
  for(int pairIndex = 0; pairIndex < pairingsSize; ++pairIndex) {
    const float indexOfPair = pairIndex;
    if(*pairIdentifierScalars->GetTuple(pairIndex) != -1) // except diagonal
      pairIdentifierScalars->SetTuple(pairIndex, &indexOfPair);
  }

  float s{0.0};

  if(!deathScalars != !birthScalars)
    return -2;
  bool is2D = !deathScalars && !birthScalars;
  bool is3D = !is2D;

  if(pairingsSize < 1 || !vertexIdentifierScalars || !nodeTypeScalars
     || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram.resize((unsigned long)pairingsSize);
  int nbNonCompact = 0;

  for(int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = (float)coords1[0];
    auto y1 = (float)coords1[1];
    auto z1 = (float)coords1[2];

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    auto x2 = (float)coords2[0];
    auto y2 = (float)coords2[1];
    auto z2 = (float)coords2[2];

    dataType value1 = (!birthScalars) ? (dataType)x1
                                      : (dataType)birthScalars->GetValue(2 * i);
    dataType value2 = (!deathScalars)
                        ? (dataType)y2
                        : (dataType)deathScalars->GetValue(2 * i + 1);

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram.at(pairIdentifier) = std::make_tuple(
        vertexId1, (ttk::CriticalType)nodeType1, vertexId2,
        (ttk::CriticalType)nodeType2, (dataType)persistence, pairType, value1,
        x1, y1, z1 + s, value2, x2, y2, z2 + s);

    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        this->printWrn(msg.str());
      }
    }
  }

  if(nbNonCompact > 0) {
    std::stringstream msg;
    msg << "Missed " << nbNonCompact << " pairs due to non-compactness."
        << std::endl;
    this->printWrn(msg.str());
  }

  // sort(diagram.begin(), diagram.end(),
  //     [](const diagramTuple &a, const diagramTuple &b) -> bool {
  //       return std::get<6>(a) < std::get<6>(b);
  //     });

  return 1;
}
