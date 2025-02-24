#include <ttkMergeTreeClustering.h>
#include <ttkMergeTreeFeatureTracking.h>
#include <ttkMergeTreeUtils.h>
#include <ttkMergeTreeVisualization.h>
#include <ttkUtils.h>

#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkStringArray.h>
#include <vtkTable.h>

using namespace ttk;
using namespace ftm;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreeFeatureTracking);

/**
 * Implement the filter constructor and destructor in the cpp file.
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
ttkMergeTreeFeatureTracking::ttkMergeTreeFeatureTracking() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(3);
}

ttkMergeTreeFeatureTracking::~ttkMergeTreeFeatureTracking() = default;

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreeFeatureTracking::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    if(port == 1)
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  } else
    return 0;

  return 1;
}

/**
 * Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( ttkAlgorithm::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMergeTreeFeatureTracking::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 or port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  else if(port == 2)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

/**
 * Pass VTK data to the base code and convert base code output to VTK
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
int ttkMergeTreeFeatureTracking::runCompute() {
  // Verify parameters
  if(not UseFieldDataParameters) {
    if(Backend == 0) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = true;
      keepSubtree_ = false;
      baseModule_ = 0;
    } else if(Backend == 1) {
      branchDecomposition_ = false;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
      baseModule_ = 0;
    } else if(Backend == 3) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
      baseModule_ = 1;
    } else if(Backend == 4) {
      branchDecomposition_ = true;
      normalizedWasserstein_ = false;
      keepSubtree_ = true;
      baseModule_ = 2;
    } else {
      baseModule_ = 0;
    }
  }
  if(baseModule_ == 0) {
    if(isPersistenceDiagram_) {
      branchDecomposition_ = true;
    }
    if(not branchDecomposition_) {
      if(normalizedWasserstein_)
        printMsg("NormalizedWasserstein is set to false since branch "
                 "decomposition is not asked.");
      normalizedWasserstein_ = false;
    }
    if(normalizedWasserstein_)
      printMsg("Computation with normalized Wasserstein.");
    else
      printMsg("Computation without normalized Wasserstein.");
    epsilonTree2_ = epsilonTree1_;
    epsilon2Tree2_ = epsilon2Tree1_;
    epsilon3Tree2_ = epsilon3Tree1_;
    printMsg("BranchDecomposition: " + std::to_string(branchDecomposition_));
    printMsg("NormalizedWasserstein: "
             + std::to_string(normalizedWasserstein_));
    printMsg("KeepSubtree: " + std::to_string(keepSubtree_));
  }
  if(baseModule_ == 1) {
    printMsg("Using Branch Mapping Distance.");
    std::string metric;
    if(branchMetric_ == 0)
      metric = "Wasserstein Distance first degree";
    else if(branchMetric_ == 1)
      metric = "Wasserstein Distance second degree";
    else if(branchMetric_ == 2)
      metric = "Persistence difference";
    else if(branchMetric_ == 3)
      metric = "Shifting cost";
    else
      return 1;
    printMsg("BranchMetric: " + metric);
  }
  if(baseModule_ == 2) {
    printMsg("Using Path Mapping Distance.");
    std::string metric;
    if(pathMetric_ == 0)
      metric = "Persistence difference";
    else
      return 1;
    printMsg("PathMetric: " + metric);
  }

  // --- Call base
  execute<float>(
    intermediateTrees, intermediateTrees2, outputMatchings, distances);
  treesNodeCorrMesh = getTreesNodeCorr();

  return 1;
}

int ttkMergeTreeFeatureTracking::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &ttkNotUsed(inputTrees2)) {
  const unsigned int numInputs = inputTrees.size();
  // --------------------------------------------------------------------------
  // --- Create output
  // --------------------------------------------------------------------------
  auto vtkOutputTrees = vtkMultiBlockDataSet::GetData(outputVector, 0);
  auto vtkOutputMatchings = vtkMultiBlockDataSet::GetData(outputVector, 1);
  auto vtkOutputCurves = vtkUnstructuredGrid::GetData(outputVector, 2);

  // ------------------------------------------------------
  // --- Trees
  // ------------------------------------------------------
  if(isPersistenceDiagram_ and not OutputSegmentation) {
    vtkOutputTrees->SetNumberOfBlocks(numInputs);
  } else {
    vtkOutputTrees->SetNumberOfBlocks((OutputSegmentation ? 3 : 2));
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockNodes
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    vtkBlockNodes->SetNumberOfBlocks(numInputs);
    vtkOutputTrees->SetBlock(0, vtkBlockNodes);
    vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockArcs
      = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    vtkBlockArcs->SetNumberOfBlocks(numInputs);
    vtkOutputTrees->SetBlock(1, vtkBlockArcs);
    if(OutputSegmentation) {
      vtkSmartPointer<vtkMultiBlockDataSet> const vtkBlockSegs
        = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      vtkBlockSegs->SetNumberOfBlocks(numInputs);
      vtkOutputTrees->SetBlock(2, vtkBlockSegs);
    }
  }

  std::vector<FTMTree_MT *> trees;
  mergeTreeToFTMTree<float>(intermediateTrees, trees);
  double prevXMax = 0;
  std::vector<std::vector<SimplexId>> nodeCorr(intermediateTrees.size());
  for(unsigned int i = 0; i < numInputs; ++i) {
    // Declare vtk objects
    vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputArc
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkDataSet *vtkOutputSegmentation{};
    if(treesSegmentation[i])
      vtkOutputSegmentation = treesSegmentation[i]->NewInstance();

    // Fill vtk objects
    ttkMergeTreeVisualization visuMaker;
    visuMaker.setPlanarLayout(PlanarLayout);
    visuMaker.setBranchDecompositionPlanarLayout(
      BranchDecompositionPlanarLayout);
    visuMaker.setBranchSpacing(BranchSpacing);
    // visuMaker.setNonImportantBranchSpacing(NonImportantBranchSpacing);
    visuMaker.setRescaleTreesIndividually(RescaleTreesIndividually);
    visuMaker.setOutputSegmentation(OutputSegmentation);
    visuMaker.setDimensionSpacing(DimensionSpacing);
    visuMaker.setDimensionToShift(DimensionToShift);
    visuMaker.setDimensionsShift(XShift, YShift, ZShift);
    visuMaker.setImportantPairs(ImportantPairs);
    visuMaker.setMaximumImportantPairs(MaximumImportantPairs);
    visuMaker.setMinimumImportantPairs(MinimumImportantPairs);
    visuMaker.setImportantPairsSpacing(ImportantPairsSpacing);
    visuMaker.setNonImportantPairsSpacing(NonImportantPairsSpacing);
    visuMaker.setNonImportantPairsProximity(NonImportantPairsProximity);
    visuMaker.setExcludeImportantPairsHigher(ExcludeImportantPairsHigher);
    visuMaker.setExcludeImportantPairsLower(ExcludeImportantPairsLower);
    visuMaker.setIsPersistenceDiagram(isPersistenceDiagram_);
    visuMaker.setTreesNodes(treesNodes);
    visuMaker.copyPointData(treesNodes[i], treesNodeCorrMesh[i]);
    visuMaker.setTreesNodeCorrMesh(treesNodeCorrMesh);
    visuMaker.setTreesSegmentation(treesSegmentation);
    visuMaker.setVtkOutputNode(vtkOutputNode);
    if(isPersistenceDiagram_)
      visuMaker.setVtkOutputArc(vtkOutputNode);
    else
      visuMaker.setVtkOutputArc(vtkOutputArc);
    visuMaker.setVtkOutputSegmentation(vtkOutputSegmentation);
    visuMaker.setPrintTreeId(i);
    visuMaker.setPrintClusterId(0);
    visuMaker.setDebugLevel(this->debugLevel_);
    visuMaker.setIsPDSadMax(mixtureCoefficient_ == 0);

    visuMaker.setShiftMode(2); // Line
    visuMaker.setPrevXMaxOffset(prevXMax);

    visuMaker.makeTreesOutput<float>(trees);
    prevXMax = visuMaker.getPrevXMax();
    nodeCorr[i] = visuMaker.getNodeCorr()[i];

    // Field data
    vtkOutputNode->GetFieldData()->ShallowCopy(treesNodes[i]->GetFieldData());
    if(not isPersistenceDiagram_)
      vtkOutputArc->GetFieldData()->ShallowCopy(treesArcs[i]->GetFieldData());
    /*if(treesSegmentation[i])
      ttkMergeTreeClustering::addFieldData(treesSegmentation[i],
      vtkOutputNode);*/
    if(OutputSegmentation)
      vtkOutputSegmentation->GetFieldData()->ShallowCopy(
        treesSegmentation[i]->GetFieldData());

    // Construct multiblock
    if(isPersistenceDiagram_ and not OutputSegmentation) {
      vtkOutputTrees->SetBlock(i, vtkOutputNode);
    } else {
      vtkMultiBlockDataSet::SafeDownCast(vtkOutputTrees->GetBlock(0))
        ->SetBlock(i, vtkOutputNode);
      if(not isPersistenceDiagram_)
        vtkMultiBlockDataSet::SafeDownCast(vtkOutputTrees->GetBlock(1))
          ->SetBlock(i, vtkOutputArc);
      if(OutputSegmentation) {
        int const segBlockID = 1 + !isPersistenceDiagram_;
        vtkMultiBlockDataSet::SafeDownCast(vtkOutputTrees->GetBlock(segBlockID))
          ->SetBlock(i, vtkOutputSegmentation);
      }
    }
  }

  // ------------------------------------------------------
  // --- Matchings
  // ------------------------------------------------------
  vtkOutputMatchings->SetNumberOfBlocks(intermediateTrees.size() - 1);
  for(size_t i = 0; i < intermediateTrees.size() - 1; ++i) {
    // Declare vtk objects
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputMatching
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNode1
      = vtkUnstructuredGrid::SafeDownCast(
        vtkMultiBlockDataSet::SafeDownCast(vtkOutputTrees->GetBlock(0))
          ->GetBlock(i));
    vtkSmartPointer<vtkUnstructuredGrid> const vtkOutputNode2
      = vtkUnstructuredGrid::SafeDownCast(
        vtkMultiBlockDataSet::SafeDownCast(vtkOutputTrees->GetBlock(0))
          ->GetBlock(i + 1));
    std::vector<std::vector<SimplexId>> nodeCorrTemp;
    nodeCorrTemp.emplace_back(nodeCorr[i]);
    nodeCorrTemp.emplace_back(nodeCorr[i + 1]);

    // Fill vtk objects
    ttkMergeTreeVisualization visuMakerMatching;
    visuMakerMatching.setVtkOutputMatching(vtkOutputMatching);
    visuMakerMatching.setOutputMatching(outputMatchings[i]);
    visuMakerMatching.setVtkOutputNode1(vtkOutputNode2);
    visuMakerMatching.setVtkOutputNode2(vtkOutputNode1);
    visuMakerMatching.setNodeCorr1(nodeCorrTemp);
    visuMakerMatching.setDebugLevel(this->debugLevel_);

    visuMakerMatching.makeMatchingOutput<float>(trees[i], trees[i + 1]);

    // Field data
    vtkNew<vtkDoubleArray> vtkDistance{};
    vtkDistance->SetName("Distance");
    vtkDistance->SetNumberOfTuples(1);
    vtkDistance->SetTuple1(0, distances[i]);
    vtkOutputMatching->GetFieldData()->AddArray(vtkDistance);

    // Construct multiblock
    vtkOutputMatchings->SetBlock(i, vtkOutputMatching);
  }

  // ------------------------------------------------------
  // --- Curves
  // ------------------------------------------------------
  vtkSmartPointer<vtkUnstructuredGrid> const vtuCurves
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(distances.size());

  for(unsigned int i = 0; i < distances.size(); ++i) {
    points->SetPoint(i, i, distances[i], 0);
    if(i != 0) {
      vtkIdType pointIds[2];
      pointIds[0] = i - 1;
      pointIds[1] = i;
      vtuCurves->InsertNextCell(VTK_LINE, 2, pointIds);
    }
  }
  vtuCurves->SetPoints(points);

  vtkOutputCurves->ShallowCopy(vtuCurves);

  return 1;
}

int ttkMergeTreeFeatureTracking::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // --- Get input object from input vector
  auto blocks = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  auto blocks2 = vtkMultiBlockDataSet::GetData(inputVector[1], 0);

  // --- Load blocks
  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> inputTrees, inputTrees2;
  loadBlocks(inputTrees, blocks);
  loadBlocks(inputTrees2, blocks2);

  bool doCompute = oldBlocks != blocks;
  std::vector<ttk::ftm::MergeTree<float>> tempTrees, tempTrees2;

  // Construct trees
  bool const useSadMaxPairs = (mixtureCoefficient_ == 0); // only for PD support
  isPersistenceDiagram_ = constructTrees<float>(
    inputTrees, (doCompute ? intermediateTrees : tempTrees), treesNodes,
    treesArcs, treesSegmentation, useSadMaxPairs);
  if(not isPersistenceDiagram_
     or (mixtureCoefficient_ != 0 and mixtureCoefficient_ != 1)) {
    auto &inputTrees2ToUse
      = (not isPersistenceDiagram_ ? inputTrees2 : inputTrees);
    constructTrees(inputTrees2ToUse,
                   (doCompute ? intermediateTrees2 : tempTrees2),
                   !useSadMaxPairs);
  }

  // --- Load field data parameters
  if(UseFieldDataParameters) {
    printMsg("Load parameters from field data.");
    std::vector<std::string> paramNames;
    getParamNames(paramNames);
    for(auto paramName : paramNames) {
      auto array = blocks->GetFieldData()->GetArray(paramName.c_str());
      if(array) {
        double const value = array->GetTuple1(0);
        setParamValueFromName(paramName, value);
        printMsg(" - " + paramName + " = " + std::to_string(value));
      } else
        printMsg(" - " + paramName + " was not found in the field data.");
    }
  }

  if(doCompute) {
    auto res = runCompute();
    if(res != 1)
      return res;
  }
  Timer t_output;
  auto res = runOutput(outputVector, inputTrees, inputTrees2);
  printMsg("Output", 1, t_output.getElapsedTime(), this->threadNumber_);
  oldBlocks = blocks;
  return res;
}
