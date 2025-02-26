/// \ingroup vtk
/// \class ttkMergeTreeFeatureTracking
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2021.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeTreeFeatureTracking
/// module.
///
/// This VTK filter uses the ttk::MergeTreeFeatureTracking module to compute
/// the distance matrix of a group of merge trees.
///
/// \param Input vtkMultiBlockDataset
/// \param Output vtkTable
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreeFeatureTracking
/// \sa ttkAlgorithm
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreeClustering/">Merge
///   Tree Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/mergeTreePGA/">Merge
///   Tree Principal Geodesic Analysis example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramPGA/">
///   Persistence Diagram Principal Geodesic Analysis example</a> \n

#pragma once

// VTK Module
#include <ttkMergeTreeFeatureTrackingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MergeTreeFeatureTracking.h>

class TTKMERGETREEFEATURETRACKING_EXPORT ttkMergeTreeFeatureTracking
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreeFeatureTracking // and we inherit from the base
                                            // class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  // Execution Options
  int Backend = 0;
  bool oldBD = branchDecomposition_;
  bool oldNW = normalizedWasserstein_;
  bool oldKS = keepSubtree_;

  bool UseFieldDataParameters = false;

  // Output options
  bool OutputTrees = true;
  bool OutputSegmentation = false;
  bool PlanarLayout = false;
  bool BranchDecompositionPlanarLayout = false;
  bool PathPlanarLayout = false;
  double BranchSpacing = 1.;
  double NonImportantBranchSpacing = 1.;
  bool RescaleTreesIndividually = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  double XShift = 1.0;
  double YShift = 0.0;
  double ZShift = 0.0;
  double ImportantPairs = 50.;
  int MaximumImportantPairs = 0;
  int MinimumImportantPairs = 0;
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 1.;
  double NonImportantPairsProximity = 0.05;
  std::string ExcludeImportantPairsLower = "";
  std::string ExcludeImportantPairsHigher = "";

  //
  vtkMultiBlockDataSet *oldBlocks = 0;
  std::vector<std::vector<int>> treesNodeCorrMesh;
  std::vector<ttk::ftm::MergeTree<float>> intermediateTrees, intermediateTrees2;
  std::vector<vtkUnstructuredGrid *> treesNodes, treesArcs;
  std::vector<vtkDataSet *> treesSegmentation;
  std::vector<
    std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>
    outputMatchings;
  std::vector<float> distances;

  void doCompute() {
    oldBlocks = 0;
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input Options
#define ttkMergeTreeFeatureTrackingSetMacro(name, type) \
  void Set##name(type v) {                              \
    name = v;                                           \
    Modified();                                         \
    doCompute();                                        \
  }

  void SetEpsilon1UseFarthestSaddle(bool epsilon1UseFarthestSaddle) {
    epsilon1UseFarthestSaddle_ = epsilon1UseFarthestSaddle;
    Modified();
    doCompute();
  }
  bool GetEpsilon1UseFarthestSaddle() {
    return epsilon1UseFarthestSaddle_;
  }

  void SetEpsilonTree1(double epsilonTree1) {
    epsilonTree1_ = epsilonTree1;
    Modified();
    doCompute();
  }
  double SetEpsilonTree1() {
    return epsilonTree1_;
  }

  void SetEpsilon2Tree1(double epsilon2Tree1) {
    epsilon2Tree1_ = epsilon2Tree1;
    Modified();
    doCompute();
  }
  double SetEpsilon2Tree1() {
    return epsilon2Tree1_;
  }

  void SetEpsilon3Tree1(double epsilon3Tree1) {
    epsilon3Tree1_ = epsilon3Tree1;
    Modified();
    doCompute();
  }
  double SetEpsilon3Tree1() {
    return epsilon3Tree1_;
  }

  void SetPersistenceThreshold(double persistenceThreshold) {
    persistenceThreshold_ = persistenceThreshold;
    Modified();
    doCompute();
  }
  double SetPersistenceThreshold() {
    return persistenceThreshold_;
  }

  void SetDeleteMultiPersPairs(bool doDelete) {
    deleteMultiPersPairs_ = doDelete;
    Modified();
    doCompute();
  }
  bool SetDeleteMultiPersPairs() {
    return deleteMultiPersPairs_;
  }

  void SetBranchMetric(int m) {
    branchMetric_ = m;
    Modified();
    doCompute();
  }

  void SetPathMetric(int m) {
    pathMetric_ = m;
    Modified();
    doCompute();
  }

  // Execution Options
  void SetBackend(int newBackend) {
    if(Backend == 2) { // Custom
      oldBD = branchDecomposition_;
      oldNW = normalizedWasserstein_;
      oldKS = keepSubtree_;
    }
    if(newBackend == 2) { // Custom
      branchDecomposition_ = oldBD;
      normalizedWasserstein_ = oldNW;
      keepSubtree_ = oldKS;
    }
    Backend = newBackend;
    Modified();
    doCompute();
  }
  vtkGetMacro(Backend, int);

  void SetAssignmentSolver(int assignmentSolver) {
    assignmentSolverID_ = assignmentSolver;
    Modified();
    doCompute();
  }
  int GetAssignmentSolver() {
    return assignmentSolverID_;
  }

  void SetBranchDecomposition(bool branchDecomposition) {
    branchDecomposition_ = branchDecomposition;
    Modified();
    doCompute();
  }
  int GetBranchDecomposition() {
    return branchDecomposition_;
  }

  void SetNormalizedWasserstein(bool normalizedWasserstein) {
    normalizedWasserstein_ = normalizedWasserstein;
    Modified();
    doCompute();
  }
  int GetNormalizedWasserstein() {
    return normalizedWasserstein_;
  }

  void SetKeepSubtree(bool keepSubtree) {
    keepSubtree_ = keepSubtree;
    Modified();
    doCompute();
  }
  int GetKeepSubtree() {
    return keepSubtree_;
  }

  void SetDistanceSquaredRoot(bool distanceSquaredRoot) {
    distanceSquaredRoot_ = distanceSquaredRoot;
    Modified();
    doCompute();
  }
  int GetDistanceSquaredRoot() {
    return distanceSquaredRoot_;
  }

  vtkSetMacro(UseFieldDataParameters, bool);
  vtkGetMacro(UseFieldDataParameters, bool);

  ttkMergeTreeFeatureTrackingSetMacro(mixtureCoefficient_, double);
  vtkGetMacro(mixtureCoefficient_, double);

  // Output Options
  vtkSetMacro(OutputTrees, bool);
  vtkGetMacro(OutputTrees, bool);

  vtkSetMacro(OutputSegmentation, bool);
  vtkGetMacro(OutputSegmentation, bool);

  vtkSetMacro(PlanarLayout, bool);
  vtkGetMacro(PlanarLayout, bool);

  vtkSetMacro(BranchDecompositionPlanarLayout, bool);
  vtkGetMacro(BranchDecompositionPlanarLayout, bool);

  vtkSetMacro(PathPlanarLayout, bool);
  vtkGetMacro(PathPlanarLayout, bool);

  vtkSetMacro(BranchSpacing, double);
  vtkGetMacro(BranchSpacing, double);

  vtkSetMacro(NonImportantBranchSpacing, double);
  vtkGetMacro(NonImportantBranchSpacing, double);

  vtkSetMacro(RescaleTreesIndividually, bool);
  vtkGetMacro(RescaleTreesIndividually, bool);

  vtkSetMacro(DimensionSpacing, double);
  vtkGetMacro(DimensionSpacing, double);

  vtkSetMacro(DimensionToShift, int);
  vtkGetMacro(DimensionToShift, int);

  vtkSetMacro(XShift, double);
  vtkGetMacro(XShift, double);

  vtkSetMacro(YShift, double);
  vtkGetMacro(YShift, double);

  vtkSetMacro(ZShift, double);
  vtkGetMacro(ZShift, double);

  vtkSetMacro(ImportantPairs, double);
  vtkGetMacro(ImportantPairs, double);

  vtkSetMacro(MaximumImportantPairs, int);
  vtkGetMacro(MaximumImportantPairs, int);

  vtkSetMacro(MinimumImportantPairs, int);
  vtkGetMacro(MinimumImportantPairs, int);

  vtkSetMacro(ImportantPairsSpacing, double);
  vtkGetMacro(ImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsSpacing, double);
  vtkGetMacro(NonImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsProximity, double);
  vtkGetMacro(NonImportantPairsProximity, double);

  vtkSetMacro(ExcludeImportantPairsLower, const std::string &);
  vtkGetMacro(ExcludeImportantPairsLower, std::string);

  vtkSetMacro(ExcludeImportantPairsHigher, const std::string &);
  vtkGetMacro(ExcludeImportantPairsHigher, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreeFeatureTracking *New();
  vtkTypeMacro(ttkMergeTreeFeatureTracking, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkMergeTreeFeatureTracking();
  ~ttkMergeTreeFeatureTracking() override;

  /**
   * Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int runCompute();

  int runOutput(
    vtkInformationVector *outputVector,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees,
    std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> &inputTrees2);
};
