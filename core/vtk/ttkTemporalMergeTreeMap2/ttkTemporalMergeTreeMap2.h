/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkTemporalMergeTreeMap2
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::TemporalMergeTreeMap2 module.
///
/// This VTK filter uses the ttk::TemporalMergeTreeMap2 module to compute an averaging of
/// the data values of an input point data array defined on the input
/// vtkDataSet.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/TemporalMergeTreeMap2/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::TemporalMergeTreeMap2
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkTemporalMergeTreeMap2Module.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>

/* Note on including VTK modules
 *
 * Each VTK module that you include a header from needs to be specified in this
 * module's vtk.module file, either in the DEPENDS or PRIVATE_DEPENDS (if the
 * header is included in the cpp file only) sections.
 *
 * In order to find the corresponding module, check its location within the VTK
 * source code. The VTK module name is composed of the path to the header. You
 * can also find the module name within the vtk.module file located in the same
 * directory as the header file.
 *
 * For example, vtkSphereSource.h is located in directory VTK/Filters/Sources/,
 * so its corresponding VTK module is called VTK::FiltersSources. In this case,
 * the vtk.module file would need to be extended to
 *
 * NAME
 *   ttkTemporalMergeTreeMap2
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <TemporalMergeTreeMap2.h>

class TTKTEMPORALMERGETREEMAP2_EXPORT ttkTemporalMergeTreeMap2
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::TemporalMergeTreeMap2 // and we inherit from the base class
{
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};
  bool useSlidingWindow = false;
  ttk::SimplexId windowSize = 5;
  ttk::SimplexId layoutMode = 2; // 1=matchings, 2=barycenter
  void dfs_linearization(
      ttk::SimplexId curr_node,
      std::vector<double> &lin,
      std::vector<ttk::SimplexId> &seg,
      std::vector<ttk::SimplexId> &bar,
      std::vector<double> &nodePositions,
      std::vector<std::vector<ttk::SimplexId>> &memiChildren,
      std::vector<std::vector<double>> &memiSegmentScalars,
      std::vector<ttk::SimplexId> &memiSizes,
      std::vector<ttk::SimplexId> &memiSegs,
      std::vector<ttk::SimplexId> &branchNodeIDs,
      std::vector<double> &memiScalars,
      std::vector<double> &memiOrdering,
      std::vector<ttk::SimplexId> prevMatching,
      std::vector<double> prevOrdering);
  void computeBaryBranchOrdering(
    vtkMultiBlockDataSet* mtmb, 
    vtkMultiBlockDataSet* members,
    std::vector<double> &ordering_branches);

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  void SetWindowSize(int s) {
    windowSize = s;
    Modified();
  }
  vtkGetMacro(windowSize, int);

  void SetLayoutMode(int m) {
    layoutMode = m;
    Modified();
  }
  vtkGetMacro(layoutMode, int);

  void SetUseSlidingWindow(bool b) {
    useSlidingWindow = b;
    Modified();
  }
  vtkGetMacro(useSlidingWindow, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkTemporalMergeTreeMap2 *New();
  vtkTypeMacro(ttkTemporalMergeTreeMap2, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkTemporalMergeTreeMap2();
  ~ttkTemporalMergeTreeMap2() override = default;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 9: Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
