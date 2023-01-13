/// \ingroup vtk
/// \class ttkGradientGraph
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date November 2022.
///
/// \brief TTK VTK-filter that wraps the ttk::GradientGraph module.
///
/// This VTK filter uses the ttk::GradientGraph module to compute
/// the Saddle Maximum and the split tree
/// using the Ascending Segmentation and the critical points
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
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
///   - standalone/GradientGraph/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::GradientGraph
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkGradientGraphModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkUnstructuredGrid.h>

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
 *   ttkGradientGraph
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <GradientGraph.h>

class TTKGRADIENTGRAPH_EXPORT ttkGradientGraph
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::GradientGraph // and we inherit from the base class
{
private:
  std::string OutputArrayName{"Pairs"};

public:
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkGradientGraph *New();
  vtkTypeMacro(ttkGradientGraph, ttkAlgorithm);

protected:
  ttkGradientGraph();
  ~ttkGradientGraph() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  template <class triangulationType = ttk::AbstractTriangulation>
  int getSkeletonArcs(vtkUnstructuredGrid *outputSkeletonArcs,
                      std::vector<std::vector<ttk::SimplexId>> &gradientGraph,
                      const triangulationType *triangulation);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
