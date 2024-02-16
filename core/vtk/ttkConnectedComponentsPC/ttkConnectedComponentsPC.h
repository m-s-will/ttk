/// \ingroup vtk
/// \class ttkConnectedComponentsPC
/// \author Robin Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK VTK-filter that wraps the ttk::ConnectedComponentsPC module.
///
/// Given an input order field, this filter computes its ascending and
/// descending segmentation by assigning every vertex to its minimum or maximum
/// in gradient or inverse gradient direction. For convienience a hash (no hash
/// collision detection) of both segmentations can be created to represent the
/// Morse-Smale segmentation.
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output output scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet) with resulting arrays attached to point data
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// This filter can be used like any other VTK filter (for instance, by using
/// the sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/ConnectedComponentsPC/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttk::ConnectedComponentsPC
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleSegmentation_at/">Morse-Smale
///   segmentation example</a> \n

#pragma once

// VTK Module
#include <ttkConnectedComponentsPCModule.h>

// ttk code includes
#include <ConnectedComponentsPC.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkPolyData;

class TTKCONNECTEDCOMPONENTSPC_EXPORT ttkConnectedComponentsPC
  : public ttkAlgorithm,
    protected ttk::ConnectedComponentsPC {

public:
  static ttkConnectedComponentsPC *New();

  vtkTypeMacro(ttkConnectedComponentsPC, ttkAlgorithm);


  vtkSetMacro(MinSize, int);
  vtkGetMacro(MinSize, int);

  vtkSetMacro(ForceFeatureScalarField, bool);
  vtkGetMacro(ForceFeatureScalarField, bool);


protected:
  ttkConnectedComponentsPC();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  int MinSize{};
  bool ForceFeatureScalarField{false};

};