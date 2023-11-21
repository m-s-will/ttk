/// \ingroup vtk
/// \class ttkTrackingGraph
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::TrackingGraph module.
///
/// This VTK filter uses the ttk::TrackingGraph module to compute an averaging
/// of the data values of an input point data array defined on the input
/// vtkDataSet.
///
#pragma once

// VTK Module
#include <ttkTrackingGraphModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <TrackingGraph.h>

class vtkPolyData;
class vtkMultiBlockDataSet;

class TTKTRACKINGGRAPH_EXPORT ttkTrackingGraph : public ttkAlgorithm,
                                                 protected ttk::TrackingGraph {
private:
public:
  static ttkTrackingGraph *New();
  vtkTypeMacro(ttkTrackingGraph, ttkAlgorithm);

protected:
  ttkTrackingGraph();
  ~ttkTrackingGraph() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int CountNodesAndEdges(int &nNodes,
                         int &nEdges,
                         std::vector<int> &nodeIdxOffsets,
                         vtkMultiBlockDataSet *similarities,
                         vtkMultiBlockDataSet *features);
  int Validate(vtkMultiBlockDataSet *similarities,
               vtkMultiBlockDataSet *features);

  int GenerateTrackingGraphFromFeatures(vtkPolyData *output,
                                        vtkMultiBlockDataSet *similarities,
                                        vtkMultiBlockDataSet *features);
  int GenerateTrackingGraphFromMatrix(vtkPolyData *output,
                                      vtkMultiBlockDataSet *similarities);
};