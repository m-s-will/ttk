#pragma once

// VTK Module
#include <ttkGaussianModeClusteringModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>
#include <vtkPolyData.h>

// TTK Base Includes
#include <GaussianModeClustering.h>

// std includes
#include <map>
#include <utility>
#include <vector>

class TTKGAUSSIANMODECLUSTERING_EXPORT ttkGaussianModeClustering
  : public ttkSimilarityAlgorithm,
    protected ttk::GaussianModeClustering {

private:
  int ClusteringType{0};
  double ScalarThreshold{0.0};

  // For clustering points
  std::vector<std::map<int, std::vector<int>>> clustersPerTimestep_;
  int nClusters0_;
  int nClusters1_;

public:
  vtkSetMacro(ClusteringType, int);
  vtkGetMacro(ClusteringType, int);

  vtkSetMacro(ScalarThreshold, double);
  vtkGetMacro(ScalarThreshold, double);

  static ttkGaussianModeClustering *New();
  vtkTypeMacro(ttkGaussianModeClustering, ttkSimilarityAlgorithm);

protected:
  ttkGaussianModeClustering();
  ~ttkGaussianModeClustering();

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;

  int AddClusterIds(vtkDataObject *inputDataObjects, const size_t t);

  int FormatClusters(vtkDataObject *inputDataObjects,
                     vtkPolyData *outputPoints,
                     const size_t t);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector);
};
