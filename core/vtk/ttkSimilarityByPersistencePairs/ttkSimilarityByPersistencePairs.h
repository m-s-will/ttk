#pragma once

// VTK Module
#include <ttkSimilarityByPersistencePairsModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <SimilarityByPersistencePairs.h>

class TTKSIMILARITYBYPERSISTENCEPAIRS_EXPORT ttkSimilarityByPersistencePairs
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByPersistencePairs {

public:
  static ttkSimilarityByPersistencePairs *New();
  vtkTypeMacro(ttkSimilarityByPersistencePairs, ttkSimilarityAlgorithm);

  vtkSetMacro(Lifting, double);
  vtkGetMacro(Lifting, double);

  vtkSetMacro(MaxJump, double);
  vtkGetMacro(MaxJump, double);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, const std::string &);
  vtkGetMacro(DistanceAlgorithm, std::string);

  template <typename dataType>
  int getDiagram(std::vector<std::tuple<int,
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
                 const int diagramNumber);

protected:
  ttkSimilarityByPersistencePairs();
  ~ttkSimilarityByPersistencePairs();

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;

private:
  // Metric weights
  double PX{0.};
  double PY{0.};
  double PZ{0.};
  double PE{1.}; // extrema
  double PS{1.}; // saddles
  // Metric config
  double MaxJump{0.1};
  double Lifting{0.0};
  std::string DistanceAlgorithm{"0"}; // distance between PPs
  std::string WassersteinMetric{"2"}; // wass vs inf (bottleneck)
  int PVAlgorithm{-1};
};
