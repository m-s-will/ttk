#pragma once

// VTK Module
#include <ttkSimilarityByGradientModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityByGradient.h>

// stl
#include <numeric>
#include <set>
#include <vector>

class TTKSIMILARITYBYGRADIENT_EXPORT ttkSimilarityByGradient
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByGradient {

private:
  bool UseNeighbourhood{false};
  int NeighbourhoodSize{0};
  float NeighbourhoodDistance{0.0};

public:
  static ttkSimilarityByGradient *New();
  vtkTypeMacro(ttkSimilarityByGradient, ttkSimilarityAlgorithm);

  vtkSetMacro(UseNeighbourhood, bool);
  vtkGetMacro(UseNeighbourhood, bool);

  vtkSetMacro(NeighbourhoodSize, int);
  vtkGetMacro(NeighbourhoodSize, int);

  vtkSetMacro(NeighbourhoodDistance, float);
  vtkGetMacro(NeighbourhoodDistance, float);

protected:
  ttkSimilarityByGradient();
  ~ttkSimilarityByGradient();

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;
};