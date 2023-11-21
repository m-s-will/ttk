#pragma once

// VTK Module
#include <ttkSimilarityByGradientV3Module.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityByGradientV3.h>

// stl
#include <vector>
#include <set>
#include <numeric>

class TTKSIMILARITYBYGRADIENTV3_EXPORT ttkSimilarityByGradientV3
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByGradientV3 {

private:
  bool UseNeighbourhood{false};
  int NeighbourhoodSize{0};
  float NeighbourhoodDistance{0.0};

public:

  static ttkSimilarityByGradientV3 *New();
  vtkTypeMacro(ttkSimilarityByGradientV3, ttkSimilarityAlgorithm);

  vtkSetMacro(UseNeighbourhood, bool);
  vtkGetMacro(UseNeighbourhood, bool);

  vtkSetMacro(NeighbourhoodSize, int);
  vtkGetMacro(NeighbourhoodSize, int);

  vtkSetMacro(NeighbourhoodDistance, float);
  vtkGetMacro(NeighbourhoodDistance, float);

protected:
  ttkSimilarityByGradientV3();
  ~ttkSimilarityByGradientV3();

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;
};