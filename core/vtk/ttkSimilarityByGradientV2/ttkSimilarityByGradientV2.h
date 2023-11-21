#pragma once

// VTK Module
#include <ttkSimilarityByGradientV2Module.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityByGradientV2.h>

// stl
#include <vector>
#include <set>
#include <numeric>

class TTKSIMILARITYBYGRADIENTV2_EXPORT ttkSimilarityByGradientV2
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByGradientV2 {

private:

public:

  static ttkSimilarityByGradientV2 *New();
  vtkTypeMacro(ttkSimilarityByGradientV2, ttkSimilarityAlgorithm);

protected:
  ttkSimilarityByGradientV2();
  ~ttkSimilarityByGradientV2();

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;
};