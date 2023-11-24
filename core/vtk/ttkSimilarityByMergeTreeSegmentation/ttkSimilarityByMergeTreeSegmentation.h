#pragma once

// VTK Module
#include <ttkSimilarityByMergeTreeSegmentationModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityByMergeTreeSegmentation.h>

class TTKSIMILARITYBYMERGETREESEGMENTATION_EXPORT
  ttkSimilarityByMergeTreeSegmentation
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityByMergeTreeSegmentation {

public:
  static ttkSimilarityByMergeTreeSegmentation *New();
  vtkTypeMacro(ttkSimilarityByMergeTreeSegmentation, ttkSimilarityAlgorithm);

protected:
  ttkSimilarityByMergeTreeSegmentation();
  ~ttkSimilarityByMergeTreeSegmentation() override;

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;
};
