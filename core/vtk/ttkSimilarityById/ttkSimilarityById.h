#pragma once

// VTK Module
#include <ttkSimilarityByIdModule.h>

// VTK Includes
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <SimilarityById.h>

class TTKSIMILARITYBYID_EXPORT ttkSimilarityById
  : public ttkSimilarityAlgorithm,
    protected ttk::SimilarityById {

private:
public:
  static ttkSimilarityById *New();
  vtkTypeMacro(ttkSimilarityById, ttkSimilarityAlgorithm);

protected:
  ttkSimilarityById();
  ~ttkSimilarityById();

  int ComputeSimilarityMatrix(vtkImageData *similarityMatrix,
                              vtkDataObject *inputDataObjects0,
                              vtkDataObject *inputDataObjects1) override;
};
