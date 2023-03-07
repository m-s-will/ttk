#pragma once

// VTK Module
#include <ttkMergeTreeModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <MergeTree.h>

class TTKMERGETREE_EXPORT ttkMergeTree : public ttkAlgorithm,
                                         public ttk::mt::MergeTree {
private:
  int Type{0};

public:
  vtkGetMacro(Type, int);
  vtkSetMacro(Type, int);

  static ttkMergeTree *New();
  vtkTypeMacro(ttkMergeTree, ttkAlgorithm);

protected:
  ttkMergeTree();
  ~ttkMergeTree() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};