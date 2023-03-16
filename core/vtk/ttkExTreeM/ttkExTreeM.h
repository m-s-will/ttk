#pragma once

// VTK Module
#include <ttkExTreeMModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <ExTreeM.h>

class vtkUnstructuredGrid;

class TTKEXTREEM_EXPORT ttkExTreeM
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ExTreeM // and we inherit from the base class
{
private:
  int Type{0};

public:
  vtkGetMacro(Type, int);
  vtkSetMacro(Type, int);

  static ttkExTreeM *New();
  vtkTypeMacro(ttkExTreeM, ttkAlgorithm);

protected:
  ttkExTreeM();
  ~ttkExTreeM() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class triangulationType = ttk::AbstractTriangulation>
  int getMergeTree(vtkUnstructuredGrid *outputSkeletonArcs,
                   std::vector<ExTreeM::Branch> &mergeTree,
                   const triangulationType *triangulation);
};
