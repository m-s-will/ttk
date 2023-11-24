/// \ingroup vtk
/// \class ttkComparingSimilarityMatrices
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date  2022-03-07.
///
/// \brief TTK VTK-filter that wraps the ttk::ComparingSimilarityMatrices
/// module.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.

#pragma once

// VTK Module
#include <ttkComparingSimilarityMatricesModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <ComparingSimilarityMatrices.h>

class TTKCOMPARINGSIMILARITYMATRICES_EXPORT ttkComparingSimilarityMatrices
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ComparingSimilarityMatrices // and we inherit from the base
                                               // class
{
private:
  // vector of event data
  std::vector<ttk::ComparingSimilarityMatrices::Events> matrixEvents;

public:
  static ttkComparingSimilarityMatrices *New();
  vtkTypeMacro(ttkComparingSimilarityMatrices, ttkAlgorithm);

protected:
  ttkComparingSimilarityMatrices();
  ~ttkComparingSimilarityMatrices() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
