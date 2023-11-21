/// \ingroup vtk
/// \class ttkMergeSimilarityMatrices
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2023-05-19.
///
/// \brief TTK VTK-filter that wraps the ttk::MergeSimilarityMatrices module.
///
///
/// \param Input vtkMultiBlockDataSet.
/// \param Output vtkMultiBlockDataSet.
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeSimilarityMatrices/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeSimilarityMatrices
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMergeSimilarityMatricesModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkSimilarityAlgorithm.h>

// TTK Base Includes
#include <MergeSimilarityMatrices.h>

// std includes
#include <set>
#include <unordered_map>
#include <vector>

class TTKMERGESIMILARITYMATRICES_EXPORT ttkMergeSimilarityMatrices
  : public ttkAlgorithm
  ,
    protected ttk::MergeSimilarityMatrices
{
private:

public:
  static ttkMergeSimilarityMatrices *New();
  vtkTypeMacro(ttkMergeSimilarityMatrices, ttkAlgorithm);

protected:
  ttkMergeSimilarityMatrices();
  ~ttkMergeSimilarityMatrices() override = default;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
