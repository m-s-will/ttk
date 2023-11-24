/// \ingroup vtk
/// \class ttkSimilarityMatrixTemporalDownsampling
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2022-03-03.
///
/// \brief TTK VTK-filter that wraps the
/// ttk::SimilarityMatrixTemporalDownsampling module.
///
/// This VTK filter uses the ttk::SimilarityMatrixTemporalDownsampling module to
/// ???
///
/// \param Input vtkMultiBlockDataSet.
/// \param Output vtkMultiBlockDataSet.
///

#pragma once

// VTK Module
#include <ttkSimilarityMatrixTemporalDownsamplingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <SimilarityMatrixTemporalDownsampling.h>

class TTKSIMILARITYMATRIXTEMPORALDOWNSAMPLING_EXPORT
  ttkSimilarityMatrixTemporalDownsampling
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::SimilarityMatrixTemporalDownsampling // and we inherit from
                                                        // the base class
{
private:
  int SamplingInterval{0};

public:
  vtkSetMacro(SamplingInterval, int);
  vtkGetMacro(SamplingInterval, int);

  static ttkSimilarityMatrixTemporalDownsampling *New();
  vtkTypeMacro(ttkSimilarityMatrixTemporalDownsampling, ttkAlgorithm);

protected:
  ttkSimilarityMatrixTemporalDownsampling();
  ~ttkSimilarityMatrixTemporalDownsampling() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
