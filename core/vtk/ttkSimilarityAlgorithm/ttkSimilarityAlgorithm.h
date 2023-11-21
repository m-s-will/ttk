/// \ingroup vtk
/// \class ttkSimilarityAlgorithm
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TODO
///
/// TODO
///
/// \sa ttk::SimilarityAlgorithm
/// \sa ttkAlgorithm

#pragma once

#include <unordered_map>

// VTK Module
#include <ttkSimilarityAlgorithmModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkImageData;
class vtkMultiBlockDataSet;
class vtkDataArray;
class vtkFieldData;

class TTKSIMILARITYALGORITHM_EXPORT ttkSimilarityAlgorithm
  : public ttkAlgorithm {

private:
  vtkSmartPointer<vtkMultiBlockDataSet> PreviousInputs;

public:
  static ttkSimilarityAlgorithm *New();
  vtkTypeMacro(ttkSimilarityAlgorithm, ttkAlgorithm);

  static std::string GetIdArrayName(vtkFieldData *fieldData);
  static int GetIndexIdMaps(vtkDataArray *&indexIdMapP,
                            vtkDataArray *&indexIdMapC,
                            vtkFieldData *fieldData);
  static int AddIndexIdMap(vtkImageData *similarityMatrix,
                           vtkDataArray *indexIdMap,
                           const bool isMapForCurrentTimestep);
  static int AddIndexIdMaps(vtkImageData *similarityMatrix,
                            vtkDataArray *indexIdMapR,
                            vtkDataArray *indexIdMapC);

  static int AddIndexIdMaps(
    vtkImageData *similarityMatrix,
    const std::unordered_map<ttk::SimplexId, ttk::SimplexId> &idIndexMapP,
    const std::unordered_map<ttk::SimplexId, ttk::SimplexId> &idIndexMapC,
    const std::string &idArrayName);

  static int
    BuildIdIndexMap(std::unordered_map<ttk::SimplexId, ttk::SimplexId> &,
                    const vtkDataArray *indexIdMap);

protected:
  ttkSimilarityAlgorithm();
  ~ttkSimilarityAlgorithm();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  virtual int
    ComputeSimilarityMatrix(vtkImageData *vtkNotUsed(similarityMatrix),
                            vtkDataObject *vtkNotUsed(inputDataObjects0),
                            vtkDataObject *vtkNotUsed(inputDataObjects1)) {
    return 0;
  };
};