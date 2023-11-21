/// \ingroup base
/// \class ttk::MergeSimilarityMatrices
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2023-05-19.
///
///
/// \b Related \b publication: \n
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The MergeSimilarityMatrices class ------.
   */
  class MergeSimilarityMatrices : virtual public Debug {

  public:
    MergeSimilarityMatrices();

    template <class dataType, typename IT>
    int initializeSimilarityMatrix(dataType *matrix, const IT nFeatures0,
                                const IT nFeatures1) const {
      
      for(IT i = 0, j = nFeatures0 * nFeatures1; i < j; i++)
        matrix[i] = 0;

      return 1;
    }

    template <class dataType, typename IT>
    int addSubMatrixToMatrix(dataType *outMatrix,
      const dataType *inSubMatrix,
      const IT *ids0,
      const IT *ids1, 
      const IT nSubFeatures0, 
      const IT nSubFeatures1,
      const IT nFeatures0,
      const std::unordered_map<IT,IT>& columnMap,
      const std::unordered_map<IT,IT>& rowMap) const {

        // Go through each row in matrix
        for(int i = 0; i < nSubFeatures1; i++) {

          // Get id of the feature repped in the row
          IT subRowId = ids1[i];
          IT indexRow = rowMap.at(subRowId);

          // Go through each column
          for (int j = 0; j < nSubFeatures0; j++) {
            // Get id of feature repped in column
            // add to row in map
            IT subColId = ids0[j];
            IT indexCol = columnMap.at(subColId);
            outMatrix[indexRow * nFeatures0 + indexCol] = inSubMatrix[i * nSubFeatures0 + j];
          }
        }


      return 1;
    }

    template <class dataType, typename IT>
    int addSimilarityMatrixData(dataType *outMatrix, 
      const dataType *inMatrix, 
      const IT nFeatures0, 
      const IT nFeatures1,
      const std::unordered_map<IT,IT>& columnMap,
      const std::unordered_map<IT,IT>& rowMap) const {

      return 1;
    }

  }; // MergeSimilarityMatrices class

} // namespace ttk
