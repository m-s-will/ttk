/// \ingroup base
/// \class ttk::SimilarityMatrixTemporalDownsampling
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2022-03-03.
///
/// \b Related \b publication: \n

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  class SimilarityMatrixTemporalDownsampling : virtual public Debug {

  public:
    SimilarityMatrixTemporalDownsampling();

    template <typename DT>
    int multiplyMatrices(unsigned char *resMatrix,
                         const DT *prevMatrix,
                         const DT *curMatrix,
                         const int *prevDims,
                         const int *curDims) {

// Set all values to zero
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < curDims[1]; i++) {
        for(int j = 0; j < prevDims[0]; j++) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif
          resMatrix[i * prevDims[0] + j] = 0;
        }
      }

// multiply matrices, remembering that t-1 is on the columns and t on the rows
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < curDims[1]; i++) { // row
        for(int j = 0; j < prevDims[0]; j++) { // column
          for(int k = 0; k < curDims[0]; k++) { // column

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
            resMatrix[i * prevDims[0] + j]
              = resMatrix[i * prevDims[0] + j]
                + curMatrix[i * curDims[0] + k]
                    * prevMatrix[k * prevDims[0] + j];
          }
        }
      }

      return 1;
    }

  }; // SimilarityMatrixTemporalDownsampling class

} // namespace ttk
