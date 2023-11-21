/// \ingroup base
/// \class ttk::SimilarityByGradientV2
/// \author Wito Engelke <wito.engelke@googlemail.com>
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date June 2020.
///
/// TODO
///
/// \b Related \b publication: \n
/// 'SimilarityByGradientV2'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <PathCompression.h>

namespace ttk {
  class SimilarityByGradientV2 : virtual public Debug {

  public:
    SimilarityByGradientV2() {
      this->setDebugMsgPrefix("SimilarityByGradientV2");

      #ifdef TTK_ENABLE_MPI
        //if (ttk::MPIsize_ == 1)
          hasMPISupport_ = true;
      #endif

    };
    ~SimilarityByGradientV2(){};

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    template <typename IT, typename TT>
    int computeMap(std::vector<IT> &matchesVals, std::set<IT> &outsideMatches, TT * triangulation, const IT *segmentation, const IT *criticalPointGlobalIds0, const IT nFeatures0) const {

      for (IT i = 0; i < nFeatures0; i++) {
        // Get the local id of the extrema and add to matches if owner of extrema
        IT localId = triangulation->getVertexLocalId(criticalPointGlobalIds0[i]);
        if (triangulation->getVertexRank(localId) == ttk::MPIrank_) {
          IT match = segmentation[localId];
          matchesVals.push_back(match);
          IT localIdMatch = triangulation->getVertexLocalId(match);

          // If the matching critical point is not owned by the rank, add it to outside matches
          if (localIdMatch < 0){
            outsideMatches.insert(match);
          }
        }

      }

      return 1;
    }

    template <typename IT, typename IF>
    int computeSimilarityMatrixMPI(int *matrix,
                                const std::vector<IT> &matchesVals,
                                const IT *criticalPointVertexIds1,
                                const IT nFeatures0,
                                const IT nFeatures1,
                                const IF indexFunction,
                                const std::string& msg) const {
      ttk::Timer timer;
      this->printMsg(msg, 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      for(int i = 0, j = nFeatures0 * nFeatures1; i < j; i++)
        matrix[i] = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nFeatures0; i++) {
        IT matchIdx = -1;
        if (i < IT(matchesVals.size()))
          matchIdx = matchesVals[i];

        // add matches to matix. Note: it is possible that a matched extrema is not in the feature list (i.e., not tracked)
        for(IT j = 0; j < nFeatures1; j++) {
          if(criticalPointVertexIds1[j] == matchIdx) {
            matrix[indexFunction(i, j, nFeatures0, nFeatures1)] = 1;
            break;
          }
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };


    template <typename IT, typename IF>
    int computeSimilarityMatrix(int *matrix,
                                const IT *segmentation,
                                const IT *criticalPointVertexIds0,
                                const IT *criticalPointVertexIds1,
                                const IT nFeatures0,
                                const IT nFeatures1,
                                const IF indexFunction,
                                const std::string& msg) const {
      ttk::Timer timer;
      this->printMsg(msg, 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      for(int i = 0, j = nFeatures0 * nFeatures1; i < j; i++)
        matrix[i] = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nFeatures0; i++) {
        const IT &seedIdx = criticalPointVertexIds0[i];
        IT matchIdx = -1;
        matchIdx = segmentation[seedIdx];

        // add matches to matix. Note: it is possible that a matched extrema is not in the feature list (i.e., not tracked)
        for(IT j = 0; j < nFeatures1; j++) {
          if(criticalPointVertexIds1[j] == matchIdx) {
            matrix[indexFunction(i, j, nFeatures0, nFeatures1)] = 1;
            break;
          }
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };

  }; // SimilarityByGradientV2 class

} // namespace ttk