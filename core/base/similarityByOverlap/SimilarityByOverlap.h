/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::SimilarityByOverlap
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %SimilarityByOverlap class that computes for
/// each vertex of a triangulation the average scalar value of itself and its
/// direct neighbors.
///
/// \b Related \b publication: \n
/// 'SimilarityByOverlap'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

#include <Debug.h>

#include <unordered_map>

namespace ttk {

  class SimilarityByOverlap : virtual public Debug {

  public:
    SimilarityByOverlap() {
      this->setDebugMsgPrefix("SimilarityByOverlap");
    };
    ~SimilarityByOverlap(){};

    template <typename DT, typename IT>
    int computeIdIndexMap(std::unordered_map<IT, IT> &idIndexMap,
                          const DT *ids,
                          const IT nIds) const {

      IT idIndex = 0;
      for(IT i = 0; i < nIds; i++) {
        auto l = static_cast<const IT>(ids[i]);
        if(l >= 0 && idIndexMap.find(l) == idIndexMap.end())
          idIndexMap.insert({l, idIndex++});
      }

      return 1;
    };

    template <typename DT, typename IT>
    int computeAdjacencyMatrix(
      int *adjacencyMatrix,
      const DT *ids0,
      const DT *ids1,
      const int nVertices,
      const std::unordered_map<IT, IT> &idIndexMap0,
      const std::unordered_map<IT, IT> &idIndexMap1) const {

      ttk::Timer timer;

      const IT nIds0 = idIndexMap0.size();
      const IT nIds1 = idIndexMap1.size();

      if(nIds0 < 1)
        return this->printWrn("Number of first ids smaller than 1.");
      if(nIds1 < 1)
        return this->printWrn("Number of second ids smaller than 1.");

      const std::string msg = "Computing Overlap " + std::to_string(nIds0) + "x"
                              + std::to_string(nIds1);
      this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      const IT nLables = nIds0 * nIds1;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT i = 0; i < nLables; i++) {
        adjacencyMatrix[i] = 0;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT i = 0; i < nVertices; i++) {
        auto l0 = static_cast<const IT>(ids0[i]);
        auto l1 = static_cast<const IT>(ids1[i]);
        if(l0 >= 0 && l1 >= 0) {
#pragma omp atomic update
          adjacencyMatrix[idIndexMap1.at(l1) * nIds0 + idIndexMap0.at(l0)]++;
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }
  };
} // namespace ttk
