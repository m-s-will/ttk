/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::SimilarityByDistance
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %SimilarityByDistance class that computes for
/// each vertex of a triangulation the average scalar value of itself and its
/// direct neighbors.
///
/// \b Related \b publication: \n
/// 'SimilarityByDistance'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

#include <Debug.h>

namespace ttk {

  class SimilarityByDistance : virtual public Debug {

  public:
    SimilarityByDistance() {
      this->setDebugMsgPrefix("SimilarityByDistance");
    };
    ~SimilarityByDistance(){};

    template <typename DT>
    int computeDistanceMatrix(DT *distanceMatrix,
                              const DT *coords0,
                              const DT *coords1,
                              const int nPoints0,
                              const int nPoints1) const {

      ttk::Timer timer;

      const std::string msg = "Computing Distance Matrix ("
                              + std::to_string(nPoints0) + "x"
                              + std::to_string(nPoints1) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < nPoints0; i++) {
        for(int j = 0; j < nPoints1; j++) {

          const int i3 = i * 3;
          const int j3 = j * 3;

          const DT dx = coords0[i3 + 0] - coords1[j3 + 0];
          const DT dy = coords0[i3 + 1] - coords1[j3 + 1];
          const DT dz = coords0[i3 + 2] - coords1[j3 + 2];

          distanceMatrix[j * nPoints0 + i]
            = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename DT>
    int normalizeDistanceMatrix(DT *distanceMatrix,
                                const int nPoints0,
                                const int nPoints1,
                                const DT maxDistance = -1) const {

      ttk::Timer timer;

      const std::string msg = "Normalizing Distance Matrix ("
                              + std::to_string(nPoints0) + "x"
                              + std::to_string(nPoints1) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const int n = nPoints0 * nPoints1;

      DT maxValue
        = (maxDistance < 0 && n > 0) ? distanceMatrix[0] : maxDistance;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
#endif
      {
        // if maxDistance < 0 then search for max value inside matrix
        if(maxDistance < 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp for reduction(max : maxValue)
#endif
          for(int i = 1; i < n; i++) {
            maxValue = std::max(maxValue, distanceMatrix[i]);
          }
        }

#ifdef TTK_ENABLE_OPENMP
#pragma omp for
#endif
        for(int i = 0; i < n; i++) {
          distanceMatrix[i] = std::max(
            std::min(static_cast<DT>(1) - distanceMatrix[i] / maxValue,
                     static_cast<DT>(1)),
            static_cast<DT>(0));
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }
  };
} // namespace ttk
