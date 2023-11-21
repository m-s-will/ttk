#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <functional>

#include <algorithm>

namespace ttk {
  class FeatureCorrespondences : virtual public Debug {
  public:
    FeatureCorrespondences() {
      this->setDebugMsgPrefix("FeatureCorrespondences");
    }

    template <typename DT, typename IDX>
    int sortAndReduceCorrespondencesPerFeature_(DT *oMatrix,
                                                const DT *iMatrix,
                                                const int iN,
                                                const int jN,
                                                const int nCandidates,
                                                const bool ascending,
                                                const IDX matrixIdx) const {

      using VT = std::pair<DT, int>;
      std::vector<VT> sortedValues;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) private(sortedValues)
#endif
      for(int i = 0; i < iN; i++) {
        // initialize vector
        sortedValues.resize(jN);

        // get values
        for(int j = 0; j < jN; j++) {
          const int idx = matrixIdx(i, j);
          auto &p = sortedValues[j];
          p.first = iMatrix[idx];
          p.second = idx;
        }

        // serial sort
        if(ascending)
          std::sort(sortedValues.begin(), sortedValues.end(), std::less<VT>());
        else
          std::sort(
            sortedValues.begin(), sortedValues.end(), std::greater<VT>());

        // copy maximum candidates
        const int cN = std::min(nCandidates, jN);
        for(int c = 0; c < cN; c++) {
          const auto &p = sortedValues[c];
          oMatrix[p.second] = p.first;
        }
      }

      return 1;
    }

    template <typename DT>
    int sortAndReduceCorrespondencesPerFeature(DT *oMatrix,
                                               const DT *iMatrix,
                                               const int dimX,
                                               const int dimY,
                                               const int nCandidates,
                                               const bool ascending) const {
      ttk::Timer timer;

      const std::string msg = "Computing " + std::to_string(nCandidates)
                              + (ascending ? " smallest" : " largest")
                              + " correspondences per feature.";
      this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      const int n = dimX * dimY;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < n; i++) {
        oMatrix[i] = 0;
      }
      this->printMsg(msg, 0.2, timer.getElapsedTime(), this->threadNumber_,
                     debug::LineMode::REPLACE);

      // forward optimization
      this->sortAndReduceCorrespondencesPerFeature_(
        oMatrix, iMatrix, dimX, dimY, nCandidates, ascending,
        [=](const int &i, const int &j) { return j * dimX + i; });
      this->printMsg(msg, 0.6, timer.getElapsedTime(), this->threadNumber_,
                     debug::LineMode::REPLACE);

      // backward optimization
      this->sortAndReduceCorrespondencesPerFeature_(
        oMatrix, iMatrix, dimY, dimX, nCandidates, ascending,
        [=](const int &i, const int &j) { return i * dimX + j; });
      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    template <typename DT, typename F>
    int mapEachElement(DT *oMatrix,
                       const DT *iMatrix,
                       const int dimX,
                       const int dimY,
                       const F maskFunction) const {
      ttk::Timer timer;

      const std::string msg = "Computing Mask";
      this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      const int n = dimX * dimY;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < n; i++) {
        oMatrix[i] = maskFunction(iMatrix[i]);
      }
      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }
  };
} // namespace ttk