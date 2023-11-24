/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::SimilarityById
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %SimilarityById class that computes for
/// each vertex of a triangulation the average scalar value of itself and its
/// direct neighbors.
///
/// \b Related \b publication: \n
/// 'SimilarityById'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

#include <Debug.h>

// std includes
#include <set>

namespace ttk {

  class SimilarityById : virtual public Debug {

  public:
    SimilarityById() {
      this->setDebugMsgPrefix("SimilarityById");
    };
    ~SimilarityById(){};

    template <typename DT>
    int computeUniqueIds(DT *uIds0,
                         DT *uIds1,
                         int &nUIds0,
                         int &nUIds1,
                         const DT *ids0,
                         const DT *ids1,
                         const int nPoints0,
                         const int nPoints1) const {

      ttk::Timer timer;

      const std::string msg = "Computing number of unique feature ids";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      std::set<DT> uniqueIds0;
      std::set<DT> uniqueIds1;

      for(int i = 0; i < nPoints0; i++) {
        uniqueIds0.insert(ids0[i]);
      }

      for(int i = 0; i < nPoints1; i++) {
        uniqueIds1.insert(ids1[i]);
      }

      // Set num values
      nUIds0 = uniqueIds0.size();
      nUIds1 = uniqueIds1.size();

      int index = 0;
      for(auto i = uniqueIds0.begin(); i != uniqueIds0.end(); i++) {
        uIds0[index] = *i;
        index++;
      }
      index = 0;
      for(auto i = uniqueIds1.begin(); i != uniqueIds1.end(); i++) {
        uIds1[index] = *i;
        index++;
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename DT>
    int computeIdentityMatrix(unsigned char *identityMatrix,
                              const DT *uIds0,
                              const DT *uIds1,
                              const int nUIds0,
                              const int nUIds1) const {

      ttk::Timer timer;

      const std::string msg = "Computing Identity Matrix ("
                              + std::to_string(nUIds0) + "x"
                              + std::to_string(nUIds1) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < nUIds0; i++) {
        for(int j = 0; j < nUIds1; j++) {

          const int id0 = uIds0[i];
          const int id1 = uIds1[j];

          identityMatrix[j * nUIds0 + i] = (id0 == id1) ? 1 : 0;
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }
  };
} // namespace ttk
