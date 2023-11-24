/// \ingroup base
/// \class ttk::SimilarityByMergeTreeSegmentation
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 08/13/2022
///
/// This module computes the spatial overlap of two merge tree segmentations
/// (MTS) and records the results in a similarity matrix. Both MTSs need to
/// be defined on the same domain, and assign to each vertex of the domain the
/// integer id of the corresponding merge tree edge. An overlap of two segments
/// (merge tree edges) imply the overlap of the contained sub- and superlevel
/// set components. For details see the related publication.
///
/// \b Related \b Publication:
/// "Dynamic Nested Tracking Graphs".
/// Jonas Lukasczyk, Christoph Garth, Gunther H. Weber, Tim Biedert, Ross
/// Maciejewski, Heike Leitte. IEEE Transactions on Visualization and Computer
/// Graphics. 2019.

#pragma once

#include <Debug.h>

namespace ttk {

  class SimilarityByMergeTreeSegmentation : virtual public Debug {

  public:
    SimilarityByMergeTreeSegmentation() {
      this->setDebugMsgPrefix("SimilarityByMergeTreeSegmentation");
    };
    ~SimilarityByMergeTreeSegmentation(){};

    template <typename IT, typename DT>
    int getBaseRepresentative(IT &base,

                              const DT baselevel,
                              const DT *scalars,
                              const IT *next) const {
      IT nextN = next[base];
      while(nextN >= 0 && scalars[nextN] > baselevel) {
        base = nextN;
        nextN = next[base];
      }

      return 1;
    }

    template <typename IT, typename DT>
    int computeSegmentationOverlap(IT *similarityMatrix,

                                   const IT *seg0,
                                   const IT *seg1,
                                   const IT nVertices,
                                   const IT *next0,
                                   const IT *next1,
                                   const DT *scalars0,
                                   const DT *scalars1,
                                   const IT nNodes0,
                                   const IT nNodes1) const {
      ttk::Timer timer;
      const std::string msg
        = "Computing Segmentation Overlap (" + std::to_string(nVertices) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      // clear similarity matrix
      const IT nNodePairs = nNodes0 * nNodes1;
      for(IT i = 0; i < nNodePairs; i++) {
        similarityMatrix[i] = 0;
      }

      // compute overlap of base edges
      std::vector<IT> baseSimilarityMatrix(nNodePairs, 0);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT v = 0; v < nVertices; v++) {
        const auto &n0 = seg0[v];
        const auto &n1 = seg1[v];

        const auto &bl = std::min(scalars0[n0], scalars1[n1]);
        IT base0 = n0;
        IT base1 = n1;
        this->getBaseRepresentative<IT, DT>(base0, bl, scalars0, next0);
        this->getBaseRepresentative<IT, DT>(base1, bl, scalars1, next1);

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
        baseSimilarityMatrix[base1 * nNodes0 + base0]++;
      }

      // propagate base similarities towards roots
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT i = 0; i < nNodes1; i++) {
        const auto offset = nNodes0 * i;

        for(IT j = 0; j < nNodes0; j++) {
          const auto &overlap = baseSimilarityMatrix[offset + j];
          if(overlap < 1)
            continue;

          IT c0 = j;
          IT c1 = i;
          while(true) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif
            similarityMatrix[c1 * nNodes0 + c0] += overlap;

            const auto &cNext0 = next0[c0];
            const auto &cNext1 = next1[c1];

            if(cNext0 < 0 && cNext1 < 0)
              break;

            if(cNext0 < 0) {
              c1 = cNext1;
            } else if(cNext1 < 0) {
              c0 = cNext0;
            } else {
              const auto &cNextScalar0 = scalars0[cNext0];
              const auto &cNextScalar1 = scalars1[cNext1];

              if(cNextScalar0 == cNextScalar1) {
                c0 = cNext0;
                c1 = cNext1;
              } else if(cNextScalar0 < cNextScalar1) {
                c1 = cNext1;
              } else {
                c0 = cNext0;
              }
            }
          };
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }
  }; // class
} // namespace ttk
