#pragma once

#include <BaseClass.h>
#include <iostream>

#include <algorithm>
#include <vector>

namespace ttk {
  template <typename scalarType>
  std::vector<std::tuple<float, ttk::SimplexId, int>> populateVector(const size_t nVerts,
                    const float *const scalars,
                    const ttk::SimplexId *const globalIds,
                    const char *const ghostCells) {
    std::vector<std::tuple<float, ttk::SimplexId, int>> outVector;
    for (int i = 0; i < nVerts; i++){
      std::cout << "Ghostcell: " << std::to_string(ghostCells[i]) << std::endl; 
      if ((int)ghostCells[i] == 0){
        float scalarValue = scalars[i];
        ttk::SimplexId globalId = globalIds[i];
        int localId = i;
        outVector.emplace_back(scalarValue, globalId, localId);
      }
    }
    return outVector;
  }

  template <typename scalarType>
  void sortVerticesDistributed(std::vector<std::tuple<scalarType, int, int>> &values) {
    std::sort(values.begin(), values.end());
  }


  /**
   * @brief Sort vertices according to scalars disambiguated by offsets
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars array of size nVerts, main vertex comparator
   * @param[in] offsets array of size nVerts, disambiguate scalars on plateaux
   * @param[out] order array of size nVerts, computed order of vertices
   * @param[in] nThreads number of parallel threads
   */
  template <typename scalarType, typename idType>
  void sortVertices(const size_t nVerts,
                    const scalarType *const scalars,
                    const idType *const offsets,
                    SimplexId *const order,
                    const int nThreads) {

    // array of pre-sorted vertices
    std::vector<SimplexId> sortedVertices(nVerts);

    TTK_FORCE_USE(nThreads);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < sortedVertices.size(); ++i) {
      sortedVertices[i] = i;
    }

    if(offsets != nullptr) {
      TTK_PSORT(
        nThreads, sortedVertices.begin(), sortedVertices.end(),
        [&](const SimplexId a, const SimplexId b) {
          return (scalars[a] < scalars[b])
                 || (scalars[a] == scalars[b] && offsets[a] < offsets[b]);
        });
    } else {
      TTK_PSORT(nThreads, sortedVertices.begin(), sortedVertices.end(),
                [&](const SimplexId a, const SimplexId b) {
                  return (scalars[a] < scalars[b])
                         || (scalars[a] == scalars[b] && a < b);
                });
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif // TTK_ENABLE_OPENMP
    for(size_t i = 0; i < sortedVertices.size(); ++i) {
      order[sortedVertices[i]] = i;
    }
  }

  /**
   * @brief Precondition an order array to be consumed by the base layer API
   *
   * @param[in] nVerts number of vertices
   * @param[in] scalars pointer to scalar field buffer of size @p nVerts
   * @param[out] order pointer to pre-allocated order buffer of size @p nVerts
   * @param[in] nThreads number of threads to be used
   */
  template <typename scalarType>
  inline void preconditionOrderArray(const size_t nVerts,
                                     const scalarType *const scalars,
                                     SimplexId *const order,
                                     const int nThreads
                                     = ttk::globalThreadNumber_) {
    ttk::sortVertices(
      nVerts, scalars, static_cast<int *>(nullptr), order, nThreads);
  }
} // namespace ttk
