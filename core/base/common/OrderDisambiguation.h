#pragma once

#include <BaseClass.h>
#include <iostream>

#include <algorithm>
#include <vector>

namespace ttk {

  struct value{
    float scalar;
    SimplexId globalId;
    int localId;
    int ordering = 0;

    value(float _scalar, SimplexId _globalId, int _localId) 
        : scalar(_scalar)
        , globalId(_globalId)
        , localId(_localId)
    {}
    
  };

  inline std::vector<value> populateVector(const size_t nVerts,
                    const float *const scalars,
                    const SimplexId *const globalIds,
                    const char *const ghostCells) {
    std::vector<value> outVector;
    for (size_t i = 0; i < nVerts; i++){
      if ((int)ghostCells[i] == 0){
        float scalarValue = scalars[i];
        SimplexId globalId = globalIds[i];
        std::cout << "GlobalId: " << std::to_string(globalId) << std::endl; 
        int localId = i;
        outVector.emplace_back(scalarValue, globalId, localId);
      }
    }
    return outVector;
  }

  inline void sortVerticesDistributed(std::vector<value> &values) {
    std::sort(values.begin(), values.end(), [](value v1, value v2){
      return (v1.scalar < v2.scalar)
                 || (v1.scalar == v2.scalar && v1.globalId < v2.globalId);
    });
  }


  // send the highest burstSize values and decrease the vector by that amount
  // check if there are actually that many elements in the vector
  inline std::vector<value> returnVectorForBurstsize(std::vector<value> &values, size_t burstSize) {
    std::vector<value> outVector;
    if(burstSize > values.size()){
      outVector = {values.begin(), values.end()};
      values.clear(); 
    } else {
      outVector = {values.end() - burstSize, values.end()};
      values.erase(values.end() - burstSize, values.end());
    }

    return outVector;
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
