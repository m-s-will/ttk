/// \ingroup base
/// \class ttk::PathCompression
/// \author Robin G. C. Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK processing package for the computation of Morse-Smale
/// segmentations using Path Compression.
///
/// Given an input order field, this class computes its ascending and descending
/// segmentation by assigning every vertex to its minimum or maximum in gradient
/// or inverse gradient direction. For convienience a hash (no hash collision
/// detection) of both segmentations can be created to represent the Morse-Smale
/// segmentation.
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttkPathCompression.cpp %for a usage example
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleSegmentation_at/">Morse-Smale
///   segmentation example</a> \n

#pragma once

// base code includes
#include <Triangulation.h>

using ttk::SimplexId;

namespace ttk {
  namespace pcp {
#ifdef TTK_ENABLE_64BIT_IDS
    /**
     * @brief Get a hash value from two keys
     *
     * @param a First Hash key
     * @param b Second hash key
     *
     * @return Hash value
     */
    constexpr unsigned long long int getHash(const unsigned long long int a,
                                             const unsigned long long int b) {
      return (a * b + (a * a) + (b * b) + (a * a * a) * (b * b * b))
             % ULONG_LONG_MAX;
    }
#else
    /**
     * @brief Get a hash value from two keys
     *
     * @param a First Hash key
     * @param b Second hash key
     *
     * @return Hash value
     */
    constexpr unsigned int getHash(const unsigned int a, const unsigned int b) {
      return (a * b + (a * a) + (b * b) + (a * a * a) * (b * b * b)) % UINT_MAX;
    }
#endif
  } // namespace pcp
  class PathCompression : public virtual Debug {
  public:
    PathCompression();

    struct globalIdOwner {
      ttk::SimplexId globalId;
      int ownerRank;
      ttk::SimplexId ascendingTarget = -1;
      ttk::SimplexId descendingTarget = -1;

      globalIdOwner(ttk::SimplexId _globalId = -1, int _owner = -1)
        : globalId(_globalId), ownerRank(_owner) {
      }

      bool operator<(const globalIdOwner &other) const {
        return globalId < other.globalId;
      }

      bool operator==(const globalIdOwner &other) const {
        return globalId == other.globalId;
      }
    };


    /** @brief Pointers to pre-allocated segmentation point data arrays */
    struct OutputSegmentation {
      SimplexId *ascending_;
      SimplexId *descending_;
      SimplexId *morseSmale_;
    };

    /**
     * Compute necessary triangulation information
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      data->preconditionVertexNeighbors();
#if TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI())
        data->preconditionDistributedVertices();
#endif
    }

    /**
     * @brief Main function for computing the Morse-Smale complex.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] outSegmentation segmentations as a struct
     * @param[in] orderArray order field
     * @param[in] triangulation triangulation
     *
     * @pre PathCompression::preconditionTriangulation must be
     * called prior to this.
     *
     * @return 0 on success
     */
    template <typename triangulationType>
    inline int execute(OutputSegmentation &outSegmentation,
                       const SimplexId *const orderArray,
                       const triangulationType &triangulation);

    /**
     * @brief Compute the ascending and descending segmentation in one run
     *
     * This function computes the ascending and descending segmentation on the
     * order field. First, the ascending and descending segmentation is set to
     * the largest and smallest neighbor of each vertex. Then, using path
     * compression, the vertices are assigned to their minimum/maximum in
     * positive/negative gradient direction.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] ascSegmentation ascending segmentation
     * @param[out] dscSegmentation descending segmentation
     * @param[in] orderArray order array
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename triangulationType>
    int computePathCompression(SimplexId *const ascSegmentation,
                               SimplexId *const dscSegmentation,
                               const SimplexId *const orderArray,
                               const triangulationType &triangulation) const;

    /**
     * @brief Compute the ascending or descending segmentation
     *
     * This function computes the ascending or descending segmentation on the
     * order field. First, the ascending or descending segmentation is set to
     * the largest or smallest neighbor of each vertex. Then, using path
     * compression, the vertices are assigned to their minimum/maximum in
     * positive/negative gradient direction.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] segmentation segmentation
     * @param[in] computeAscending compute the ascending or descending
     * segmentation
     * @param[in] orderArray order array
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename triangulationType>
    int computePathCompressionSingle(
      SimplexId *const segmentation,
      const bool computeAscending,
      const SimplexId *const orderArray,
      const triangulationType &triangulation) const;

    /**
     * Enable/Disable computation of the geometrical embedding of the
     * manifolds of the critical points.
     */
    inline void setComputeSegmentation(const bool doAscending,
                                       const bool doDescending,
                                       const bool doMorseSmale) {
      this->ComputeAscendingSegmentation = doAscending;
      this->ComputeDescendingSegmentation = doDescending;
      this->ComputeMSSegmentationHash = doMorseSmale;
    }
    /**
     * @brief Computes a MS segmentation hash
     *
     * Computes a hash from the ascending and descending segmentation as keys.
     * The function does not check for hash conflicts.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] morseSmaleSegmentation
     * @param[out] ascSegmentation ascending segmentation
     * @param[out] dscSegmentation descending segmentation
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename triangulationType>
    int computeMSHash(SimplexId *const morseSmaleSegmentation,
                      const SimplexId *const ascSegmentation,
                      const SimplexId *const dscSegmentation,
                      const triangulationType &triangulation) const;

  protected:
    // Compute ascending segmentation?
    bool ComputeAscendingSegmentation{true};

    // Compute descending segmentation?
    bool ComputeDescendingSegmentation{true};

    // Compute Morse-Smale segmentation hash?
    bool ComputeMSSegmentationHash{true};
  };
} // namespace ttk

template <typename triangulationType>
int ttk::PathCompression::execute(OutputSegmentation &outSegmentation,
                                  const SimplexId *const orderArray,
                                  const triangulationType &triangulation) {
  if(orderArray == nullptr)
    return this->printErr("Input offset field pointer is null.");

  Timer t;

  this->printMsg("Start computing segmentations", 0.0, t.getElapsedTime(),
                 this->threadNumber_);

  if((ComputeAscendingSegmentation && ComputeDescendingSegmentation)
     || ComputeMSSegmentationHash) {
    computePathCompression(outSegmentation.ascending_,
                           outSegmentation.descending_, orderArray,
                           triangulation);
  } else if(ComputeAscendingSegmentation) {
    computePathCompressionSingle(
      outSegmentation.ascending_, true, orderArray, triangulation);
  } else if(ComputeDescendingSegmentation) {
    computePathCompressionSingle(
      outSegmentation.descending_, false, orderArray, triangulation);
  }

  if(ComputeMSSegmentationHash) {
    computeMSHash(outSegmentation.morseSmale_, outSegmentation.ascending_,
                  outSegmentation.descending_, triangulation);
  }

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompression(
  SimplexId *const ascSegmentation,
  SimplexId *const dscSegmentation,
  const SimplexId *const orderArray,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;
  std::vector<globalIdOwner> foreignVertices;

  TTK_FORCE_USE(foreignVertices);


  const SimplexId nVertices = triangulation.getNumberOfVertices();
  std::vector<SimplexId> lActiveVertices;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_) private(lActiveVertices)
  {
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  lActiveVertices.reserve(nVertices);
#endif // TTK_ENABLE_OPENMP
    // find the largest and smallest neighbor for each vertex
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId{0};
      SimplexId const numNeighbors = triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &dmi = dscSegmentation[i];
      dmi = i;

      bool hasSmallerNeighbor = false;
      SimplexId &ami = ascSegmentation[i];
      ami = i;

#ifdef TTK_ENABLE_MPI
        if(ttk::isRunningWithMPI()) {
          if(triangulation.getVertexRank(i) == ttk::MPIrank_) {
            // check all neighbors
            for(SimplexId n = 0; n < numNeighbors; n++) {
              triangulation.getVertexNeighbor(i, n, neighborId);

              if(orderArray[neighborId] < orderArray[ami]) {
                ami = neighborId;
                hasSmallerNeighbor = true;
              } else if(orderArray[neighborId] > orderArray[dmi]) {
                dmi = neighborId;
                hasLargerNeighbor = true;
              }
            }
          } else {
            globalIdOwner GIO = {triangulation.getVertexGlobalId(i),
                                 triangulation.getVertexRank(i)};
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
            foreignVertices.push_back(GIO);
          }
        } else {
          // check all neighbors
          for(SimplexId n = 0; n < numNeighbors; n++) {
            triangulation.getVertexNeighbor(i, n, neighborId);

            if(orderArray[neighborId] < orderArray[ami]) {
              ami = neighborId;
              hasSmallerNeighbor = true;
            } else if(orderArray[neighborId] > orderArray[dmi]) {
              dmi = neighborId;
              hasLargerNeighbor = true;
            }
          }
        }
#else
      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(orderArray[neighborId] < orderArray[ami]) {
          ami = neighborId;
          hasSmallerNeighbor = true;
        } else if(orderArray[neighborId] > orderArray[dmi]) {
          dmi = neighborId;
          hasLargerNeighbor = true;
        }
      }
#endif
      if(hasLargerNeighbor || hasSmallerNeighbor) {
        lActiveVertices.push_back(i);
      }
    }

    size_t lnActiveVertices = lActiveVertices.size();
    size_t currentIndex = 0;

    // compress paths until no changes occur
    while(lnActiveVertices > 0) {
      for(size_t i = 0; i < lnActiveVertices; i++) {
        SimplexId const &v = lActiveVertices[i];
        SimplexId &vDsc = dscSegmentation[v];
        SimplexId &vAsc = ascSegmentation[v];

// compress paths
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif // TTK_ENABLE_OPENMP
        vDsc = dscSegmentation[vDsc];

#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif // TTK_ENABLE_OPENMP
        vAsc = ascSegmentation[vAsc];

        // check if fully compressed
        if(vDsc != dscSegmentation[vDsc] || vAsc != ascSegmentation[vAsc]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP
#ifdef TTK_ENABLE_MPI
  // now we need to transform local ids into global ids to correctly work
  // over all ranks
  if(ttk::isRunningWithMPI()) {
    this->printMsg("using mpi");
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(ttk::SimplexId i = 0; i < nVertices; i++) {
        dscManifold[i] = triangulation.getVertexGlobalId(dscManifold[i]);
        ascManifold[i] = triangulation.getVertexGlobalId(ascManifold[i]);
    }
    this->printMsg(
      "Rank " + std::to_string(ttk::MPIrank_) + ", finished own values", 1,
      localTimer.getElapsedTime());

    // now we need to request the values we still need from other ranks
    // R0 builds up out transferance map over ranks
    MPI_Barrier(ttk::MPIcomm_);
    // MPI_Reduce to get the size of everything we need
    int localSize = foreignVertices.size();
    this->printMsg("Localsize: " + std::to_string(localSize));
    int totalSize;
    MPI_Reduce(&localSize, &totalSize, 1, MPI_INT, MPI_SUM, 0, ttk::MPIcomm_);
    MPI_Bcast(&totalSize, 1, MPI_INT, 0, ttk::MPIcomm_);
    this->printMsg("Rank " + std::to_string(ttk::MPIrank_)
                   + " got the totalsize " + std::to_string(totalSize));
    std::vector<globalIdOwner> edgesWithTargets(totalSize);
    std::vector<int> sizes(ttk::MPIsize_);
    std::vector<int> displacements(ttk::MPIsize_);
    std::vector<globalIdOwner> sendValues;
    int receivedSize;
    std::vector<globalIdOwner> receivedIds;
    std::vector<globalIdOwner> allValuesFromRanks;
    if(ttk::MPIrank_ == 0) {
      std::vector<globalIdOwner> edges(totalSize);
      // construct the set with the globalids not owned by R0
      // first we use MPI_Gather to get the size of each rank to populate
      // the displacements
      MPI_Gather(
        &localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, ttk::MPIcomm_);
      for(int i = 0; i < ttk::MPIsize_; i++) {
        sizes[i] = sizes[i] * sizeof(globalIdOwner);
      }
      displacements[0] = 0;
      // build our displacements
      for(int i = 1; i < ttk::MPIsize_; i++) {
        displacements[i] = displacements[i - 1] + sizes[i - 1];
      }

      // then we use MPI_Gatherv to get the data that is needed
      // R0 gets all the needed global ids and their original owners,
      // which rank actually needs those gids doesn't matter, because in
      // the end we send the finished edges back
      MPI_Gatherv(foreignVertices.data(),
                  foreignVertices.size() * sizeof(globalIdOwner), MPI_CHAR,
                  edges.data(), sizes.data(), displacements.data(), MPI_CHAR, 0,
                  ttk::MPIcomm_);
      this->printMsg("R0 received " + std::to_string(edges.size())
                     + " ids which are needed");

      // we have all the gids which are needed by _some_ rank and the
      // owner of them, now we have to request from the owners to which
      // vertices these gids are pointing to build our map
      std::vector<std::vector<globalIdOwner>> valuesFromRanks;
      valuesFromRanks.resize(ttk::MPIsize_);
      for(globalIdOwner currentId : edges) {
        valuesFromRanks[currentId.ownerRank].push_back(currentId);
      }
      this->printMsg("R0 reordered Ids");

      for(int i = 0; i < ttk::MPIsize_; i++) {
        sizes[i] = valuesFromRanks[i].size() * sizeof(globalIdOwner);
      }
      displacements[0] = 0;
      // build our displacements
      for(int i = 1; i < ttk::MPIsize_; i++) {
        displacements[i] = displacements[i - 1] + sizes[i - 1];
      }
      // we turn our vector of vectors into a 1D vector to send it via
      // scatter
      for(auto &&v : valuesFromRanks) {
        allValuesFromRanks.insert(allValuesFromRanks.end(), v.begin(), v.end());
      }
    } else { // the other ranks
      // first send the number of ids this rank needs to the root, then
      // the ids themselves the NULL attributes are only relevant for the
      // root rank
      MPI_Gather(&localSize, 1, MPI_INT, NULL, 0, MPI_INT, 0, ttk::MPIcomm_);
      MPI_Gatherv(foreignVertices.data(),
                  foreignVertices.size() * sizeof(globalIdOwner), MPI_CHAR,
                  NULL, NULL, NULL, MPI_CHAR, 0, ttk::MPIcomm_);
    }

    // we need to receive the results to which the gids are pointing from
    // the ranks and build our map
    // we broadcast sizes and displacements and gather all the results
    // with allgatherv
    MPI_Bcast(sizes.data(), ttk::MPIsize_, MPI_INT, 0, ttk::MPIcomm_);
    MPI_Bcast(displacements.data(), ttk::MPIsize_, MPI_INT, 0, ttk::MPIcomm_);
    this->printMsg("R" + std::to_string(ttk::MPIrank_)
                   + " received sizes and displacements from R0");

    receivedSize = sizes[ttk::MPIrank_] / sizeof(globalIdOwner);
    this->printMsg("R" + std::to_string(ttk::MPIrank_) + " owns "
                   + std::to_string(receivedSize) + " ids");

    // and then the actual gids
    receivedIds.resize(receivedSize);
    MPI_Scatterv(allValuesFromRanks.data(), sizes.data(), displacements.data(),
                 MPI_CHAR, receivedIds.data(),
                 receivedSize * sizeof(globalIdOwner), MPI_CHAR, 0,
                 ttk::MPIcomm_);
    this->printMsg("R" + std::to_string(ttk::MPIrank_) + " received ids");

    // now we need to find to where the gids point and send the values
    // back to R0
    sendValues.resize(receivedSize);
    for(ttk::SimplexId i = 0; i < receivedSize; i++) {
      globalIdOwner currentVal = receivedIds[i];
      ttk::SimplexId lId = triangulation.getVertexLocalId(currentVal.globalId);
      currentVal.ascendingTarget = ascManifold[lId];
      currentVal.descendingTarget = dscManifold[lId];
      sendValues[i] = currentVal;
    }
    this->printMsg("R" + std::to_string(ttk::MPIrank_)
                   + " is done with their values");

    MPI_Allgatherv(sendValues.data(), sendValues.size() * sizeof(globalIdOwner),
                   MPI_CHAR, edgesWithTargets.data(), sizes.data(),
                   displacements.data(), MPI_CHAR, ttk::MPIcomm_);
    this->printMsg("R" + std::to_string(ttk::MPIrank_) + " got the results");

    // now each rank has a vector consisting of gIds, the ranks to which
    // they belong and the ascending / descending target we have all the
    // information on R0 which we need to resolve any manifolds stretching
    // over multiple ranks
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToAscendingMap;
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToDescendingMap;

    for(size_t i = 0; i < edgesWithTargets.size(); i++) {
      globalIdOwner currentVal = edgesWithTargets[i];
      gIdToAscendingMap.insert(
        std::make_pair(currentVal.globalId, currentVal.ascendingTarget));
      gIdToDescendingMap.insert(
        std::make_pair(currentVal.globalId, currentVal.descendingTarget));
    }

    // now we need to check for graphs in the map and iteratively compress
    // them
    bool changed = true;
    while(changed) {
      changed = false;
      for(auto &it : gIdToAscendingMap) {
        if(gIdToAscendingMap.count(it.second) && (it.first != it.second)
          && (gIdToAscendingMap[it.first] != gIdToAscendingMap[it.second])) {
          gIdToAscendingMap[it.first] = gIdToAscendingMap[it.second];
          changed = true;
        }
      }

      for(auto &it : gIdToDescendingMap) {
        if(gIdToDescendingMap.count(it.second) && (it.first != it.second)
          && (gIdToDescendingMap[it.first] != gIdToDescendingMap[it.second])) {
          gIdToDescendingMap[it.first] = gIdToDescendingMap[it.second];
          changed = true;
        }
      }
    }
    // now each rank simply needs to walk over their vertices and replace
    // ones aiming to ghostcells with the correct ones from the map
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(ttk::SimplexId i = 0; i < nVertices; i++) {
      ttk::SimplexId gid = triangulation.getVertexGlobalId(i);
      ttk::SimplexId descVal = dscManifold[i];
      ttk::SimplexId ascVal = ascManifold[i];
      if(gIdToDescendingMap.count(descVal)) {
        dscManifold[i] = gIdToDescendingMap[descVal];
      } else if(gIdToDescendingMap.count(gid)) {
        dscManifold[i] = gIdToDescendingMap[gid];
      }
      if(gIdToAscendingMap.count(ascVal)) {
        ascManifold[i] = gIdToAscendingMap[ascVal];
      } else if(gIdToAscendingMap.count(gid)) {
        ascManifold[i] = gIdToAscendingMap[gid];
      }
    }
  } else {
    this->printMsg("not using mpi");
  }
#endif // TTK_ENABLE_MPI


  this->printMsg("Asc. and Desc. segmentation computed", 1.0,
                 localTimer.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompressionSingle(
  SimplexId *const segmentation,
  const bool computeAscending,
  const SimplexId *const orderArray,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;
  std::vector<globalIdOwner> foreignVertices;

  TTK_FORCE_USE(foreignVertices);

  const SimplexId nVertices = triangulation.getNumberOfVertices();
  std::vector<SimplexId> lActiveVertices;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_)
  {
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

#pragma omp for schedule(static)
#else // TTK_ENABLE_OPENMP
  lActiveVertices.reserve(nVertices);
#endif // TTK_ENABLE_OPENMP
    // find the largest neighbor for each vertex
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId{0};
      SimplexId const numNeighbors = triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &mi = segmentation[i];
      mi = i;

      #ifdef TTK_ENABLE_MPI
        if(ttk::isRunningWithMPI()) {
          if(triangulation.getVertexRank(i) == ttk::MPIrank_) {
            // check all neighbors
            for(SimplexId n = 0; n < numNeighbors; n++) {
              triangulation.getVertexNeighbor(i, n, neighborId);

              if(computeAscending) {
                if(orderArray[neighborId] < orderArray[mi]) {
                  mi = neighborId;
                  hasLargerNeighbor = true;
                }
              } else {
                if(orderArray[neighborId] > orderArray[mi]) {
                  mi = neighborId;
                  hasLargerNeighbor = true;
                }
              }
            }
          } else {
            globalIdOwner GIO = {triangulation.getVertexGlobalId(i),
                                 triangulation.getVertexRank(i)};
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
            foreignVertices.push_back(GIO);
          }
        } else {
          // check all neighbors
          for(SimplexId n = 0; n < numNeighbors; n++) {
            triangulation.getVertexNeighbor(i, n, neighborId);

            if(computeAscending) {
              if(orderArray[neighborId] < orderArray[mi]) {
                mi = neighborId;
                hasLargerNeighbor = true;
              }
            } else {
              if(orderArray[neighborId] > orderArray[mi]) {
                mi = neighborId;
                hasLargerNeighbor = true;
              }
            }
          }
        }
#else
      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(computeAscending) {
          if(orderArray[neighborId] < orderArray[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        } else {
          if(orderArray[neighborId] > orderArray[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        }
      }
#endif

      if(hasLargerNeighbor) {
        lActiveVertices.push_back(i);
      }
    }

    size_t lnActiveVertices = lActiveVertices.size();
    size_t currentIndex = 0;

    // compress paths until no changes occur
    while(lnActiveVertices > 0) {
      for(size_t i = 0; i < lnActiveVertices; i++) {
        SimplexId const &v = lActiveVertices[i];
        SimplexId &vMan = segmentation[v];

// compress paths
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic read
#endif // TTK_ENABLE_OPENMP
        vMan = segmentation[vMan];

        // check if fully compressed
        if(vMan != segmentation[vMan]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
#ifdef TTK_ENABLE_OPENMP
  }
#endif // TTK_ENABLE_OPENMP

#ifdef TTK_ENABLE_MPI
// now we need to transform local ids into global ids to correctly work
// over all ranks
if(ttk::isRunningWithMPI()) {
  // this->printMsg("using mpi for global share");
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < nVertices; i++) {
    manifold[i] = triangulation.getVertexGlobalId(manifold[i]);
  }
  this->printMsg(
    "Rank " + std::to_string(ttk::MPIrank_) + ", finished own values", 1,
    localTimer.getElapsedTime());

  // now we need to request the values we still need from other ranks
  // R0 builds up out transferance map over ranks
  MPI_Barrier(ttk::MPIcomm_);
  // MPI_Reduce to get the size of everything we need
  int localSize = foreignVertices.size();
  this->printMsg("Localsize: " + std::to_string(localSize));
  int totalSize;
  MPI_Reduce(&localSize, &totalSize, 1, MPI_INT, MPI_SUM, 0, ttk::MPIcomm_);
  MPI_Bcast(&totalSize, 1, MPI_INT, 0, ttk::MPIcomm_);
  this->printMsg("Rank " + std::to_string(ttk::MPIrank_) + " got the totalsize "
                 + std::to_string(totalSize));
  std::vector<globalIdOwner> edgesWithTargets(totalSize);
  std::vector<int> sizes(ttk::MPIsize_);
  std::vector<int> displacements(ttk::MPIsize_);
  std::vector<globalIdOwner> sendValues;
  int receivedSize;
  std::vector<globalIdOwner> receivedIds;
  std::vector<globalIdOwner> allValuesFromRanks;
  if(ttk::MPIrank_ == 0) {
    std::vector<globalIdOwner> edges(totalSize);
    // construct the set with the globalids not owned by R0

    // first we use MPI_Gather to get the size of each rank to populate
    // the displacements
    MPI_Gather(
      &localSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, ttk::MPIcomm_);
    for(int i = 0; i < ttk::MPIsize_; i++) {
      sizes[i] = sizes[i] * sizeof(globalIdOwner);
    }
    displacements[0] = 0;
    // build our displacements
    for(int i = 1; i < ttk::MPIsize_; i++) {
      displacements[i] = displacements[i - 1] + sizes[i - 1];
    }

    // then we use MPI_Gatherv to get the data that is needed
    // R0 gets all the needed global ids and their original owners,
    // which rank actually needs those gids doesn't matter, because in
    // the end we send the finished edges back
    MPI_Gatherv(foreignVertices.data(),
                foreignVertices.size() * sizeof(globalIdOwner), MPI_CHAR,
                edges.data(), sizes.data(), displacements.data(), MPI_CHAR, 0,
                ttk::MPIcomm_);
    this->printMsg("R0 received " + std::to_string(edges.size())
                   + " ids which are needed");

    // we have all the gids which are needed by _some_ rank and the
    // owner of them, now we have to request from the owners to which
    // vertices these gids are pointing to build our map
    std::vector<std::vector<globalIdOwner>> valuesFromRanks;
    valuesFromRanks.resize(ttk::MPIsize_);
    for(globalIdOwner currentId : edges) {
      valuesFromRanks[currentId.ownerRank].push_back(currentId);
    }
    this->printMsg("R0 reordered Ids");

    for(int i = 0; i < ttk::MPIsize_; i++) {
      sizes[i] = valuesFromRanks[i].size() * sizeof(globalIdOwner);
    }
    displacements[0] = 0;
    // build our displacements
    for(int i = 1; i < ttk::MPIsize_; i++) {
      displacements[i] = displacements[i - 1] + sizes[i - 1];
    }
    // we turn our vector of vectors into a 1D vector to send it via
    // scatter
    for(auto &&v : valuesFromRanks) {
      allValuesFromRanks.insert(allValuesFromRanks.end(), v.begin(), v.end());
    }
  } else { // the other ranks
    // first send the number of ids this rank needs to the root, then
    // the ids themselves the NULL attributes are only relevant for the
    // root rank
    MPI_Gather(&localSize, 1, MPI_INT, NULL, 0, MPI_INT, 0, ttk::MPIcomm_);
    MPI_Gatherv(foreignVertices.data(),
                foreignVertices.size() * sizeof(globalIdOwner), MPI_CHAR, NULL,
                NULL, NULL, MPI_CHAR, 0, ttk::MPIcomm_);
  }

  // we need to receive the results to which the gids are pointing from
  // the ranks and build our map
  // we broadcast sizes and displacements and gather all the results
  // with allgatherv
  MPI_Bcast(sizes.data(), ttk::MPIsize_, MPI_INT, 0, ttk::MPIcomm_);
  MPI_Bcast(displacements.data(), ttk::MPIsize_, MPI_INT, 0, ttk::MPIcomm_);
  this->printMsg("R" + std::to_string(ttk::MPIrank_)
                 + " received sizes and displacements from R0");

  receivedSize = sizes[ttk::MPIrank_] / sizeof(globalIdOwner);
  this->printMsg("R" + std::to_string(ttk::MPIrank_) + " owns "
                 + std::to_string(receivedSize) + " ids");

  // and then the actual gids
  receivedIds.resize(receivedSize);
  MPI_Scatterv(allValuesFromRanks.data(), sizes.data(), displacements.data(),
               MPI_CHAR, receivedIds.data(),
               receivedSize * sizeof(globalIdOwner), MPI_CHAR, 0,
               ttk::MPIcomm_);
  this->printMsg("R" + std::to_string(ttk::MPIrank_) + " received ids");

  // now we need to find to where the gids point and send the values
  // back to R0
  sendValues.resize(receivedSize);
  for(ttk::SimplexId i = 0; i < receivedSize; i++) {
    globalIdOwner currentVal = receivedIds[i];
    ttk::SimplexId lId = triangulation.getVertexLocalId(currentVal.globalId);
    if(computeAscending) {
      currentVal.ascendingTarget = manifold[lId];
    } else {
      currentVal.descendingTarget = manifold[lId];
    }
    sendValues[i] = currentVal;
  }
  this->printMsg("R" + std::to_string(ttk::MPIrank_)
                 + " is done with their values");

  MPI_Allgatherv(sendValues.data(), sendValues.size() * sizeof(globalIdOwner),
                 MPI_CHAR, edgesWithTargets.data(), sizes.data(),
                 displacements.data(), MPI_CHAR, ttk::MPIcomm_);
  this->printMsg("R" + std::to_string(ttk::MPIrank_) + " got the results");

  // now each rank has a vector consisting of gIds, the ranks to which
  // they belong and the ascending / descending target we have all the
  // information on R0 which we need to resolve any manifolds stretching
  // over multiple ranks
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToManifoldMap;

  for(size_t i = 0; i < edgesWithTargets.size(); i++) {
    globalIdOwner currentVal = edgesWithTargets[i];
    if(computeAscending) {
      gIdToManifoldMap.insert(
        std::make_pair(currentVal.globalId, currentVal.ascendingTarget));
    } else {
      gIdToManifoldMap.insert(
        std::make_pair(currentVal.globalId, currentVal.descendingTarget));
    }
  }

  // now we need to check for graphs in the map and iteratively compress
  // them
  bool changed = true;
  while(changed) {
    changed = false;
    for(auto &it : gIdToManifoldMap) {
      if(gIdToManifoldMap.count(it.second) && (it.first != it.second)
         && (gIdToManifoldMap[it.first] != gIdToManifoldMap[it.second])) {
        gIdToManifoldMap[it.first] = gIdToManifoldMap[it.second];
        changed = true;
      }
    }
  }
  // now each rank simply needs to walk over their vertices and replace
  // ones aiming to ghostcells with the correct ones from the map
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < nVertices; i++) {
    ttk::SimplexId gid = triangulation.getVertexGlobalId(i);
    ttk::SimplexId val = manifold[i];
    if(gIdToManifoldMap.count(val)) {
      manifold[i] = gIdToManifoldMap[val];
    } else if(gIdToManifoldMap.count(gid)) {
      manifold[i] = gIdToManifoldMap[gid];
    }
  }
}
#endif // TTK_ENABLE_MPI


  if(computeAscending) {
    this->printMsg("Ascending segmentation computed", 1.0,
                   localTimer.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  } else {
    this->printMsg("Descending segmentation computed", 1.0,
                   localTimer.getElapsedTime(), this->threadNumber_,
                   debug::LineMode::NEW, debug::Priority::DETAIL);
  }

  return 0; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computeMSHash(
  SimplexId *const morseSmaleSegmentation,
  const SimplexId *const ascSegmentation,
  const SimplexId *const dscSegmentation,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;

  const size_t nVerts = triangulation.getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(static) num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleSegmentation[i]
      = pcp::getHash(ascSegmentation[i], dscSegmentation[i]);
  }

  this->printMsg("Morse-Smale segmentation hash computed", 1.0,
                 localTimer.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  return 0;
}