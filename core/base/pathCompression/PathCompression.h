/// \ingroup base
/// \class ttk::PathCompression
/// \author Robin Maack <maack@rptu.de>
/// \date January 2023.
///
/// \brief TTK processing package for the computation of Morse-Smale
/// segmentations using Path Compression.

#pragma once

// base code includes
#include <Triangulation.h>

using ttk::SimplexId;

#ifdef TTK_ENABLE_64BIT_IDS
constexpr unsigned long long int hash_max = ULONG_LONG_MAX;

constexpr unsigned long long int getHash(
  const unsigned long long int a, const unsigned long long int b)  {
  return (a*b + (a*a) + (b*b) + (a*a*a)*(b*b*b)) % hash_max;
	//return std::rotl(a,1) ^ b;
}
#else
constexpr unsigned int hash_max = UINT_MAX;

constexpr unsigned int getHash(
  const unsigned int a, const unsigned int b) {
  return (a*b + (a*a) + (b*b) + (a*a*a)*(b*b*b)) % hash_max;
  //return std::rotl(a,1) ^ b;
}
#endif

namespace ttk {
  class PathCompression : public virtual Debug {
  public:
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

    PathCompression();

    /** @brief Pointers to pre-allocated segmentation point data arrays */
    struct OutputManifold {
      SimplexId *ascending_;
      SimplexId *descending_;
      SimplexId *morseSmale_;
    };

    /**
     * Main function for computing the Morse-Smale complex.
     *
     * @pre PathCompression::preconditionTriangulation must be
     * called prior to this.
     */
    template <typename triangulationType>
    inline int execute(OutputManifold &outManifold,
                       const SimplexId *const offsets,
                       const triangulationType &triangulation,
                       const ttk::SimplexId *globalIds = nullptr);

    /**
     * Enable/Disable computation of the geometrical embedding of the
     * manifolds of the critical points.
     */
    inline void setComputeSegmentation(const bool doAscending,
                                       const bool doDescending,
                                       const bool doMorseSmale) {
      this->ComputeAscendingSegmentation = doAscending;
      this->ComputeDescendingSegmentation = doDescending;
      this->ComputeFinalSegmentation = doMorseSmale;
    }

    /**
     * Set the input triangulation and preprocess the needed
     * mesh traversal queries.
     */
    inline void preconditionTriangulation(AbstractTriangulation *const data) {
      data->preconditionVertexNeighbors();
#if TTK_ENABLE_MPI
      if(ttk::isRunningWithMPI())
        data->preconditionDistributedVertices();
#endif
    }

    template <typename triangulationType>
    int computePathCompression(SimplexId *const ascManifold,
                               SimplexId *const dscManifold,
                               const SimplexId *const orderArr,
                               const triangulationType &triangulation,
                               const ttk::SimplexId *globalIds = nullptr) const;

    template <typename triangulationType>
    int computePathCompressionSingle(SimplexId *const manifold,
                                     const bool computeAscending,
                                     const SimplexId *const orderArr,
                                     const triangulationType &triangulation,
                                     const ttk::SimplexId *globalIds
                                     = nullptr) const;

    template <typename triangulationType>
    int computeFinalSegmentation(
      SimplexId *const morseSmaleManifold,
      const SimplexId *const ascManifold,
      const SimplexId *const desManifold,
      const triangulationType &triangulation) const;

    bool ComputeDescendingSeparatrices2{false};
    bool ComputeAscendingSegmentation{true};
    bool ComputeDescendingSegmentation{true};
    bool ComputeFinalSegmentation{true};

  };
} // namespace ttk

// ---------------- //
//  Execute method  //
// ---------------- //

template <typename triangulationType>
int ttk::PathCompression::execute(OutputManifold &outManifold,
                                  const SimplexId *const orderArray,
                                  const triangulationType &triangulation,
                                  const ttk::SimplexId *globalIds) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(orderArray == nullptr) {
    this->printErr("Input offset field pointer is null.");
    return -1;
  }
#endif

  if((ComputeAscendingSegmentation && ComputeDescendingSegmentation) ||
    ComputeFinalSegmentation) {
    computePathCompression(outManifold.ascending_, outManifold.descending_,
                           orderArray, triangulation, globalIds);
  } else if(ComputeAscendingSegmentation) {
    computePathCompressionSingle(
      outManifold.ascending_, true, orderArray, triangulation, globalIds);
  } else if(ComputeDescendingSegmentation) {
    computePathCompressionSingle(
      outManifold.descending_, false, orderArray, triangulation, globalIds);
  }

  if(ComputeFinalSegmentation) {
    computeFinalSegmentation(outManifold.morseSmale_, outManifold.ascending_,
      outManifold.descending_, triangulation);
  }

  return 0;
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompression(
  SimplexId *const ascManifold,
  SimplexId *const dscManifold,
  const SimplexId *const orderArr,
  const triangulationType &triangulation,
  const ttk::SimplexId *globalIds) const {

  ttk::Timer localTimer;
  bool useMPI = false;
  std::vector<globalIdOwner> foreignVertices;

  TTK_FORCE_USE(useMPI);
  TTK_FORCE_USE(foreignVertices);
  TTK_FORCE_USE(globalIds);
#if TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI() && globalIds != nullptr) {
    useMPI = true;
  }
#endif
  //this->printWrn(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  const std::string msg = "Computing Asc + Desc Segmentation";
  this->printMsg(msg,0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation.getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP
  this->printMsg("in 1 thread mode");
  SimplexId nActiveVertices;
  std::vector<SimplexId> activeVertices;
  activeVertices.reserve(nVertices);
  // find maxima and intialize vector of not fully compressed vertices
  for(SimplexId i = 0; i < nVertices; i++) {
    SimplexId neighborId;
    const SimplexId numNeighbors =
      triangulation.getVertexNeighborNumber(i);

    bool hasLargerNeighbor = false;
    SimplexId &dmi = dscManifold[i];
    dmi = i;

    bool hasSmallerNeighbor = false;
    SimplexId &ami = ascManifold[i];
    ami = i;

    // check all neighbors
    for(SimplexId n = 0; n < numNeighbors; n++) {
      triangulation.getVertexNeighbor(i, n, neighborId);

      if(orderArr[neighborId] < orderArr[ami]) {
        ami = neighborId;
        hasSmallerNeighbor = true;
      } else if(orderArr[neighborId] > orderArr[dmi]) {
        dmi = neighborId;
        hasLargerNeighbor = true;
      }
    }

    if(hasLargerNeighbor || hasSmallerNeighbor) {
      activeVertices.push_back(i);
    }
  }

  nActiveVertices = activeVertices.size();
  size_t currentIndex = 0;

  // compress paths until no changes occur
  while(nActiveVertices > 0) {
    for(SimplexId i = 0; i < nActiveVertices; i++) {
      SimplexId &v = activeVertices[i];
      SimplexId &vDsc = dscManifold[v];
      SimplexId &vAsc = ascManifold[v];

      // compress paths
      vDsc = dscManifold[vDsc];
      vAsc = ascManifold[vAsc];

      // check if not fully compressed
      if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
        activeVertices[currentIndex++] = v;
      }
    }

    nActiveVertices = currentIndex;
    currentIndex = 0;
  }

#ifdef TTK_ENABLE_OPENMP
} else {
  this->printMsg("in multithread mode");

#pragma omp parallel num_threads(threadNumber_)
  {
    std::vector<SimplexId> lActiveVertices;
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

    // find the biggest neighbor for each vertex
    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &dmi = dscManifold[i];
      dmi = i;

      bool hasSmallerNeighbor = false;
      SimplexId &ami = ascManifold[i];
      ami = i;

      // check all neighbors
#ifdef TTK_ENABLE_MPI
      if(useMPI) {
        if(triangulation.getVertexRank(i) == ttk::MPIrank_) {
          for(SimplexId n = 0; n < numNeighbors; n++) {
            triangulation.getVertexNeighbor(i, n, neighborId);

            if(orderArr[neighborId] < orderArr[ami]) {
              ami = neighborId;
              hasSmallerNeighbor = true;
            } else if(orderArr[neighborId] > orderArr[dmi]) {
              dmi = neighborId;
              hasLargerNeighbor = true;
            }
          }
        } else {
          globalIdOwner GIO = {globalIds[i], triangulation.getVertexRank(i)};
          foreignVertices.push_back(GIO);
        }
      } else {
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);

          if(orderArr[neighborId] < orderArr[ami]) {
            ami = neighborId;
            hasSmallerNeighbor = true;
          } else if(orderArr[neighborId] > orderArr[dmi]) {
            dmi = neighborId;
            hasLargerNeighbor = true;
          }
        }
      }
#else
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(orderArr[neighborId] < orderArr[ami]) {
          ami = neighborId;
          hasSmallerNeighbor = true;
        } else if(orderArr[neighborId] > orderArr[dmi]) {
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
        SimplexId &v = lActiveVertices[i];
        SimplexId &vDsc = dscManifold[v];
        SimplexId &vAsc = ascManifold[v];

        // compress paths
        #pragma omp atomic read
        vDsc = dscManifold[vDsc];

        #pragma omp atomic read
        vAsc = ascManifold[vAsc];

        // check if fully compressed
        if(vDsc != dscManifold[vDsc] || vAsc != ascManifold[vAsc]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
  }
}
#endif // TTK_ENABLE_OPENMP

#ifdef TTK_ENABLE_MPI
// now we need to transform local ids into global ids to correctly work
// over all ranks
if(useMPI) {
  this->printMsg("using mpi");
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(ttk::SimplexId i = 0; i < nVertices; i++) {
    dscManifold[i] = globalIds[dscManifold[i]];
    ascManifold[i] = globalIds[ascManifold[i]];
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
    ttk::SimplexId gid = globalIds[i];
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
}
#endif // TTK_ENABLE_MPI
this->printMsg(msg, 1, localTimer.getElapsedTime(), this->threadNumber_);

return 1; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computePathCompressionSingle(
  SimplexId *const manifold,
  const bool computeAscending,
  const SimplexId *const orderArr,
  const triangulationType &triangulation,
  const ttk::SimplexId *globalIds) const {

  ttk::Timer localTimer;
  bool useMPI = false;
  std::vector<globalIdOwner> foreignVertices;

  TTK_FORCE_USE(useMPI);
  TTK_FORCE_USE(foreignVertices);
  TTK_FORCE_USE(globalIds);

#if TTK_ENABLE_MPI
  if(ttk::isRunningWithMPI() && globalIds != nullptr) {
    this->printMsg("Using mpi");
    useMPI = true;
  }
#endif
  //this->printWrn(ttk::debug::Separator::L1);
  // print the progress of the current subprocedure (currently 0%)
  const std::string msg = "Computing "+std::string(computeAscending? "Ascending" : "Descending")+" Segmentation";
  this->printMsg(msg,0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE);

  /* compute the Descending Maifold iterating over each vertex, searching
   * the biggest neighbor and compressing its path to its maximum */
  const SimplexId nVertices = triangulation.getNumberOfVertices();

#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  SimplexId nActiveVertices;
  std::vector<SimplexId> activeVertices;
  activeVertices.reserve(nVertices);
  // find maxima and intialize vector of not fully compressed vertices
  for(SimplexId i = 0; i < nVertices; i++) {
    SimplexId neighborId;
    const SimplexId numNeighbors =
      triangulation.getVertexNeighborNumber(i);

    bool hasLargerNeighbor = false;
    SimplexId &mi = manifold[i];
    mi = i;

    // check all neighbors
    for(SimplexId n = 0; n < numNeighbors; n++) {
      triangulation.getVertexNeighbor(i, n, neighborId);

      if(computeAscending) {
        if(orderArr[neighborId] < orderArr[mi]) {
          mi = neighborId;
          hasLargerNeighbor = true;
        }
      } else {
        if(orderArr[neighborId] > orderArr[mi]) {
          mi = neighborId;
          hasLargerNeighbor = true;
        }
      }
    }

    if(hasLargerNeighbor) {
      activeVertices.push_back(i);
    }
  }

  nActiveVertices = activeVertices.size();
  size_t currentIndex = 0;

  // compress paths until no changes occur
  while(nActiveVertices > 0) {
    for(SimplexId i = 0; i < nActiveVertices; i++) {
      SimplexId &v = activeVertices[i];
      SimplexId &vMan = manifold[v];

      // compress paths
      vMan = manifold[vMan];

      // check if not fully compressed
      if(vMan != manifold[vMan]) {
        activeVertices[currentIndex++] = v;
      }
    }

    nActiveVertices = currentIndex;
    currentIndex = 0;
  }

#ifdef TTK_ENABLE_OPENMP
} else {

  #pragma omp parallel num_threads(threadNumber_)
  {
    std::vector<SimplexId> lActiveVertices;
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

    // find the biggest neighbor for each vertex
    #pragma omp for schedule(static)
    for(SimplexId i = 0; i < nVertices; i++) {
      SimplexId neighborId;
      SimplexId numNeighbors =
        triangulation.getVertexNeighborNumber(i);

      bool hasLargerNeighbor = false;
      SimplexId &mi = manifold[i];
      mi = i;

      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);

        if(computeAscending) {
          if(orderArr[neighborId] < orderArr[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        } else {
          if(orderArr[neighborId] > orderArr[mi]) {
            mi = neighborId;
            hasLargerNeighbor = true;
          }
        }
      }

      if(hasLargerNeighbor) {
        lActiveVertices.push_back(i);
      }
    }

    size_t lnActiveVertices = lActiveVertices.size();
    size_t currentIndex = 0;

    // compress paths until no changes occur
    while(lnActiveVertices > 0) {
      for(size_t i = 0; i < lnActiveVertices; i++) {
        SimplexId &v = lActiveVertices[i];
        SimplexId &vMan = manifold[v];

        // compress paths
        #pragma omp atomic read
        vMan = manifold[vMan];

        // check if fully compressed
        if(vMan != manifold[vMan]) {
          lActiveVertices[currentIndex++] = v;
        }
      }

      lnActiveVertices = currentIndex;
      currentIndex = 0;
    }
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg(msg,1,localTimer.getElapsedTime(),this->threadNumber_);

  return 1; // return success
}

template <typename triangulationType>
int ttk::PathCompression::computeFinalSegmentation(
  SimplexId *const morseSmaleManifold,
  const SimplexId *const ascManifold,
  const SimplexId *const dscManifold,
  const triangulationType &triangulation) const {
    ttk::Timer localTimer;

  this->printMsg("Computing MSC Manifold",
                 0, localTimer.getElapsedTime(), this->threadNumber_,
                 ttk::debug::LineMode::REPLACE);

  const size_t nVerts = triangulation.getNumberOfVertices();


#ifdef TTK_ENABLE_OPENMP
if(threadNumber_ == 1) {
#endif // TTK_ENABLE_OPENMP

  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = getHash(ascManifold[i], dscManifold[i]);
  }

#ifdef TTK_ENABLE_OPENMP
} else {

  #pragma omp parallel for schedule(static) num_threads(threadNumber_)
  for(size_t i = 0; i < nVerts; ++i) {
    morseSmaleManifold[i] = getHash(ascManifold[i], dscManifold[i]);
  }
}
#endif // TTK_ENABLE_OPENMP

  this->printMsg("Computed MSC Manifold",
                 1, localTimer.getElapsedTime(), this->threadNumber_);

  return 0;
}
