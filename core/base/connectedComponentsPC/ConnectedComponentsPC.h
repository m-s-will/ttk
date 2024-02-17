/// \ingroup base
/// \class ttk::ConnectedComponentsPC
/// \author Robin G. C. Maack <maack@rptu.de>
/// \date May 2023.
///
/// \brief TTK processing package for the computation of Morse-Smale
/// segmentations using Path Compression.
///
/// Given an input scalar field, this class computes its ascending and
/// descending segmentation by assigning every vertex to its minimum or maximum
/// in gradient or inverse gradient direction. For convienience a hash (no hash
/// collision detection) of both segmentations can be created to represent the
/// Morse-Smale segmentation.
///
/// \b Related \b publication \n
/// "Parallel Computation of Piecewise Linear Morse-Smale Segmentations" \n
/// Robin G. C. Maack, Jonas Lukasczyk, Julien Tierny, Hans Hagen,
/// Ross Maciejewski, Christoph Garth \n
/// IEEE Transactions on Visualization and Computer Graphics \n
///
/// \sa ttkConnectedComponentsPC.cpp %for a usage example
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
  class ConnectedComponentsPC : public virtual Debug {
  public:
    ConnectedComponentsPC();

    struct globalIdOwner {
      ttk::SimplexId globalId;
      int ownerRank;
      ttk::SimplexId target = -1;

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
     * @param[in] scalarArray scalar field
     * @param[in] triangulation triangulation
     *
     * @pre ConnectedComponentsPC::preconditionTriangulation must be
     * called prior to this.
     *
     * @return 0 on success
     */
    template <typename dataType, typename triangulationType>
    inline int execute(SimplexId *segmentation_,
                       const SimplexId minSize,
                       const dataType *const scalarArray,
                       const triangulationType &triangulation);

    /**
     * @brief Compute the ascending or descending segmentation
     *
     * This function computes the ascending or descending segmentation on the
     * scalar field. First, the ascending or descending segmentation is set to
     * the largest or smallest neighbor of each vertex. Then, using path
     * compression, the vertices are assigned to their minimum/maximum in
     * positive/negative gradient direction.
     *
     * @tparam triangulationType type of triangulation
     * @param[out] segmentation segmentation
     * @param[in] scalarArray scalar array
     * @param[in] triangulation triangulation
     * @return 0 on success
     */
    template <typename dataType, typename triangulationType>
    int computeConnectedComponentsPC(
      SimplexId *const segmentation,
      const SimplexId minSize,
      const dataType *const scalarArray,
      const triangulationType &triangulation) const;
  };
} // namespace ttk

template <typename dataType, typename triangulationType>
int ttk::ConnectedComponentsPC::execute(
  SimplexId *segmentation_,
  const SimplexId minSize,
  const dataType *const scalarArray,
  const triangulationType &triangulation) {
  if(scalarArray == nullptr)
    return this->printWrn("FeatureMask is null, taking everthing as foreground.");

  Timer t;

  this->printMsg("Start computing segmentations", 0.0, t.getElapsedTime(),
                 this->threadNumber_);

  computeConnectedComponentsPC(
    segmentation_, minSize, scalarArray, triangulation);

  this->printMsg("Data-set ("
                   + std::to_string(triangulation.getNumberOfVertices())
                   + " points) processed",
                 1.0, t.getElapsedTime(), this->threadNumber_);

  return 0;
}

template <typename dataType, typename triangulationType>
int ttk::ConnectedComponentsPC::computeConnectedComponentsPC(
  SimplexId *const segmentation,
  const SimplexId minSize,
  const dataType *const scalarArray,
  const triangulationType &triangulation) const {

  ttk::Timer localTimer;
  std::vector<globalIdOwner> foreignVertices;

  TTK_FORCE_USE(foreignVertices);

  const SimplexId nVertices = triangulation.getNumberOfVertices();
  std::vector<int> featureMask(nVertices, 1);

  if (scalarArray != nullptr){
    for (SimplexId i = 0; i < nVertices; i++){
      featureMask[i] = scalarArray[i] != 0;
    }
  }

  std::vector<SimplexId> lActiveVertices;
  this->printMsg("Starting to compute active vertices");
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(threadNumber_) private(lActiveVertices)
  {
    lActiveVertices.reserve(std::ceil(nVertices / threadNumber_));

#pragma omp for
#else // TTK_ENABLE_OPENMP
  lActiveVertices.reserve(nVertices);
#endif // TTK_ENABLE_OPENMP
    // find the largest neighbor for each vertex
    for(SimplexId i = 0; i < nVertices; i++) {
      bool hasLargerNeighbor = false;
      if(featureMask[i] == 0) {
        segmentation[i] = -1;
      } else {
        SimplexId neighborId{0};
        SimplexId const numNeighbors = triangulation.getVertexNeighborNumber(i);

        SimplexId &mi = segmentation[i];
        // i is the local id in the rank
        mi = i;
#ifdef TTK_ENABLE_MPI
        if(ttk::isRunningWithMPI()) {
          if(triangulation.getVertexRank(i) == ttk::MPIrank_) {
            for(SimplexId n = 0; n < numNeighbors; n++) {
              triangulation.getVertexNeighbor(i, n, neighborId);
              if(featureMask[neighborId] == 1 && neighborId > mi) {
                mi = neighborId;
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
          for(SimplexId n = 0; n < numNeighbors; n++) {
            triangulation.getVertexNeighbor(i, n, neighborId);
            if(featureMask[neighborId] == 1 && neighborId > mi) {
              mi = neighborId;
              hasLargerNeighbor = true;
            }
          }
        }
#else
      // check all neighbors
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(featureMask[neighborId] == 1 && neighborId > mi) {
          mi = neighborId;
          hasLargerNeighbor = true;
        }
      }
#endif
      }
      if(hasLargerNeighbor) {
        //this->printMsg("Vertex " + std::to_string(i) + " has larger neighbor",
        //               1, localTimer.getElapsedTime());
        lActiveVertices.push_back(i);
      }
    }
    this->printMsg("Finished computing active vertices, starting first PC");
    size_t lnActiveVertices = lActiveVertices.size();
    size_t currentIndex = 0;
    // this->printMsg("Starting compressing paths for thread");
    //  compress paths until no changes occur
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

this->printMsg("Finished first PC, starting second vertex calculation");

//  second pathcompression to merge multiple segments in one connected component
#ifdef TTK_ENABLE_OPENMP
#pragma omp for schedule(static)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVertices; i++) {
    SimplexId neighborId{0};
    SimplexId const numNeighbors = triangulation.getVertexNeighborNumber(i);

    SimplexId &mi = segmentation[i];

#ifdef TTK_ENABLE_MPI
    if(ttk::isRunningWithMPI()) {
      if(triangulation.getVertexRank(i) == ttk::MPIrank_) {
        // check all neighbors
        if(mi != -1) {
          for(SimplexId n = 0; n < numNeighbors; n++) {
            triangulation.getVertexNeighbor(i, n, neighborId);
            if(segmentation[neighborId] > mi) {
              mi = segmentation[neighborId];
            }
          }
        }
      } else {
        globalIdOwner GIO = {
          triangulation.getVertexGlobalId(i), triangulation.getVertexRank(i)};
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
        foreignVertices.push_back(GIO);
      }
    } else {
      // check all neighbors
      if(mi != -1) {
        for(SimplexId n = 0; n < numNeighbors; n++) {
          triangulation.getVertexNeighbor(i, n, neighborId);
          if(segmentation[neighborId] > mi) {
            mi = segmentation[neighborId];
          }
        }
      }
    }
#else
    // check all neighbors
    if(mi != -1) {
      for(SimplexId n = 0; n < numNeighbors; n++) {
        triangulation.getVertexNeighbor(i, n, neighborId);
        if(segmentation[neighborId] > mi) {
          mi = segmentation[neighborId];
        }
      }
    }
#endif
  }

  this->printMsg(
    "Finished second vertex calculation, starting final compression");
  std::vector<ttk::SimplexId> currentSeg(nVertices);
  int step = 0;
  bool same = false;
  while(!same) {
    same = true;
    if(step % 2 == 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(ttk::SimplexId i = 0; i < nVertices; i++) {
        if(featureMask[i] == 1) {
          ttk::SimplexId nextSeg = segmentation[segmentation[i]];
          if(nextSeg != currentSeg[i]) {
            currentSeg[i] = nextSeg;
            same = false;
          }
        }
      }
    } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(ttk::SimplexId i = 0; i < nVertices; i++) {
        if (featureMask[i] == 1) {
          ttk::SimplexId nextSeg = currentSeg[currentSeg[i]];
          if(nextSeg != segmentation[i]) {
            segmentation[i] = nextSeg;
            same = false;
          }
        }
      }
    }
    step++;
  }
  this->printMsg("Finished final compression");

  // we need to compress everything pointing to mi to the new segmentation
  // this->printMsg("Starting compressing paths for thread");
  // compress paths until no changes occur

  // remove all the components with less than x vertices
  this->printMsg("Starting to remove small components");
  std::vector<ttk::SimplexId> componentSizes(nVertices, 0);

  #pragma omp parallel
  {
    std::vector<ttk::SimplexId> componentSizes_priv(nVertices, 0);
    for(ttk::SimplexId i = 0; i < nVertices; i++) {
      if(featureMask[i] == 1) {
        componentSizes_priv[segmentation[i]]++;
      }
    }
    #pragma omp critical
    for (ttk::SimplexId i = 0; i < nVertices; i++)
    {
      componentSizes[i] += componentSizes_priv[i];
    }
  }

  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for  num_threads(this->threadNumber_)
  #endif
  for(ttk::SimplexId i = 0; i < nVertices; i++) {
    if(featureMask[i] == 1) {
      if(componentSizes[segmentation[i]] < minSize) {
        segmentation[i] = -1;
      }
    }
  }
  this->printMsg("Finished removing small components");


  // this->printMsg("Finished compressing paths for thread");
#ifdef TTK_ENABLE_MPI
  // now we need to transform local ids into global ids to correctly work
  // over all ranks
  if(ttk::isRunningWithMPI()) {
    this->printMsg("using mpi");
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
    for(ttk::SimplexId i = 0; i < nVertices; i++) {
      if(segmentation[i] != -1)
        segmentation[i] = triangulation.getVertexGlobalId(segmentation[i]);
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
      currentVal.target = segmentation[lId];
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
    std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToSegmentationMap;

    for(size_t i = 0; i < edgesWithTargets.size(); i++) {
      globalIdOwner currentVal = edgesWithTargets[i];
      gIdToSegmentationMap.insert(
        std::make_pair(currentVal.globalId, currentVal.target));
    }

    // now we need to check for graphs in the map and iteratively compress
    // them
    bool changed = true;
    while(changed) {
      changed = false;
      for(auto &it : gIdToSegmentationMap) {
        if(gIdToSegmentationMap.count(it.second) && (it.first != it.second)
           && (gIdToSegmentationMap[it.first]
               != gIdToSegmentationMap[it.second])) {
          gIdToSegmentationMap[it.first] = gIdToSegmentationMap[it.second];
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
      ttk::SimplexId target = segmentation[i];
      if (target != -1){
        if(gIdToSegmentationMap.count(target)) {
          segmentation[i] = gIdToSegmentationMap[target];
        } else if(gIdToSegmentationMap.count(gid)) {
          segmentation[i] = gIdToSegmentationMap[gid];
        }
      }
    }
  } else {
    this->printMsg("not using mpi");
  }
#endif // TTK_ENABLE_MPI

  this->printMsg("Segmentation computed", 1.0, localTimer.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0; // return success
}
