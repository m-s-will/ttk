/// \ingroup base
/// \class ttk::PathCompressionDistributedTest
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %PathCompressionDistributedTest class that computes
/// for each vertex of a triangulation the vertices to which ascending and
/// descending manifolds it belongs.
///
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <limits.h>
#include <map>
#include <stdint.h>
#include <unordered_map>

namespace ttk {

  /**
   * The PathCompressionDistributedTest class provides methods to compute for
   * each vertex of a triangulation the vertices to which ascending and
   * descending manifolds it belongs.
   */
  class PathCompressionDistributedTest : virtual public Debug {

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

    PathCompressionDistributedTest();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    std::vector<ttk::SimplexId>
      compressArray(const std::vector<ttk::SimplexId> &input) const {
      ttk::Timer compressTimer;

      std::vector<ttk::SimplexId> output(input.size());
      std::unordered_map<ttk::SimplexId, ttk::SimplexId> uniquesMap;
      ttk::SimplexId counter = 0;
      // assemble the output by creating a map of unique values while running
      // over the array
      for(size_t i = 0; i < input.size(); i++) {
        if(uniquesMap.find(input[i]) == uniquesMap.end()) {
          uniquesMap[input[i]] = counter;
          counter++;
        }
        output[i] = uniquesMap[input[i]];
      }
      this->printMsg("#Unique Segmentations: " + std::to_string(counter), 1,
                     compressTimer.getElapsedTime());
      return output;
    }

    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeCompression(ttk::SimplexId *descendingManifold,
                           ttk::SimplexId *ascendingManifold,
                           const dataType *inputData,
                           const int *rankArray,
                           const ttk::SimplexId *globalIds,
                           const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Minima
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing compression",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_);

        ttk::SimplexId nVertices = triangulation->getNumberOfVertices();

        std::vector<ttk::SimplexId> previousDesc(nVertices);
        std::vector<ttk::SimplexId> currentDesc(nVertices);
        std::vector<ttk::SimplexId> previousAsc(nVertices);
        std::vector<ttk::SimplexId> currentAsc(nVertices);
        int numProcs;
        int rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int intTag = 100;
        int structTag = 101;
        this->printMsg("Initializing MPI done for rank "
                       + std::to_string(rank));

        std::vector<globalIdOwner> foreignVertices;
        std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdTolIdMap;
        // for the first step we initialize each vertex with the id of their
        // largest / smallest neighbor. Afterwards we only compare the arrays
        for(ttk::SimplexId i = 0; i < nVertices; i++) {
          gIdTolIdMap.insert(std::make_pair(globalIds[i], i));
          int nNeighbors = triangulation->getVertexNeighborNumber(i);
          ttk::SimplexId neighborId;
          dataType smallest = inputData[i];
          dataType largest = inputData[i];
          // if there is no larger / smaller neighbor, the vertex points to
          // itself and is therefore a maximum / minimum we do not need to check
          // for equality, because we use the order array
          previousDesc[i] = i;
          previousAsc[i] = i;
          // if the vertex belongs to ourselves, we don't need to strictly point
          // to ourselves, but to the largest neighbor
          if(rankArray[i] == rank) {
            for(int j = 0; j < nNeighbors; j++) {
              triangulation->getVertexNeighbor(i, j, neighborId);
              // and for the largest neighbor to get to the ascending manifold
              if(inputData[neighborId] > largest) {
                previousAsc[i] = neighborId;
                largest = inputData[neighborId];
              }
              // we're checking for the smallest neighbor to get the descending
              // manifold
              if(inputData[neighborId] < smallest) {
                previousDesc[i] = neighborId;
                smallest = inputData[neighborId];
              }
            }
          } else {
            globalIdOwner GIO = {globalIds[i], rankArray[i]};
            foreignVertices.push_back(GIO);
          }
        }

        this->printMsg("Initializing values done for rank "
                       + std::to_string(rank));
        // now we swap between the two arrays until nothing changes anymore ergo
        // all paths are finished
        int step = 0;
        bool same = false;
        while(!same) {
          this->printMsg("Rank " + std::to_string(rank) + ", running Step "
                         + std::to_string(step));
          same = true;
          if(step % 2 == 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for(ttk::SimplexId i = 0; i < nVertices; i++) {
              int nextDesc = previousDesc[previousDesc[i]];
              int nextAsc = previousAsc[previousAsc[i]];
              if(nextDesc != currentDesc[i]) {
                currentDesc[i] = nextDesc;
                same = false;
              }
              if(nextAsc != currentAsc[i]) {
                currentAsc[i] = nextAsc;
                same = false;
              }
            }
          } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for(ttk::SimplexId i = 0; i < nVertices; i++) {
              int nextDesc = currentDesc[currentDesc[i]];
              int nextAsc = currentAsc[currentAsc[i]];
              if(nextDesc != previousDesc[i]) {
                previousDesc[i] = nextDesc;
                same = false;
              }
              if(nextAsc != previousAsc[i]) {
                previousAsc[i] = nextAsc;
                same = false;
              }
            }
          }

          step++;
        }

// now we need to transform local ids into global ids to correctly work over all
// ranks
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(ttk::SimplexId i = 0; i < nVertices; i++) {
          currentDesc[i] = globalIds[currentDesc[i]];
          currentAsc[i] = globalIds[currentAsc[i]];
        }
        this->printMsg("Rank " + std::to_string(rank)
                         + ", finished own values in Step "
                         + std::to_string(step),
                       1, localTimer.getElapsedTime());

        // now we need to request the values we still need from other ranks
        // R0 builds up out transferance map over ranks
        std::vector<globalIdOwner> edgesWithTargets;
        MPI_Barrier(MPI_COMM_WORLD);
        // MPI_Reduce to get the size of everything we need
        int localSize = foreignVertices.size();
        this->printMsg("Localsize: " + std::to_string(localSize));
        int totalSize;
        MPI_Reduce(&localSize, &totalSize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&totalSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        this->printMsg("Rank " + std::to_string(rank) + " got the totalsize " + std::to_string(totalSize));

        std::vector<globalIdOwner> edges;
        if(rank == 0) {
          edges.resize(totalSize);
          // construct the set with the globalids not owned by R0
          int displacements[numProcs];

          // first we use MPI_Gather to get the size of each rank to populate the displacements
          int sizes[numProcs];
          MPI_Gather(&localSize, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
          for (int i = 0; i < numProcs; i++){
            sizes[i] = sizes[i] * sizeof(globalIdOwner);

          }
          displacements[0] = 0;
          // build our displacements
          for (int i = 1; i < numProcs; i++){
            displacements[i] = displacements[i-1] + sizes[i-1];
          }

          // then we use MPI_Gatherv to get the data that is needed
          // R0 gets all the needed global ids and their original owners,
          // which rank actually needs those gids doesn't matter, because in the
          // end we send the finished edges back
          MPI_Gatherv(foreignVertices.data(),
                      foreignVertices.size() * sizeof(globalIdOwner),
                      MPI_CHAR,
                      edges.data(),
                      sizes,
                      displacements,
                      MPI_CHAR,
                      0,
                      MPI_COMM_WORLD);
          this->printMsg("R0 received " + std::to_string(edges.size())
                         + " ids which are needed");

          // we have all the gids which are needed by _some_ rank and the owner
          // of them, now we have to request from the owners to which vertices
          // these gids are pointing to build our map
          std::vector<std::vector<globalIdOwner>> valuesFromRanks;
          valuesFromRanks.resize(numProcs);
          for(globalIdOwner currentId : edges) {
            // das ist ein tragendes print
            //this->printMsg("Gid " + std::to_string(currentId.globalId) + " belongs to rank " + std::to_string(currentId.ownerRank));
            valuesFromRanks[currentId.ownerRank].push_back(currentId);
          }
          this->printMsg("R0 reordered Ids");

          // for r = 0, we don't need to send anything, we process everything
          // locally
          std::vector<globalIdOwner> fromRank0 = valuesFromRanks[0];
          std::vector<globalIdOwner> edgesForR0(fromRank0.size());
          for(size_t i = 0; i < fromRank0.size(); i++) {
            globalIdOwner currentVal = fromRank0[i];
            ttk::SimplexId gId = currentVal.globalId;
            ttk::SimplexId lId = gIdTolIdMap[gId];
            ttk::SimplexId ascendingTarget = currentAsc[lId];
            ttk::SimplexId descendingTarget = currentDesc[lId];
            currentVal.ascendingTarget = ascendingTarget;
            currentVal.descendingTarget = descendingTarget;
            edgesForR0[i] = currentVal;
          }

          edgesWithTargets.assign(edgesForR0.begin(), edgesForR0.end());
          this->printMsg("R0 started building edgesWithTargets locally");

          // for r = 1, .., numProcs, we need to send and receive data. We first
          // send everything to the ranks and then receive afterwards
          for(int r = 1; r < numProcs; r++) {
            std::vector<globalIdOwner> fromThisRank = valuesFromRanks[r];
            int nValues = fromThisRank.size();
            MPI_Send(&nValues, 1, MPI_INT, r, intTag, MPI_COMM_WORLD);
            if(nValues > 0) {
              MPI_Send(fromThisRank.data(), nValues * sizeof(globalIdOwner), MPI_CHAR, r, structTag,
                       MPI_COMM_WORLD);
            }
          }
          this->printMsg("R0 sent the needed ids to their owners");

          // we need to receive the results to which the gids are pointing from
          // the ranks and build our map
          for(int r = 1; r < numProcs; r++) {
            ttk::SimplexId nValues = valuesFromRanks[r].size();
            std::vector<globalIdOwner> receivedIds;
            receivedIds.resize(nValues, {0, 0});
            if(nValues > 0) {
              MPI_Status status;
              MPI_Recv(receivedIds.data(), nValues * sizeof(globalIdOwner), MPI_CHAR, r, structTag,
                       MPI_COMM_WORLD, &status);
              int amount;
              MPI_Get_count(&status, MPI_CHAR, &amount);
              this->printMsg("R0 received " + std::to_string(amount / sizeof(globalIdOwner)) + " Ids with their targets from R" + std::to_string(r));

              edgesWithTargets.insert(
                edgesWithTargets.end(), receivedIds.begin(), receivedIds.end());
            }
          }
          this->printMsg(
            "R0 received Ids with their targets from the owners");
        } else { // the other ranks
          // first send the number of ids this rank needs to the root, then the ids themselves
          // the NULL attributes are only relevant for the root rank
          MPI_Gather(&localSize, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Gatherv(foreignVertices.data(),
                      foreignVertices.size() * sizeof(globalIdOwner),
                      MPI_CHAR,
                      NULL,
                      NULL,
                      NULL,
                      MPI_CHAR,
                      0,
                      MPI_COMM_WORLD);

          // we receive a variable amount of values from R0
          int receivedSize;
          MPI_Recv(&receivedSize, 1, MPI_INT, 0, intTag, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          if(receivedSize > 0) {
            std::vector<globalIdOwner> receivedIds;
            receivedIds.resize(receivedSize, {0, 0});
            MPI_Recv(receivedIds.data(), receivedSize * sizeof(globalIdOwner), MPI_CHAR, 0, structTag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            this->printMsg("Rank " + std::to_string(rank) + " is the owner of "
                           + std::to_string(receivedSize)
                           + " ids and received them from R0");

            // now we need to find to where these gids point and send the values
            // back to R0
            std::vector<globalIdOwner> sendValues(receivedSize);

            for(ttk::SimplexId i = 0; i < receivedSize; i++) {
              globalIdOwner currentVal = receivedIds[i];
              ttk::SimplexId gId = currentVal.globalId;
              ttk::SimplexId lId = gIdTolIdMap[gId];
              ttk::SimplexId ascendingTarget = currentAsc[lId];
              ttk::SimplexId descendingTarget = currentDesc[lId];
              currentVal.ascendingTarget = ascendingTarget;
              currentVal.descendingTarget = descendingTarget;
              sendValues[i] = currentVal;
            }
            this->printMsg("Rank " + std::to_string(rank)
                           + " wants to sent " + std::to_string(receivedSize) + " owned ids with their targets to R0");

            MPI_Ssend(sendValues.data(), receivedSize * sizeof(globalIdOwner), MPI_CHAR, 0, structTag,
                     MPI_COMM_WORLD);
            this->printMsg("Rank " + std::to_string(rank)
                           + " sent owned ids with their targets to R0");
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        edgesWithTargets.resize(totalSize);
        // all ranks receive the values from R0

        MPI_Bcast(edgesWithTargets.data(), totalSize * sizeof(globalIdOwner), MPI_CHAR, 0,
                  MPI_COMM_WORLD);
        this->printMsg("Rank " + std::to_string(rank) + " got the results");

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
               && (gIdToAscendingMap[it.first]
                   != gIdToAscendingMap[it.second])) {
              gIdToAscendingMap[it.first] = gIdToAscendingMap[it.second];
              changed = true;
            }
          }

          for(auto &it : gIdToDescendingMap) {
            if(gIdToDescendingMap.count(it.second) && (it.first != it.second)
               && (gIdToDescendingMap[it.first]
                   != gIdToDescendingMap[it.second])) {
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
          ttk::SimplexId descVal = currentDesc[i];
          if(gIdToDescendingMap.count(descVal))
            currentDesc[i] = gIdToDescendingMap[descVal];
          ttk::SimplexId ascVal = currentAsc[i];
          if(gIdToAscendingMap.count(ascVal))
            currentAsc[i] = gIdToAscendingMap[ascVal];
        }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(ttk::SimplexId i = 0; i < nVertices; i++) {
          descendingManifold[i] = currentDesc[i];
          ascendingManifold[i] = currentAsc[i];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Compression",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  }; // PathCompressionDistributedTest class

} // namespace ttk
