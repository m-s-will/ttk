/// \ingroup base
/// \class ttk::PathCompressionDistributedTest
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %PathCompressionDistributedTest class that computes for each vertex of a
/// triangulation the vertices to which ascending and descending manifolds it belongs.
///
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <map>
#include <unordered_map>
#include <mpi.h>
#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "size_t size not found"
#endif


namespace ttk {


  /**
   * The PathCompressionDistributedTest class provides methods to compute for each vertex of a
   * triangulation the vertices to which ascending and descending manifolds it belongs.
   */
  class PathCompressionDistributedTest : virtual public Debug {

  public:

    struct globalIdOwner {
      ttk::SimplexId globalId;
      int ownerRank;
      ttk::SimplexId ascendingTarget = -1;
      ttk::SimplexId descendingTarget = -1;


      globalIdOwner(ttk::SimplexId _globalId = -1, int _owner = -1)
        : globalId(_globalId), ownerRank(_owner) {}

      bool operator<(const globalIdOwner& other) const {
        return globalId < other.globalId;
      }

      bool operator==(const globalIdOwner& other) const{
        return globalId == other.globalId;
      }
    };



    PathCompressionDistributedTest();
    
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    std::vector<int> compressArray(const std::vector<int>& input) const {
      ttk::Timer compressTimer;
    
      std::vector<int> output(input.size());
      std::unordered_map<int, int> uniquesMap;
      int counter = 0;
      // assemble the output by creating a map of unique values while running over the array
      for (size_t i = 0; i < input.size(); i++){
        if ( uniquesMap.find(input[i]) == uniquesMap.end() ){
          uniquesMap[input[i]] = counter;
          counter++;
        }
        output[i] = uniquesMap[input[i]];
      }
      this->printMsg("#Unique Segmentations: " + std::to_string(counter), 1, compressTimer.getElapsedTime()); 
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

        size_t nVertices = triangulation->getNumberOfVertices();
        
        std::vector<ttk::SimplexId> previousDesc(nVertices);
        std::vector<ttk::SimplexId> currentDesc(nVertices);
        std::vector<ttk::SimplexId> previousAsc(nVertices);
        std::vector<ttk::SimplexId> currentAsc(nVertices);
        int numProcs;
        int rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Datatype mpi_values;
        const int nitems = 4;
        int blocklengths[4] = {1, 1, 1, 1};

        #ifdef TTK_ENABLE_64BIT_IDS
          MPI_Datatype MIT = MPI_LONG_LONG_INT;
        #else
          MPI_Datatype MIT = MPI_INT;
        #endif
        MPI_Datatype types[4] = {MIT, MPI_INT, MIT, MIT};
        MPI_Aint offsets[4];
        offsets[0] = offsetof(globalIdOwner, globalId);
        offsets[1] = offsetof(globalIdOwner, ownerRank);
        offsets[2] = offsetof(globalIdOwner, ascendingTarget);
        offsets[3] = offsetof(globalIdOwner, descendingTarget);
        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_values);
        MPI_Type_commit(&mpi_values);
        int intTag = 100;
        int structTag = 101;
        this->printMsg("Initializing MPI done for rank " + std::to_string(rank));

        std::vector<globalIdOwner> foreignVertices;
        std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdTolIdMap;
        // for the first step we initialize each vertex with the id of their largest / smallest neighbor. Afterwards we only compare the arrays
        for(size_t i = 0; i < nVertices; i++) {
          gIdTolIdMap.insert(std::make_pair(globalIds[i],i));
          int nNeighbors = triangulation->getVertexNeighborNumber(i);
          // local id in this rank
          ttk::SimplexId localNeighborId;
          // global id over all ranks
          ttk::SimplexId neighborId;
          float smallest = inputData[i];
          float largest = inputData[i];
          // if there is no larger / smaller neighbor, the vertex points to itself and is therefore a maximum / minimum
          // we do not need to check for equality, because we use the order array
          previousDesc[i] = globalIds[i];
          previousAsc[i] = globalIds[i];
          // if the vertex belongs to ourselves, we don't need to strictly point to ourselves, but to the largest neighbor
          if (rankArray[i] == rank) {
            for(int j = 0; j < nNeighbors; j++) {
              triangulation->getVertexNeighbor(i, j, localNeighborId);
              neighborId = globalIds[localNeighborId];
              // and for the largest neighbor to get to the ascending manifold
              if (inputData[neighborId] > largest){
                previousAsc[i] = neighborId;
                largest = inputData[neighborId];      
              }
              // we're checking for the smallest neighbor to get the descending manifold
              if (inputData[neighborId] <= smallest){
                previousDesc[i] = neighborId;
                smallest = inputData[neighborId];
              }
            }
          } else {
            globalIdOwner GIO = {globalIds[i], rankArray[i]};
            foreignVertices.push_back(GIO);

          }
        }

        this->printMsg("Initializing values done for rank " + std::to_string(rank));
        // now we swap between the two arrays until nothing changes anymore ergo all paths are finished
        int step = 0;
        bool same = false;
        int nextDesc = 0;
        int nextAsc = 0;
        while (!same){
          this->printMsg("Rank " + std::to_string(rank) + ", running Step "+std::to_string(step));
          same = true;
          if (step % 2 == 0){
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for (size_t i = 0; i < nVertices; i++){
              nextDesc = previousDesc[gIdTolIdMap[previousDesc[i]]];
              nextAsc = previousAsc[gIdTolIdMap[previousAsc[i]]];
              if (nextDesc != currentDesc[i]){currentDesc[i] = nextDesc; same = false;}
              if (nextAsc != currentAsc[i]){currentAsc[i] = nextAsc; same = false;}
            }
          } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for (size_t i = 0; i < nVertices; i++){
              nextDesc = currentDesc[gIdTolIdMap[currentDesc[i]]];
              nextAsc = currentAsc[gIdTolIdMap[currentAsc[i]]];
              if (nextDesc != previousDesc[i]){previousDesc[i] = nextDesc; same = false;}
              if (nextAsc != previousAsc[i]){previousAsc[i] = nextAsc; same = false;}
            }
          }
          
          step++;
        }

        this->printMsg("Rank " + std::to_string(rank) + ", finished own values in Step "+std::to_string(step), 1, localTimer.getElapsedTime());

        // now we need to request the values we still need from other ranks
        // rank 0 builds up out transferance map over ranks
        std::vector<globalIdOwner> edgesWithTargets;
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
          // construct the set with the globalids not owned by rank 0
          std::set<globalIdOwner> edges(foreignVertices.begin(), foreignVertices.end());
          // receive the data from all ranks
          // rank 0 gets all the needed global ids and their original owners, which rank actually needs those gids doesn't matter,
          // because in the end we send the finished edges back
          for (int r = 0; r < numProcs; r++){
            if (r != 0){
              size_t receivedSize;
              MPI_Recv(&receivedSize, 1, my_MPI_SIZE_T, r, intTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              if (receivedSize > 0){
                std::vector<globalIdOwner> receivedIds(receivedSize);
                MPI_Recv(receivedIds.data(), receivedSize, mpi_values, r, structTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::copy(receivedIds.begin(), receivedIds.end(), std::inserter(edges, edges.end()));
              }
            }
          }
          this->printMsg("Rank 0 received Ids which are needed");

          // we have all the gids which are needed by _some_ rank and the owner of them,
          // now we have to request from the owners to which vertices these gids are pointing to build our map
          std::vector<std::vector<globalIdOwner>> valuesFromRanks;
          valuesFromRanks.resize(numProcs);
          auto edgesIt = edges.begin();
          while (edgesIt != edges.end()){
            globalIdOwner currentId = *edgesIt;
            int owner = currentId.ownerRank;
            valuesFromRanks[owner].push_back(currentId);
            edgesIt++;
          }
          this->printMsg("Rank 0 reordered Ids");

          // for r = 0, we don't need to send anything, we process everything locally
          std::vector<globalIdOwner> fromRank0 = valuesFromRanks[0];
          std::vector<globalIdOwner> edgesForR0(fromRank0.size());
          for (size_t i = 0; i < fromRank0.size(); i++){
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
          this->printMsg("Rank 0 started building edgesWithTargets locally");

          // for r = 1, .., numProcs, we need to send and receive data. We first send everything to the ranks and then receive afterwards
          for (int r = 1; r < numProcs; r++){
            std::vector<globalIdOwner> fromThisRank = valuesFromRanks[r];
            size_t nValues = fromThisRank.size();
            MPI_Send(&nValues, 1, my_MPI_SIZE_T, r, intTag, MPI_COMM_WORLD);
            if (nValues > 0){
              MPI_Send(fromThisRank.data(), nValues, mpi_values, r, structTag, MPI_COMM_WORLD);
            }
          }
          this->printMsg("Rank 0 sent the needed ids to their owners");

          // we need to receive the results to which the gids are pointing from the ranks and build our map
          for (int r = 1; r < numProcs; r++){
            std::vector<globalIdOwner> fromThisRank = valuesFromRanks[r];
            size_t nValues = fromThisRank.size();
            std::vector<globalIdOwner> receivedIds;
            receivedIds.resize(nValues, {0,0});
            if (nValues > 0){
              MPI_Recv(receivedIds.data(), nValues, mpi_values, r, structTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              edgesWithTargets.insert(edgesWithTargets.end(), receivedIds.begin(), receivedIds.end());
            }
          }
          this->printMsg("Rank 0 received Ids with their targets from the owners");

          size_t totalSize = edgesWithTargets.size();
          MPI_Bcast(&totalSize, 1, my_MPI_SIZE_T, 0, MPI_COMM_WORLD);
          MPI_Bcast(edgesWithTargets.data(), totalSize, mpi_values, 0, MPI_COMM_WORLD);
          this->printMsg("Rank 0 broadcasted the result");

        } else {  // the other ranks
          size_t nValues = foreignVertices.size();
          MPI_Send(&nValues, 1, my_MPI_SIZE_T, 0, intTag, MPI_COMM_WORLD);
          if (nValues > 0){
            MPI_Send(foreignVertices.data(), nValues, mpi_values, 0, structTag, MPI_COMM_WORLD);
          }
          this->printMsg("Rank " + std::to_string(rank) + " sent needed ids to Rank 0");

          // we receive a variable amount of values from R0
          size_t receivedSize;
          MPI_Recv(&receivedSize, 1, my_MPI_SIZE_T, 0, intTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          if (receivedSize > 0){
            std::vector<globalIdOwner> receivedIds(receivedSize);
            MPI_Recv(receivedIds.data(), receivedSize, mpi_values, 0, structTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            this->printMsg("Rank " + std::to_string(rank) + " received needed ids of which itself is the owner from Rank 0");

            // now we need to find to where these gids point and send the values back to R0
            std::vector<globalIdOwner> sendValues(receivedSize);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for (size_t i = 0; i < receivedSize; i++){
              globalIdOwner currentVal = receivedIds[i];
              ttk::SimplexId gId = currentVal.globalId;
              ttk::SimplexId lId = gIdTolIdMap[gId];
              ttk::SimplexId ascendingTarget = currentAsc[lId];
              ttk::SimplexId descendingTarget = currentDesc[lId];
              currentVal.ascendingTarget = ascendingTarget;
              currentVal.descendingTarget = descendingTarget;
              sendValues[i] = currentVal;
            }

            MPI_Send(sendValues.data(), receivedSize, mpi_values, 0, structTag, MPI_COMM_WORLD);
            this->printMsg("Rank " + std::to_string(rank) + " sent owned ids with their targets to Rank 0");
            // receive the values from R0
            size_t totalSize;
            MPI_Bcast(&totalSize, 1, my_MPI_SIZE_T, 0, MPI_COMM_WORLD);
            edgesWithTargets.resize(totalSize);
            MPI_Bcast(edgesWithTargets.data(), totalSize, mpi_values, 0, MPI_COMM_WORLD);
            this->printMsg("Rank " + std::to_string(rank) + " got the results");

          }
        }


        // now each rank has a vector consisting of gIds, the ranks to which they belong and the ascending / descending target
        // we have all the information on R0 which we need to resolve any manifolds stretching over multiple ranks
        std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToAscendingMap;
        std::unordered_map<ttk::SimplexId, ttk::SimplexId> gIdToDescendingMap;

        for (size_t i = 0; i < edgesWithTargets.size(); i++){
          globalIdOwner currentVal = edgesWithTargets[i];
          gIdToAscendingMap.insert(std::make_pair(currentVal.globalId, currentVal.ascendingTarget));
          gIdToDescendingMap.insert(std::make_pair(currentVal.globalId, currentVal.descendingTarget));
        }

        // now we need to check for graphs in the map and iteratively compress them
        bool changed = true;
        while (changed){
          changed = false;
          for (auto& it: gIdToAscendingMap) {
            if (gIdToAscendingMap.count(it.second) && (it.first != it.second) && (gIdToAscendingMap[it.first] != gIdToAscendingMap[it.second])){
              gIdToAscendingMap[it.first] = gIdToAscendingMap[it.second];
              changed = true;
            }
          }

          for (auto& it: gIdToDescendingMap) {
            if (gIdToDescendingMap.count(it.second) && (it.first != it.second) && (gIdToDescendingMap[it.first] != gIdToDescendingMap[it.second])){
              gIdToDescendingMap[it.first] = gIdToDescendingMap[it.second];
              changed = true;
            }
          }
        }
        // now each rank simply needs to walk over foreignVertices and replace the values at these gIds with the ones from the map

        for (globalIdOwner val : foreignVertices){
          ttk::SimplexId gId = val.globalId;
          ttk::SimplexId lId = gIdTolIdMap[gId];

          currentAsc[lId] = gIdToAscendingMap[gId];
          currentDesc[lId] = gIdToDescendingMap[gId];
        }


        // compress the arrays into the ranges of 0 - #segmentation areas
        //currentDesc = this->compressArray(currentDesc);
        //currentAsc = this->compressArray(currentAsc);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for (size_t i = 0; i < nVertices; i++){
            descendingManifold[i] = currentDesc[i];
            ascendingManifold[i] = currentAsc[i];
        }
        
        
        
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
