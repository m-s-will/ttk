///
/// \ingroup base
/// \class ttk::PairExtrema
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date November 2022.
///
/// This module defines the %PairExtrema class that computes the Saddle Maximum and the split tree
/// using the ascending Segmentation and the critical points
///
/// \b Related \b publication \n
/// "A Hybrid Parallel Algorithm for Computing and Tracking Level Set Topology" \n
/// Maadasmamy et al. \n
/// 2012 19th International Conference on High Performance Computing

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <chrono>
#include <limits.h>
#include <set>
#include <numeric>


#define duration(a) \
  std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

namespace ttk {

  /**
   * The PairExtrema class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class PairExtrema : virtual public Debug {

  public:
    PairExtrema();

    struct Branch {
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
        vertices; // order, globalId, first pair is the maximum
      Branch *parentBranch = nullptr;
    };

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      if(triangulation)
        triangulation->preconditionVertexNeighbors();
      return 0;
    }

    int constructPersistencePairs(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &pairs,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &maximaTriplets,
      std::vector<std::array<ttk::SimplexId, 15>> &saddleTriplets,
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      ttk::SimplexId globalMin) {
      int step = 0;
      bool changed = true;
      // std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
      // saddleMaxPairs{};
      std::vector<ttk::SimplexId> maximumPointer(maximaLocalToGlobal.size());
      std::iota(std::begin(maximumPointer), std::end(maximumPointer), 0);
      // ttk::SimplexId nrOfPairs = 0;
      ttk::SimplexId globalMax = maximaLocalToGlobal.size() - 1;
      while(changed) {
        ttk::Timer stepTimer;
        changed = false;
        this->printMsg(ttk::debug::Separator::L2, ttk::debug::Priority::DETAIL);
        this->printMsg("Running step " + std::to_string(step),
                       //  + ", Pairs: " + std::to_string(nrOfPairs),
                       0, 0, ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

        ttk::Timer buildTimer;

        // make localtoglobal smaller by deleting used maxima?
        std::vector<ttk::SimplexId> largestSaddlesForMax(
          maximaLocalToGlobal.size(), saddleTriplets.size());
        std::vector<omp_lock_t> maximaLocks(maximaLocalToGlobal.size());
        for(size_t i = 0; i < maximaLocks.size(); i++) {
          omp_init_lock(&maximaLocks[i]);
        }
        // this->printMsg("starting with parallel step");
#pragma omp parallel num_threads(this->threadNumber_)
        {
#pragma omp for schedule(guided)
          for(ttk::SimplexId i = 0; i < (ttk::SimplexId)saddleTriplets.size(); i++) {
            auto &triplet = saddleTriplets[i];
            ttk::SimplexId temp;
            for(int p = 0; p < triplet[14]; p++) {
              // TODO if OMP 5.1 is widespread: use omp atomic compare
              const auto &max = triplet[p];
              // this->printMsg("max: " + std::to_string(max));
              if(max != globalMax) {
#pragma omp atomic read
                temp = largestSaddlesForMax[max];
                if(i < temp) {
                  omp_set_lock(&maximaLocks[max]);
                  // save only maximum saddle
                  largestSaddlesForMax[max]
                    = std::min(i, largestSaddlesForMax[max]);
                  omp_unset_lock(&maximaLocks[max]);
                }
              }
            }
          }
          //#pragma omp single
          // this->printMsg("Finished building largestSaddlesForMax");

#pragma omp single
          this->printMsg("Finished building largestSaddlesForMax", 0.2,
                         buildTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                         ttk::debug::Priority::DETAIL);
          ttk::Timer pairTimer;

          std::vector<SimplexId> lActiveMaxima;
          lActiveMaxima.reserve(maximaLocalToGlobal.size());
          //#pragma omp single
          // this->printMsg("starting finding pairs");

#pragma omp for schedule(guided)
          for(size_t i = 0; i < largestSaddlesForMax.size(); i++) {
            ttk::SimplexId maximum = i;
            auto largestSaddle = largestSaddlesForMax[maximum];
            if(largestSaddle < (ttk::SimplexId)saddleTriplets.size()) {
              const auto &triplet = saddleTriplets[largestSaddle];
              if(triplet[0]
                 == maximum) { // smalles maximum reachable from the saddle
                changed = true;
                //#pragma omp atomic
                // nrOfPairs++;
                pairs[maximum]
                  = std::make_pair(saddlesLocalToGlobal[largestSaddle],
                                   maximaLocalToGlobal[maximum]);
                auto largestMax = triplet[triplet[14] - 1];
                maximumPointer[maximum] = largestMax;
                maximaTriplets[maximum] = (std::make_pair(largestSaddle, largestMax));
                lActiveMaxima.push_back(maximum);
              }
            }
          }
          //#pragma omp single
          // this->printMsg("Finished finding pairs");

#pragma omp single
          this->printMsg("Finished finding pairs and swapping pointers", 0.5,
                         pairTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                         ttk::debug::Priority::DETAIL);

          ttk::Timer compressTimer;
          // use pathcompression on the maximumPointer
          size_t lnActiveMaxima = lActiveMaxima.size();
          size_t currentIndex = 0;

          while(lnActiveMaxima > 0) {
            for(size_t i = 0; i < lnActiveMaxima; i++) {
              ttk::SimplexId &v = lActiveMaxima[i];
              ttk::SimplexId &nextPointer = maximumPointer[v];

// compress paths
#pragma omp atomic read
              nextPointer = maximumPointer[nextPointer];
              if(nextPointer != maximumPointer[nextPointer]) {
                lActiveMaxima[currentIndex++] = v;
              }
            }
            lnActiveMaxima = currentIndex;
            currentIndex = 0;
          }
#pragma omp single
          this->printMsg("Did pathcompression on maximumpointer", 0.8,
                         compressTimer.getElapsedTime(),
                         ttk::debug::LineMode::NEW,
                         ttk::debug::Priority::DETAIL);
          ttk::Timer replaceTimer;
          // replace the values with their maximumPointers, delete the saddle if
          // max and min are the same
          //#pragma omp single
          // this->printMsg("starting replacing saddletriplets");

#pragma omp for schedule(guided)
          for(size_t i = 0; i < saddleTriplets.size(); i++) {
            auto &triplet = saddleTriplets[i];
            for(int r = 0; r < triplet[14]; r++) {
              triplet[r] = maximumPointer[triplet[r]];
            }
            sortAndRemoveUniques(triplet);
            if(triplet[14] == 1) {
              triplet[14] = 0;
            }
          }
          //#pragma omp single
          // this->printMsg("finished replacing saddletriplets");

#pragma omp single
          this->printMsg("Replaced values and deleted unnecessary saddles.", 1,
                         replaceTimer.getElapsedTime(),
                         ttk::debug::LineMode::NEW,
                         ttk::debug::Priority::DETAIL);

        } // end parallel

        this->printMsg("Finished step " + std::to_string(step), 1,
                       stepTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);
        step++;
      }
      // the global max is always in the last position of the
      // maximaLocalToGlobal vector and needs to connect with the global
      // minimum
      pairs[pairs.size() - 1] = std::make_pair(
        globalMin, maximaLocalToGlobal[maximaLocalToGlobal.size() - 1]);
      maximaTriplets[maximaTriplets.size() - 1] = std::make_pair(
        globalMin,
        maximaTriplets.size() - 1); // maximaTriplets[globalmax] = (globalmin, globalmax);

      return 1;
    }

    void sortAndRemoveUniques(std::array<ttk::SimplexId, 15> &triplet) {
      std::sort(triplet.begin(), triplet.begin() + triplet[14]);
      int tempPointer = 1;
      for(int p = 1; p < triplet[14]; p++) {
        if(triplet[p - 1]
           != triplet[p]) { // if we have a new value, we step ahead
          triplet[tempPointer] = triplet[p];
          tempPointer++;
        }
      }
      triplet[14] = tempPointer;
    }

    int constructMergeTree(
      std::vector<Branch> &branches,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
        &maximaTriplets, // maximaTriplets[max] = (saddle, largestMax)
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      const ttk::SimplexId *order) {

      maximaTriplets[maximaTriplets.size() - 1].second = maximaTriplets.size() - 1;
      // compress the maximaTriplets
      bool same = false;
      ttk::Timer compressTimer;

      while(!same) {
        same = true;
        for(size_t i = 0; i < maximaTriplets.size(); i++) {
          auto nextBranch = maximaTriplets[maximaTriplets[i].second];
          if(nextBranch.second != maximaTriplets[i].second
             && nextBranch.first
                  < maximaTriplets[i].first) { // we need to follow along larger
                                         // saddles to the largestMax
            maximaTriplets[i].second = nextBranch.second;
            same = false;
          }
        }
      }
      this->printMsg("Compressed maximaTriplets", 0, compressTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);
      ttk::Timer maximaVectorTimer;

      for(ttk::SimplexId b = 0; b < (ttk::SimplexId)maximaTriplets.size(); b++) {
        auto &triplet = maximaTriplets[b]; /// maximaTriplets[max] = (saddle, biggerMax wo
                                     /// angedockt wird)
        auto &branch = branches[b]; /// maximaTriplets[max] = (saddle, biggerMax
                                         /// wo angedockt wird)
        auto branchMaxId = maximaLocalToGlobal[b];
        branch.vertices.emplace_back(order[branchMaxId], branchMaxId);

        auto parent = triplet.second;
        if(parent != b) {
          auto saddle = saddlesLocalToGlobal[triplet.first];
          auto orderForSaddle = order[saddle];
          branch.vertices.emplace_back(orderForSaddle, saddle);
          branch.parentBranch
            = &branches[parent]; // branches[mainBranch].parentBranch
                                      // = branches[mainBranch];
          branches[parent].vertices.emplace_back(orderForSaddle, saddle);
        } else {
          branch.vertices.emplace_back(-1, triplet.first); // triplet.first == globalMin
        }
      }
      this->printMsg(
        "Built up maxima vectors", 0, maximaVectorTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

      ttk::Timer mgTimer;
#pragma omp parallel for num_threads(this->threadNumber_)
      for(size_t i = 0; i < branches.size(); i++) {
        auto vect = &branches[i].vertices;
        std::sort(vect->begin(), vect->end(), std::greater<>());
      }
      this->printMsg("Built up mergetree", 0, mgTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

      return 1;
    }
    template <typename triangulationType>
    int constructSegmentation(
      ttk::SimplexId *segmentation,
      const std::vector<Branch> &branches,
      const ttk::SimplexId *order,
      ttk::SimplexId *descendingManifold,
      ttk::SimplexId *tempArray,
      const triangulationType *triangulation) {

      ttk::Timer segmentationTimer;
      auto nVertices = triangulation->getNumberOfVertices();
      const auto &trunkSaddle
        = branches[branches.size() - 1]
            .vertices[branches[branches.size() - 1].vertices.size() - 2];
#pragma omp parallel for num_threads(this->threadNumber_)
      for(ttk::SimplexId i = 0; i < nVertices; i++) {
        auto orderForVertex = order[i];
        if(orderForVertex <= trunkSaddle.first){
          segmentation[i] = trunkSaddle.second;
          continue;
        }
        auto maximum
          = tempArray[descendingManifold[i]]; // global maximum id, anywhere in
                                              // 0 - nVertices
        // auto localMax = tempArray[maximum]; // local maximum id from 0 -
        // nMaxima
        Branch const *cBranch = &branches[maximum];
        auto lowestOrder = (*(cBranch->vertices.rbegin())).first;
        while(lowestOrder >= orderForVertex) { // finding the branch on which we are
          cBranch = cBranch->parentBranch;
          lowestOrder = (*(cBranch->vertices.rbegin())).first;
        }
        auto vect = &cBranch->vertices;
        auto lower = std::lower_bound(
          vect->rbegin(), vect->rend(), std::make_pair(orderForVertex, i));
        segmentation[i] = (*lower).second;
      }
      this->printMsg(
        "Finished phase 3 of segmentation: ", 1, segmentationTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

      return 1;
    }

    template <typename triangulationType>
    int
      findAscPaths(std::vector<std::array<ttk::SimplexId, 15>> &saddleTriplets,
                   std::vector<ttk::SimplexId> &maximaLocalToGlobal,
                   std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
                   ttk::SimplexId *maxima,
                   ttk::SimplexId *saddles,
                   const ttk::SimplexId *order,
                   ttk::SimplexId *descendingManifold,
                   ttk::SimplexId *tempArray,
                   std::vector<ttk::SimplexId> &saveGlobalIds,
                   const triangulationType *triangulation,
                   ttk::SimplexId nMaxima,
                   ttk::SimplexId nSaddles) {
      // construct the maximumLists for each saddle, the maxima which can be reached from this saddle
      // and the pathLists for each maximum, the saddles points which can reach this maximum
      ttk::Timer buildTimer;
      // sort the maxima and saddles by their order, the maxima ascending, the
      // saddles descending
      TTK_PSORT(this->threadNumber_, maxima, maxima + nMaxima,
                [&](ttk::SimplexId p1, ttk::SimplexId p2) {
                  return (order[p1] < order[p2]);
                });

      TTK_PSORT(this->threadNumber_, saddles, saddles + nSaddles,
                [&](ttk::SimplexId p1, ttk::SimplexId p2) {
                  return (order[p1] > order[p2]);
                });

#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp for nowait
        for(ttk::SimplexId i = 0; i < nMaxima; i++) {
          maximaLocalToGlobal[i] = maxima[i];
          saveGlobalIds[i] = tempArray[maxima[i]];
          tempArray[maxima[i]] = i;
        }

#pragma omp for nowait
        for(ttk::SimplexId i = 0; i < nSaddles; i++) {
          saddlesLocalToGlobal[i] = saddles[i];
        }
      }
      this->printMsg("Finished sorting and building the normalization arrays",
                     0, buildTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                     ttk::debug::Priority::DETAIL);
      ttk::Timer tripletTimer;

#pragma omp parallel for num_threads(this->threadNumber_)
      for(ttk::SimplexId i = 0; i < nSaddles; i++) {
        const auto &gId = saddles[i];
        const auto &nNeighbors = triangulation->getVertexNeighborNumber(gId);
        auto &triplet = saddleTriplets[i];
        auto &thisOrder = order[gId];
        ttk::SimplexId neighborId = 0;
        triplet[14] = 0;
        for(int j = 0; j < nNeighbors; j++) {
          triangulation->getVertexNeighbor(gId, j, neighborId);
          //  get the manifold result for this neighbor
          if(order[neighborId] > thisOrder) {
            // this->printMsg("maximum " +
            // std::to_string(descendingManifold[neighborId]) + " reachable from
            // neighbor " + std::to_string(neighborId) + " for saddle " +
            // std::to_string(gId));
            triplet[triplet[14]] = tempArray[descendingManifold[neighborId]];
            triplet[14]++;
          }
        }
        sortAndRemoveUniques(triplet);
      }

      this->printMsg("Finished building the saddleTriplets", 0,
                     tripletTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                     ttk::debug::Priority::DETAIL);
      return 1;
    }

    template <class triangulationType>
    int computePairs(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
      std::vector<Branch> &branches,
      ttk::SimplexId *segmentation,
      const ttk::SimplexId *minimaIds,
      ttk::SimplexId *saddle2Ids,
      ttk::SimplexId *maximaIds,
      ttk::SimplexId *descendingManifold,
      ttk::SimplexId *tempArray,
      const ttk::SimplexId *order,
      const triangulationType *triangulation,
      ttk::SimplexId &nMinima,
      ttk::SimplexId &nSaddle2,
      ttk::SimplexId &nMaxima) {

      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
      // print input parameters in table format
      this->printMsg(
        {{"#Threads", std::to_string(this->threadNumber_)},
         {"#Vertices", std::to_string(triangulation->getNumberOfVertices())}});
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing extremum pairs",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_);

        ttk::Timer preTimer;
        this->printMsg("Allocating memory", 0, ttk::debug::LineMode::REPLACE);
        persistencePairs.resize(nMaxima);
        branches.resize(nMaxima);
        std::vector<ttk::SimplexId> saveGlobalIds(nMaxima);
        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> maximaTriplets(
          nMaxima);
        std::vector<std::array<ttk::SimplexId, 15>> saddleTriplets(nSaddle2);
        std::vector<ttk::SimplexId> maximaLocalToGlobal(nMaxima);
        std::vector<ttk::SimplexId> saddlesLocalToGlobal(nSaddle2);
        this->printMsg("Allocating memory", 1, preTimer.getElapsedTime(),
                       this->threadNumber_);

        ttk::Timer ascTimer;
        ttk::SimplexId globalMin = 0;
#pragma omp parallel for num_threads(this->threadNumber_)
        for(ttk::SimplexId min = 0; min < nMinima; min++) {
          if(order[minimaIds[min]] == 0) {
            globalMin = minimaIds[min];
          }
        }
        this->printMsg(
          "Starting with findAscPaths", 0, ttk::debug::LineMode::REPLACE);
        findAscPaths<triangulationType>(
          saddleTriplets, maximaLocalToGlobal, saddlesLocalToGlobal, maximaIds,
          saddle2Ids, order, descendingManifold, tempArray, saveGlobalIds,
          triangulation, nMaxima, nSaddle2);
        this->printMsg("Finished with findAscPaths", 1, ascTimer.getElapsedTime(), this->threadNumber_);

        ttk::Timer pairTimer;
        this->printMsg(
          "Starting with PersistencePairs", 0, ttk::debug::LineMode::REPLACE);
        constructPersistencePairs(persistencePairs, maximaTriplets,
                                  saddleTriplets, maximaLocalToGlobal,
                                  saddlesLocalToGlobal, globalMin);
        this->printMsg("Finished with PersistencePairs", 1, pairTimer.getElapsedTime(), this->threadNumber_);

        ttk::Timer mergeTreeTimer;
        this->printMsg("Starting with mergetree computation", 0, ttk::debug::LineMode::REPLACE);
        constructMergeTree(branches, maximaTriplets,
                           maximaLocalToGlobal, saddlesLocalToGlobal, order);
        this->printMsg("Finished with mergetree computation", 1, mergeTreeTimer.getElapsedTime(), this->threadNumber_);

        ttk::Timer segmentationTimer;
        this->printMsg("Starting with mergetree segmentation", 0, ttk::debug::LineMode::REPLACE);
        constructSegmentation<triangulationType>(
          segmentation, branches,
          order, descendingManifold,
          tempArray, triangulation);
        this->printMsg("Finished mergetree segmentation", 1,
                       segmentationTimer.getElapsedTime(), this->threadNumber_);
//  print the progress of the current subprocedure with elapsed time
// transform the manifold back
#pragma omp parallel for num_threads(this->threadNumber_)
        for(ttk::SimplexId i = 0; i < nMaxima; i++) {
          tempArray[maximaIds[i]] = saveGlobalIds[i];
        }

        this->printMsg("Computing extremum pairs",
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

  }; // PairExtrema class

} // namespace ttk
