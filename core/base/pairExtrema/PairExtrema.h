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
//#include <valgrind/callgrind.h>
//#include <gperftools/profiler.h>

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
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &branches,
      std::vector<std::set<ttk::SimplexId>> &triplets,
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
          maximaLocalToGlobal.size(), triplets.size());
        std::vector<omp_lock_t> maximaLocks(maximaLocalToGlobal.size());
        for(size_t i = 0; i < maximaLocks.size(); i++) {
          omp_init_lock(&maximaLocks[i]);
        }

#pragma omp parallel num_threads(this->threadNumber_)
        {
#pragma omp for schedule(guided)
          for(ttk::SimplexId i = 0; i < (ttk::SimplexId)triplets.size(); i++) {
            auto &maxList = triplets.at(i);
            ttk::SimplexId temp;
            for(auto &max : maxList) {
              // TODO if OMP 5.1 is widespread: use omp atomic compare
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
#pragma omp single
          this->printMsg("Finished building largestSaddlesForMax", 0.2,
                         buildTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                         ttk::debug::Priority::DETAIL);
          ttk::Timer pairTimer;

          std::vector<SimplexId> lActiveMaxima;
          lActiveMaxima.reserve(maximaLocalToGlobal.size());

#pragma omp for schedule(guided)
          for(size_t i = 0; i < largestSaddlesForMax.size(); i++) {
            ttk::SimplexId maximum = i;
            auto largestSaddle = largestSaddlesForMax[maximum];
            if(largestSaddle < (ttk::SimplexId)triplets.size()) {
              if(*(triplets.at(largestSaddle).begin()) == maximum) {
                changed = true;
                //#pragma omp atomic
                // nrOfPairs++;
                pairs[maximum]
                  = std::make_pair(saddlesLocalToGlobal[largestSaddle],
                                   maximaLocalToGlobal[maximum]);
                auto largestMax = *(
                  triplets.at(largestSaddle)
                    .rbegin()); // largest maximum reachable from paired saddle
                maximumPointer[maximum] = largestMax;
                branches[maximum] = (std::make_pair(largestSaddle, largestMax));
                lActiveMaxima.push_back(maximum);
              }
            }
          }
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
#pragma omp for schedule(guided)
          for(size_t i = 0; i < triplets.size(); i++) {
            std::set<ttk::SimplexId> newSet{};
            for(auto &toReplace : triplets[i]) {
              newSet.emplace(maximumPointer[toReplace]); // expensive
            }
            if(newSet.size() == 1) {
              triplets[i].clear();
            } else {
              triplets[i] = newSet;
            }
          }
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
      branches[branches.size() - 1] = std::make_pair(
        globalMin,
        branches.size() - 1); // branches[globalmax] = (globalmin, globalmax);

      return 1;
    }

    int constructMergeTree(
      std::vector<Branch> &maximaVectors,
      std::vector<std::vector<ttk::SimplexId>> &maximaOrders,
      std::vector<Branch> &mergeTree,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
        &branches, // branches[max] = (saddle, largestMax)
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      const ttk::SimplexId *order,
      const ttk::SimplexId &globalMin) {

      branches[branches.size() - 1].second = branches.size() - 1;
      // compress the branches
      bool same = false;
      ttk::Timer compressTimer;

      while(!same) {
        same = true;
        for(size_t i = 0; i < branches.size(); i++) {
          auto nextBranch = branches[branches[i].second];
          if(nextBranch.second != branches[i].second
             && nextBranch.first
                  < branches[i].first) { // we need to follow along larger
                                         // saddles to the largestMax
            branches[i].second = nextBranch.second;
            same = false;
          }
        }
      }
      this->printMsg("Compressed branches", 0, compressTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);
      ttk::Timer maximaVectorTimer;

      for(ttk::SimplexId b = 0; b < branches.size(); b++) {
        auto &triplet = branches[b]; /// branches[max] = (saddle, biggerMax wo
                                     /// angedockt wird)
        auto &branch = maximaVectors[b]; /// branches[max] = (saddle, biggerMax
                                         /// wo angedockt wird)
        auto branchMaxId = maximaLocalToGlobal[b];
        branch.vertices.emplace_back(order[branchMaxId], branchMaxId);

        auto parent = triplet.second;
        if(parent != b) {
          auto saddle = saddlesLocalToGlobal[triplet.first];
          auto orderForSaddle = order[saddle];
          branch.vertices.emplace_back(orderForSaddle, saddle);
          branch.parentBranch
            = &maximaVectors[parent]; // maximaVectors[mainBranch].parentBranch
                                      // = maximaVectors[mainBranch];
          maximaVectors[parent].vertices.emplace_back(orderForSaddle, saddle);
        } else {
          auto globalMin = triplet.first;
          branch.vertices.emplace_back(-1, globalMin);
        }
      }
      this->printMsg(
        "Built up maxima vectors", 0, maximaVectorTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

      ttk::SimplexId minSaddle = -1;
      ttk::Timer mgTimer;
      /*#pragma omp parallel for schedule(dynamic, 4)
         num_threads(this->threadNumber_) for(size_t i = 0; i <
         maximaVectors.size(); i++) {
              // we want descending saddles and the lower the saddle id value,
         the
              // higher the saddle scalar
              auto vect = &maximaVectors[i];
              std::sort(maximaVectors[i].begin(), maximaVectors[i].end());
              maximaOrders[i].resize(vect->size());
              for(size_t j = 0; j < vect->size(); j++) {
                maximaOrders[i][j] = order[saddlesLocalToGlobal[vect->at(j)]];
              }
            }*/
#pragma omp parallel for num_threads(this->threadNumber_)
      for(size_t i = 0; i < maximaVectors.size(); i++) {
        auto vect = &maximaVectors[i].vertices;
        std::sort(vect->begin(), vect->end(), std::greater<>());
      }
      this->printMsg("Built up mergetree", 0, mgTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

      return 1;
    }
    template <typename triangulationType>
    int constructSegmentation(
      ttk::SimplexId *segmentation,
      const std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &branches,
      const std::vector<Branch> &mergeTree,
      const std::vector<Branch> &maximaVectors,
      const std::vector<std::vector<ttk::SimplexId>> &maximaOrders,
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      const ttk::SimplexId *order,
      const ttk::SimplexId *descendingManifold,
      ttk::SimplexId *tempArray,
      const triangulationType *triangulation) {

      ttk::SimplexId globalMax = maximaLocalToGlobal.size() - 1;

      // auto minSaddle = mergeTree[mergeTree.size() - 1].first;
      // ttk::SimplexId minOrder = order[minSaddle];

      ttk::Timer phase3Timer;
      auto nVertices = triangulation->getNumberOfVertices();
      // std::vector<ttk::SimplexId>
      // orderInverse(triangulation->getNumberOfVertices()); phase 3: all
      // unmarked vertices after flood fill need tree traversal

      ttk::Timer trunkTimer;
      // maximaVect[0] (5, 3, 2, 1)
      // maximaOrders[0] ( order(saddle5), order(saddle3))
#pragma omp parallel for num_threads(this->threadNumber_)
      for(ttk::SimplexId i = 0; i < nVertices; i++) {
        auto maximum = descendingManifold[i]; // global maximum id, anywhere in
                                              // 0 - nVertices
        auto orderForVertex = order[i];
        auto localMax = tempArray[maximum]; // local maximum id from 0 - nMaxima
        Branch const *cBranch = &maximaVectors[localMax];
        auto lowestOrder = (*(cBranch->vertices.rbegin())).first;
        while(lowestOrder
              >= orderForVertex) { // finding the branch on which we are
          // this->printMsg(std::to_string(lowestOrder));
          // this->printMsg("x");
          cBranch = cBranch->parentBranch;
          // this->printMsg("y");
          lowestOrder = (*(cBranch->vertices.rbegin())).first;
          // this->printMsg("z");
        }
        // auto orders = &maximaOrders[prev];
        auto vect = &cBranch->vertices;
        auto lower = std::lower_bound(
          vect->rbegin(), vect->rend(), std::make_pair(orderForVertex, i));
        // auto saddleId = vect->at((lower) - orders->rbegin());
        segmentation[i] = (*lower).second;
      }
      this->printMsg(
        "Finished phase 3 of segmentation: ", 1, trunkTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

      return 1;
    }

    template <typename triangulationType>
    int findAscPaths(std::vector<std::set<ttk::SimplexId>> &triplets,
                     std::vector<ttk::SimplexId> &maximaLocalToGlobal,
                     std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
                     ttk::SimplexId *maxima,
                     ttk::SimplexId *saddles,
                     const ttk::SimplexId *order,
                     const ttk::SimplexId *descendingManifold,
                     ttk::SimplexId *tempArray,
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
        for(size_t i = 0; i < nMaxima; i++) {
          maximaLocalToGlobal[i] = maxima[i];
          tempArray[maxima[i]] = i;
          // descendingManifold[maxima[i]] = i
          // go once over descendingManifold and transfrom from global to local
        }

#pragma omp for nowait
        for(size_t i = 0; i < nSaddles; i++) {
          saddlesLocalToGlobal[i] = saddles[i];
        }
      }
      this->printMsg("Finished sorting and building the normalization arrays",
                     0, buildTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                     ttk::debug::Priority::DETAIL);
      ttk::Timer tripletTimer;

#pragma omp parallel for num_threads(this->threadNumber_)
      for(size_t i = 0; i < nSaddles; i++) {
        const auto &gId = saddles[i];
        const auto &nNeighbors = triangulation->getVertexNeighborNumber(gId);
        auto triplet = &triplets[i];
        ttk::SimplexId neighborId;
        for(int j = 0; j < nNeighbors; j++) {
          triangulation->getVertexNeighbor(gId, j, neighborId);
          //  get the manifold result for this neighbor
          //  problematic if manifold is dense and not sparse, because we need
          //  the id of the point to which it is ascending, not the id of the
          //  segmentation
          if(order[neighborId] > order[gId]) {
            triplet->emplace(tempArray[descendingManifold[neighborId]]);
          }
        }
      }

      this->printMsg("Finished building the triplets", 0,
                     tripletTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                     ttk::debug::Priority::DETAIL);
      return 1;
    }

    template <class triangulationType = ttk::AbstractTriangulation>
    int computePairs(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
      std::vector<Branch> &mergeTree,
      ttk::SimplexId *segmentation,
      const ttk::SimplexId *minimaIds,
      ttk::SimplexId *saddle2Ids,
      ttk::SimplexId *maximaIds,
      const ttk::SimplexId *descendingManifold,
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
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
        {"#Minima", std::to_string(nMinima)}
      });
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

        // std::vector<ttk::SimplexId> saddles(saddle2Ids, saddle2Ids +
        // nSaddle2); std::vector<ttk::SimplexId> maxima(maximaIds, maximaIds +
        // nMaxima);
        ttk::Timer preTimer;
        this->printMsg("Allocating memory", 0, ttk::debug::LineMode::REPLACE);
        persistencePairs.resize(nMaxima);
        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> branches(
          nMaxima);
        std::vector<Branch> maximaVectors(nMaxima);
        std::vector<std::vector<ttk::SimplexId>> maximaOrders(nMaxima);
        std::vector<std::set<ttk::SimplexId>> triplets(nSaddle2);
        std::vector<ttk::SimplexId> maximaLocalToGlobal(nMaxima);
        std::vector<ttk::SimplexId> saddlesLocalToGlobal(nSaddle2);
        this->printMsg("Allocating memory", 1, preTimer.getElapsedTime(),
                       this->threadNumber_);

        // CALLGRIND_START_INSTRUMENTATION;
        // ProfilerStart("pairextrema.prof");
        ttk::Timer ascTimer;
        this->printMsg(
          "Starting with findAscPaths", 0, ttk::debug::LineMode::REPLACE);
        findAscPaths<triangulationType>(
          triplets, maximaLocalToGlobal, saddlesLocalToGlobal, maximaIds,
          saddle2Ids, order, descendingManifold, tempArray, triangulation,
          nMaxima, nSaddle2);
        this->printMsg("Finished with findAscPaths", 1, ascTimer.getElapsedTime(), this->threadNumber_);

        ttk::Timer pairTimer;
        this->printMsg(
          "Starting with PersistencePairs", 0, ttk::debug::LineMode::REPLACE);
        constructPersistencePairs(persistencePairs, branches, triplets,
                                  maximaLocalToGlobal, saddlesLocalToGlobal,
                                  minimaIds[0]);
        this->printMsg("Finished with PersistencePairs", 1, pairTimer.getElapsedTime(), this->threadNumber_);

        ttk::Timer mergeTreeTimer;
        this->printMsg("Starting with mergetree computation", 0, ttk::debug::LineMode::REPLACE);
        constructMergeTree(maximaVectors, maximaOrders, mergeTree, branches,
                           maximaLocalToGlobal, saddlesLocalToGlobal, order,
                           minimaIds[0]);
        this->printMsg("Finished with mergetree computation", 1, mergeTreeTimer.getElapsedTime(), this->threadNumber_);

        ttk::Timer segmentationTimer;
        this->printMsg("Starting with mergetree segmentation", 0, ttk::debug::LineMode::REPLACE);
        constructSegmentation<triangulationType>(
          segmentation, branches, mergeTree, maximaVectors, maximaOrders,
          maximaLocalToGlobal, saddlesLocalToGlobal, order, descendingManifold,
          tempArray, triangulation);
        this->printMsg("Finished mergetree segmentation", 1, segmentationTimer.getElapsedTime(), this->threadNumber_);
        // ProfilerStop();
        // CALLGRIND_STOP_INSTRUMENTATION;
        // CALLGRIND_DUMP_STATS;
        //  print the progress of the current subprocedure with elapsed time

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
