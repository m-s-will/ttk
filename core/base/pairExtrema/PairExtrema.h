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
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal) {
      int step = 0;
      bool changed = true;
      // std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
      // saddleMaxPairs{};
      std::vector<ttk::SimplexId> maximumPointer(maximaLocalToGlobal.size());
      std::iota(std::begin(maximumPointer), std::end(maximumPointer), 0);
      ttk::SimplexId nrOfPairs = 0;
      ttk::SimplexId globalMax = maximaLocalToGlobal.size() - 1;
      while(changed) {
        ttk::Timer stepTimer;
        changed = false;
        this->printMsg(ttk::debug::Separator::L2, ttk::debug::Priority::DETAIL);
        this->printMsg("Running step " + std::to_string(step)
                         + ", Pairs: " + std::to_string(nrOfPairs),
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

#pragma omp parallel for schedule(guided) num_threads(this->threadNumber_)
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

        // locks werden teuer wenn es nur noch wenige Maxima gibt
        // eher alle eigenes largestSaddlesForMax und dann reduction?

        this->printMsg("Finished building largestSaddlesForMax", 0.2,
                       buildTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);
        ttk::Timer pairTimer;

#pragma omp parallel for schedule(guided) num_threads(this->threadNumber_)
        for(size_t i = 0; i < largestSaddlesForMax.size(); i++) {
          ttk::SimplexId maximum = i;
          auto largestSaddle = largestSaddlesForMax[maximum];
          if(largestSaddle < (ttk::SimplexId)triplets.size()) {
            if(*(triplets.at(largestSaddle).begin()) == maximum) {
              changed = true;
#pragma omp atomic
              nrOfPairs++;
              pairs[maximum]
                = std::make_pair(saddlesLocalToGlobal[largestSaddle],
                                 maximaLocalToGlobal[maximum]);
              auto largestMax = *(
                triplets.at(largestSaddle)
                  .rbegin()); // largest maximum reachable from paired saddle
              maximumPointer[maximum] = largestMax;
              branches[maximum] = (std::make_pair(largestSaddle, largestMax));
            }
          }
        }

        this->printMsg("Finished finding pairs and swapping pointers", 0.5,
                       pairTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

        ttk::Timer compressTimer;
        // use pathcompression on the maximumPointer
        int compressionStep = 0;
        bool same = false;
        std::vector<ttk::SimplexId> nextMaximumPointer(
          maximaLocalToGlobal.size());
        while(!same) {
          same = true;
          if(compressionStep % 2 == 0) {
#pragma omp parallel for num_threads(this->threadNumber_)
            for(size_t i = 0; i < maximumPointer.size(); i++) {
              int nextPointer = maximumPointer[maximumPointer[i]];
              if(nextPointer != nextMaximumPointer[i]) {
                nextMaximumPointer[i] = nextPointer;
                same = false;
              }
            }
          } else {
#pragma omp parallel for num_threads(this->threadNumber_)
            for(size_t i = 0; i < maximumPointer.size(); i++) {
              int nextPointer = nextMaximumPointer[nextMaximumPointer[i]];
              if(nextPointer != maximumPointer[i]) {
                maximumPointer[i] = nextPointer;
                same = false;
              }
            }
          }
          compressionStep++;
        }
        this->printMsg("Did pathcompression on maximumpointer in "
                         + std::to_string(compressionStep) + " steps",
                       0.8, compressTimer.getElapsedTime(),
                       ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);
        ttk::Timer replaceTimer;
// replace the values with their maximumPointers, delete the saddle if max and
// min are the same
#pragma omp parallel num_threads(this->threadNumber_)
        {
#pragma omp for schedule(guided) nowait
          for(size_t i = 0; i < triplets.size(); i++) {
            std::set<ttk::SimplexId> newSet{};
            for(auto &toReplace : triplets[i]) {
              newSet.emplace(maximumPointer[toReplace]);
            }
            if(newSet.size() == 1) {
              triplets[i].clear();
              // maybe delete completely? and shift localtoglobal
            } else {
              triplets[i] = newSet;
            }
          }
        }

        this->printMsg("Replaced values and deleted unnecessary saddles.", 1,
                       replaceTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);

        this->printMsg("Finished step " + std::to_string(step), 1,
                       stepTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                       ttk::debug::Priority::DETAIL);
        step++;
      }
      return 1;
    }

    int constructMergeTree(
      std::vector<std::vector<ttk::SimplexId>> &maximaVectors,
      std::vector<std::vector<ttk::SimplexId>> &maximaOrders,
      std::vector<ttk::SimplexId> &searchTree,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &mergeTree,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
        &branches, // branches[max] = (saddle, largestMax)
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      const ttk::SimplexId *order,
      const ttk::SimplexId &globalMin) {

      searchTree.resize(maximaLocalToGlobal.size()
                        + saddlesLocalToGlobal.size());
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
      this->printMsg("Compressed branches", 0, compressTimer.getElapsedTime());
      ttk::Timer maximaVectorTimer;
      for(size_t i = 0; i < branches.size() - 1; i++) {
        auto branch = branches[i];
        auto max = i;
        auto saddle = branch.first;
        auto largestMax = branch.second;
        maximaVectors[max].emplace_back(saddle);
        maximaVectors[largestMax].emplace_back(saddle);
      }
      this->printMsg(
        "Built up maxima vectors", 0, maximaVectorTimer.getElapsedTime());

      ttk::SimplexId minSaddle = -1;
      // omp_lock_t minSaddleLock;
      // omp_init_lock(&minSaddleLock);
      ttk::Timer mgTimer;
      //#pragma omp parallel for schedule(guided)
      //num_threads(this->threadNumber_)
      for(size_t i = 0; i < maximaVectors.size(); i++) {
        auto vect = maximaVectors[i];
        std::sort(vect.begin(),
                  vect.end()); // we want descending saddles and the lower the
                               // saddle id value, the higher the saddle scalar
        minSaddle = std::max(minSaddle, vect[vect.size() - 1]);
        mergeTree.emplace_back(std::make_pair(
          maximaLocalToGlobal[i], saddlesLocalToGlobal[vect[0]]));
        searchTree[i] = vect[0];
        for(size_t j = 0; j < vect.size() - 1; j++) {
          if(vect[j] != vect[j + 1]) {
            mergeTree.emplace_back(
              std::make_pair(saddlesLocalToGlobal[vect[j]],
                             saddlesLocalToGlobal[vect[j + 1]]));
            searchTree[maximaLocalToGlobal.size() + vect[j]] = vect[j + 1];
          }
        }
      }
      this->printMsg("Built up mergetree", 0, mgTimer.getElapsedTime());

#pragma omp parallel for schedule(dynamic, 4) num_threads(this->threadNumber_)
      for(size_t i = 0; i < maximaVectors.size(); i++) {
        auto vect = maximaVectors[i];
        maximaOrders[i].resize(vect.size());
        ;
        for(size_t j = 0; j < vect.size(); j++) {
          maximaOrders[i][j] = order[saddlesLocalToGlobal[vect[j]]];
        }
      }
      // connect minimal saddle and global min
      mergeTree.emplace_back(
        std::make_pair(saddlesLocalToGlobal[minSaddle], globalMin));
      searchTree[maximaLocalToGlobal.size() + minSaddle] = globalMin;
      return 1;
    }
    template <typename triangulationType>
    int constructSegmentation(
      ttk::SimplexId *segmentation,
      std::vector<ttk::SimplexId> &searchTree,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &mergeTree,
      const std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &branches,
      const std::vector<std::vector<ttk::SimplexId>> &maximaVectors,
      const std::vector<std::vector<ttk::SimplexId>> &maximaOrders,
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      const ttk::SimplexId *order,
      const ttk::SimplexId *ascendingManifold,
      const ttk::SimplexId *tempArray,
      const triangulationType *triangulation) {

      // everwhere, where is no -1, is definitely correct, some vertices are not
      // reached correctly
      auto minSaddle = mergeTree[mergeTree.size() - 1].first;

      std::vector<int> marked(triangulation->getNumberOfVertices(), 0);
      // this->printMsg("#Branches in the MergeTree: " +
      // std::to_string(mergeTree.size()));
      ttk::SimplexId minOrder = order[minSaddle];
      ttk::Timer phase1Timer;
// phase 1: everything below the minimum saddle belongs to it
#pragma omp parallel for num_threads(this->threadNumber_)
      for(size_t i = 0; i < triangulation->getNumberOfVertices(); i++) {
        if(order[i] <= minOrder) {
          segmentation[i] = minSaddle;
          marked[i] = 1;
        }
      }
      auto nrPhase1 = std::count(marked.begin(), marked.end(), 1);
      this->printMsg("Finished phase 1 of segmentation with "
                       + std::to_string(nrPhase1) + " marked vertices.",
                     0.33, phase1Timer.getElapsedTime());
      ttk::Timer phase2Timer;
// phase 2: flood fill
#pragma omp parallel for schedule(dynamic, 4) num_threads(this->threadNumber_)
      for(size_t i = 0; i < mergeTree.size(); i++) {
        auto branch = mergeTree[i];
        auto mark = branch.first;
        // this->printMsg("Working on branch " + std::to_string(branch.first) +
        // "-" + std::to_string(branch.second));
        auto top = order[mark];
        auto bottom = order[branch.second];
        // this->printMsg("Top " + std::to_string(top) + ", bottom " +
        // std::to_string(bottom));
        std::deque<ttk::SimplexId> queue{mark};
        marked[mark] = 1;
        while(!queue.empty()) {
          auto element = queue.front();
          queue.pop_front();
          segmentation[element] = mark;
          // this->printMsg("Id " + std::to_string(element) + " belongs to
          // branch " + std::to_string(mark));
          auto nNeighbors = triangulation->getVertexNeighborNumber(element);
          ttk::SimplexId neighborId;
          for(int j = 0; j < nNeighbors; j++) {
            triangulation->getVertexNeighbor(element, j, neighborId);
            if(bottom < order[neighborId] && order[neighborId] < top
               && marked[neighborId] == 0) {
              queue.push_back(neighborId);
              marked[neighborId] = 1;
            }
          }
        }
      }
      auto nrPhase2 = std::count(marked.begin(), marked.end(), 1) - nrPhase1;
      this->printMsg("Finished phase 2 (flood fill) of segmentation with "
                       + std::to_string(nrPhase2) + " marked vertices.",
                     0.66, phase2Timer.getElapsedTime());

      this->printMsg("dealing with "
                     + std::to_string(triangulation->getNumberOfVertices()
                                      - nrPhase1 - nrPhase2)
                     + " remaining vertices");
      ttk::Timer phase3Timer;
// phase 3: all unmarked vertices after flood fill need tree traversal
#pragma omp parallel for schedule(dynamic, 4) num_threads(this->threadNumber_)
      for(size_t i = 0; i < triangulation->getNumberOfVertices(); i++) {
        if(marked[i] == 0){
          auto maximum
            = ascendingManifold[i]; // global, we need local from tempArray
          auto thisOrder = order[i];
          auto prev = tempArray[maximum];
          auto vect = maximaVectors[prev];
          auto orders = maximaOrders[prev];
          bool useSearchTree = false;
          int maxVectorsteps = 0;
          while(*(orders.rbegin()) > thisOrder
                && branches[prev].second != prev) {
            prev = branches[prev].second;
            vect = maximaVectors[prev];
            orders = maximaOrders[prev];
          }
          /*if(branches[prev].second == prev)
          {
            this->printMsg("No further branching for order " +
          std::to_string(thisOrder) + ", max: " +
          std::to_string(maximaLocalToGlobal[prev]) + ", highest order: " +
          std::to_string(orders[0]) + ", lowest order: " +
          std::to_string(*(orders.rbegin())));
          }*/
          auto lower
            = std::lower_bound(orders.begin(), orders.end(), thisOrder);
          // this->printMsg("orders.size: " + std::to_string(orders.size()) + ",
          // thisOrder: " + std::to_string(thisOrder) + ", lowerOrder: " +
          // std::to_string(*lower) + ", lowerOrder-1: " +
          // std::to_string(*(lower-1)));
          prev = vect[std::distance(orders.begin(), lower) - 1];

          segmentation[i] = saddlesLocalToGlobal[prev];
          marked[i] = 1;
        }
      }
      this->printMsg(
        "Finished phase 3  of segmentation", 1, phase3Timer.getElapsedTime());

      return 1;
    }

    template <typename triangulationType>
    int findAscPaths(
      std::vector<std::set<ttk::SimplexId>> &triplets,
      std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &maxima,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &saddles,
      const ttk::SimplexId *order,
      const ttk::SimplexId *ascendingManifold,
      ttk::SimplexId *tempArray,
      const triangulationType *triangulation) {
      // construct the maximumLists for each saddle, the maxima which can be reached from this saddle
      // and the pathLists for each maximum, the saddles points which can reach this maximum
      ttk::Timer buildTimer;
      // sort the maxima and saddles by their order, the maxima ascending, the
      // saddles descending
      TTK_PSORT(this->threadNumber_, maxima.begin(), maxima.end(),
                [](std::pair<ttk::SimplexId, ttk::SimplexId> p1,
                   std::pair<ttk::SimplexId, ttk::SimplexId> p2) {
                  return (p1.second < p2.second);
                });

      TTK_PSORT(this->threadNumber_, saddles.begin(), saddles.end(),
                [](std::pair<ttk::SimplexId, ttk::SimplexId> p1,
                   std::pair<ttk::SimplexId, ttk::SimplexId> p2) {
                  return (p1.second > p2.second);
                });

#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp for nowait
        for(size_t i = 0; i < maxima.size(); i++) {
          maximaLocalToGlobal[i] = maxima[i].first;
          tempArray[maxima[i].first] = i;
        }

#pragma omp for nowait
        for(size_t i = 0; i < saddles.size(); i++) {
          saddlesLocalToGlobal[i] = saddles[i].first;
        }
      }
      this->printMsg("Finished sorting and building the normalization arrays",
                     0, buildTimer.getElapsedTime(), ttk::debug::LineMode::NEW,
                     ttk::debug::Priority::DETAIL);

      ttk::Timer tripletTimer;

#pragma omp parallel for schedule(guided) num_threads(this->threadNumber_)
      for(size_t i = 0; i < saddles.size(); i++) {
        auto gId = saddles[i].first;
        auto nNeighbors = triangulation->getVertexNeighborNumber(gId);
        ttk::SimplexId neighborId;
        for(int j = 0; j < nNeighbors; j++) {
          triangulation->getVertexNeighbor(gId, j, neighborId);
          // get the manifold result for this neighbor
          // problematic if manifold is dense and not sparse, because we need
          // the id of the point to which it is ascending, not the id of the
          // segmentation
          if(order[neighborId] > saddles[i].second) {
            triplets[i].emplace(tempArray[ascendingManifold[neighborId]]);
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
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &mergeTree,
      ttk::SimplexId *segmentation,
      const ttk::SimplexId *minimaIds,
      const ttk::SimplexId *saddle2Ids,
      const ttk::SimplexId *maximaIds,
      const ttk::SimplexId *ascendingManifold,
      ttk::SimplexId *tempArray,
      const ttk::SimplexId *order,
      const triangulationType *triangulation,
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

        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> saddles;
        saddles.resize(nSaddle2);
        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> maxima;
        maxima.resize(nMaxima);
        ttk::Timer preTimer;
        // ttk::SimplexId maxTime = 0;
        // auto t0 = timeNow();
#pragma omp parallel for num_threads(this->threadNumber_)
        for(ttk::SimplexId i = 0; i < nMaxima; i++) {
          auto gId = maximaIds[i];
          auto thisOrder = order[gId];
          maxima[i] = (std::make_pair(gId, thisOrder));
        }

#pragma omp parallel for num_threads(this->threadNumber_)
        for(ttk::SimplexId i = 0; i < nSaddle2; i++) {
          auto gId = saddle2Ids[i];
          auto thisOrder = order[gId];
          saddles[i] = std::make_pair(gId, thisOrder);
        }

        persistencePairs.resize(maxima.size());
        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> branches;
        std::vector<ttk::SimplexId> searchTree;
        branches.resize(maxima.size());
        std::vector<std::vector<ttk::SimplexId>> maximaVectors(maxima.size());
        std::vector<std::vector<ttk::SimplexId>> maximaOrders(maxima.size());

        std::vector<std::set<ttk::SimplexId>> triplets(saddles.size());
        std::vector<ttk::SimplexId> maximaLocalToGlobal(maxima.size());
        std::vector<ttk::SimplexId> saddlesLocalToGlobal(saddles.size());
        this->printMsg(
          "Finished with Preprocessing, starting with findAscPaths", 0,
          preTimer.getElapsedTime());

        ttk::Timer ascTimer;
        findAscPaths<triangulationType>(
          triplets, maximaLocalToGlobal, saddlesLocalToGlobal, maxima, saddles,
          order, ascendingManifold, tempArray, triangulation);

        this->printMsg(
          "Finished with findAscPaths, starting with PersistencePairs", 0,
          ascTimer.getElapsedTime());
        ttk::Timer pairTimer;
        constructPersistencePairs(persistencePairs, branches, triplets,
                                  maximaLocalToGlobal, saddlesLocalToGlobal);
        this->printMsg(
          "Finished constructing persistence pairs, starting with mergetree", 0,
          pairTimer.getElapsedTime());

        ttk::Timer mergeTreeTimer;
        constructMergeTree(maximaVectors, maximaOrders, searchTree, mergeTree,
                           branches, maximaLocalToGlobal, saddlesLocalToGlobal,
                           order, minimaIds[0]);
        this->printMsg("Finished constructing mergetree", 0,
                       mergeTreeTimer.getElapsedTime());
        ttk::Timer segmentationTimer;
        // constructSegmentation<triangulationType>(
        //   segmentation, searchTree, mergeTree, branches, maximaVectors,
        //   maximaOrders, maximaLocalToGlobal, saddlesLocalToGlobal, order,
        //   ascendingManifold, tempArray, triangulation);
        this->printMsg("Finished mergetree segmentation", 0,
                       segmentationTimer.getElapsedTime());

        // the global max is always in the last position of the
        // maximaLocalToGlobal vector and needs to connect with the global
        // minimum
        persistencePairs[persistencePairs.size() - 1]
          = std::make_pair(minimaIds[0], maxima[maxima.size() - 1].first);
        // print the progress of the current subprocedure with elapsed time
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
