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
      std::vector<std::vector<ttk::SimplexId>> &pairs,
      std::vector<std::set<ttk::SimplexId>> &triplets,
      std::vector<std::set<ttk::SimplexId>> &largestSaddlesForMax,
      std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      const ttk::SimplexId *order) {
      int step = 0;
      bool changed = true;
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> saddleMaxPairs{};
      // std::vector<ttk::SimplexId> maximumPointer(maximaLocalToGlobal.size());
      // std::iota(std::begin(maximumPointer), std::end(maximumPointer), 0);
      std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>>
        mergeTree{}; // tuples are saddle, smaller maximum, larger maximum to
                     // which it merges
      std::set<ttk::SimplexId> saddlesToDelete{};
      saddleMaxPairs.reserve(maximaLocalToGlobal.size());
      this->printMsg("Nr of Maxima: "
                     + std::to_string(maximaLocalToGlobal.size()));
      while(changed) {
        ttk::Timer stepTimer;
        changed = false;
        saddleMaxPairs.clear();
        saddlesToDelete.clear();
        this->printMsg("============================================");
        this->printMsg("Running step " + std::to_string(step)
                         + ", Pairs size: " + std::to_string(pairs.size()),
                       0, 0);
        ttk::Timer pairTimer;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
#endif
        {

#ifdef TTK_ENABLE_OPENMP
          std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
            saddleMaxPairs_priv{};
          std::vector<std::vector<ttk::SimplexId>> pairs_priv{};
#pragma omp for schedule(guided) nowait
#endif
          for(size_t i = 0; i < largestSaddlesForMax.size(); i++) {
            ttk::SimplexId maximum = i; // only the maximum id
            // this->printMsg("maximum: " + std::to_string(maximum));
            ttk::SimplexId largestSaddle
              = *(largestSaddlesForMax.at(maximum).begin());
            // this->printMsg("largestSaddle: " +
            // std::to_string(largestSaddle));
            if(largestSaddle != -1) {
              if(triplets.at(largestSaddle).size() > 0) {
                ttk::SimplexId maximumForSaddle
                  = *(triplets.at(largestSaddle)
                        .begin()); // the lowest maximum id and val for
                                   // the current saddle
                // we want to check if smallest maximum per saddle and largest
                // saddle per maximum match
                if(maximumForSaddle == maximum) {
                  changed = true;
#ifdef TTK_ENABLE_OPENMP
                  {
                    pairs_priv.emplace_back(
                      std::initializer_list<ttk::SimplexId>{
                        saddlesLocalToGlobal[largestSaddle],
                        order[saddlesLocalToGlobal[largestSaddle]],
                        maximaLocalToGlobal[maximumForSaddle],
                        order[maximaLocalToGlobal[maximumForSaddle]]});
                    saddleMaxPairs_priv.push_back(
                      std::make_pair(largestSaddle, maximum));
                  }
#else
                  {
                    pairs.emplace_back(std::initializer_list<ttk::SimplexId>{
                      saddlesLocalToGlobal[largestSaddle],
                      order[saddlesLocalToGlobal[largestSaddle]],
                      maximaLocalToGlobal[maximumForSaddle],
                      order[maximaLocalToGlobal[maximumForSaddle]]});
                    saddleMaxPairs.push_back(
                      std::make_pair(largestSaddle, maximum));
                  }
#endif
                }
              }
            }
          }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(deletions)
          {
            saddleMaxPairs.insert(saddleMaxPairs.end(),
                                  saddleMaxPairs_priv.begin(),
                                  saddleMaxPairs_priv.end());
            pairs.insert(pairs.end(), pairs_priv.begin(), pairs_priv.end());
          }
#endif
        }
        this->printMsg(
          "Finished finding pairs", 0.33, pairTimer.getElapsedTime());
        ttk::Timer swappingTimer;

        for(size_t i = 0; i < saddleMaxPairs.size(); i++) {
          auto pair = saddleMaxPairs[i];
          auto saddle = pair.first;
          auto max = pair.second;
          auto &maxSet
            = triplets.at(saddle); // set containing all the maxima reachable
                                   // from the saddle to be removed
          auto largestMax = *(maxSet.rbegin());
          // maximumPointer[max] = largestMax;
          mergeTree.push_back(std::make_tuple(saddle, max, largestMax));
          auto &saddleSet = largestSaddlesForMax.at(max);
          maxSet.erase(maxSet.begin());
          saddleSet.erase(saddleSet.begin());
          if(maxSet.size() < 2)
            saddlesToDelete.insert(saddle);
          // now we need to swap the pointers, all saddles which could reach
          // maximum pair.first, now can reach all the other maxima of saddle
          // we can find this by just checking largestSaddlesForMax for the
          // maximum
          for(auto &saddleId : saddleSet) {
            auto &triplet = triplets.at(saddleId);
            triplet.erase(max);
            // if the maximum can be reached from other saddles, we
            // need to replace it by the other maxima from the
            // chosen saddle
            triplet.insert(maxSet.begin(), maxSet.end());
            // triplet.insert(largestMax);
            if(triplet.size() < 2)
              saddlesToDelete.insert(saddleId);
          }
          for(auto &maxId : maxSet) {
            auto &saddleList = largestSaddlesForMax.at(maxId);
            // if the maximum can be reached from other saddles, we
            // need to replace it by the other maxima from the
            // chosen saddle
            saddleList.insert(saddleSet.begin(), saddleSet.end());
          }
        }

        this->printMsg("Finished swapping pointers for "
                         + std::to_string(saddleMaxPairs.size()) + " pairs.",
                       0.66, swappingTimer.getElapsedTime());
        ttk::Timer delTimer;

        for(auto saddle : saddlesToDelete) {
          auto &maxList = triplets.at(saddle);
          for(auto &m : maxList) {
            auto &saddleSet = largestSaddlesForMax.at(m);
            saddleSet.erase(saddle);
          }
          triplets.at(saddle).clear();
        }

        this->printMsg("Finished deletion of "
                         + std::to_string(saddlesToDelete.size())
                         + " irrelevant saddles.",
                       1, delTimer.getElapsedTime());
        this->printMsg("Finished step " + std::to_string(step), 1,
                       stepTimer.getElapsedTime());
        step++;
      }
      return 1;
    }

    template <typename triangulationType>
    int findAscPaths(
      std::vector<std::set<ttk::SimplexId>> &triplets,
      std::vector<std::set<ttk::SimplexId>> &largestSaddlesForMax,
      std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      std::vector<ttk::SimplexId> &saddlesLocalToGlobal,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &maxima,
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &saddles,
      const ttk::SimplexId *order,
      const ttk::SimplexId *ascendingManifold,
      const triangulationType *triangulation) {
      // construct the maximumLists for each saddle, the maxima which can be reached from this saddle
      // and the pathLists for each maximum, the saddles points which can reach this maximum
      ttk::Timer sortTimer;
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

      std::unordered_map<ttk::SimplexId, ttk::SimplexId> maximaGlobalToLocal{};
      std::unordered_map<ttk::SimplexId, ttk::SimplexId> saddlesGlobalToLocal{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
      {
        std::unordered_map<ttk::SimplexId, ttk::SimplexId>
          maximaGlobalToLocal_priv{};
        std::unordered_map<ttk::SimplexId, ttk::SimplexId>
          saddlesGlobalToLocal_priv{};
#pragma omp for nowait
        for(size_t i = 0; i < maxima.size(); i++) {
          maximaLocalToGlobal[i] = maxima[i].first;
          maximaGlobalToLocal_priv[maxima[i].first] = i;
        }
#else
      for(size_t i = 0; i < maxima.size(); i++) {
        maximaLocalToGlobal[i] = maxima[i].first;
        maximaGlobalToLocal[maxima[i].first] = i;
      }
#endif

#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
        for(size_t i = 0; i < saddles.size(); i++) {
          saddlesLocalToGlobal[i] = saddles[i].first;
          saddlesGlobalToLocal_priv[saddles[i].first] = i;
        }
#else
      for(size_t i = 0; i < saddles.size(); i++) {
        saddlesLocalToGlobal[i] = saddles[i].first;
        saddlesGlobalToLocal[saddles[i].first] = i;
      }
#endif

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
        {
          maximaGlobalToLocal.insert(
            maximaGlobalToLocal_priv.begin(), maximaGlobalToLocal_priv.end());
          saddlesGlobalToLocal.insert(
            saddlesGlobalToLocal_priv.begin(), saddlesGlobalToLocal_priv.end());
        }
      }
#endif
      //#pragma omp parallel num_threads(this->threadNumber_)
      {
        ttk::SimplexId neighborId, gId, thisOrder, maximum;
        int nNeighbors;
        //#pragma omp for nowait
        for(size_t i = 0; i < saddles.size(); i++) {
          gId = saddles[i].first;
          thisOrder = saddles[i].second;
          nNeighbors = triangulation->getVertexNeighborNumber(gId);
          for(int j = 0; j < nNeighbors; j++) {
            triangulation->getVertexNeighbor(gId, j, neighborId);
            // get the manifold result for this neighbor
            // problematic if manifold is dense and not sparse, because we need
            // the id of the point to which it is ascending, not the id of the
            // segmentation
            if(order[neighborId] > thisOrder) {
              maximum = ascendingManifold[neighborId];
              auto localMax = maximaGlobalToLocal.at(maximum);
              // this->printMsg("Saddle: " + std::to_string(gId) + " local: " +
              // std::to_string(saddlesGlobalToLocal.at(gId)) + ", max: " +
              // std::to_string(maximum) + " local: " +
              // std::to_string(localMax)); #pragma omp critical
              largestSaddlesForMax[localMax].emplace(i);
              triplets[i].emplace(localMax);
            }
          }
        }
      }
      return 1;
    }

    template <class triangulationType = ttk::AbstractTriangulation>
    int computePairs(std::vector<std::vector<ttk::SimplexId>> &persistencePairs,
                     const char *inputCriticalPoints,
                     const ttk::SimplexId *ascendingManifold,
                     const ttk::SimplexId *order,
                     const ttk::SimplexId *criticalOrder,
                     const triangulationType *triangulation,
                     const ttk::SimplexId *criticalGlobalIds,
                     ttk::SimplexId &nCriticalPoints) {

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
        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> maxima;
        ttk::SimplexId globalMin = -2;
        ttk::SimplexId globalMax = -2;
        ttk::SimplexId highestOrder = 0;
        const int dimension = triangulation->getDimensionality();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
        {
          ttk::SimplexId highestOrder_private = 0;
          ttk::SimplexId globalMax_private = -1;
          std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> saddles_priv;
          std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> maxima_priv;

#endif
          ttk::SimplexId gId = -1;
          ttk::SimplexId thisOrder = -1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp for schedule(guided) nowait
#endif
          for(ttk::SimplexId i = 0; i < nCriticalPoints; i++) {
            gId = criticalGlobalIds[i];
            // extract the global minimum id
            thisOrder = criticalOrder[i];
            if(thisOrder == 0) {
              globalMin = gId;
            }
#ifdef TTK_ENABLE_OPENMP
            else if(thisOrder > highestOrder_private) {
              highestOrder_private = thisOrder;
              globalMax_private = gId;
            }
#else
          else if(thisOrder > highestOrder) {
            highestOrder = thisOrder;
            globalMax = gId;
          }
#endif
            if(inputCriticalPoints[i] == 3) { // maxima
#ifdef TTK_ENABLE_OPENMP
              maxima_priv.push_back(std::make_pair(gId, thisOrder));
#else
            maxima.push_back(std::make_pair(gId, thisOrder));
#endif
            } else if(inputCriticalPoints[i] == (dimension - 1)
                      || inputCriticalPoints[i] == 4) { // saddles
              ttk::SimplexId neighborId;
              int nNeighbors = triangulation->getVertexNeighborNumber(gId);
              std::set<ttk::SimplexId> reachableMaxima;
              for(int j = 0; j < nNeighbors; j++) {
                triangulation->getVertexNeighbor(gId, j, neighborId);
                // only check the upper link
                if(order[neighborId] > thisOrder) {
                  reachableMaxima.insert(ascendingManifold[neighborId]);
                }
              }
              // we only care about vertices with an upperlink of at least 2
              if(reachableMaxima.size() > 1) {
#ifdef TTK_ENABLE_OPENMP
                saddles_priv.push_back(std::make_pair(gId, thisOrder));
#else
              saddles.push_back(std::make_pair(gId, thisOrder));
#endif
              }
            }
          }
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(maxima)
          {
            saddles.insert(
              saddles.end(), saddles_priv.begin(), saddles_priv.end());
            maxima.insert(maxima.end(), maxima_priv.begin(), maxima_priv.end());
            if(highestOrder_private > highestOrder) {
              highestOrder = highestOrder_private;
              globalMax = globalMax_private;
            }
          }
#endif
        }

        std::vector<std::set<ttk::SimplexId>> triplets(saddles.size());
        std::vector<std::set<ttk::SimplexId>> largestSaddlesForMax(
          maxima.size());
        std::vector<ttk::SimplexId> maximaLocalToGlobal(maxima.size());
        std::vector<ttk::SimplexId> saddlesLocalToGlobal(saddles.size());
        ttk::Timer ascTimer;
        findAscPaths<triangulationType>(
          triplets, largestSaddlesForMax, maximaLocalToGlobal,
          saddlesLocalToGlobal, maxima, saddles, order, ascendingManifold,
          triangulation);
        auto it = std::find(
          maximaLocalToGlobal.begin(), maximaLocalToGlobal.end(), globalMax);
        auto pos = it - maximaLocalToGlobal.begin();
        largestSaddlesForMax[pos] = {-1};
        this->printMsg(
          "Finished with Preprocessing, starting with PersistencePairs", 0,
          ascTimer.getElapsedTime());
        constructPersistencePairs(persistencePairs, triplets,
                                  largestSaddlesForMax, maximaLocalToGlobal,
                                  saddlesLocalToGlobal, order);
        persistencePairs.emplace_back(std::initializer_list<ttk::SimplexId>{
          globalMin, 0, globalMax, highestOrder});
        if(persistencePairs.size() != maxima.size()) {
          this->printMsg("ALARM! MAXIMA SIND NICHT GEPAIRT!!!!");
        }
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
