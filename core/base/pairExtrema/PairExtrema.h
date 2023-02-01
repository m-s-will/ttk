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
      std::vector<ttk::SimplexId> maximumPointer(maximaLocalToGlobal.size());
      std::iota(std::begin(maximumPointer), std::end(maximumPointer), 0);
      std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>>
        mergeTree{}; // tuples are saddle, smaller maximum, larger maximum to
                     // which it merges
      std::set<ttk::SimplexId> saddlesToDelete{};
      std::vector<size_t> tripletSize(triplets.size());
      std::iota(std::begin(tripletSize), std::end(tripletSize), 0);
      for(size_t i = 0; i < triplets.size(); i++){
        tripletSize[i] = triplets[i].size();
      }
      saddleMaxPairs.reserve(maximaLocalToGlobal.size());
      this->printMsg("Nr of Maxima: "
                     + std::to_string(maximaLocalToGlobal.size()));

      std::vector<bool> maximaToRemove(largestSaddlesForMax.size(), false);

      while(changed && (pairs.size() < (maximaLocalToGlobal.size() - 1))) {
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
            ttk::SimplexId largestSaddle = -1;
            if (!maximaToRemove.at(maximum)){
              largestSaddle = *(largestSaddlesForMax.at(maximum).begin());
            }
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
                  maximaToRemove[maximum] = true;
#ifdef TTK_ENABLE_OPENMP
                  {
                    pairs_priv.emplace_back(
                      std::initializer_list<ttk::SimplexId>{
                        saddlesLocalToGlobal[largestSaddle],
                        order[saddlesLocalToGlobal[largestSaddle]],
                        maximaLocalToGlobal[maximum],
                        order[maximaLocalToGlobal[maximum]]});
                    saddleMaxPairs_priv.push_back(
                      std::make_pair(largestSaddle, maximum));
                  }
#else
                  {
                    pairs.emplace_back(std::initializer_list<ttk::SimplexId>{
                      saddlesLocalToGlobal[largestSaddle],
                      order[saddlesLocalToGlobal[largestSaddle]],
                      maximaLocalToGlobal[maximum],
                      order[maximaLocalToGlobal[maximum]]});
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
        std::pair<ttk::SimplexId, ttk::SimplexId> pair;
        ttk::SimplexId saddle, max, largestMax;
        for(size_t i = 0; i < saddleMaxPairs.size(); i++) {
          pair = saddleMaxPairs[i];
          saddle = pair.first;
          max = pair.second;
          auto &maxSet
            = triplets.at(saddle); // set containing all the maxima reachable
                                  // from the saddle to be removed
          largestMax = *(maxSet.rbegin());
          maximumPointer[max] = largestMax;
          mergeTree.push_back(std::make_tuple(saddle, max, largestMax));
          auto &saddleSet = largestSaddlesForMax.at(max);
          maxSet.erase(maxSet.begin());
          saddleSet.erase(saddleSet.begin());
          //auto &thisSize = tripletSize.at(saddle);
          //thisSize--;
          //if (thisSize < 2)
          //  saddlesToDelete.insert(saddle);
          if(maxSet.size() < 2)
            saddlesToDelete.insert(saddle);
          // now we need to swap the pointers, all saddles which could reach
          // maximum pair.first, now can reach all the other maxima of saddle
          // we can find this by just checking largestSaddlesForMax for the
          // maximum
          for(auto &saddleId : saddleSet) {
            auto &triplet = triplets.at(saddleId);
            //auto &size = tripletSize.at(saddleId);
            triplet.erase(max);
            // if the maximum can be reached from other saddles, we
            // need to replace it by the other maxima from the
            // chosen saddle
            triplet.insert(maxSet.begin(), maxSet.end());
            // triplet.insert(largestMax);
            //size--;
            //if (size < 2)
            //  saddlesToDelete.insert(saddleId);
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


        /*
        this->printMsg("Triplets beforehand:");
        for (size_t i = 0; i < triplets.size(); i++){
          this->printMsg("Saddle " + std::to_string(i));
          this->printMsg("Maxima:");
          for (auto &t: triplets[i]){
            this->printMsg(std::to_string(t));
          }
        }
        for (size_t i = 0; i < largestSaddlesForMax.size(); i++){
          this->printMsg("Maximum " + std::to_string(i));
          this->printMsg("Saddles:");
          for (auto &l : largestSaddlesForMax[i]){
            this->printMsg(std::to_string(l));
          }
        }


        //collecting everything to do in parallel
        std::vector<std::set<ttk::SimplexId>> toInsertForSaddles(triplets.size());
        std::set<ttk::SimplexId> toRemoveForSaddles{};
        std::vector<std::set<ttk::SimplexId>> toInsertForMaxima(largestSaddlesForMax.size());
        #pragma omp parallel
        {
          std::set<ttk::SimplexId> saddlesToDelete_priv{};
          std::vector<std::set<ttk::SimplexId>> toInsertForSaddles_priv(triplets.size());
          std::set<ttk::SimplexId> toRemoveForSaddles_priv{};
          std::vector<std::set<ttk::SimplexId>> toInsertForMaxima_priv(largestSaddlesForMax.size());

          #pragma omp for schedule(guided)
          for(size_t i = 0; i < saddleMaxPairs.size(); i++) {
            std::pair<ttk::SimplexId, ttk::SimplexId> pair = saddleMaxPairs[i];
            ttk::SimplexId saddle = pair.first;
            ttk::SimplexId max = pair.second;
            //this->printMsg("Working on saddle-max pair " + std::to_string(saddle) + " " + std::to_string(max));
            auto &maxSet = triplets.at(saddle);
            //ttk::SimplexId largestMax = *(maxSet.rbegin());
            //mergeTree.push_back(std::make_tuple(saddle, max, largestMax));
            auto &saddleSet = largestSaddlesForMax.at(max);
            maxSet.erase(maxSet.begin());
            saddleSet.erase(saddleSet.begin());
            if(maxSet.size() < 2)
              saddlesToDelete_priv.insert(saddle);
            for(auto &saddleId : saddleSet) {
              toRemoveForSaddles_priv.insert(max);
              toInsertForSaddles_priv[saddleId].insert(maxSet.begin(), maxSet.end());
              if (triplets.at(saddleId).size() < 2){
                //this->printMsg("Removing saddle " + std::to_string(saddlesLocalToGlobal[i]) + " after removing paired maxima " + std::to_string(i));
                saddlesToDelete_priv.insert(saddleId);
              }
              //this->printMsg("saddleid " + std::to_string(saddleId) + " has to remove max " + std::to_string(max));
            }
            for(auto &maxId : maxSet) {
              toInsertForMaxima_priv[maxId].insert(saddleSet.begin(), saddleSet.end());
            }
          }

          #pragma omp critical
          {
            saddlesToDelete.insert(saddlesToDelete_priv.begin(), saddlesToDelete_priv.end());
            for (size_t i = 0; i < toInsertForSaddles.size(); i++)
              toInsertForSaddles[i].insert(toInsertForSaddles_priv[i].begin(), toInsertForSaddles_priv[i].end());
            for (size_t i = 0; i < toInsertForMaxima.size(); i++)
              toInsertForMaxima[i].insert(toInsertForMaxima_priv[i].begin(), toInsertForMaxima_priv[i].end());
            //for (size_t i = 0; i < toRemoveForSaddles.size(); i++)
            //  toRemoveForSaddles[i].insert(toRemoveForSaddles_priv[i].begin(), toRemoveForSaddles_priv[i].end());
              toRemoveForSaddles.insert(toRemoveForSaddles_priv.begin(), toRemoveForSaddles_priv.end());
          }
        }

        // do everything in parallel
        #pragma omp parallel
        {
          std::set<ttk::SimplexId> saddlesToDelete_priv{};

          #pragma omp for schedule(guided) nowait
          for(size_t i = 0; i < toInsertForMaxima.size(); i++){
            if (!maximaToRemove[i])
            {
              //this->printMsg("Inserting " + std::to_string(toInsertForMaxima[i].size()) + " values for max " + std::to_string(i));
              largestSaddlesForMax.at(i).insert(toInsertForMaxima[i].begin(), toInsertForMaxima[i].end());
            } else {
              //this->printMsg("Clearing values for max " + std::to_string(i));
              largestSaddlesForMax.at(i).clear();
            }
          }

          #pragma omp for schedule(guided)
          for(size_t i = 0; i < toInsertForSaddles.size(); i++){
            //this->printMsg("Inserting " + std::to_string(toInsertForSaddles[i].size()) + " values for saddle " + std::to_string(i));
            triplets.at(i).insert(toInsertForSaddles[i].begin(), toInsertForSaddles[i].end());
          }

          #pragma omp for schedule(guided) nowait
          for(size_t i = 0; i < triplets.size(); i++){
            //this->printMsg("Erasing " + std::to_string(toRemoveForSaddles.size()) + " values for saddle " + std::to_string(i));

            auto v = triplets.at(i);
            if (v.size() > 0){
              std::set<ttk::SimplexId> out;
              std::set_difference(std::begin(v), std::end(v),
                          std::begin(toRemoveForSaddles), std::end(toRemoveForSaddles),
                          std::inserter(out, out.end()));
              triplets.at(i) = out;
            }
          }

          #pragma omp critical
          {
            saddlesToDelete.insert(saddlesToDelete_priv.begin(), saddlesToDelete_priv.end());
          }
        }
        */
        this->printMsg("Finished swapping pointers for "
                         + std::to_string(saddleMaxPairs.size()) + " pairs.",
                       0.66, swappingTimer.getElapsedTime());
        ttk::Timer delTimer;

        for(auto &saddleToDel : saddlesToDelete) {
          //this->printMsg("Removing saddle " + std::to_string(saddleToDel));
          auto &maxList = triplets.at(saddleToDel);
          for(auto &m : maxList) {
            auto &saddleSet = largestSaddlesForMax.at(m);
            if (saddleSet.size() > 0)
              saddleSet.erase(saddleToDel);
          }
          triplets.at(saddleToDel).clear();
        }

        this->printMsg("Finished deletion of "
                         + std::to_string(saddlesToDelete.size())
                         + " irrelevant saddles.",
                       1, delTimer.getElapsedTime());
        this->printMsg("Finished step " + std::to_string(step), 1,
                       stepTimer.getElapsedTime());
        step++;
        /*
        for (size_t i = 0; i < triplets.size(); i++){
          this->printMsg("Saddle " + std::to_string(i));
          this->printMsg("Maxima:");
          for (auto &t: triplets[i]){
            this->printMsg(std::to_string(t));
          }
        }
        for (size_t i = 0; i < largestSaddlesForMax.size(); i++){
          this->printMsg("Maximum " + std::to_string(i) + " " + std::to_string(maximaLocalToGlobal[i]));
          this->printMsg("Saddles:");
          for (auto &l : largestSaddlesForMax[i]){
            this->printMsg(std::to_string(l) + " " + std::to_string(saddlesLocalToGlobal[l]));
          }
        } */

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

    ttk::SimplexId neighborId, gId, thisOrder, maximum;
    int nNeighbors;
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
          // std::to_string(localMax));
          triplets[i].emplace(localMax);
          largestSaddlesForMax[localMax].emplace(i);
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
