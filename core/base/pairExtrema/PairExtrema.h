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
      std::vector<std::set<ttk::SimplexId>> &triplets,
      const std::vector<ttk::SimplexId> &maximaLocalToGlobal,
      const std::vector<ttk::SimplexId> &saddlesLocalToGlobal) {
      int step = 0;
      bool changed = true;
      // std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
      // saddleMaxPairs{};
      std::vector<ttk::SimplexId> maximumPointer(maximaLocalToGlobal.size());
      std::iota(std::begin(maximumPointer), std::end(maximumPointer), 0);
      std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>>
        mergeTree(maximaLocalToGlobal.size());
      ttk::SimplexId nrOfPairs = 0;
      while(changed) {
        ttk::Timer stepTimer;
        changed = false;
        this->printMsg(ttk::debug::Separator::L2);
        this->printMsg("Running step " + std::to_string(step)
                         + ", Pairs: " + std::to_string(nrOfPairs),
                       0, 0);

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

        // locks werden teuer wenn es nur noch wenige Maxima gibt
        // eher alle eigenes largestSaddlesForMax und dann reduction?

        this->printMsg("Finished building largestSaddlesForMax", 0.2,
                       buildTimer.getElapsedTime());
        ttk::Timer pairTimer;

#pragma omp parallel for schedule(guided) num_threads(this->threadNumber_)
        for(size_t i = 0; i < largestSaddlesForMax.size() - 1; i++) {
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
              mergeTree[maximum]
                = (std::make_tuple(largestSaddle, maximum, largestMax));
            }
          }
        }

        this->printMsg("Finished finding pairs and swapping pointers", 0.5,
                       pairTimer.getElapsedTime());

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
                       0.8, compressTimer.getElapsedTime());
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
                       replaceTimer.getElapsedTime());

        this->printMsg("Finished step " + std::to_string(step), 1,
                       stepTimer.getElapsedTime());
        step++;
      }
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

      // std::unordered_map<ttk::SimplexId, ttk::SimplexId>
      // maximaGlobalToLocal{}; maximaGlobalToLocal.reserve(maxima.size());
#pragma omp parallel num_threads(this->threadNumber_)
      {
        // std::unordered_map<ttk::SimplexId, ttk::SimplexId>
        //   maximaGlobalToLocal_priv{};
#pragma omp for nowait
        for(size_t i = 0; i < maxima.size(); i++) {
          maximaLocalToGlobal[i] = maxima[i].first;
          // maximaGlobalToLocal_priv[maxima[i].first] = i;
          tempArray[maxima[i].first] = i;
        }

#pragma omp for nowait
        for(size_t i = 0; i < saddles.size(); i++) {
          saddlesLocalToGlobal[i] = saddles[i].first;
        }
      }
      /*
      #pragma omp critical
              {
                maximaGlobalToLocal.insert(
                  maximaGlobalToLocal_priv.begin(),
      maximaGlobalToLocal_priv.end());
              }
            }*/
      this->printMsg("Finished sorting and building the normalization arrays",
                     0, buildTimer.getElapsedTime());

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
            // triplets[i].emplace(maximaGlobalToLocal[ascendingManifold[neighborId]]);
            triplets[i].emplace(tempArray[ascendingManifold[neighborId]]);
          }
        }
      }
      this->printMsg(
        "Finished building the triplets", 0, tripletTimer.getElapsedTime());
      return 1;
    }

    template <class triangulationType = ttk::AbstractTriangulation>
    int computePairs(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
      const ttk::SimplexId *minimaIds,
      const ttk::SimplexId *saddle1Ids,
      const ttk::SimplexId *saddle2Ids,
      const ttk::SimplexId *maximaIds,
      const ttk::SimplexId *ascendingManifold,
      ttk::SimplexId *tempArray,
      const ttk::SimplexId *order,
      const triangulationType *triangulation,
      ttk::SimplexId &nMinima,
      ttk::SimplexId &nSaddle1,
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
        ttk::SimplexId globalMin = -2;
        const int dimension = triangulation->getDimensionality();
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

        // this->printMsg(std::to_string(maxTime / 1000));
        persistencePairs.resize(maxima.size());
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
        constructPersistencePairs(persistencePairs, triplets,
                                  maximaLocalToGlobal, saddlesLocalToGlobal);
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
