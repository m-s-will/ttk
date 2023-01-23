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
#include <limits.h>
#include <Triangulation.h>

namespace ttk {
  struct PairCmp {
    bool
      operator()(const std::pair<ttk::SimplexId, ttk::SimplexId> &lhs,
                 const std::pair<ttk::SimplexId, ttk::SimplexId> &rhs) const {
      return lhs.second < rhs.second;
    }
  };

  struct InversePairCmp {
    bool
      operator()(const std::pair<ttk::SimplexId, ttk::SimplexId> &lhs,
                 const std::pair<ttk::SimplexId, ttk::SimplexId> &rhs) const {
      return lhs.second > rhs.second;
    }
  };

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
      std::unordered_map<
        ttk::SimplexId,
        std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>> &triplets,
      std::unordered_map<ttk::SimplexId,
                         std::set<std::pair<ttk::SimplexId, ttk::SimplexId>,
                                  InversePairCmp>> &largestSaddlesForMax,
      const ttk::SimplexId *order) {
      int step = 0;
      bool changed = true;
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> saddleMaxPairs{};
      std::vector<ttk::SimplexId> saddlesToDelete{};
      saddleMaxPairs.reserve(triplets.size());
      saddlesToDelete.reserve(triplets.size());
      while(triplets.size() > 0 && changed) {
        changed = false;
        saddleMaxPairs.clear();
        saddlesToDelete.clear();
        this->printMsg("============================================");
        this->printMsg("Running step " + std::to_string(step)
                       + ", #triplets: " + std::to_string(triplets.size())
                       + ", Pairs size: " + std::to_string(pairs.size()));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
#endif
        {
          std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
            saddleMaxPairs_priv{};
          std::vector<ttk::SimplexId> saddlesToDelete_priv{};
          std::vector<std::vector<ttk::SimplexId>> pairs_priv{};
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif
          for(size_t i = 0; i < triplets.size(); i++) {
            auto it = triplets.begin();
            advance(it, i);
            if(it->second.size() < 2) {
              changed = true;
#ifdef TTK_ENABLE_OPENMP
              saddlesToDelete_priv.push_back(it->first);
#else
              saddlesToDelete.push_back(it->first);
#endif
            } else {
              ttk::SimplexId saddle = it->first;
              std::pair<ttk::SimplexId, ttk::SimplexId> pair
                = *(it->second.begin()); // the lowest maximum id and val for
                                         // the current saddle
              ttk::SimplexId lowId = pair.first; // only the maximum id
              std::pair<ttk::SimplexId, ttk::SimplexId> largestSaddlePair
                = *(largestSaddlesForMax.at(lowId).begin());

              // we want to check if smallest maximum per saddle and largest
              // saddle per maximum match
              if(largestSaddlePair.first == saddle) {
                changed = true;
#ifdef TTK_ENABLE_OPENMP
                {
                  pairs_priv.emplace_back(std::initializer_list<ttk::SimplexId>{
                    saddle, largestSaddlePair.second, lowId, pair.second});
                  saddleMaxPairs_priv.push_back(std::make_pair(saddle, lowId));
                }
#else
                {
                  pairs.emplace_back(std::initializer_list<ttk::SimplexId>{
                    saddle, largestSaddlePair.second, lowId, pair.second});
                  saddleMaxPairs.push_back(std::make_pair(saddle, lowId));
                }
#endif
              }
            }
          }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(deletions)
          {
            saddleMaxPairs.insert(saddleMaxPairs.end(),
                                  saddleMaxPairs_priv.begin(),
                                  saddleMaxPairs_priv.end());
            saddlesToDelete.insert(saddlesToDelete.end(),
                                   saddlesToDelete_priv.begin(),
                                   saddlesToDelete_priv.end());
            pairs.insert(pairs.end(), pairs_priv.begin(), pairs_priv.end());
          }
#endif
        }
        this->printMsg("Finished finding pairs for step " + std::to_string(step)
                       + ", starting swapping of pointers for "
                       + std::to_string(saddleMaxPairs.size()) + " pairs.");

        // now we need to swap pointers
        //#ifdef TTK_ENABLE_OPENMP
        //#pragma omp parallel for num_threads(this->threadNumber_)
        //#endif
        for(size_t i = 0; i < saddleMaxPairs.size(); i++) {
          auto pair = saddleMaxPairs[i];
          auto saddle = pair.first;
          auto max = pair.second;
          auto &maxSet
            = triplets.at(saddle); // set containing all the maxima reachable
                                   // from the saddle to be removed
          auto maxPair = *(
            maxSet.begin()); // pair contains the maximum with which this
                             // saddle is paired and the order of the maximum
          auto &saddleSet = largestSaddlesForMax.at(max);
          // auto largestSaddlePair = *(saddleSet.begin());
          maxSet.erase(maxSet.begin());

          // now we need to swap the pointers, all saddles which could reach
          // maximum pair.first, now can reach all the other maxima of saddle
          // we can find this by just checking largestSaddlesForMax for the
          // maximum
          for(auto &val : saddleSet) {
            auto saddleId = val.first;
            auto &triplet = triplets.at(saddleId);
            //#ifdef TTK_ENABLE_OPENMP
            //#pragma omp critical
            //#endif
            {
              triplet.erase(maxPair);
              // if the maximum can be reached from other saddles, we
              // need to replace it by the other maxima from the
              // chosen saddle
              triplet.insert(maxSet.begin(), maxSet.end());
            }
          }

          for(auto &val : maxSet) {
            auto maxId = val.first;
            auto &saddleList = largestSaddlesForMax.at(maxId);
            // if the maximum can be reached from other saddles, we
            // need to replace it by the other maxima from the
            // chosen saddle
            //#ifdef TTK_ENABLE_OPENMP
            //#pragma omp critical
            //#endif
            saddleList.insert(saddleSet.begin(), saddleSet.end());
          }
          //#ifdef TTK_ENABLE_OPENMP
          //#pragma omp critical
          //#endif
          largestSaddlesForMax.erase(max);
        }

        this->printMsg("Finished removing pairs for step "
                       + std::to_string(step) + ", starting deletion for "
                       + std::to_string(saddlesToDelete.size())
                       + " irrelevant saddles.");
        //#ifdef TTK_ENABLE_OPENMP
        //#pragma omp parallel for num_threads(this->threadNumber_)
        //#endif
        for(size_t i = 0; i < saddlesToDelete.size(); i++) {
          auto saddle = saddlesToDelete[i];
          std::pair<ttk::SimplexId, ttk::SimplexId> saddlePair
            = std::make_pair(saddle, order[saddle]);
          auto &maxList = triplets.at(saddle);
          for(auto &m : maxList) {
            auto &saddleList = largestSaddlesForMax.at(m.first);
            //#ifdef TTK_ENABLE_OPENMP
            //#pragma omp critical
            //#endif
            saddleList.erase(saddlePair);
          }
          //#ifdef TTK_ENABLE_OPENMP
          //#pragma omp critical
          //#endif
          triplets.erase(saddle);
        }
        step++;
      }
      return 1;
    }

    template <typename triangulationType>
    int findAscPaths(
      std::unordered_map<
        ttk::SimplexId,
        std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>> &triplets,
      std::unordered_map<ttk::SimplexId,
                         std::set<std::pair<ttk::SimplexId, ttk::SimplexId>,
                                  InversePairCmp>> &largestSaddlesForMax,
      const std::vector<ttk::SimplexId> saddles,
      const ttk::SimplexId *order,
      const ttk::SimplexId *ascendingManifold,
      const triangulationType *triangulation) {
      // construct the maximumLists for each saddle, the maxima which can be reached from this saddle
      // and the pathLists for each maximum, the saddles points which can reach this maximum
      ttk::SimplexId neighborId;
      for (auto const& saddle: saddles){
        ttk::SimplexId gId = saddle;
        ttk::SimplexId thisOrder = order[gId];
        int nNeighbors = triangulation->getVertexNeighborNumber(gId);
        for(int j = 0; j < nNeighbors; j++) {
          triangulation->getVertexNeighbor(gId, j, neighborId);
          // get the manifold result for this neighbor
          // problematic if manifold is dense and not sparse, because we need
          // the id of the point to which it is ascending, not the id of the segmentation
          if(order[neighborId] > thisOrder) {
            ttk::SimplexId maximum = ascendingManifold[neighborId];
            if(largestSaddlesForMax.find(maximum)
               != largestSaddlesForMax.end()) {
              // always save the largest saddles reachable from the maximum
              // separately
              largestSaddlesForMax.at(maximum).emplace(gId, thisOrder);
            } else {
              largestSaddlesForMax[maximum] = {std::make_pair(gId, thisOrder)};
            }
            if(triplets.find(gId) != triplets.end()) {
              triplets.at(gId).emplace(maximum, order[maximum]);
            } else {
              triplets[gId] = {std::make_pair(maximum, order[maximum])};
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

        std::vector<ttk::SimplexId> saddles;
        ttk::SimplexId globalMin = -2;
        ttk::SimplexId globalMax = -2;
        ttk::SimplexId highestOrder = 0;
        const int dimension = triangulation->getDimensionality();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
        {
          ttk::SimplexId highestOrder_private = 0;
          ttk::SimplexId globalMax_private = -1;
          std::vector<ttk::SimplexId> saddles_priv;

#endif
          ttk::SimplexId gId = -1;
          ttk::SimplexId thisOrder = -1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp for nowait
#endif
          for(ttk::SimplexId i = 0; i < nCriticalPoints; i++) {
            gId = criticalGlobalIds[i];
            // this->printMsg("i: " + std::to_string(i) + ", gid: " +
            // std::to_string(gId) + ", thread: " +
            // std::to_string(omp_get_thread_num()));

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
            if(inputCriticalPoints[i] == (dimension - 1)
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
                saddles_priv.push_back(gId);
#else
              saddles.push_back(gId);
#endif
              }
            }
          }
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical(maxima)
          {
            saddles.insert(
              saddles.end(), saddles_priv.begin(), saddles_priv.end());
            if(highestOrder_private > highestOrder) {
              highestOrder = highestOrder_private;
              globalMax = globalMax_private;
            }
          }
#endif
        }

        std::unordered_map<
          ttk::SimplexId,
          std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
          triplets{};
        std::unordered_map<
          ttk::SimplexId,
          std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, InversePairCmp>>
          largestSaddlesForMax{};
        this->printMsg("Starting with AscPaths, saddles length: "
                       + std::to_string(saddles.size()));
        findAscPaths<triangulationType>(triplets, largestSaddlesForMax, saddles,
                                        order, ascendingManifold,
                                        triangulation);

        largestSaddlesForMax[globalMax] = {std::make_pair(globalMin, 0)};
        this->printMsg(
          "Finished with Preprocessing, starting with PersistencePairs");
        constructPersistencePairs(
          persistencePairs, triplets, largestSaddlesForMax, order);
        persistencePairs.emplace_back(std::initializer_list<ttk::SimplexId>{
          globalMin, 0, globalMax, highestOrder});

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
