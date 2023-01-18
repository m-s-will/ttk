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
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
        &triplets,
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>,
                        InversePairCmp>> &largestSaddlesForMax,
      const ttk::SimplexId *order) {
      int step = 0;
      bool changed = true;
      std::vector<ttk::SimplexId> saddlesToRemove{};
      std::vector<ttk::SimplexId> maxToRemove{};
      saddlesToRemove.reserve(triplets.size());
      maxToRemove.reserve(triplets.size());
      while(triplets.size() > 0 && changed) {
        changed = false;
        saddlesToRemove.clear();
        maxToRemove.clear();
        this->printMsg("============================================");
        this->printMsg("Running step " + std::to_string(step)
                       + ", #triplets: " + std::to_string(triplets.size())
                       + ", Pairs size: " + std::to_string(pairs.size()));
        for(auto it = triplets.begin(); it != triplets.end(); ++it) {
          if(it->second.size() < 2) {
            changed = true;
            saddlesToRemove.push_back(it->first);
          } else {
            ttk::SimplexId saddle = it->first;
            std::pair<ttk::SimplexId, ttk::SimplexId> pair
              = *(it->second.begin()); // the lowest maximum id and val for the
                                       // current saddle
            ttk::SimplexId lowId = pair.first; // only the maximum id
            std::pair<ttk::SimplexId, ttk::SimplexId> largestSaddlePair
              = *(largestSaddlesForMax.at(lowId).begin());

            std::string text = "Triplet with Saddle " + std::to_string(saddle)
                               + ", maxima ids and values: ";
            this->printMsg("Connected Maxima:");
            std::string maximText;
            for(auto &curpair : it->second) {
              text += " " + std::to_string(curpair.first) + " "
                      + std::to_string(curpair.second) + ",";
              auto saddlesForMax = largestSaddlesForMax.at(curpair.first);
              maximText = "Maximum " + std::to_string(curpair.first)
                          + ", saddle ids and values: ";
              for(auto &cursad : saddlesForMax) {
                maximText += " " + std::to_string(cursad.first) + " "
                             + std::to_string(cursad.second) + ",";
              }
              this->printMsg(maximText);
            }
            this->printMsg(text);
            this->printMsg("================================");

            // we want to check if smallest maximum per saddle and largest
            // saddle per maximum match
            if(largestSaddlePair.first == saddle) {
              changed = true;
              pairs.emplace_back(std::initializer_list<ttk::SimplexId>{
                saddle, largestSaddlePair.second, lowId, pair.second});
              saddlesToRemove.push_back(saddle);
              maxToRemove.push_back(lowId);
            }
          }
        }
        this->printMsg("Finished pairing for step " + std::to_string(step)
                       + ", starting removal for "
                       + std::to_string(saddlesToRemove.size()) + " saddles.");
        int saddleStep = 0;
        // now we need to swap pointers

        for(auto saddle : saddlesToRemove) {
          if(saddleStep % 1000 == 0)
            this->printMsg(
              "Removing values for saddle " + std::to_string(saddle),
              (float)saddleStep / saddlesToRemove.size());
          saddleStep++;
          auto saddleSet
            = triplets.at(saddle); // set containing all the maxima reachable
                                   // from the saddle to be removed
          auto pair = *(
            saddleSet.begin()); // pair contains the maximum with which this
                                // saddle is paired and the order of the maximum
          saddleSet.erase(saddleSet.begin());
          if(saddleSet.size() != 0) {
            // now we need to swap the pointers, all saddles which could reach
            // maximum pair.first, now can reach all the other maxima of saddle
            // we can find this by just checking largestSaddleForMax for the
            // maximum
            int didWeDo = 0;
            auto saddlesForThisMax = largestSaddlesForMax.at(pair.first);
            std::vector<ttk::SimplexId> saddlesChanged{};
            for(auto &triplet : triplets) {
              if(triplet.second.erase(pair)
                 > 0) { // if the maximum can be reached from other saddles, we
                        // need to replace it by the other maxima from the
                        // chosen saddle
                triplet.second.insert(saddleSet.begin(), saddleSet.end());
                didWeDo++;
                saddlesChanged.push_back(triplet.first);
              }
            }
            if(true) {
              this->printMsg("Saddle " + std::to_string(saddle)
                             + ", max das erased wird: "
                             + std::to_string(pair.first));
              this->printMsg("#Triplets we had to change: "
                             + std::to_string(didWeDo));
              for(auto &val : saddlesChanged) {
                this->printMsg(std::to_string(val));
              }
              this->printMsg("#Saddles for this max: "
                             + std::to_string(saddlesForThisMax.size()));
              for(auto &val : saddlesForThisMax) {
                this->printMsg(std::to_string(val.first));
              }
            }
          } else {
            triplets.erase(saddle);
            std::pair<ttk::SimplexId, ttk::SimplexId> saddlePair
              = std::make_pair(saddle, order[saddle]);
            for(auto &saddleList : largestSaddlesForMax) {
              saddleList.second.erase(saddlePair);
            }
          }
        }

        this->printMsg("Finished removing saddles for step "
                       + std::to_string(step) + ", starting removal for "
                       + std::to_string(maxToRemove.size()) + " maxima.");
        int maxStep = 0;
        for(auto max : maxToRemove) {
          if(maxStep % 1000 == 0)
            this->printMsg("Removing max  " + std::to_string(max),
                           (float)maxStep / maxToRemove.size());
          maxStep++;
          auto maxSet = largestSaddlesForMax.at(max);
          auto largestSaddlePair = *(maxSet.begin());
          // we also need to update largestSaddlesForMax for all maxima
          // having had this saddle as their largest
          // maxSet.erase(maxSet.begin());
          for(auto &saddleList : largestSaddlesForMax) {
            if(saddleList.second.find(largestSaddlePair)
               != saddleList.second.end()) { // if the saddle can be reached
                                             // from this maximum, we
              // need to replace it by the other saddles from the
              // chosen maximum
              saddleList.second.insert(maxSet.begin(), maxSet.end());
            }
          }
          largestSaddlesForMax.erase(max);
        }
        step++;
      }
      return 1;
    }

    template <typename triangulationType>
    int findAscPaths(
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
        &triplets,
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>,
                        InversePairCmp>> &largestSaddlesForMax,
      const std::map<ttk::SimplexId, ttk::SimplexId> maxima,
      const std::map<ttk::SimplexId, ttk::SimplexId> saddles,
      const ttk::SimplexId *order,
      const ttk::SimplexId *ascendingManifold,
      const triangulationType *triangulation) {
      // construct the maximumLists for each saddle, the maxima which can be reached from this saddle
      // and the pathLists for each maximum, the saddles points which can reach this maximum
      ttk::SimplexId neighborId;
      for (auto const& saddle: saddles){
        ttk::SimplexId gId = saddle.first;
        ttk::SimplexId thisOrder = saddle.second;
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
              triplets.at(gId).emplace(maximum, maxima.at(maximum));
            } else {
              triplets[gId] = {std::make_pair(maximum, maxima.at(maximum))};
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

        std::map<ttk::SimplexId, ttk::SimplexId> maxima;
        std::map<ttk::SimplexId, ttk::SimplexId> saddles;
        ttk::SimplexId globalMin = -2;
        ttk::SimplexId globalMax = -2;
        ttk::SimplexId highestOrder = 0;
        ttk::SimplexId gId = -1;
        ttk::SimplexId thisOrder = -1;
        const int dimension = triangulation->getDimensionality();

        for (ttk::SimplexId i = 0; i < nCriticalPoints; i++){
          gId = criticalGlobalIds[i];

          // extract the global minimum id
          thisOrder = criticalOrder[i];
          if(thisOrder == 0) {
            globalMin = gId;
          } else if(thisOrder > highestOrder) {
            highestOrder = thisOrder;
            globalMax = gId;
          }
          if(inputCriticalPoints[i] == 3) { // maxima
            maxima.emplace(gId, thisOrder);
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
              saddles.emplace(gId, thisOrder);
            }
          }
        }
        std::map<ttk::SimplexId,
                 std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
          triplets{};
        std::map<
          ttk::SimplexId,
          std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, InversePairCmp>>
          largestSaddlesForMax{};
        findAscPaths<triangulationType>(triplets, largestSaddlesForMax, maxima,
                                        saddles, order, ascendingManifold,
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
