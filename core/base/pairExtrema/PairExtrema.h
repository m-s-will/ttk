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

    int constructJoinTree(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &joinTree,
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &pathLists,
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &maximumLists,
      const std::map<ttk::SimplexId, ttk::SimplexId> maxima,
      const std::map<ttk::SimplexId, ttk::SimplexId> saddles,
      const ttk::SimplexId globalMin) {
      // l1-3
      for (auto const& saddle: saddles){
        pathLists[saddle.first] = std::set<ttk::SimplexId>();
      }
      // l4+5
      std::vector<ttk::SimplexId> growingNodes;
      for(auto const& maximum: maxima){
        growingNodes.push_back(maximum.first);
      }
      std::vector<ttk::SimplexId> newNodes;
      //l6
      while (!growingNodes.empty()){
        //l7-14
        for(auto const &ci : growingNodes) {
          auto Li = pathLists[ci];

          // get the critical point with maximum function value in Li
          ttk::SimplexId ck = globalMin;
          ttk::SimplexId val = std::numeric_limits<ttk::SimplexId>::lowest();
          for (auto const& potential: Li){
            const auto saddleval = saddles.at(potential);
            //this->printMsg("SaddleVal: " + std::to_string(saddleval));
            if (saddleval > val){
              ck = potential;
              val = saddleval;
            }
          }
          joinTree.emplace_back(ck, ci);
          // TODO: union of those elements, later on check UF and don't walk
          // over whole joinTree the representative element is the largest
          // maximum
          Li.erase(ck);
          pathLists[ci] = Li;
          auto mk = maximumLists[ck];
          ttk::SimplexId mkSize = mk.size();
          // get number of join tree arcs incident on ck
          ttk::SimplexId joinArcs = 0;
          for (auto const& arc: joinTree){
            if((arc.first == ck) || (arc.second == ck))
              joinArcs++;
          }
          if (joinArcs == mkSize){
            newNodes.emplace_back(ck);
          }
        }
        // l15-20
        //  apply union find
        for (auto const& ck: newNodes){
          for (auto const& arc: joinTree){
            if(arc.first == ck) {
              auto ci = arc.second;
              auto Lk = pathLists[ck];
              auto Li = pathLists[ci];
              Lk.insert(Li.begin(), Li.end());
              pathLists[ck] = Lk;
              for (auto const& keylists: maximumLists){
                auto Mj = keylists.second;
                std::set<ttk::SimplexId>::iterator search = Mj.find(ci);
                if (search != Mj.end()){
                  Mj.erase(search);
                  Mj.insert(ck);
                  maximumLists[keylists.first] = Mj;
                }
              }
            }
          }
        }
        //l21+22
        growingNodes.swap(newNodes);
        newNodes.clear();
      }
      return 1;
    }

    int constructJoinTreeV2(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &joinTree,
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
        &triplets,
      std::map<ttk::SimplexId, ttk::SimplexId> &largestSaddleForMax,
      const std::map<ttk::SimplexId, ttk::SimplexId> saddles,
      const ttk::SimplexId globalMin) {
      int step = 0;
      ttk::SimplexId smallestId = -1;
      ttk::SimplexId smallestVal = std::numeric_limits<ttk::SimplexId>::max();
      while(triplets.size() > 1) {
        this->printMsg("============================================");
        this->printMsg("Running step " + std::to_string(step)
                       + ", #triplets: " + std::to_string(triplets.size())
                       + ", Jointree size: " + std::to_string(joinTree.size()));
        step++;
        for(auto it = triplets.begin(); it != triplets.end();) {
          if(it->second.size() == 0) {
            this->printMsg("Removing saddle " + std::to_string(it->first));
            it = triplets.erase(it);
          } else {
            ttk::SimplexId saddle = it->first;
            std::pair<ttk::SimplexId, ttk::SimplexId> pair
              = *(it->second.begin());
            ttk::SimplexId lowId = pair.first;
            // we want to check if smallest maximum per saddle and largest
            // saddle per maximum match
            if(largestSaddleForMax[lowId] == saddle) {
              if(saddles.at(saddle) < smallestVal) {
                smallestVal = saddles.at(saddle);
                smallestId = saddle;
              }
              // we have a match
              this->printMsg("Saddle-lowest Max: " + std::to_string(saddle)
                             + "-" + std::to_string(lowId));
              joinTree.emplace_back(saddle, lowId);
              it->second.erase(it->second.begin());
              // now we need to swap the pointers
              for(auto &triplet : triplets) {
                if(triplet.second.erase(pair) > 0) {
                  this->printMsg("Erasing maximum " + std::to_string(lowId)
                                 + " for saddle "
                                 + std::to_string(triplet.first));
                  // triplet.second.insert(it->second.begin(),
                  // it->second.end());
                  //   if we had to erase something in this triplet, we have to
                  //   connect the saddles
                  if(saddle != triplet.first) {
                    this->printMsg("Adding edge between "
                                   + std::to_string(saddle) + " and "
                                   + std::to_string(triplet.first));
                    joinTree.emplace_back(saddle, triplet.first);
                    if(saddles.at(triplet.first) < smallestVal) {
                      smallestVal = saddles.at(triplet.first);
                      smallestId = triplet.first;
                    }
                  }
                }
              }
            }
            ++it;
          }
        }
      }
      // the final saddle is connected to the global minimum
      joinTree.emplace_back(smallestId, globalMin);
      return 1;
    }

    int constructPersistencePairs(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &pairs,
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
        &triplets,
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>,
                        InversePairCmp>> &largestSaddlesForMax) {
      int step = 0;
      bool changed = true;
      while(triplets.size() > 1 && changed) {
        changed = false;
        this->printMsg("============================================");
        this->printMsg("Running step " + std::to_string(step)
                       + ", #triplets: " + std::to_string(triplets.size())
                       + ", Pairs size: " + std::to_string(pairs.size()));
        step++;
        for(auto it = triplets.begin(); it != triplets.end();) {
          if(it->second.size() == 0) {
            changed = true;
            it = triplets.erase(it);
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
              pairs.emplace_back(saddle, lowId);
              it->second.erase(it->second.begin());
              // now we need to swap the pointers
              for(auto &triplet : triplets) {
                if(triplet.second.erase(pair)
                   > 0) { // if the maximum can be reached from this saddle, we
                          // need to replace it by the other maxima from the
                          // chosen saddle
                  triplet.second.insert(it->second.begin(), it->second.end());
                }
              }
              // we also need to update largestSaddlesForMax for all maxima
              // having had this saddle as their largest
              largestSaddlesForMax.at(lowId).erase(
                largestSaddlesForMax.at(lowId).begin());
              for(auto &saddleList : largestSaddlesForMax) {
                if(saddleList.second.erase(largestSaddlePair)
                   > 0) { // if the saddle can be reached from this maximum, we
                          // need to replace it by the other saddles from the
                          // chosen maximum
                  saddleList.second.insert(
                    largestSaddlesForMax.at(lowId).begin(),
                    largestSaddlesForMax.at(lowId).end());
                }
              }
              // for the persistence pairs, each saddle is only in one pair
              it = triplets.erase(it);
            } else {
              ++it;
            }
          }
        }
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
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &pathLists,
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &maximumLists,
      const std::map<ttk::SimplexId, ttk::SimplexId> maxima,
      const std::map<ttk::SimplexId, ttk::SimplexId> saddles,
      const ttk::SimplexId *ascendingManifold,
      const triangulationType *triangulation) {
      // construct the maximumLists for each saddle, the maxima which can be reached from this saddle
      // and the pathLists for each maximum, the saddles points which can reach this maximum
      ttk::SimplexId neighborId;
      for (auto const& saddle: saddles){
        ttk::SimplexId gId = saddle.first;
        int nNeighbors = triangulation->getVertexNeighborNumber(gId);
        for(int j = 0; j < nNeighbors; j++) {
          triangulation->getVertexNeighbor(gId, j, neighborId);
          // get the manifold result for this neighbor
          // problematic if manifold is dense and not sparse, because we need
          // the id of the point to which it is ascending, not the id of the segmentation
          ttk::SimplexId maximum = ascendingManifold[neighborId];
          if (pathLists.find(maximum) != pathLists.end()){
            pathLists.at(maximum).insert(gId);
            // always save the largest saddles reachable from the maximum
            // separately
            largestSaddlesForMax.at(maximum).emplace(
              std::make_pair(gId, saddles.at(gId)));
          } else {
            pathLists[maximum] = {gId};
            largestSaddlesForMax[maximum]
              = {std::make_pair(gId, saddles.at(gId))};
          }
          if (maximumLists.find(gId) != maximumLists.end()){
            maximumLists.at(gId).insert(maximum);
            triplets.at(gId).emplace(
              std::make_pair(maximum, maxima.at(maximum)));
          } else {
            maximumLists[gId] = {maximum};
            triplets[gId] = {std::make_pair(maximum, maxima.at(maximum))};
          }
        }
      }
      return 1;
    }

    template <class triangulationType = ttk::AbstractTriangulation>
    int computePairs(std::vector<std::vector<ttk::SimplexId>> &joinTree,
                     const char *inputCriticalPoints,
                     const ttk::SimplexId *ascendingManifold,
                     const ttk::SimplexId *order,
                     const triangulationType *triangulation,
                     const ttk::SimplexId *criticalGlobalIds,
                     const ttk::SimplexId *allGlobalIds,
                     ttk::SimplexId &nCriticalPoints,
                     ttk::SimplexId &nPoints) {

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

        // std::map<ttk::SimplexId, ttk::SimplexId> maxima;
        // std::map<ttk::SimplexId, ttk::SimplexId> saddles;
        ttk::SimplexId globalMin = -2;
        ttk::SimplexId globalMax = -2;
        ttk::SimplexId highestOrder = 0;
        ttk::SimplexId gId = -1;
        const int dimension = triangulation->getDimensionality();
        std::map<ttk::SimplexId, ttk::SimplexId> orderMap;
        for(ttk::SimplexId i = 0; i < nPoints; i++) {
          orderMap[allGlobalIds[i]] = order[i];
        }
        for (ttk::SimplexId i = 0; i < nCriticalPoints; i++){
          gId = criticalGlobalIds[i];

          // extract the global minimum id
          if(orderMap.at(gId) == 0) {
            globalMin = gId;
          } else if(orderMap.at(gId) > highestOrder) {
            highestOrder = orderMap.at(gId);
            globalMax = gId;
          }

          if(inputCriticalPoints[i] == (dimension - 1)
             || inputCriticalPoints[i] == 4) {
            ttk::SimplexId neighborId;
            int nNeighbors = triangulation->getVertexNeighborNumber(gId);
            std::set<ttk::SimplexId> reachableMaxima;
            for(int j = 0; j < nNeighbors; j++) {
              triangulation->getVertexNeighbor(gId, j, neighborId);
              // only check the upper link
              if(orderMap.at(neighborId) > orderMap.at(gId)) {
                reachableMaxima.insert(ascendingManifold[neighborId]);
              }
            }
            // we only care about vertices with an upperlink of at least 2
            if(reachableMaxima.size() > 1) {
              for(auto &maximum : reachableMaxima) {
                std::vector<ttk::SimplexId> emplacer{
                  gId, orderMap.at(gId), maximum, orderMap.at(maximum)};
                joinTree.emplace_back(emplacer);
              }
            }
          }
        }
        std::vector<ttk::SimplexId> emplacer{
          globalMin, 0, globalMax, highestOrder};
        joinTree.emplace_back(emplacer);
        /*std::map<ttk::SimplexId, std::set<ttk::SimplexId>> pathLists{};
        std::map<ttk::SimplexId, std::set<ttk::SimplexId>> maximumLists{};
        std::map<ttk::SimplexId,
                 std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
          triplets{};
        std::map<
          ttk::SimplexId,
          std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, InversePairCmp>>
          largestSaddlesForMax{};*/
        // findAscPaths<triangulationType>(
        //   triplets, largestSaddlesForMax, pathLists, maximumLists, maxima,
        //   saddles, ascendingManifold, triangulation);
        // largestSaddlesForMax.at(globalMax) = {std::make_pair(globalMin, 0)};
        // this->printMsg("Finished with AscPaths, starting with JoinTree");
        // std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> joinTreeV2{};
        //  constructJoinTree(
        //   joinTreeV2, pathLists, maximumLists, maxima, saddles, globalMin);
        //   constructJoinTreeV2(
        //     joinTree, triplets, largestSaddleForMax, saddles, globalMin);
        // constructPersistencePairs(joinTree, triplets, largestSaddlesForMax);

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
