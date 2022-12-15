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
      const ttk::SimplexId globalMin) {
      int step = 0;
      while(triplets.size() > 1) {
        this->printMsg("============================================");
        this->printMsg("Running step " + std::to_string(step)
                       + ", #triplets: " + std::to_string(triplets.size())
                       + ", Jointree size: " + std::to_string(joinTree.size()));
        step++;
        for(auto it = triplets.begin(); it != triplets.end();) {
          if(it->second.size() <= 1) {
            this->printMsg("Removing saddle " + std::to_string(it->first));
            it = triplets.erase(it);
          } else {
            ttk::SimplexId saddle = it->first;
            ttk::SimplexId lowId = it->second.begin()->first;
            this->printMsg("Saddle-lowest Max: " + std::to_string(saddle) + "-"
                           + std::to_string(lowId));
            ++it;
            // we want to check if smallest maximum per saddle and largest
            // saddle per maximum match
            if(largestSaddleForMax[lowId] == saddle) {
              // we have a match
              joinTree.emplace_back(saddle, lowId);
              // now we need to swap the pointers
              for(auto &triplet : triplets) {
                for(auto vectit = triplet.second.begin();
                    vectit != triplet.second.end();) {
                  if(vectit->first == lowId) {
                    this->printMsg("Erasing maximum " + std::to_string(lowId)
                                   + " for saddle "
                                   + std::to_string(triplet.first));
                    vectit = triplet.second.erase(vectit);
                    // if we had to erase something in this triplet, we have to
                    // connect the saddles
                    if(saddle != triplet.first)
                      joinTree.emplace_back(saddle, triplet.first);
                  } else {
                    ++vectit;
                  }
                }
              }
            }
          }
        }
      }
      // the final saddle is connected to the global minimum
      this->printMsg("Final saddle: " + std::to_string(triplets.begin()->first)
                     + ", size: "
                     + std::to_string(triplets.begin()->second.size()));
      // joinTree.emplace_back(triplets.begin()->first, globalMin);
      return 1;
    }

    template <typename triangulationType>
    int findAscPaths(
      std::map<ttk::SimplexId,
               std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
        &triplets,
      std::map<ttk::SimplexId, ttk::SimplexId> &largestSaddleForMax,
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
            // always save the largest saddle reachable from the maximum
            // separately
            if(saddles.at(largestSaddleForMax.at(maximum)) < saddles.at(gId)) {
              largestSaddleForMax[maximum] = gId;
            }
          } else {
            this->printMsg("Created Pathlist for maximum "
                           + std::to_string(maximum));
            pathLists[maximum] = {gId};
            largestSaddleForMax[maximum] = gId;
          }
          if (maximumLists.find(gId) != maximumLists.end()){
            maximumLists.at(gId).insert(maximum);
            triplets.at(gId).emplace(
              std::make_pair(maximum, maxima.at(maximum)));
          } else {
            this->printMsg("Created Maximumlist for saddle "
                           + std::to_string(gId));
            maximumLists[gId] = {maximum};
            triplets[gId] = {std::make_pair(maximum, maxima.at(maximum))};
          }
        }
      }
      return 1;
    }

    template <class triangulationType = ttk::AbstractTriangulation>
    int computePairs(
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &joinTree,
      const char *inputCriticalPoints,
      const ttk::SimplexId *ascendingManifold,
      const ttk::SimplexId *order,
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
        const int dimension = triangulation->getDimensionality();
        // first extract the maxima and the saddles we want
        for (ttk::SimplexId i = 0; i < nCriticalPoints; i++){
          // extract the global minimum id
          if(order[i] == 0) {
            globalMin = criticalGlobalIds[i];
          }
          if (inputCriticalPoints[i] == 3){
            this->printMsg("Added " + std::to_string(criticalGlobalIds[i])
                           + " to maxima.");
            maxima.emplace(criticalGlobalIds[i], order[i]);
          } else if(inputCriticalPoints[i] == (dimension - 1)) {
            ttk::SimplexId gId = criticalGlobalIds[i];
            ttk::SimplexId neighborId;
            int nNeighbors = triangulation->getVertexNeighborNumber(gId);
            std::set<ttk::SimplexId> reachableMaxima;
            reachableMaxima.insert(ascendingManifold[gId]);
            for(int j = 0; j < nNeighbors; j++) {
              triangulation->getVertexNeighbor(gId, j, neighborId);
              reachableMaxima.insert(ascendingManifold[neighborId]);
            }
            // we only care about vertices with an upperlink of at least 2
            if(reachableMaxima.size() > 1) {
              this->printMsg("Added " + std::to_string(gId) + " to saddles.");
              this->printMsg("#reachable Maxima "
                             + std::to_string(reachableMaxima.size()));
              saddles.emplace(gId, order[i]);
            }
          }
        }

        std::map<ttk::SimplexId, std::set<ttk::SimplexId>> pathLists{};
        std::map<ttk::SimplexId, std::set<ttk::SimplexId>> maximumLists{};
        std::map<ttk::SimplexId,
                 std::set<std::pair<ttk::SimplexId, ttk::SimplexId>, PairCmp>>
          triplets{};
        std::map<ttk::SimplexId, ttk::SimplexId> largestSaddleForMax{};
        findAscPaths<triangulationType>(
          triplets, largestSaddleForMax, pathLists, maximumLists, maxima,
          saddles, ascendingManifold, triangulation);
        this->printMsg("Finished with AscPaths, starting with JoinTree");
        std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> joinTreeV2{};
        constructJoinTree(
          joinTreeV2, pathLists, maximumLists, maxima, saddles, globalMin);
        constructJoinTreeV2(joinTree, triplets, largestSaddleForMax, globalMin);

        this->printMsg("====================================");
        for(auto pair : joinTree) {
          this->printMsg(std::to_string(pair.first) + " - "
                         + std::to_string(pair.second));
        }
        this->printMsg("====================================");
        for(auto pair : joinTreeV2) {
          this->printMsg(std::to_string(pair.first) + " - "
                         + std::to_string(pair.second));
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
