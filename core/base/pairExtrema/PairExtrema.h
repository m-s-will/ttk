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

    template <typename dataType>
    int constructJoinTree(
      std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId>> &joinTree,
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &pathLists,
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &maximumLists,
      const std::map<ttk::SimplexId, dataType> maxima,
      const std::map<ttk::SimplexId, dataType> saddles,
      const dataType globalMin) {
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
        for (auto const& ci: growingNodes){
          //this->printMsg("ci: " + std::to_string(ci));
          auto Li = pathLists[ci];
          //this->printMsg("Li size: " + std::to_string(Li.size()));

          // get the critical point with maximum function value in Li
          ttk::SimplexId ck = globalMin;
          dataType val = std::numeric_limits<dataType>::lowest();
          for (auto const& potential: Li){
            const auto saddleval = saddles.at(potential);
            //this->printMsg("SaddleVal: " + std::to_string(saddleval));
            if (saddleval > val){
              ck = potential;
              val = saddleval;
            }
          }
          joinTree.emplace_back(ck, ci);
          Li.erase(ck);
          pathLists[ci] = Li;
          auto mk = maximumLists[ck];
          ttk::SimplexId mkSize = mk.size();
          // get number of join tree arcs incident on ck
          ttk::SimplexId joinArcs = 0;
          for (auto const& arc: joinTree){
            if ((std::get<0>(arc) == ck) || (std::get<1>(arc) == ck))
              joinArcs++;
          }
          if (joinArcs == mkSize){
            newNodes.emplace_back(ck);
          }
        }
        //l15-20
        for (auto const& ck: newNodes){
          for (auto const& arc: joinTree){
            if (std::get<0>(arc) == ck){
              auto ci = std::get<1>(arc);
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

    template <typename dataType, typename triangulationType>
    int findAscPaths(
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &pathLists,
      std::map<ttk::SimplexId, std::set<ttk::SimplexId>> &maximumLists,
      const std::map<ttk::SimplexId, dataType> saddles,
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
          this->printMsg("Neighbor " + std::to_string(neighborId) + " of point "
                         + std::to_string(gId) + " reaching maximum "
                         + std::to_string(maximum));
          if (pathLists.find(maximum) != pathLists.end()){
            pathLists.at(maximum).insert(gId);
          } else {
            this->printMsg("Created Pathlist for maximum "
                           + std::to_string(maximum));
            pathLists[maximum] = {gId};
          }
          if (maximumLists.find(gId) != maximumLists.end()){
            maximumLists.at(gId).insert(maximum);
          } else {
            this->printMsg("Created Maximumlist for saddle "
                           + std::to_string(gId));
            maximumLists[gId] = {maximum};
          }
        }

      }
      return 1;
    }
    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computePairs(
      std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId>> &joinTree,
      const char *inputCriticalPoints,
      const ttk::SimplexId *ascendingManifold,
      const dataType *order,
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

        std::map<ttk::SimplexId, dataType> maxima;
        std::map<ttk::SimplexId, dataType> saddles;
        dataType globalMin = -2;
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
            //} else if (inputCriticalPoints[i] == 2){
            // we need to throw out genus changing saddles
            ttk::SimplexId gId = criticalGlobalIds[i];
            ttk::SimplexId neighborId;
            int nNeighbors = triangulation->getVertexNeighborNumber(gId);
            std::vector<ttk::SimplexId> upperLink;
            std::vector<ttk::SimplexId> lowerLink;
            std::set<ttk::SimplexId> reachableMaxima;
            dataType critOrder = order[gId];
            reachableMaxima.insert(ascendingManifold[gId]);
            for(int j = 0; j < nNeighbors; j++) {
              triangulation->getVertexNeighbor(gId, j, neighborId);
              reachableMaxima.insert(ascendingManifold[neighborId]);
              if(order[neighborId] < critOrder)
                lowerLink.emplace_back(neighborId);
              if(order[neighborId] > critOrder)
                upperLink.emplace_back(neighborId);
            }
            // TODO: for join tree?? for Split tree invert upper and lower link
            // and check for reachable minima
            if(reachableMaxima.size() > 1
               && !(upperLink.size() == 0 && lowerLink.size() > 0)) {
              this->printMsg("#Upperlink: " + std::to_string(upperLink.size())
                             + ", #lowerlink: "
                             + std::to_string(lowerLink.size()));
              this->printMsg("Added " + std::to_string(gId) + " to saddles.");
              saddles.emplace(gId, order[i]);
            }
          }
        }

        std::map<ttk::SimplexId, std::set<ttk::SimplexId>> pathLists{};
        std::map<ttk::SimplexId, std::set<ttk::SimplexId>> maximumLists{};

        findAscPaths<dataType, triangulationType>(pathLists, maximumLists, saddles, ascendingManifold, triangulation);

        constructJoinTree<dataType>(
          joinTree, pathLists, maximumLists, maxima, saddles, globalMin);

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
