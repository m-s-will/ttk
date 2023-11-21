/// \ingroup base
/// \class ttk::BranchDecomposition
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 29.07.2021
///
/// This module assigns to each vertex of a tracking graph a branch id based on
/// a given attribute. First, all birth nodes are assigned a unique branch id,
/// and then the algorithm iterates over every vertex in order of time and then
/// either inherits the branch id of its largest predecessor (but only if the
/// current vertex is also the largest successor of this predecessor), or the
/// vertex gets assinged a new unique branch id.
///
/// \b Related \b publications: \n
/// 'Nested Tracking Graphs'
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///
/// 'Dynamic Nested Tracking Graphs'
/// Jonas Lukasczyk, Christoph Garth, Gunther H. Weber, Tim Biedert, Ross
/// Maciejewski, and Heike Leitte. IEEE Transactions on Visualization and
/// Computer Graphics, Vol. 26, No. 1, 2020.

#pragma once

// ttk common includes
#include <Debug.h>
#include <TrackingGraph.h>
#include <Triangulation.h>

namespace ttk {

  class BranchDecomposition : virtual public Debug {

  public:
    BranchDecomposition() {
      this->setDebugMsgPrefix("BranchDecomposition");
    };

    template <typename IT, typename DT>
    int computeBranchDecompositionByAttribute(int *branchIdPoints,
                                              int *branchIdEdges,
                                              ttk::TrackingGraph &trackingGraph,
                                              const IT *time,
                                              const DT *attribute,
                                              const int attributeAssociation) {
      ttk::Timer globalTimer;
      const int nNodes = trackingGraph.inEdges.size();
      const int nEdges = trackingGraph.numberOfEdges;

      std::vector<int> nodesSortedByTime(nNodes);
      // sort all nodes by time in ascending order
      {
        ttk::Timer timer;
        const std::string msg
          = "Sorting Nodes by Time (#" + std::to_string(nNodes) + ")";
        this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

        for(int i = 0; i < nNodes; i++)
          nodesSortedByTime[i] = i;

        std::sort(nodesSortedByTime.begin(), nodesSortedByTime.end(),
                  [=](int a, int b) { return time[a] < time[b]; });

        this->printMsg(msg, 1, timer.getElapsedTime(), 1);
      }

      // Sorting Edges by attribute
      {
        ttk::Timer timer;
        const std::string msg = "Sorting Edges by Attribute";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

        const auto compareAttributeU = [=](const ttk::TrackingGraph::Edge &a,
                                           const ttk::TrackingGraph::Edge &b) {
          return attribute[a.u] > attribute[b.u];
        };
        const auto compareAttributeV = [=](const ttk::TrackingGraph::Edge &a,
                                           const ttk::TrackingGraph::Edge &b) {
          return attribute[a.v] > attribute[b.v];
        };
        const auto compareAttributeE = [=](const ttk::TrackingGraph::Edge &a,
                                           const ttk::TrackingGraph::Edge &b) {
          return attribute[a.e] > attribute[b.e];
        };

        if(attributeAssociation == 0) {
#pragma omp parallel for num_threads(this->threadNumber_)
          for(int i = 0; i < nNodes; i++) {
            auto &outEdges = trackingGraph.outEdges[i];
            auto &inEdges = trackingGraph.inEdges[i];
            std::sort(outEdges.begin(), outEdges.end(), compareAttributeV);
            std::sort(inEdges.begin(), inEdges.end(), compareAttributeU);
          }
        } else {
#pragma omp parallel for num_threads(this->threadNumber_)
          for(int i = 0; i < nNodes; i++) {
            auto &outEdges = trackingGraph.outEdges[i];
            auto &inEdges = trackingGraph.inEdges[i];
            std::sort(outEdges.begin(), outEdges.end(), compareAttributeE);
            std::sort(inEdges.begin(), inEdges.end(), compareAttributeE);
          }
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      // Initializing Branches
      {
        ttk::Timer timer;
        const std::string msg = "Initializing Branches";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

#pragma omp parallel num_threads(this->threadNumber_)
        {
#pragma omp for
          for(int i = 0; i < nNodes; i++) {
            branchIdPoints[i] = trackingGraph.inEdges[i].size() < 1 ? i : -1;
          }

#pragma omp for
          for(int i = 0; i < nEdges; i++) {
            branchIdEdges[i] = -1;
          }
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      // Propagating Branches
      {
        ttk::Timer timer;
        const std::string msg = "Propagating Branches";
        this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

        // propagate branch id along points
        for(int i = 0; i < nNodes; i++) {
          const auto &v = nodesSortedByTime[i];
          if(branchIdPoints[v] != -1)
            continue;

          branchIdPoints[v] = v;

          // propagate id only if max incoming edge is max outgoing edge
          if(trackingGraph.inEdges.size() >= 1) {
            const auto &maxIncomingEdge = trackingGraph.inEdges[v][0];
            if(trackingGraph.outEdges[maxIncomingEdge.u][0].v == v) {
              branchIdPoints[v] = branchIdPoints[maxIncomingEdge.u];
              branchIdEdges[maxIncomingEdge.e] = branchIdPoints[v];
            }
          }
        }

        for(int i = 0; i < nNodes; i++) {
          const auto &v = nodesSortedByTime[i];
          const auto &inEdges = trackingGraph.inEdges[v];
          if(inEdges.size() == 1)
            branchIdEdges[inEdges[0].e] = branchIdPoints[inEdges[0].v];
          const auto &outEdges = trackingGraph.outEdges[v];
          if(outEdges.size() == 1)
            branchIdEdges[outEdges[0].e] = branchIdPoints[outEdges[0].u];
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      this->printMsg(ttk::debug::Separator::L2);
      this->printMsg("Complete", 1, globalTimer.getElapsedTime());
      this->printMsg(ttk::debug::Separator::L1);

      return 1;
    }
  }; // BranchDecomposition class

} // namespace ttk