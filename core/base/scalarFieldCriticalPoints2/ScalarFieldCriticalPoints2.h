#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>


namespace ttk {

  class ScalarFieldCriticalPoints2 : virtual public Debug {

  public:
    ScalarFieldCriticalPoints2(){
      this->setDebugMsgPrefix("ScalarFieldCriticalPoints2");
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    int computeNumberOfReachableExtrema(
      const std::vector<SimplexId>& linkVertices,
      const int nLinkVertices,
      const SimplexId* manifold
    ) const {
      if(nLinkVertices<1)
        return 0;

      int numberOfReachableExtrema = 0;
      SimplexId extremumID = -1;
      for(int i=0; i<nLinkVertices; i++){
        const auto& extremumID_= manifold[linkVertices[i]];
        if(extremumID!=extremumID_){
          if(numberOfReachableExtrema>0){
            return 2;
          } else {
            numberOfReachableExtrema=1;
            extremumID=extremumID_;
          }
        }
      }

      return numberOfReachableExtrema;
    }

    template <typename TT = ttk::AbstractTriangulation>
    int computeNumberOfLinkComponents(
      std::vector<SimplexId>& linkVertices,
      const int nLinkVertices,
      const TT* triangulation
    ) const {

      // compute map
      std::unordered_map<SimplexId,SimplexId> linkVerticesMap;
      for(int i=0; i<nLinkVertices; i++){
        const SimplexId v = linkVertices[i];
        linkVerticesMap.insert({v,v});
      }

      // compute link edges
      for(int i=0; i<nLinkVertices; i++){
        const SimplexId vId = linkVertices[i];

        const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(vId);
        for(SimplexId n=0; n<nNeighbors; n++){
          SimplexId uId = -1;
          triangulation->getVertexNeighbor(vId, n, uId);

          // only consider edges in one direction
          if(uId<vId)
            continue;

          // only consider edges that are part of the link
          auto u = linkVerticesMap.find(uId);
          if(u == linkVerticesMap.end())
            continue;

          auto v = linkVerticesMap.find(vId);

          // find
          while(u->first!=u->second){
            u = linkVerticesMap.find(u->second);
          }
          while(v->first!=v->second){
            v = linkVerticesMap.find(v->second);
          }

          // union
          u->second = v->first;
        }
      }

      // count components
      int nComponents = 0;
      for(auto kv : linkVerticesMap)
        if(kv.first==kv.second)
          nComponents++;

      return nComponents;
    }

    template <typename TT = ttk::AbstractTriangulation>
    int computeCritialPoints(
      std::vector<std::vector<std::vector<SimplexId>>>& cp,
      const SimplexId* order,
      const SimplexId* ascManifold,
      const SimplexId* desManifold,
      const TT* triangulation
    ) const {
      ttk::Timer timer;

      const std::string msg{"Computing Critical Points"};
      this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const SimplexId nVertices = triangulation->getNumberOfVertices();

      #pragma omp parallel num_threads(this->threadNumber_)
      {
        const int threadId = omp_get_thread_num();
        const int nThreads = omp_get_num_threads();

        #pragma omp single
        {
          for(int i=0; i<4; i++)
            cp[i].resize(nThreads);
        }

        auto& cp0 = cp[0][threadId];
        auto& cp1 = cp[1][threadId];
        auto& cp2 = cp[2][threadId];
        auto& cp3 = cp[3][threadId];

        std::vector<SimplexId> lowerLinkVertices(32); // room for max 32 vertices
        std::vector<SimplexId> upperLinkVertices(32); // room for max 32 vertices
        int nLowerLinkVertices=0;
        int nUpperLinkVertices=0;

        #pragma omp for
        for(SimplexId v=0; v<nVertices; v++){

          // std::unordered_map<SimplexId,SimplexId> upperLinkVertices;
          // std::unordered_map<SimplexId,SimplexId> lowerLinkVertices;

          // compute lower and upper link vertices
          nLowerLinkVertices=0;
          nUpperLinkVertices=0;

          const SimplexId orderV = order[v];
          const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(SimplexId n=0; n<nNeighbors; n++){
            SimplexId u = -1;
            triangulation->getVertexNeighbor(v, n, u);
            const SimplexId orderN = order[u];
            if(orderV<orderN){
              upperLinkVertices[nUpperLinkVertices++]=u;
            } else {
              lowerLinkVertices[nLowerLinkVertices++]=u;
            }
          }

          // unsigned char type = 1;
          // check if min
          if(nLowerLinkVertices<1)
            cp0.emplace_back(v);
            // type *= 2;

          // check if max
          if(nUpperLinkVertices<1)
            cp3.emplace_back(v);
            // type *= 3;

          // check if lower and upper link lead to more than one extremum
          const int numberOfReachableMinima = this->computeNumberOfReachableExtrema(
            lowerLinkVertices,
            nLowerLinkVertices,
            ascManifold
          );
          if(numberOfReachableMinima>1 && this->computeNumberOfLinkComponents(lowerLinkVertices,nLowerLinkVertices,triangulation)>1)
            cp1.emplace_back(v);
            // type *= 5;

          const int numberOfReachableMaxima = this->computeNumberOfReachableExtrema(
            upperLinkVertices,
            nUpperLinkVertices,
            desManifold
          );
          if(numberOfReachableMaxima>1 && this->computeNumberOfLinkComponents(upperLinkVertices,nUpperLinkVertices,triangulation)>1)
            cp2.emplace_back(v);
            // type *= 7;

          // if(type>1)
          //   criticalPointsOfThread.emplace_back(v,type);
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    // template <typename CT, typename TT = ttk::AbstractTriangulation>
    // int computeCellAndPointArray(
    //   float* coords,
    //   SimplexId* ids,
    //   CT* offsets,
    //   CT* connectivity,
    //   // SimplexId* outputOrder,
    //   // const SimplexId* inputOrder,
    //   const std::vector<std::vector<CriticalPoint>>& criticalPoints,
    //   const TT* triangulation
    // ) const {
    //   const size_t nThreads = criticalPoints.size();

    //   size_t offset3 = 0;
    //   size_t offset = 0;
    //   for(size_t i=0; i<nThreads; i++){
    //     const auto& cp_ = criticalPoints[i];
    //     const size_t n = cp_.size();
    //     for(size_t j=0; j<n; j++){
    //       const auto& cp__ = cp_[j];
    //       triangulation->getVertexPoint(cp__.idx, *(coords+offset3), *(coords+offset3+1), *(coords+offset3+2));
    //       offset3 += 3;

    //       ids[offset] = cp__.idx;
    //       offsets[offset] = offset;
    //       // outputOrder[offset] = inputOrder[cp__.idx];
    //       connectivity[offset] = offset;
    //       offset++;
    //     }
    //   }
    //   offsets[offset] = offset;

    //   return 1;
    // }

    int computeIdArray(
      SimplexId* ids,
      const std::vector<std::vector<SimplexId>>& criticalPoints
    ) const {
      const size_t nThreads = criticalPoints.size();

      size_t offset = 0;
      for(size_t t=0; t<nThreads; t++){
        const auto& cp_ = criticalPoints[t];
        const size_t n = cp_.size();
        for(size_t j=0; j<n; j++){
          ids[offset++] = cp_[j];
        }
      }

      return 1;
    }
  }; // ScalarFieldCriticalPoints2 class
} // namespace ttk
