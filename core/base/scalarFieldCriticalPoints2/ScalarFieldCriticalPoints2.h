#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>


namespace ttk {

  class ScalarFieldCriticalPoints2 : virtual public Debug {

  public:
    struct CriticalPoint {
      SimplexId idx;
      unsigned char reachableExtremaCount;
      CriticalPoint(){}
      CriticalPoint(SimplexId idx_, unsigned char reachableExtremaCount_): idx(idx_),reachableExtremaCount(reachableExtremaCount_){}
    };

    ScalarFieldCriticalPoints2(){
      this->setDebugMsgPrefix("ScalarFieldCriticalPoints2");
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    int computeNumberOfReachableExtrema(
      const std::unordered_map<SimplexId,SimplexId>& linkVertices,
      const SimplexId* manifold
    ) const {
      if(linkVertices.size()<1)
        return 0;

      int numberOfReachableExtrema = 0;
      SimplexId extremumID = -1;
      for(auto kv : linkVertices) {
        const auto& extremumID_= manifold[kv.first];
        if(extremumID!=extremumID_){
          if(numberOfReachableExtrema>0){
            return 2;
          } else {
            numberOfReachableExtrema++;
            extremumID=extremumID_;
          }
        }
      }

      return numberOfReachableExtrema;
    }

    template <typename TT = ttk::AbstractTriangulation>
    int computeNumberOfLinkComponents(
      std::unordered_map<SimplexId,SimplexId>& linkVertices,
      const TT* triangulation
    ) const {
      // compute link edges
      std::vector<std::pair<SimplexId,SimplexId>> linkEdges;
      for(auto kv : linkVertices) {
        const SimplexId v = kv.first;
        const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(v);
        for(SimplexId n=0; n<nNeighbors; n++){
          SimplexId u = -1;
          triangulation->getVertexNeighbor(v, n, u);

          // only add edges in one direction
          if(u<v)
            continue;

          auto search = linkVertices.find(u);
          if(search != linkVertices.end()){
            linkEdges.emplace_back(v,u);
          }
        }
      }

      // compute union find
      const SimplexId nEdges = linkEdges.size();
      for(SimplexId e=0; e<nEdges; e++){
        const auto& edge = linkEdges[e];
        auto u = linkVertices.find(edge.first);
        auto v = linkVertices.find(edge.second);

        // find
        while(u->first!=u->second){
          u = linkVertices.find(u->second);
        }
        while(v->first!=v->second){
          v = linkVertices.find(v->second);
        }

        // union
        u->second = v->first;
        // v->second = u->first;
      }

      // count components
      int nComponents = 0;
      for(auto kv : linkVertices)
        if(kv.first==kv.second)
          nComponents++;

      return nComponents;
    }

    template <typename TT = ttk::AbstractTriangulation>
    int computeCritialPoints(
      std::vector<std::vector<std::vector<CriticalPoint>>>& cp,
      const SimplexId* order,
      const SimplexId* ascManifold,
      const SimplexId* desManifold,
      const TT* triangulation
    ) const {
      ttk::Timer timer;

      const std::string msg{"Computing Critical Points"};
      this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const SimplexId nVertices = triangulation->getNumberOfVertices();

      std::vector<std::vector<std::vector<CriticalPoint>>> cp_(4);
      for(int i=0; i<4; i++)
        cp_.resize(this->threadNumber_);

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

        #pragma omp for
        for(SimplexId v=0; v<nVertices; v++){
          // compute lower and upper link vertices
          std::unordered_map<SimplexId,SimplexId> upperLinkVertices;
          std::unordered_map<SimplexId,SimplexId> lowerLinkVertices;
          const SimplexId orderV = order[v];
          const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(SimplexId n=0; n<nNeighbors; n++){
            SimplexId u = -1;
            triangulation->getVertexNeighbor(v, n, u);
            const SimplexId orderN = order[u];
            if(orderV<orderN){
              upperLinkVertices.insert({u,u});
            } else {
              lowerLinkVertices.insert({u,u});
            }
          }

          // unsigned char type = 1;
          // check if min
          if(lowerLinkVertices.size()<1)
            cp0.emplace_back(v,0);
            // type *= 2;

          // check if max
          if(upperLinkVertices.size()<1)
            cp3.emplace_back(v,0);
            // type *= 3;

          // check if lower and upper link lead to more than one extremum
          const int numberOfReachableMinima = this->computeNumberOfReachableExtrema(
            lowerLinkVertices,
            ascManifold
          );
          if(numberOfReachableMinima>1 && this->computeNumberOfLinkComponents(lowerLinkVertices,triangulation)>1)
            cp1.emplace_back(v,numberOfReachableMinima);
            // type *= 5;

          const int numberOfReachableMaxima = this->computeNumberOfReachableExtrema(
            upperLinkVertices,
            desManifold
          );
          if(numberOfReachableMaxima>1 && this->computeNumberOfLinkComponents(upperLinkVertices,triangulation)>1)
            cp2.emplace_back(v,numberOfReachableMaxima);
            // type *= 7;

          // if(type>1)
          //   criticalPointsOfThread.emplace_back(v,type);
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename CT, typename TT = ttk::AbstractTriangulation>
    int computeCellAndPointArray(
      float* coords,
      SimplexId* ids,
      CT* offsets,
      CT* connectivity,
      // SimplexId* outputOrder,
      // const SimplexId* inputOrder,
      const std::vector<std::vector<CriticalPoint>>& criticalPoints,
      const TT* triangulation
    ) const {
      const size_t nThreads = criticalPoints.size();

      size_t offset3 = 0;
      size_t offset = 0;
      for(size_t i=0; i<nThreads; i++){
        const auto& cp_ = criticalPoints[i];
        const size_t n = cp_.size();
        for(size_t j=0; j<n; j++){
          const auto& cp__ = cp_[j];
          triangulation->getVertexPoint(cp__.idx, *(coords+offset3), *(coords+offset3+1), *(coords+offset3+2));
          offset3 += 3;

          ids[offset] = cp__.idx;
          offsets[offset] = offset;
          // outputOrder[offset] = inputOrder[cp__.idx];
          connectivity[offset] = offset;
          offset++;
        }
      }
      offsets[offset] = offset;

      return 1;
    }

  }; // ScalarFieldCriticalPoints2 class
} // namespace ttk
