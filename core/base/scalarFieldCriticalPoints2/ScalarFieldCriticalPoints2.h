#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>


namespace ttk {

  class ScalarFieldCriticalPoints2 : virtual public Debug {

  public:
    struct CriticalPoint {
      SimplexId idx;
      unsigned char type;
      CriticalPoint(){}
      CriticalPoint(SimplexId idx_, unsigned char type_): idx(idx_),type(type_){}
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
      std::vector<CriticalPoint>& criticalPoints,
      const SimplexId* order,
      const SimplexId* ascManifold,
      const SimplexId* desManifold,
      const TT* triangulation
    ) const {
      ttk::Timer timer;

      const std::string msg{"Computing Critical Points"};
      this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const SimplexId nVertices = triangulation->getNumberOfVertices();

      std::vector<std::vector<CriticalPoint>> criticalPointsPerThread(this->threadNumber_);
      #pragma omp parallel num_threads(this->threadNumber_)
      {
        int threadId = omp_get_thread_num();
        auto& criticalPointsOfThread = criticalPointsPerThread[threadId];

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

          unsigned char type = 1;
          // check if min
          if(lowerLinkVertices.size()<1)
            type *= 2;

          // check if max
          if(upperLinkVertices.size()<1)
            type *= 3;

          // check if lower and upper link lead to more than one extremum
          const int numberOfReachableMinima = this->computeNumberOfReachableExtrema(
            lowerLinkVertices,
            ascManifold
          );
          if(numberOfReachableMinima>1 && this->computeNumberOfLinkComponents(lowerLinkVertices,triangulation)>1)
            type *= 5;

          const int numberOfReachableMaxima = this->computeNumberOfReachableExtrema(
            upperLinkVertices,
            desManifold
          );
          if(numberOfReachableMaxima>1 && this->computeNumberOfLinkComponents(upperLinkVertices,triangulation)>1)
            type *= 7;

          if(type>1)
            criticalPointsOfThread.emplace_back(v,type);
        }

        #pragma omp barrier

        #pragma omp single
        {
          SimplexId nCriticalPoints = 0;
          for(int t=0; t<this->threadNumber_; t++)
            nCriticalPoints += criticalPointsPerThread[t].size();
          criticalPoints.resize(nCriticalPoints);

          this->printMsg(msg, 0.8, timer.getElapsedTime(), this->threadNumber_, ttk::debug::LineMode::REPLACE);
        }

        // compute offset
        SimplexId offset = 0;
        for(int t=0; t<threadId; t++){
          offset += criticalPointsPerThread[t].size();
        }

        // write to output
        SimplexId nCriticalPointsOfThread = criticalPointsOfThread.size();
        for(int c=0; c<nCriticalPointsOfThread; c++){
          criticalPoints[offset++] = criticalPointsOfThread[c];
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    template <typename CT, typename TT = ttk::AbstractTriangulation>
    int computeCellAndPointArray(
      std::vector<CriticalPoint>& criticalPoints,
      float* coords,
      unsigned char* type,
      CT* offsets,
      CT* connectivity,
      const TT* triangulation
    ) const {
      ttk::Timer timer;

      const size_t nCriticalPoints = criticalPoints.size();

      const std::string msg{"Computing Output ("+std::to_string(nCriticalPoints)+")"};
      this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);


      #pragma omp parallel for num_threads(this->threadNumber_)
      for(size_t i = 0; i < nCriticalPoints; i++) {
        const auto &cp = criticalPoints[i];
        const size_t offset = i*3;
        triangulation->getVertexPoint(cp.idx, *(coords+offset), *(coords+offset+1), *(coords+offset+2));
        type[i] = cp.type;
        offsets[i] = i;
        connectivity[i] = i;
      }
      offsets[nCriticalPoints] = nCriticalPoints;

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

  }; // ScalarFieldCriticalPoints2 class
} // namespace ttk
