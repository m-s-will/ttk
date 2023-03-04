#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <math.h>

#include "lut.cpp"

#include <bitset>
// #include <fstream>

namespace ttk {

  class ScalarFieldCriticalPoints2 : virtual public Debug {

  public:
    ScalarFieldCriticalPoints2(){
      this->setDebugMsgPrefix("ScalarFieldCriticalPoints2");
    }

    template <typename TT = ttk::AbstractTriangulation>
    int computeLookUpTable(
      const TT* triangulation
    ) const{

      const SimplexId nVertices = triangulation->getNumberOfVertices();

      constexpr int lutSize = pow(2,14);
      std::array<unsigned char,lutSize> lut;

      for(int v=0; v<nVertices; v++){
        if(triangulation->getVertexNeighborNumber(v)==14){

          #pragma omp parallel for
          for(long unsigned int i=0; i<lutSize; i++){
            std::bitset<14> binary(i);
            // std::cout<<binary[0]<<" "<<binary<<std::endl;

            std::array<SimplexId,32> linkVertices;
            int nLinkVertices = 0;

            SimplexId u = -1;
            for(SimplexId n=0; n<14; n++){
              triangulation->getVertexNeighbor(v, n, u);

              if(binary[n]==1)
                linkVertices[nLinkVertices++]=u;
            }

            const int nComponents = this->computeNumberOfLinkComponents(
              linkVertices,
              nLinkVertices,
              triangulation
            );

            lut[i] = nComponents < 2;
          }

          // std::cout<<"writing"<<std::endl;

          // std::ofstream lutFile;
          // lutFile.open("/home/jones/external/projects/ttk/core/base/scalarFieldCriticalPoints2/lut.cpp");
          // for(auto x: lut){
          //   lutFile << std::to_string(x)<<',';
          // }
          // lutFile.close();

          break;
        }
      }
      return 1;
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    int computeNumberOfReachableExtrema(
      const std::array<SimplexId,32>& linkVertices,
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
      std::array<SimplexId,32>& linkVertices,
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
      std::array<std::vector<std::vector<SimplexId>>,4>& cp,
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

        std::array<SimplexId,32> lowerLinkVertices; // room for max 32 vertices
        std::array<SimplexId,32> upperLinkVertices; // room for max 32 vertices
        int nLowerLinkVertices=0;
        int nUpperLinkVertices=0;

        #pragma omp for
        for(SimplexId v=0; v<nVertices; v++){

          const SimplexId orderV = order[v];
          const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(v);

          if(nNeighbors==14){
            std::bitset<14> upperLinkKey;
            std::bitset<14> lowerLinkKey;

            SimplexId minId = -1;
            int nMin = 0;
            SimplexId maxId = -1;
            int nMax = 0;

            for(SimplexId n=0; n<14; n++){
              SimplexId u = -1;
              triangulation->getVertexNeighbor(v, n, u);
              const SimplexId orderN = order[u];
              if(orderV<orderN){
                upperLinkKey[n]=1;
                if(maxId != desManifold[u]){
                  maxId = desManifold[u];
                  nMax++;
                }
              } else {
                lowerLinkKey[n]=1;
                if(minId != ascManifold[u]){
                  minId = ascManifold[u];
                  nMin++;
                }
              }
            }

            if(nMin==0)
              cp0.emplace_back(v);

            if(nMax==0)
              cp3.emplace_back(v);

            if(nMin>1 && lut[lowerLinkKey.to_ulong()]==0){
              cp1.emplace_back(v);
            }
            if(nMax>1 && lut[upperLinkKey.to_ulong()]==0){
              cp2.emplace_back(v);
            }
          } else {
            // compute lower and upper link vertices
            nLowerLinkVertices=0;
            nUpperLinkVertices=0;

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

            // check if min
            if(nLowerLinkVertices<1)
              cp0.emplace_back(v);

            // check if max
            if(nUpperLinkVertices<1)
              cp3.emplace_back(v);

            // check if lower and upper link lead to more than one extremum
            const int numberOfReachableMinima = this->computeNumberOfReachableExtrema(
              lowerLinkVertices,
              nLowerLinkVertices,
              ascManifold
            );
            if(numberOfReachableMinima>1 && this->computeNumberOfLinkComponents(lowerLinkVertices,nLowerLinkVertices,triangulation)>1)
              cp1.emplace_back(v);

            const int numberOfReachableMaxima = this->computeNumberOfReachableExtrema(
              upperLinkVertices,
              nUpperLinkVertices,
              desManifold
            );
            if(numberOfReachableMaxima>1 && this->computeNumberOfLinkComponents(upperLinkVertices,nUpperLinkVertices,triangulation)>1)
              cp2.emplace_back(v);
          }
        }
      }

      // std::cout<<n0<<" "<<n1<<" "<<n2<<" "<<n3<<std::endl;

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
