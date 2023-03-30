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

    struct CP {
      SimplexId id;
      std::vector<SimplexId> neighbors;
      CP(const SimplexId& _) : id(_){
      }
    };

    ScalarFieldCriticalPoints2(){
      this->setDebugMsgPrefix("ScalarFieldCriticalPoints2");
    }

    void toXYZ(ttk::SimplexId* xyz, const ttk::SimplexId idx) const {
      const auto d = 3 * 3;
      xyz[2] = idx / d;
      const auto idx2 = idx - (xyz[2] * d);
      xyz[1] = idx2 / 3;
      xyz[0] = idx2 % 3;
    };

    template <typename TT = ttk::AbstractTriangulation>
    int computeLookUpTable(
      const TT* triangulation
    ) const{

      const SimplexId nVertices = triangulation->getNumberOfVertices();

      //constexpr int lutSize = pow(2,14);
      constexpr int lutSize = 16384;
      std::array<unsigned char,lutSize> lut;

      std::array<int,64> lut2;
      for(int i=0; i<64; i++){
        lut2[i] = 0;
      }

      const std::array<ttk::SimplexId,64> nNeighborsLUT{
      14,10,10,0,10,6,8,0,10,8,6,0,0,0,0,0,10,6,8,0,8,4,7,0,6,4,4,0,0,0,0,0,10,8,6,0,6,4,4,0,8,7,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
      };
      constexpr int offsetsLutSize = 64*14*3;
      std::array<ttk::SimplexId,offsetsLutSize> offsetsLUT;
      for(int i=0; i<offsetsLutSize; i++)
        offsetsLUT[i] = 0;

      for(int v=0; v<nVertices; v++){
        std::cout<<v<<std::endl;
        std::cout<<" n: "<<triangulation->getVertexNeighborNumber(v)<<std::endl;

        ttk::SimplexId xyz[3];
        this->toXYZ(xyz,v);
        // build case key
        int key =
          (xyz[0]==0?1:xyz[0]==2?2:0)+
          (xyz[1]==0?4:xyz[1]==2?8:0)+
          (xyz[2]==0?16:xyz[2]==2?32:0)
        ;
        lut2[key] = triangulation->getVertexNeighborNumber(v);

        std::cout<<" [";
        int offsetCursor = key*3*14;
        for(SimplexId n=0; n<triangulation->getVertexNeighborNumber(v); n++){
          SimplexId u = -1;
          ttk::SimplexId nxyz[3];
          triangulation->getVertexNeighbor(v, n, u);
          this->toXYZ(nxyz,u);

          offsetsLUT[offsetCursor]=nxyz[0]-xyz[0];
          offsetsLUT[offsetCursor+1]=nxyz[1]-xyz[1];
          offsetsLUT[offsetCursor+2]=nxyz[2]-xyz[2];

          offsetCursor+=3;

          std::cout<<" "<<u<<"("<<nxyz[0]<<","<<nxyz[1]<<","<<nxyz[2]<<")";
        }
        std::cout<<" ]"<<std::endl;
      }

      std::cout<<"{"<<std::endl;
     for(int i=0; i<64; i++){
        std::cout<<lut2[i]<<",";
      }
      std::cout<<std::endl<<"}"<<std::endl;
      std::cout<<"{"<<std::endl;
     for(int i=0; i<offsetsLutSize; i+=3){
        std::cout<<offsetsLUT[i]<<","<<offsetsLUT[i+1]<<","<<offsetsLUT[i+2]<<",";
      }
      std::cout<<std::endl<<"}"<<std::endl;

      // for(int v=0; v<nVertices; v++){
      //   if(triangulation->getVertexNeighborNumber(v)==14){

      //     #pragma omp parallel for
      //     for(long unsigned int i=0; i<lutSize; i++){
      //       std::bitset<14> binary(i);
      //       // std::cout<<binary[0]<<" "<<binary<<std::endl;

      //       std::array<SimplexId,32> linkVertices;
      //       int nLinkVertices = 0;

      //       SimplexId u = -1;
      //       for(SimplexId n=0; n<14; n++){
      //         triangulation->getVertexNeighbor(v, n, u);

      //         if(binary[n]==1)
      //           linkVertices[nLinkVertices++]=u;
      //       }

      //       const int nComponents = this->computeNumberOfLinkComponents(
      //         linkVertices,
      //         nLinkVertices,
      //         triangulation
      //       );

      //       lut[i] = nComponents < 2;
      //     }

      //     // std::cout<<"writing"<<std::endl;

      //     // std::ofstream lutFile;
      //     // lutFile.open("/home/jones/external/projects/ttk/core/base/scalarFieldCriticalPoints2/lut.cpp");
      //     // for(auto x: lut){
      //     //   lutFile << std::to_string(x)<<',';
      //     // }
      //     // lutFile.close();

      //     break;
      //   }
      // }
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
      if(nLinkVertices<2)
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
      const std::array<SimplexId,32>& linkVertices,
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

        // general case
        std::array<SimplexId,32> lowerMem; // room for max 32 vertices
        std::array<SimplexId,32> upperMem; // room for max 32 vertices
        std::array<SimplexId,32> lowerMem2; // used by general case
        std::array<SimplexId,32> upperMem2; // used by general case

        #pragma omp for
        for(SimplexId v=0; v<nVertices; v++){
          const SimplexId orderV = order[v];
          const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(v);

          if(nNeighbors==14){
            std::bitset<14> upperLinkKey;
            std::bitset<14> lowerLinkKey;

            lowerMem[0] = -1;
            upperMem[0] = -1;
            int lowerCursor=0;
            int upperCursor=0;

            for(SimplexId n=0; n<14; n++){
              SimplexId u{0};
              triangulation->getVertexNeighbor(v, n, u);

              const SimplexId orderN = order[u];
              if(orderV<orderN){
                upperLinkKey[n]=1;
                if(upperMem[upperCursor] != desManifold[u]){
                  upperMem[++upperCursor] = desManifold[u];
                }
              } else {
                lowerLinkKey[n]=1;
                if(lowerMem[lowerCursor] != ascManifold[u]){
                  lowerMem[++lowerCursor] = ascManifold[u];
                }
              }
            }

            if(lowerCursor==0)
              cp0.emplace_back(v);

            if(upperCursor==0)
              cp3.emplace_back(v);

            // if lowerCursor or upperCursor >1 then more than one extremum reachable
            if(lowerCursor>1 && lut[lowerLinkKey.to_ulong()]==0)
              cp1.emplace_back(v);

            if(upperCursor>1 && lut[upperLinkKey.to_ulong()]==0)
              cp2.emplace_back(v);
          } else {
            lowerMem[0] = -1;
            upperMem[0] = -1;
            int lowerCursor=0;
            int upperCursor=0;
            int lowerCursor2=0;
            int upperCursor2=0;

            for(SimplexId n=0; n<nNeighbors; n++){
              SimplexId u{0};
              triangulation->getVertexNeighbor(v, n, u);

              const SimplexId orderN = order[u];
              if(orderV<orderN){
                upperMem2[upperCursor2++] = u;
                if(upperMem[upperCursor] != desManifold[u]){
                  upperMem[++upperCursor] = desManifold[u];
                }
              } else {
                lowerMem2[lowerCursor2++] = u;
                if(lowerMem[lowerCursor] != ascManifold[u]){
                  lowerMem[++lowerCursor] = ascManifold[u];
                }
              }
            }

            if(lowerCursor==0)
              cp0.emplace_back(v);

            if(upperCursor==0)
              cp3.emplace_back(v);

            if(lowerCursor>1 && this->computeNumberOfLinkComponents(lowerMem2,lowerCursor2,triangulation)>1)
              cp1.emplace_back(v);

            if(upperCursor>1 && this->computeNumberOfLinkComponents(upperMem2,upperCursor2,triangulation)>1)
              cp2.emplace_back(v);
          }
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    int processSaddle(
      std::vector<CP>& cps,
      std::array<SimplexId,32>& neighbors,
      const SimplexId v,
      const int nNeighbors
    ) const {
      cps.emplace_back(v);
      auto& cp = cps.back();
      std::sort(
        neighbors.begin()+1, // skip -1 at first position
        neighbors.begin()+nNeighbors+1
      );
      cp.neighbors.reserve(nNeighbors+1);
      cp.neighbors.emplace_back(neighbors[1]);
      for(int i=2; i<=nNeighbors; i++){
        if(cp.neighbors.back()!=neighbors[i]){
          cp.neighbors.emplace_back( neighbors[i] );
        }
      }
      return 1;
    };


    template <typename TT = ttk::AbstractTriangulation>
    int computeCritialPoints(
      std::array<std::vector<std::vector<CP>>,4>& cp,
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

        // general case
        std::array<SimplexId,32> lowerMem; // room for max 32 vertices
        std::array<SimplexId,32> upperMem; // room for max 32 vertices
        std::array<SimplexId,32> lowerMem2; // used by general case
        std::array<SimplexId,32> upperMem2; // used by general case

        #pragma omp for
        for(SimplexId v=0; v<nVertices; v++){
          const SimplexId orderV = order[v];
          const SimplexId nNeighbors = triangulation->getVertexNeighborNumber(v);

          if(nNeighbors==14){
            std::bitset<14> upperLinkKey;
            std::bitset<14> lowerLinkKey;

            lowerMem[0] = -1;
            upperMem[0] = -1;
            int lowerCursor=0;
            int upperCursor=0;

            for(SimplexId n=0; n<14; n++){
              SimplexId u{0};
              triangulation->getVertexNeighbor(v, n, u);

              const SimplexId orderN = order[u];
              if(orderV<orderN){
                upperLinkKey[n]=1;
                if(upperMem[upperCursor] != desManifold[u]){
                  upperMem[++upperCursor] = desManifold[u];
                }
              } else {
                lowerLinkKey[n]=1;
                if(lowerMem[lowerCursor] != ascManifold[u]){
                  lowerMem[++lowerCursor] = ascManifold[u];
                }
              }
            }

            if(lowerCursor==0)
              cp0.emplace_back(v);

            if(upperCursor==0)
              cp3.emplace_back(v);

            if(lowerCursor>1 && lut[lowerLinkKey.to_ulong()]==0){
              this->processSaddle(cp1,lowerMem,v,lowerCursor);
            }

            if(upperCursor>1 && lut[upperLinkKey.to_ulong()]==0){
              this->processSaddle(cp2,upperMem,v,upperCursor);
            }
          } else {
            lowerMem[0] = -1;
            upperMem[0] = -1;
            int lowerCursor=0;
            int upperCursor=0;
            int lowerCursor2=0;
            int upperCursor2=0;

            for(SimplexId n=0; n<nNeighbors; n++){
              SimplexId u{0};
              triangulation->getVertexNeighbor(v, n, u);

              const SimplexId orderN = order[u];
              if(orderV<orderN){
                upperMem2[upperCursor2++] = u;
                if(upperMem[upperCursor] != desManifold[u]){
                  upperMem[++upperCursor] = desManifold[u];
                }
              } else {
                lowerMem2[lowerCursor2++] = u;
                if(lowerMem[lowerCursor] != ascManifold[u]){
                  lowerMem[++lowerCursor] = ascManifold[u];
                }
              }
            }

            if(lowerCursor==0)
              cp0.emplace_back(v);

            if(upperCursor==0)
              cp3.emplace_back(v);

            if(lowerCursor>1 && this->computeNumberOfLinkComponents(lowerMem2,lowerCursor2,triangulation)>1){
              this->processSaddle(cp1,lowerMem,v,lowerCursor);
            }

            if(upperCursor>1 && this->computeNumberOfLinkComponents(upperMem2,upperCursor2,triangulation)>1){
              this->processSaddle(cp2,upperMem,v,upperCursor);
            }
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

    int mergeCriticalPointVectors(
      std::array<std::vector<SimplexId>,4>& criticalPoints,
      const std::array<std::vector<std::vector<SimplexId>>,4>& criticalPoints_
    ) const {
      ttk::Timer timer;
      const std::string msg{"Merging Critical Point Vectors"};
      this->printMsg(msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      #pragma omp parallel for schedule(static,1) num_threads(this->threadNumber_)
      for(int b=0; b<4; b++){
        const auto& cp = criticalPoints_[b];
        const size_t nVectors = cp.size();
        size_t nCriticalPoints = 0;
        for(size_t i=0; i<nVectors; i++)
          nCriticalPoints += cp[i].size();

        criticalPoints[b].resize(nCriticalPoints);

        this->computeIdArray(
          criticalPoints[b].data(),
          cp
        );
      }
      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      this->printMsg({
        {"#Minima",std::to_string(criticalPoints[0].size())},
        {"#1-Saddle",std::to_string(criticalPoints[1].size())},
        {"#2-Saddle",std::to_string(criticalPoints[2].size())},
        {"#Maxima",std::to_string(criticalPoints[3].size())}
      });
      return 1;
    }

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
