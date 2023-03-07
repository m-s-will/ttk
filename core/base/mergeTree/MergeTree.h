#pragma once

// ttk common includes
#include <Debug.h>
#include <Propagation.h>
#include <Triangulation.h>

#if(defined(__GNUC__) && !defined(__clang__))
#include <parallel/algorithm>
#endif

// for numerical perturbation
#include <boost/math/special_functions/next.hpp>

namespace ttk::mt {

  class MergeTree : virtual public Debug {

  public:
    MergeTree() {
      this->setDebugMsgPrefix("MergeTree");
    };
    ~MergeTree(){};

    int PreconditionTriangulation(ttk::Triangulation *triangulation) const {
      triangulation->preconditionVertexNeighbors();
      return 1;
    }

    template <typename DT, typename IT>
    int computeOrderArray(IT *orderArray,

                          const IT &nVertices,
                          const DT *rank1,
                          const IT *rank2 = nullptr) const {
      ttk::Timer timer;

      this->printMsg("Computing Order Array", 0, timer.getElapsedTime(),
                     this->threadNumber_, debug::LineMode::REPLACE);

      // init tuples
      std::vector<std::tuple<DT, IT, IT>> indices(nVertices);
      if(rank2 == nullptr) {
#pragma omp parallel for num_threads(this->threadNumber_)
        for(IT i = 0; i < nVertices; i++) {
          auto &t = indices[i];
          std::get<0>(t) = rank1[i];
          std::get<1>(t) = i;
          std::get<2>(t) = i;
        }
      } else {
#pragma omp parallel for num_threads(this->threadNumber_)
        for(IT i = 0; i < nVertices; i++) {
          auto &t = indices[i];
          std::get<0>(t) = rank1[i];
          std::get<1>(t) = rank2[i];
          std::get<2>(t) = i;
        }
      }

      this->printMsg("Computing Order Array", 0.2, timer.getElapsedTime(),
                     this->threadNumber_, debug::LineMode::REPLACE);

#if TTK_ENABLE_OPENMP && !defined __clang__
      omp_set_num_threads(this->threadNumber_);
      __gnu_parallel::sort(indices.begin(), indices.end());
      omp_set_num_threads(1);
#else
      this->printWrn("Caution, outside GCC, sequential sort");
      std::sort(indices.begin(), indices.end());
#endif

      this->printMsg("Computing Order Array", 0.8, timer.getElapsedTime(),
                     this->threadNumber_, debug::LineMode::REPLACE);

#pragma omp parallel for num_threads(this->threadNumber_)
      for(IT i = 0; i < nVertices; i++)
        orderArray[std::get<2>(indices[i])] = i;

      this->printMsg("Computing Order Array", 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    }

    /// TODO
    /// TODO
    /// TODO
    template <typename IT, typename TT>
    int initializePropagations(std::vector<Propagation<IT>> &propagations,
                               IT *temp,

                               const TT *triangulation,
                               const IT *orderArray) const {

      ttk::Timer timer;
      this->printMsg("Initialize Propagations", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      const IT nVertices = triangulation->getNumberOfVertices();
      IT writeIndex = 0;

// find discareded maxima
#pragma omp parallel for num_threads(this->threadNumber_)
      for(IT v = 0; v < nVertices; v++) {
        // check if v has larger neighbors
        bool hasLargerNeighbor = false;
        const IT &vOrder = orderArray[v];
        const IT nNeighbors = triangulation->getVertexNeighborNumber(v);
        for(IT n = 0; n < nNeighbors; n++) {
          IT u;
          triangulation->getVertexNeighbor(v, n, u);
          if(vOrder < orderArray[u]) {
            hasLargerNeighbor = true;
            break;
          }
        }

        // if v has larger neighbors then v can not be maximum
        if(hasLargerNeighbor)
          continue;

        // get local write index for this thread
        IT localWriteIndex = 0;
#pragma omp atomic capture
        localWriteIndex = writeIndex++;

        // write maximum index
        temp[localWriteIndex] = -v;
      }

      this->printMsg("Initialize Propagations (" + std::to_string(writeIndex)
                       + "|" + std::to_string(nVertices) + ")",
                     0.5, timer.getElapsedTime(), this->threadNumber_,
                     debug::LineMode::REPLACE);

      propagations.resize(writeIndex);
      for(IT i = 0; i < writeIndex; i++) {
        auto &p = propagations[i];
        p.criticalPoints.push_back(-temp[i]);
        p.branchId = i;
      }

      this->printMsg("Initialize Propagations (" + std::to_string(writeIndex)
                       + "|" + std::to_string(nVertices) + ")",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    /// TODO
    /// TODO
    /// TODO
    template <typename DT, typename IT, typename TT>
    int finalizePropagations(
      std::vector<const Propagation<IT> *> &sortedPropagations,
      IT *segmentationIds,
      std::vector<Propagation<IT>> &propagations,

      const TT *triangulation,
      const DT *scalarArray) const {

      ttk::Timer timer;
      this->printMsg("Finalizing Propagations", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      const IT nPropagations = propagations.size();
      const IT nVertices = triangulation->getNumberOfVertices();

      sortedPropagations.resize(nPropagations, nullptr);
      std::vector<std::tuple<DT, IT, IT>> persistencePairs(nPropagations);

      for(IT i = 0; i < nPropagations; i++) {
        const auto &prop = propagations[i];
        auto &pair = persistencePairs[i];
        const auto &s0 = scalarArray[prop.criticalPoints.front()];
        const auto &s1 = scalarArray[prop.criticalPoints.back()];

        std::get<0>(pair) = s0 < s1 ? s0 - s1 : s1 - s0;
        std::get<1>(pair) = prop.criticalPoints.front();
        std::get<2>(pair) = i;
      }

      std::sort(persistencePairs.begin(), persistencePairs.end());

      for(IT i = 0; i < nPropagations; i++) {
        auto &pair = persistencePairs[i];
        auto &p = propagations[std::get<2>(pair)];
        sortedPropagations[i] = &p;
        p.branchId = i;
      }

      for(IT i = 0; i < nVertices; i++) {
        auto &segId = segmentationIds[i];
        segId = propagations[segId].branchId;
      }

      this->printMsg("Finalizing Propagations (" + std::to_string(nPropagations)
                       + "|" + std::to_string(nVertices) + ")",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename IT, typename TT>
    int computeDynamicPropagation(Propagation<IT> &propagation,
                                  IT *segmentationIds,
                                  Propagation<IT> **propagationMask,
                                  IT *queueMask,

                                  const TT *triangulation,
                                  const IT *orderArray,
                                  const IT &nActivePropagations) const {

      // pointer used to compare against representative
      auto *currentPropagation = &propagation;

      // frequently used propagation members
      IT segmentationId = currentPropagation->branchId;
      auto *queue = &currentPropagation->queue;

      // add extremumIndex to queue
      {
        const auto &extremumIndex = currentPropagation->criticalPoints[0];
        queue->emplace(orderArray[extremumIndex], extremumIndex);
        queueMask[extremumIndex] = segmentationId;
      }

      IT counter = 0;

      // room for twenty neighbors
      std::vector<IT> neighbors;
      neighbors.reserve(40);

      // grow region until prop reaches a saddle and then decide if prop should
      // continue
      IT v = -1;
      while(!queue->empty()) {
        v = std::get<1>(queue->top());
        queue->pop();

        // continue if this thread has already seen this vertex
        if(propagationMask[v])
          continue;

        // get neighbors
        const IT nNeighbors = triangulation->getVertexNeighborNumber(v);
        neighbors.resize(nNeighbors, 0);

        for(IT n = 0; n < nNeighbors; n++)
          triangulation->getVertexNeighbor(v, n, neighbors[n]);

        // add neighbors to queue AND check if v is a saddle
        bool isSaddle = false;
        IT numberOfLargerNeighbors = 0;
        IT numberOfLargerNeighborsThisThreadVisited = 0;
        const IT &orderV = orderArray[v];

        for(IT n = 0; n < nNeighbors; n++) {
          const IT &u = neighbors[n];
          const IT &orderU = orderArray[u];

          // if larger neighbor
          if(orderU > orderV) {
            numberOfLargerNeighbors++;

            auto uPropagation = propagationMask[u];
            if(!uPropagation || currentPropagation != uPropagation->find())
              isSaddle = true;
            else
              numberOfLargerNeighborsThisThreadVisited++;
          } else if(queueMask[u] != segmentationId) {
            queue->emplace(orderU, u);
            queueMask[u] = segmentationId;
          }
        }

        // if v is a saddle check if current thread is the last visitor
        if(isSaddle) {

          currentPropagation->criticalPoints.push_back(v);

          IT numberOfRegisteredLargerVertices = 0;
#pragma omp atomic capture
          {
            segmentationIds[v] -= numberOfLargerNeighborsThisThreadVisited;
            numberOfRegisteredLargerVertices = segmentationIds[v];
          }

          // if this thread did not register the last remaining larger vertices
          // then terminate propagation
          if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors - 1)
            return 1;

          // merge propagations
          for(IT n = 0; n < nNeighbors; n++) {
            currentPropagation = Propagation<IT>::unify(
              currentPropagation, propagationMask[neighbors[n]], orderArray);
          }

          queue = &currentPropagation->queue;
          segmentationId = currentPropagation->branchId;
        }

        // mark vertex as visited and continue
        propagationMask[v] = currentPropagation;
        segmentationIds[v] = segmentationId;

        if(counter++ > 100) {
          counter = 0;

          IT nActivePropagations_;

#pragma omp atomic read
          nActivePropagations_ = nActivePropagations;

          if(nActivePropagations_ == 1) {
            currentPropagation->interrupted = true;
            return 1;
          }
        }
      }

      // add global minimum (if not already added as saddle)
      if(currentPropagation->criticalPoints.back() != v)
        currentPropagation->criticalPoints.push_back(v);

      return 1;
    }

    template <typename IT, typename TT>
    int computeDynamicPropagations(IT *segmentationIds,
                                   std::vector<Propagation<IT>> &propagations,
                                   Propagation<IT> **propagationMask,
                                   IT *queueMask,

                                   const TT *triangulation,
                                   const IT *orderArray) const {
      ttk::Timer timer;

      int status = 1;
      const IT nPropagations = propagations.size();
      IT nActivePropagations = nPropagations;

      this->printMsg("Computing Dynamic Propagations ("
                       + std::to_string(nPropagations) + ")",
                     0, 0, this->threadNumber_, debug::LineMode::REPLACE);

// compute propagations
#pragma omp parallel for schedule(dynamic, 1) num_threads(this->threadNumber_)
      for(IT p = 0; p < nPropagations; p++) {
        int localStatus = this->computeDynamicPropagation<IT, TT>(
          propagations[p], segmentationIds, propagationMask, queueMask,

          triangulation, orderArray, nActivePropagations);
        if(!localStatus)
          status = 0;

#pragma omp atomic update
        nActivePropagations--;
      }
      if(!status)
        return 0;

      this->printMsg("Computing Dynamic Propagations ("
                       + std::to_string(nPropagations) + ")",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename IT, typename TT>
    int computeTrunk(IT *segmentationIds,
                     std::vector<Propagation<IT>> &propagations,
                     Propagation<IT> **propagationMask,
                     IT *mask,

                     const TT *triangulation,
                     const IT *orderArray) const {
      ttk::Timer timer;
      this->printMsg(
        "Computing trunk", 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      const IT nVertices = triangulation->getNumberOfVertices();

// reconstruct sorted vertices
#pragma omp parallel for num_threads(this->threadNumber_)
      for(IT v = 0; v < nVertices; v++)
        mask[orderArray[v]] = v;

      // find last unterminated propagation
      Propagation<IT> *currentPropagation = nullptr;
      IT segmentationId = -1;
      IT trunkIndex = nVertices - 1;

      const IT nPropagations = propagations.size();
      for(IT p = 0; p < nPropagations; p++) {
        if(propagations[p].interrupted) {
          currentPropagation = &propagations[p];
          segmentationId = currentPropagation->branchId;
          break;
        }
      }

      if(currentPropagation == nullptr) {
        this->printWrn("Empty Trunk.");
        return 1;
      }

      IT nTrunkVertices = 0;
      for(; trunkIndex >= 0; trunkIndex--) {
        const IT &v = mask[trunkIndex];

        if(propagationMask[v] != nullptr)
          continue;

        nTrunkVertices++;

        // if saddle
        if(segmentationIds[v] < -1) {

          currentPropagation->criticalPoints.push_back(v);

          // merge propagations
          const IT nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(IT n = 0; n < nNeighbors; n++) {
            IT u;
            triangulation->getVertexNeighbor(v, n, u);

            currentPropagation = Propagation<IT>::unify(
              currentPropagation, propagationMask[u], orderArray,
              false // not necessary to merge queues
            );
          }

          segmentationId = currentPropagation->branchId;
        }

        propagationMask[v] = currentPropagation;
        segmentationIds[v] = segmentationId;
      }

      // add global minimum (if not already added as saddle)
      if(currentPropagation->criticalPoints.back() != mask[0])
        currentPropagation->criticalPoints.push_back(mask[0]);

      std::stringstream vFraction;
      vFraction << std::fixed << std::setprecision(2)
                << ((float)nTrunkVertices / (float)nVertices);

      this->printMsg("Computing trunk (" + std::to_string(nTrunkVertices) + "|"
                       + vFraction.str() + ")",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    /// TODO
    /// TODO
    /// TODO
    template <typename IT>
    int allocateMemory(std::vector<Propagation<IT> *> &propagationMask,
                       std::vector<IT> &queueMask,
                       IT *segmentationIds,

                       const IT &nVertices) const {
      ttk::Timer timer;

      // =============================================================
      // allocate and init memory
      // =============================================================
      this->printMsg("Allocating Memory", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      propagationMask.resize(nVertices, nullptr);
      queueMask.resize(nVertices, -1);

#pragma omp parallel for num_threads(this->threadNumber_)
      for(IT i = 0; i < nVertices; i++)
        segmentationIds[i] = -1;

      this->printMsg(
        "Allocating Memory", 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename IT>
    int invertArray(std::vector<IT> &invertedArray,
                    const IT *orderArray,
                    const IT &nVertices) const {
      ttk::Timer timer;

      this->printMsg("Inverting Order Array", 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      invertedArray.resize(nVertices);
      const IT nVerticesM1 = nVertices - 1;
      for(IT i = 0; i < nVertices; i++)
        invertedArray[i] = nVerticesM1 - orderArray[i];

      this->printMsg("Inverting Order Array", 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    }

    /// TODO
    /// TODO
    /// TODO
    template <typename IT, typename TT>
    int computeMergeTreeSegmentation(IT *segmentationIds,
                                     std::vector<Propagation<IT>> &propagations,

                                     const TT *triangulation,
                                     const IT *orderArray,
                                     const int &type
                                     = 0 // 0: split tree, 1: join tree
    ) const {
      ttk::Timer timer;

      this->printMsg(debug::Separator::L2);

      const IT nVertices = triangulation->getNumberOfVertices();

      std::vector<IT> invertedArray;
      if(type != 0 && !invertArray<IT>(invertedArray, orderArray, nVertices))
        return 0;

      const IT *orderArray_ = type == 0 ? orderArray : invertedArray.data();

      // =============================================================
      // allocate and init memory
      // =============================================================
      std::vector<Propagation<IT> *> propagationMask;
      std::vector<IT> queueMask;
      if(!this->allocateMemory<IT>(propagationMask, queueMask, segmentationIds,

                                   nVertices))
        return 0;

      // =============================================================
      // initialize propagations
      // =============================================================
      if(!this->initializePropagations<IT, TT>(propagations, queueMask.data(),

                                               triangulation, orderArray_))
        return 0;

      // =============================================================
      // execute propagations
      // =============================================================
      if(!this->computeDynamicPropagations<IT, TT>(
           segmentationIds, propagations, propagationMask.data(),
           queueMask.data(),

           triangulation, orderArray_))
        return 0;

      // =============================================================
      // compute trunk
      // =============================================================
      if(!this->computeTrunk<IT, TT>(segmentationIds, propagations,
                                     propagationMask.data(), queueMask.data(),

                                     triangulation, orderArray_))
        return 0;

      this->printMsg(debug::Separator::L2);
      this->printMsg(
        "Complete", 1, timer.getElapsedTime(), this->threadNumber_);
      this->printMsg(debug::Separator::L2);

      return 1;
    }
  };
} // namespace ttk::mt