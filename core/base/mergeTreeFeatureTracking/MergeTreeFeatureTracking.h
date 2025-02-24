/// \ingroup base
/// \class ttk::MergeTreeFeatureTracking
/// \author XXX
/// \date 2025.
///
/// This VTK filter uses the ttk::MergeTreeFeatureTracking module to compute
/// TODO
///

#pragma once

// ttk common includes
#include <Debug.h>

#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <MergeTreeBase.h>
#include <MergeTreeDistance.h>

namespace ttk {

  /**
   * The MergeTreeFeatureTracking class provides methods to compute TODO
   */
  class MergeTreeFeatureTracking : virtual public Debug,
                                   virtual public MergeTreeBase {
  protected:
    double minMaxPairWeight_ = 1.0;
    int baseModule_ = 0;
    int branchMetric_ = 0;
    int pathMetric_ = 0;

  public:
    MergeTreeFeatureTracking() {
      this->setDebugMsgPrefix(
        "MergeTreeFeatureTracking"); // inherited from Debug: prefix will be
                                     // printed at the
      // beginning of every msg
    }
    ~MergeTreeFeatureTracking() override = default;

    /**
     * Implementation of the algorithm.
     */
    template <class dataType>
    void execute(
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<ftm::MergeTree<dataType>> &trees2,
      std::vector<
        std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>
        &outputMatchings,
      std::vector<dataType> &distances) {
      Timer t_preprocessing;
      treesNodeCorr_.resize(trees.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < trees.size(); ++i) {
        preprocessingPipeline<dataType>(
          trees[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
          baseModule_ == 0 ? branchDecomposition_ : false, useMinMaxPair_, true,
          treesNodeCorr_[i]);
      }
      printTreesStats(trees);
      if(trees2.size() != 0) {
        std::vector<std::vector<int>> trees2NodeCorr(trees2.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
        for(unsigned int i = 0; i < trees.size(); ++i) {
          preprocessingPipeline<dataType>(
            trees2[i], epsilonTree2_, epsilon2Tree2_, epsilon3Tree2_,
            baseModule_ == 0 ? branchDecomposition_ : false, useMinMaxPair_,
            true, trees2NodeCorr[i]);
        }
        printTreesStats(trees2);
      }
      printMsg("Preprocessing", 1, t_preprocessing.getElapsedTime(),
               this->threadNumber_);

      Timer t_total;
      executePara<dataType>(trees, outputMatchings, distances);
      // TODO double input
      /*if(trees2.size() != 0) {
        useDoubleInput_ = true;
        std::vector<std::vector<double>> distanceMatrix2(
          trees2.size(), std::vector<double>(trees2.size()));
        executePara<dataType>(trees2, distanceMatrix2, false);
      }*/

      for(unsigned int i = 0; i < trees.size(); ++i)
        postprocessingPipeline<dataType>(&(trees[i].tree));
      if(branchDecomposition_)
        for(unsigned int i = 0; i < outputMatchings.size(); ++i)
          convertBranchDecompositionMatching<dataType>(
            &(trees[i].tree), &(trees[i + 1].tree), outputMatchings[i]);
      printMsg("Total", 1, t_total.getElapsedTime(), this->threadNumber_);
    }

    template <class dataType>
    void executePara(
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<
        std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>
        &outputMatchings,
      std::vector<dataType> &distances,
      bool isFirstInput = true) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel num_threads(this->threadNumber_)
      {
#pragma omp single nowait
#endif
        executeParaImpl<dataType>(
          trees, outputMatchings, distances, isFirstInput);
#ifdef TTK_ENABLE_OPENMP
#pragma omp taskwait
      } // pragma omp parallel
#endif
    }

    template <class dataType>
    void executeParaImpl(
      std::vector<ftm::MergeTree<dataType>> &trees,
      std::vector<
        std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>>
        &outputMatchings,
      std::vector<dataType> &distances,
      bool isFirstInput = true) {
      outputMatchings.resize(trees.size() - 1);
      distances.resize(trees.size() - 1);
      for(unsigned int ind = 0; ind < trees.size() - 1; ++ind) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp task firstprivate(ind) UNTIED() \
  shared(outputMatchings, trees, distances)
        {
#endif
          unsigned int i = ind;
          unsigned int j = ind + 1;
          // Execute
          if(baseModule_ == 0) {
            MergeTreeDistance mergeTreeDistance;
            mergeTreeDistance.setAssignmentSolver(assignmentSolverID_);
            mergeTreeDistance.setEpsilonTree1(epsilonTree1_);
            mergeTreeDistance.setEpsilonTree2(epsilonTree2_);
            mergeTreeDistance.setEpsilon2Tree1(epsilon2Tree1_);
            mergeTreeDistance.setEpsilon2Tree2(epsilon2Tree2_);
            mergeTreeDistance.setEpsilon3Tree1(epsilon3Tree1_);
            mergeTreeDistance.setEpsilon3Tree2(epsilon3Tree2_);
            mergeTreeDistance.setBranchDecomposition(branchDecomposition_);
            mergeTreeDistance.setParallelize(parallelize_);
            mergeTreeDistance.setPersistenceThreshold(persistenceThreshold_);
            mergeTreeDistance.setDebugLevel(std::min(debugLevel_, 2));
            mergeTreeDistance.setThreadNumber(this->threadNumber_);
            mergeTreeDistance.setNormalizedWasserstein(normalizedWasserstein_);
            mergeTreeDistance.setKeepSubtree(keepSubtree_);
            mergeTreeDistance.setDistanceSquaredRoot(distanceSquaredRoot_);
            mergeTreeDistance.setUseMinMaxPair(useMinMaxPair_);
            mergeTreeDistance.setMinMaxPairWeight(minMaxPairWeight_);
            mergeTreeDistance.setPreprocess(false);
            mergeTreeDistance.setSaveTree(false);
            mergeTreeDistance.setCleanTree(true);
            mergeTreeDistance.setIsCalled(true);
            mergeTreeDistance.setPostprocess(false);
            mergeTreeDistance.setIsPersistenceDiagram(isPersistenceDiagram_);
            if(useDoubleInput_) {
              double const weight = mixDistancesMinMaxPairWeight(isFirstInput);
              mergeTreeDistance.setMinMaxPairWeight(weight);
              mergeTreeDistance.setDistanceSquaredRoot(true);
            }
            distances[ind] = mergeTreeDistance.execute<dataType>(
              trees[i], trees[j], outputMatchings[ind]);
          }
#ifdef TTK_ENABLE_OPENMP
        } // end task
#endif
      } // end for i
    }

  }; // MergeTreeFeatureTracking class

} // namespace ttk
