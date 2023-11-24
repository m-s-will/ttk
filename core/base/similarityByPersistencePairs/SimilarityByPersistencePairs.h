///
/// \ingroup base
/// \class ttk::SimilarityByPersistencePairs
/// \author Maxime Soler
/// \date 09.06.2021
///
/// \brief Computes a similarity matrix from persistence diagram matchings.
///
/// This module defines the SimilarityByPersistencePairs class that computes
/// the correspondance matrix from the output of a Wasserstein-based matching
/// between persistence diagrams. 1 row = 1 persistence pair in diagram 1,
/// 1 column = 1 persistence pair in diagram 2,
/// value = 0 => no matching
/// value > 0 => matching.
///
/// Maxime Soler, Jonas Lukasczyk, Julien Tierny, 2020.
///

#pragma once

#include <BottleneckDistance.h>
#include <Debug.h>

namespace ttk {

  class SimilarityByPersistencePairs : virtual public Debug {

  public:
    SimilarityByPersistencePairs() {
      this->setDebugMsgPrefix("SimilarityByPersistencePairs");
    };
    ~SimilarityByPersistencePairs(){};

    template <class DT>
    int computeDistanceMatrix(
      std::vector<std::tuple<int,
                             ttk::CriticalType,
                             int,
                             ttk::CriticalType,
                             DT,
                             int,
                             DT,
                             float,
                             float,
                             float,
                             DT,
                             float,
                             float,
                             float>> &CTDiagram0,
      std::vector<std::tuple<int,
                             ttk::CriticalType,
                             int,
                             ttk::CriticalType,
                             DT,
                             int,
                             DT,
                             float,
                             float,
                             float,
                             DT,
                             float,
                             float,
                             float>> &CTDiagram1,
      std::vector<std::tuple<int, int, double>> &matchings,
      double px,
      double py,
      double pz,
      double ps,
      double pe,
      const std::string algorithm,
      const std::string wasserstein,
      const int pvAlgorithm,
      const double maxJump) const {
      ttk::Timer timer;

      // const std::string msg = "Computing Distance Matrix ("
      //                        + std::to_string(nPoints0) + "x"
      //                        + std::to_string(nPoints1) + ")";
      // this->printMsg(
      //  msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      ttk::BottleneckDistance bottleneckDistance_;
      bottleneckDistance_.setPersistencePercentThreshold(0);
      // (tolerance can be set earlier, when features are defined)
      bottleneckDistance_.setPX(px);
      bottleneckDistance_.setPY(py);
      bottleneckDistance_.setPZ(pz);
      bottleneckDistance_.setPS(ps);
      bottleneckDistance_.setPE(pe);
      // bottleneckDistance_.setPercentMaxJump(maxJump);
      bottleneckDistance_.setAlgorithm(algorithm);
      bottleneckDistance_.setPVAlgorithm(pvAlgorithm);
      bottleneckDistance_.setWasserstein(wasserstein);

      // int status = bottleneckDistance_.computeBottleneck(
      //   CTDiagram0,
      //   CTDiagram1,
      //   matchings
      // );

      // bottleneckDistance_.setCTDiagram1(&CTDiagram0);
      // bottleneckDistance_.setCTDiagram2(&CTDiagram1);
      // bottleneckDistance_.setOutputMatchings(&matchings);
      // int status = bottleneckDistance_.execute<double>(false);
      // if(status < 0)
      //   return -1;

      return 1;
    }
  };
} // namespace ttk
