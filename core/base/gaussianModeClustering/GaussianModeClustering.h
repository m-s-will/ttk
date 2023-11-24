/// \ingroup base
/// \class ttk::GaussianModeClustering
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date The Date Here.
///
/// \b Related \b publication: \n
///

#pragma once

#include <Debug.h>

// std includes
#include <map>
#include <math.h>
#include <set>
#include <utility>
#include <vector>

namespace ttk {

  class GaussianModeClustering : virtual public Debug {

  public:
    GaussianModeClustering() {
      this->setDebugMsgPrefix("GaussianModeClustering");
    };
    ~GaussianModeClustering(){};

    // Class used to determine which vector fields exist and can be used
    enum class ClusteringType {
      UmbrellaClustering,
      ThresholdedUmbrellaClustering,
      Unimodality
    };

    void setClusteringType(const int ct) {
      switch(ct) {
        case 0: {
          clusteringType_ = ClusteringType::UmbrellaClustering;
          break;
        }
        case 1: {
          clusteringType_ = ClusteringType::ThresholdedUmbrellaClustering;
          break;
        }
        case 2: {
          clusteringType_ = ClusteringType::Unimodality;
          break;
        }
      }
    }

    void setThreshold(const double t) {
      threshold_ = t;
    }

    // ax^2 + bx + c = 0
    template <typename DT>
    int rootsQuadratic(const DT a, const DT b, const DT c, DT &t) const {
      DT discriminant = b * b - 4 * a * c;
      DT r1;
      DT r2;

      if(discriminant > 0) {
        // this->printMsg("Diff");
        r1 = (-1 * b + std::sqrt(discriminant)) / (2 * a);
        r2 = (-1 * b - std::sqrt(discriminant)) / (2 * a);
        if(r1 > 0 && r1 < 1)
          t = r1;
        else
          t = r2;

        return 1;
      } else if(discriminant == 0) {
        // this->printMsg("Equal");
        if(a != 0)
          r1 = r2 = (-1 * b) / (2 * a);
        else
          r1 = r2 = 0;

        t = r1;
        return 1;
      } else {
        t = -1;
        return 0;
      }
    }

    template <typename DT>
    bool umbrellaCondition(const DT iCoords[3],
                           const DT jCoords[3],
                           const double pws[2],
                           const double pcs[2],
                           double &maxVal) const {

      const DT dx = iCoords[0] - jCoords[0];
      const DT dy = iCoords[1] - jCoords[1];
      const DT dz = iCoords[2] - jCoords[2];

      const double sq = dx * dx + dy * dy + dz * dz;

      // Evaluate js value at i
      const double umbrellaVal = pws[1] * std::exp(-0.5 * (sq) / pcs[1]);

      if(umbrellaVal > maxVal) {
        maxVal = umbrellaVal;
        return true;
      } else
        return false;
    }

    template <typename DT>
    bool thresholdedUmbrellaCondition(const DT iCoords[3],
                                      const DT jCoords[3],
                                      const double pws[2],
                                      const double pcs[2],
                                      double &maxVal) const {

      const DT dx = jCoords[0] - iCoords[0];
      const DT dy = jCoords[1] - iCoords[1];
      const DT dz = jCoords[2] - iCoords[2];

      const double distance = dx * dx + dy * dy + dz * dz;

      // Thresholded gaussians with same variance in all dimensions end up as
      // circles/spheres
      const double ri = std::log(threshold_ / pws[0]) * (-2.0 * pcs[0]);
      const double rj = std::log(threshold_ / pws[1]) * (-2.0 * pcs[1]);
      if(ri < 0 || rj < 0)
        return false;

      // no overlap
      if(std::sqrt(ri) + std::sqrt(rj) < std::sqrt(distance))
        return false;

      // they are connected
      if(pws[1] > maxVal) {
        maxVal = pws[1];
        return true;
      } else
        return false;
    }

    template <typename DT>
    bool unimodalityCondition(const DT iCoords[3],
                              const DT jCoords[3],
                              const double pws[2],
                              const double pcs[2],
                              double &maxVal) const {

      const DT dx = iCoords[0] - jCoords[0];
      const DT dy = iCoords[1] - jCoords[1];
      const DT dz = iCoords[2] - jCoords[2];

      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

      // if i == j
      if(distance == 0 && pws[0] > maxVal)
        return true;

      // Check if points are close enough
      if(std::abs(dx) < 3 * std::sqrt(pcs[1])
         && std::abs(dy) < 3 * std::sqrt(pcs[1])
         && std::abs(dz) < 3 * std::sqrt(pcs[1])) {

        if(distance <= 2 * std::min(std::sqrt(pcs[0]), std::sqrt(pcs[1]))
           && pws[1] > pws[0]) {
          maxVal = pws[1];
          return true;
        }
      }

      return false;
    }

    template <typename DT>
    bool intersectionCondition(const DT iCoords[3],
                               const DT jCoords[3],
                               const double pws[2],
                               const double pcs[2],
                               double &maxVal) const {

      const DT dx = iCoords[0] - jCoords[0];
      const DT dy = iCoords[1] - jCoords[1];
      const DT dz = iCoords[2] - jCoords[2];

      const double sq = dx * dx + dy * dy + dz * dz;

      double a = ((-0.5 / pcs[0]) - (-0.5 / pcs[1])) * sq;
      double b = 2 * (-0.5 / pcs[1]) * sq;
      double c
        = (-1 * (-0.5 / pcs[1]) * sq) + std::log(pws[0]) - std::log(pws[1]);
      double t;
      int status = 0;
      status = this->rootsQuadratic<double>(a, b, c, t);
      double ti = 1 - t; // for i
      double tj = t; // for j

      const double sqTi
        = (ti * dx) * (ti * dx) + (ti * dy) * (ti * dy) + (ti * dz) * (ti * dz);
      const double sqTj
        = (tj * dx) * (tj * dx) + (tj * dy) * (tj * dy) + (tj * dz) * (tj * dz);

      // Evaluate the summed value at intersection between i and j
      const double intersectionVal
        = pws[0] * std::exp(-0.5 * (sqTj) / pcs[0])
          + pws[1] * std::exp(-0.5 * (sqTi) / pcs[1]);

      if(intersectionVal > pws[1] && pws[1] > maxVal) {
        maxVal = pws[1];
        return true;
      }

      return false;
    }

    template <typename DT>
    int computeClusters(std::map<int, std::vector<int>> &pointClusters,
                        int &nClusters,
                        const DT *coords,
                        const double *pws,
                        const double *pcs,
                        const int nPoints) const {

      ttk::Timer timer;

      const std::string msg = "Computing Clusters";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      // umbrella indices list
      std::vector<int> inCluster(nPoints);
      std::vector<double> clusterMaxVals(nPoints);
      for(int i = 0; i < nPoints; i++) {
        inCluster[i] = -1;
        clusterMaxVals[i] = 0.0;
      }

      for(int i = 0; i < nPoints; i++) {

        for(int j = 0; j < nPoints; j++) {
          const int i3 = i * 3;
          const int j3 = j * 3;

          const DT iCoords[3]
            = {coords[i3 + 0], coords[i3 + 1], coords[i3 + 2]};
          const DT jCoords[3]
            = {coords[j3 + 0], coords[j3 + 1], coords[j3 + 2]};

          const double pw[2] = {pws[i], pws[j]};
          const double pc[2] = {pcs[i], pcs[j]};

          // Evaluate condition to see if i belongs to cluster of j
          switch(clusteringType_) {
            case ClusteringType::UmbrellaClustering: {
              if(umbrellaCondition<DT>(
                   iCoords, jCoords, pw, pc, clusterMaxVals[i]))
                inCluster[i] = j;

              break;
            }
            case ClusteringType::ThresholdedUmbrellaClustering: {
              if(thresholdedUmbrellaCondition<DT>(
                   iCoords, jCoords, pw, pc, clusterMaxVals[i])) {
                if(inCluster[i] > -1 && inCluster[i] != i && i != j) {

                  inCluster[inCluster[i]] = j;
                  clusterMaxVals[inCluster[i]] = clusterMaxVals[i];
                }
                inCluster[i] = j;
              }

              break;
            }
            case ClusteringType::Unimodality: {
              if(unimodalityCondition<DT>(
                   iCoords, jCoords, pw, pc, clusterMaxVals[i]))
                inCluster[i] = j;

              break;
            }
          }
        }
      }

      // initialize umbrellas with point representatives
      for(int i = 0; i < nPoints; i++) {
        if(inCluster[i] == i) {
          pointClusters[i] = std::vector<int>(0);
          pointClusters[i].push_back(i);
          nClusters++;
        } else if(inCluster[i] == -1) {
          pointClusters[i] = std::vector<int>(0);
          pointClusters[i].push_back(-1);
        }
      }

      // loop through all points
      for(int i = 0; i < nPoints; i++) {
        // if a point is not within its own umbrella
        if(inCluster[i] != i && inCluster[i] != -1) {
          int j = i;

          // loop until you find the root representative
          while(j != inCluster[j]) {
            j = inCluster[j];
          }

          pointClusters[j].push_back(i);
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

    template <typename DT>
    int computeMatrix(unsigned char *matrix,
                      std::map<int, std::vector<int>> &clusters0,
                      std::map<int, std::vector<int>> &clusters1,
                      const DT *ids0,
                      const DT *ids1,
                      const int nClusters0,
                      const int nClusters1) const {

      ttk::Timer timer;

      const std::string msg = "Computing Matrix (" + std::to_string(nClusters0)
                              + "x" + std::to_string(nClusters1) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      int i = 0;
      for(const auto &c0 : clusters0) {
        if(c0.second[0] == -1) {
          continue;
        }
        int j = 0;
        for(const auto &c1 : clusters1) {
          if(c1.second[0] == -1)
            continue;
          int exist = 0;
          for(long unsigned int p0 = 0; p0 < c0.second.size(); p0++) {
            for(long unsigned int p1 = 0; p1 < c1.second.size(); p1++) {
              if(ids0[c0.second[p0]] == ids1[c1.second[p1]]) {
                exist = 1;
              }
            }
          }
          matrix[j * nClusters0 + i] = exist;
          ++j;
        }
        ++i;
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

  protected:
    ClusteringType clusteringType_{};
    double threshold_{};
  };
} // namespace ttk
