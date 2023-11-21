/// \ingroup base
/// \class ttk::SimilarityByGradient
/// \author Wito Engelke <wito.engelke@googlemail.com>
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date June 2020.
///
/// TODO
///
/// \b Related \b publication: \n
/// 'SimilarityByGradient'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <PathCompression.h>

#include <set>

namespace ttk {
  class SimilarityByGradient : virtual public Debug {

  public:
    SimilarityByGradient() {
      this->setDebugMsgPrefix("SimilarityByGradient");

      #ifdef TTK_ENABLE_MPI
        //if (ttk::MPIsize_ == 1)
          hasMPISupport_ = true;
      #endif

    };
    ~SimilarityByGradient(){};

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    template <typename IT, typename TT>
    int computeMap(std::vector<IT> &matchesVals, std::set<IT> &outsideMatches, TT * triangulation, const IT *segmentation, const IT *criticalPointGlobalIds0, const IT nFeatures0) const {

      for (IT i = 0; i < nFeatures0; i++) {
        // Get the local id of the extrema and add to matches if owner of extrema
        IT localId = triangulation->getVertexLocalId(criticalPointGlobalIds0[i]);
        if (triangulation->getVertexRank(localId) == ttk::MPIrank_) {
          IT match = segmentation[localId];
          matchesVals.push_back(match);
          IT localIdMatch = triangulation->getVertexLocalId(match);

          // If the matching critical point is not owned by the rank, add it to outside matches
          if (localIdMatch < 0){
            outsideMatches.insert(match);
          }
        }

      }

      return 1;
    };

    template <typename IT, typename TT>
    int gatherNeighbourhood(std::set<IT>& toTraverse, const TT * triangulation) const {
      std::set<IT> used;

      for (int s = 0; s < 1; s++) {
        std::vector<IT> toAdd(0);
        for (auto& vertex : toTraverse) {
          //this->printMsg(std::to_string(vertex));
          auto res = used.emplace(vertex);
          if (res.second) {
            // add neighbours to toTraverse
            const IT nNeighbors = triangulation->getVertexNeighborNumber(vertex);
            for (int n = 0; n < nNeighbors; n++) {
              IT u;
              triangulation->getVertexNeighbor(vertex, n, u);
              toAdd.push_back(u);
            }
          }
        }
        for (auto& vertex : toAdd) {
          toTraverse.insert(vertex);
        }
      }
      return 1;
    };

    template <typename IT, typename TT>
    int computeMapCombinatorial(std::vector<std::map<IT,int>> &matchMap, TT * triangulation, const IT *segmentation, const IT *criticalPointVertexIds0, const IT nFeatures0, const int neighborhoodSize) const {

      for (IT i = 0; i < nFeatures0; i++) {
        std::set<IT> toTraverse;
        std::set<IT> used;

        // add feature to toTraverse
        toTraverse.insert(criticalPointVertexIds0[i]);

        for (int s = 0; s <= neighborhoodSize; s++) {
          // save vertices to erase afterwards
          std::vector<IT> toErase(0);
          std::vector<IT> toAdd(0);
          for (auto& vertex : toTraverse) {
            //this->printMsg(std::to_string(vertex));
            auto res = used.emplace(vertex);
            if (res.second) {
              IT match = segmentation[vertex];
              matchMap[i][match] += 1;

              // add neighbours to toTraverse
              const IT nNeighbors = triangulation->getVertexNeighborNumber(vertex);
              for (int n = 0; n < nNeighbors; n++) {
                IT u;
                triangulation->getVertexNeighbor(vertex, n, u);
                toAdd.push_back(u);
              }
            }
            toErase.push_back(vertex);
          }

          for (auto& vertex : toErase) {
            toTraverse.erase(vertex);
          }
          for (auto& vertex : toAdd) {
            toTraverse.insert(vertex);
          }
        }
      }

      return 1;
    };

    template <typename IT, typename TT>
    int computeMapDistance(std::vector<std::map<IT,int>> &matchMap, TT * triangulation, const IT *segmentation, const IT *criticalPointVertexIds0, const IT nFeatures0, const int neighborhoodSize, const float neighbourhoodDistance) const {
      for (IT i = 0; i < nFeatures0; i++) {
        std::set<IT> toTraverse;
        std::set<IT> used;
        float pFeature[3] = {0.0, 0.0, 0.0};
        triangulation->getVertexPoint(criticalPointVertexIds0[i], pFeature[0], pFeature[1], pFeature[2]);
        // std::string msg = "Feature pos: (" + std::to_string(pFeature[0]) + ", " + std::to_string(pFeature[1]) + ", " + std::to_string(pFeature[2]) + ")\n";
        // add feature to toTraverse
        toTraverse.insert(criticalPointVertexIds0[i]);

        for (int s = 0; s <= neighborhoodSize; s++) {
          // save vertices to erase afterwards
          std::vector<IT> toErase(0);
          std::vector<IT> toAdd(0);
          for (auto& vertex : toTraverse) {
            //this->printMsg(std::to_string(vertex));
            auto res = used.emplace(vertex);
            if (res.second ) {
              IT match = segmentation[vertex];
              matchMap[i][match] += 1;

              if ( s != neighborhoodSize) {
                // add neighbours to toTraverse unless we are at the last iteration
                const IT nNeighbors = triangulation->getVertexNeighborNumber(vertex);
                for (int n = 0; n < nNeighbors; n++) {
                  IT u;
                  triangulation->getVertexNeighbor(vertex, n, u);

                  // Check if neighbor is within euclidean grid distance
                  //triangulation->vertexToPosition(u, p);
                  float p[3] = {0.0, 0.0, 0.0};
                  triangulation->getVertexPoint(u, p[0], p[1], p[2]);
                  //msg += "\t p: (" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) + ")\n";

                  float dist = ttk::Geometry::distance<float>(pFeature, p);

                  if (dist < neighbourhoodDistance)
                    toAdd.push_back(u);
                }
              }
            }
            toErase.push_back(vertex);
          }

          for (auto& vertex : toErase) {
            toTraverse.erase(vertex);
          }
          for (auto& vertex : toAdd) {
            toTraverse.insert(vertex);
          }
        }

          //this->printMsg(msg);
      }

      return 1;
    };

    template <typename IT>
    int computeManifoldSize(std::map<IT, int>& sizeMap, const IT *segmentation, const IT nVertices) const {

      for (IT i = 0; i < nVertices; i++) {
        IT manifold = segmentation[i];
        auto res = sizeMap.find(manifold);
        if (res != sizeMap.end()) {
          sizeMap[manifold] += 1;
        }
        else {
          sizeMap.insert({{manifold, 1}});
        }
      }
      return 1;
    }

    template <typename IT>
    int computeMapManifoldOverlap(std::map<IT,std::map<IT,int>> &matchMap, const IT *segmentation0, const IT *segmentation1, const IT nVertices) const {
      // go through all vertices and add to map
      for (IT i = 0; i < nVertices; i++) {
        IT manifold0 = segmentation0[i];
        IT manifold1 = segmentation1[i];
        // Check if map already contains manifold
        auto res0 = matchMap.find(manifold0);
        if (res0 != matchMap.end()) {
          // Check if manifolds own match map contains the matching manifold
          auto res1 = matchMap[manifold0].find(manifold1);
          if (res1 != matchMap[manifold0].end()) {
            matchMap[manifold0][manifold1] += 1;
          }
          else {
            // Insert matching manifold if not matched with before
            matchMap[manifold0].insert({{manifold1, 1}});
          }
        }
        else {
          // Insert new found manifold
          matchMap.insert({{manifold0, std::map<IT,int>()}});

          // Insert matching manifold if not matched with before
          matchMap[manifold0].insert({{manifold1, 1}});
        }
      }
      return 1;
    };

    template <typename IT, typename IF>
    int computeSimilarityMatrixMPI(int *matrix,
                                const std::vector<IT> &matchesVals,
                                const IT *criticalPointVertexIds1,
                                const IT nFeatures0,
                                const IT nFeatures1,
                                const IF indexFunction,
                                const std::string& msg) const {
      ttk::Timer timer;
      this->printMsg(msg, 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      for(int i = 0, j = nFeatures0 * nFeatures1; i < j; i++)
        matrix[i] = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nFeatures0; i++) {
        IT matchIdx = -1;
        if (i < IT(matchesVals.size()))
          matchIdx = matchesVals[i];

        // add matches to matix. Note: it is possible that a matched extrema is not in the feature list (i.e., not tracked)
        for(IT j = 0; j < nFeatures1; j++) {
          if(criticalPointVertexIds1[j] == matchIdx) {
            matrix[indexFunction(i, j, nFeatures0, nFeatures1)] = 1;
            break;
          }
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };


    template <typename IT, typename IF>
    int computeSimilarityMatrix(int *matrix,
                                const std::vector<std::map<IT,int>> &matchMap,
                                const IT nFeatures0,
                                const IT nFeatures1,
                                const std::unordered_map<IT,IT>& map1,
                                const IF indexFunction,
                                const std::string& msg) const {
      ttk::Timer timer;
      this->printMsg(msg, 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      for(int i = 0, j = nFeatures0 * nFeatures1; i < j; i++)
        matrix[i] = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nFeatures0; i++) {
        auto matches = matchMap[i];

        for (auto& m : matches) {
          auto res = map1.find(m.first);
          if (res != map1.end()) {
            matrix[indexFunction(i, res->second, nFeatures0, nFeatures1)] = m.second;
          }
        }

      }

      this->printMsg(msg, 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };

    template <typename IT, typename IF>
    int computeSimilarityManifoldMatrix(double *matrix,
                                const std::map<IT, std::map<IT,int>> &matchMap,
                                const std::map<IT, int> &sizeMap0,
                                const IT* criticalPointVertexIds0,
                                const IT nFeatures0,
                                const IT nFeatures1,
                                const std::unordered_map<IT,IT>& map1,
                                const IF indexFunction,
                                const std::string& msg) const {
      ttk::Timer timer;
      this->printMsg(msg, 0, 0, this->threadNumber_,
                     debug::LineMode::REPLACE);

      for(int i = 0, j = nFeatures0 * nFeatures1; i < j; i++)
        matrix[i] = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(IT i = 0; i < nFeatures0; i++) {
        IT vertex = criticalPointVertexIds0[i];
        auto innerMap = matchMap.at(vertex);
        int size0 = sizeMap0.at(vertex);

        for (auto& m : innerMap) {
          auto res = map1.find(m.first);
          if (res != map1.end()) {
            matrix[indexFunction(i, res->second, nFeatures0, nFeatures1)] = double(m.second)/ (size0);
          }
        }

      }

      this->printMsg(msg, 1, timer.getElapsedTime(),
                     this->threadNumber_);

      return 1;
    };



  }; // SimilarityByGradient class

} // namespace ttk