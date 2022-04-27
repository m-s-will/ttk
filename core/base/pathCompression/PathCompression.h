/// \ingroup base
/// \class ttk::PathCompression
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %PathCompression class that computes for each vertex
/// of a triangulation the vertices to which ascending and descending manifolds
/// it belongs.
///
/// \b Related \b publication: \n
/// 'PathCompression'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <map>
#include <unordered_map>

namespace ttk {

  /**
   * The PathCompression class provides methods to compute for each vertex of a
   * triangulation the vertices to which ascending and descending manifolds it
   * belongs.
   */
  class PathCompression : virtual public Debug {

  public:
    PathCompression();

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    std::vector<int> compressArray(const std::vector<int> &input) const {
      ttk::Timer compressTimer;

      std::vector<int> output(input.size());
      std::unordered_map<int, int> uniquesMap;
      int counter = 0;
      // assemble the output by creating a map of unique values while running
      // over the array
      for(size_t i = 0; i < input.size(); i++) {
        if(uniquesMap.find(input[i]) == uniquesMap.end()) {
          uniquesMap[input[i]] = counter;
          counter++;
        }
        output[i] = uniquesMap[input[i]];
      }
      this->printMsg(
        "#Unique Segmentations with new method: " + std::to_string(counter), 1,
        compressTimer.getElapsedTime());
      return output;
    }

    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeCompression(int *descendingManifold,
                           int *ascendingManifold,
                           const dataType *inputData,
                           const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Minima
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing compression",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_);

        int nVertices = triangulation->getNumberOfVertices();

        std::vector<int> previousDesc(nVertices);
        std::vector<int> currentDesc(nVertices);
        std::vector<int> previousAsc(nVertices);
        std::vector<int> currentAsc(nVertices);

        // for the first step we initialize each vertex with the id of their
        // largest / smallest neighbor. Afterwards we only compare the arrays
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(int i = 0; i < nVertices; i++) {
          int nNeighbors = triangulation->getVertexNeighborNumber(i);
          ttk::SimplexId neighborId;
          float smallest = inputData[i];
          float largest = inputData[i];
          // if there is no larger / smaller neighbor, the vertex points to
          // itself and is therefore a maximum / minimum we do not need to check
          // for equality, because we use the order array
          previousDesc[i] = i;
          previousAsc[i] = i;
          for(int j = 0; j < nNeighbors; j++) {
            triangulation->getVertexNeighbor(i, j, neighborId);

            // and for the largest neighbor to get to the ascending manifold
            if(inputData[neighborId] > largest) {
              previousAsc[i] = neighborId;
              largest = inputData[neighborId];
            }
            // we're checking for the smallest neighbor to get the descending
            // manifold
            if(inputData[neighborId] <= smallest) {
              previousDesc[i] = neighborId;
              smallest = inputData[neighborId];
            }
          }
        }
        // now we swap between the two arrays until nothing changes anymore ergo
        // all paths are finished
        int step = 0;
        bool same = false;
        int nextDesc = 0;
        int nextAsc = 0;
        while(!same) {
          this->printMsg("Running Step " + std::to_string(step));
          same = true;
          if(step % 2 == 0) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for(int i = 0; i < nVertices; i++) {
              // this->printMsg("Previous zu current");
              nextDesc = previousDesc[previousDesc[i]];
              nextAsc = previousAsc[previousAsc[i]];
              if(nextDesc != currentDesc[i]) {
                currentDesc[i] = nextDesc;
                same = false;
              }
              if(nextAsc != currentAsc[i]) {
                currentAsc[i] = nextAsc;
                same = false;
              }
            }
          } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
            for(int i = 0; i < nVertices; i++) {
              // this->printMsg("Current zu Previous");
              nextDesc = currentDesc[currentDesc[i]];
              nextAsc = currentAsc[currentAsc[i]];
              if(nextDesc != previousDesc[i]) {
                previousDesc[i] = nextDesc;
                same = false;
              }
              if(nextAsc != previousAsc[i]) {
                previousAsc[i] = nextAsc;
                same = false;
              }
            }
          }

          step++;
        }

        this->printMsg("Finished in Step " + std::to_string(step), 1,
                       localTimer.getElapsedTime());

        // compress the arrays into the ranges of 0 - #segmentation areas
        // currentDesc = this->compressArray(currentDesc);
        // currentAsc = this->compressArray(currentAsc);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(int i = 0; i < nVertices; i++) {
          descendingManifold[i] = currentDesc[i];
          ascendingManifold[i] = currentAsc[i];
        }

        this->printMsg("Erster Wert in fertigem Ascending "
                       + std::to_string(ascendingManifold[0]));
        this->printMsg("Erster Wert in fertigem Descending "
                       + std::to_string(descendingManifold[0]));
        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Compression",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

  }; // PathCompression class

} // namespace ttk
