/// \ingroup base
/// \class ttk::MorseSegmentation
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \date 2022.
///
/// This module defines the %MorseSegmentation class that compute the morse smale segmentation given an ascending and descending manifold
///
/// \b Related \b publication: \n
/// 'MorseSegmentation'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <unordered_map>

namespace ttk {

  /**
   * The MorseSegmentation class provides methods to compute the morse smale segmentation given an ascending and descending manifold
   */
  class MorseSegmentation : virtual public Debug {

  public:
    MorseSegmentation();


    std::vector<int> compressArray(const std::vector<int>& input) const {
      ttk::Timer compressTimer;
    
      std::vector<int> output(input.size());
      std::unordered_map<int, int> uniquesMap;
      int counter = 0;
      // assemble the output by creating a map of unique values while running over the array
      for (size_t i = 0; i < input.size(); i++){
        if ( uniquesMap.find(input[i]) == uniquesMap.end() ){
          uniquesMap[input[i]] = counter;
          counter++;
        }
        output[i] = uniquesMap[input[i]];
      }
      this->printMsg("#Unique Segmentations: " + std::to_string(counter), 1, compressTimer.getElapsedTime()); 
      return output;
    }

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }


    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeAverages(dataType *outputData,
                        const dataType *ascendingManifold,
                        const dataType *descendingManifold,
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
      // Compute Vertex Averages
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Segmentation",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

        // compute the average of each vertex in parallel
        int nVertices = triangulation->getNumberOfVertices();
        std::vector<int> segmentation(nVertices);    
        int numberOfMaxima = 0;
        // get the number of unique values of descendingManifold by just getting the maximum value
        for (int i = 0; i < nVertices; i++) {
          if (descendingManifold[i] > numberOfMaxima){
            numberOfMaxima = descendingManifold[i];
          }
        }

        // generate unique values for unique descending / ascending combinations by using the number of distinct values in descendingManifold
        // and multiplying it by the current ascending value. This ensures that adding the descending value doesn't overflow into a new "bucket"
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(int i = 0; i < nVertices; i++) {
          segmentation[i] = descendingManifold[i] + ascendingManifold[i] * numberOfMaxima;
        }

        segmentation = compressArray(segmentation);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for (int i = 0; i < nVertices; i++){
            outputData[i] = segmentation[i];
        }

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Segmentation",
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

  }; // MorseSegmentation class

} // namespace ttk
