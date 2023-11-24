/// \ingroup base
/// \class ttk::ComparingSimilarityMatrices
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2022-03-07.
///
/// \b Related \b publication: \n

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <unordered_map>

namespace ttk {

  class ComparingSimilarityMatrices : virtual public Debug {

  public:
    ComparingSimilarityMatrices();

    struct Events {
      // {alg, gt, correct}
      int continuations[3]{0, 0, 0};
      int deaths[3]{0, 0, 0};
      int splits[3]{0, 0, 0};
      int births[3]{0, 0, 0};
      int merges[3]{0, 0, 0};
    };

    template <typename DT0, typename DT1>
    int compareEvents(Events &curEvents,
                      const DT0 *algMatrix,
                      const DT0 *gtMatrix,
                      const DT1 *indexIdMapAlg0,
                      const DT1 *indexIdMapAlg1,
                      const DT1 *indexIdMapGT0,
                      const DT1 *indexIdMapGT1,
                      const long long int dims[2]) {

      std::map<int, std::vector<int>> rowFillAlg;
      std::map<int, std::vector<int>> columnFillAlg;
      std::map<int, std::vector<int>> rowFillGT;
      std::map<int, std::vector<int>> columnFillGT;

      for(long long int r = 0; r < dims[1]; r++) {
        rowFillAlg[indexIdMapAlg1[r]] = std::vector<int>(0);
        rowFillGT[indexIdMapGT1[r]] = std::vector<int>(0);
      }
      for(long long int c = 0; c < dims[0]; c++) {
        columnFillAlg[indexIdMapAlg0[c]] = std::vector<int>(0);
        columnFillGT[indexIdMapGT0[c]] = std::vector<int>(0);
      }

      //
      for(long long int r = 0; r < dims[1]; r++) {
        for(long long int c = 0; c < dims[0]; c++) {

          int valAlg = algMatrix[r * dims[0] + c];
          int valGT = gtMatrix[r * dims[0] + c];
          if(valAlg == 1) {
            rowFillAlg[indexIdMapAlg1[r]].push_back(indexIdMapAlg0[c]);
            columnFillAlg[indexIdMapAlg0[c]].push_back(indexIdMapAlg1[r]);
          }
          if(valGT == 1) {
            rowFillGT[indexIdMapGT1[r]].push_back(indexIdMapGT0[c]);
            columnFillGT[indexIdMapGT0[c]].push_back(indexIdMapGT1[r]);
          }
        }
      }

      // Count and compare events in columns: deaths, continuations, splits
      for(long long int c = 0; c < dims[0]; c++) {
        int nOnesAlg = columnFillAlg[indexIdMapAlg0[c]].size();
        int nOnesGT = columnFillGT[indexIdMapGT0[c]].size();

        // Count events possible in columns
        if(nOnesAlg == 0)
          ++curEvents.deaths[0];
        else if(nOnesAlg == 1)
          ++curEvents.continuations[0];
        else if(nOnesAlg >= 2) {
          ++curEvents.splits[0];

          // Even if the feature splits into two or more, it still continues
          ++curEvents.continuations[0];
        }

        if(nOnesGT == 0)
          ++curEvents.deaths[1];
        else if(nOnesGT == 1)
          ++curEvents.continuations[1];
        else if(nOnesGT >= 2) {
          ++curEvents.splits[1];
          // Even if the feature splits into two or more, it still continues
          ++curEvents.continuations[1];
        }
        // Check correctness
        if(nOnesAlg == nOnesGT && nOnesAlg == 0) {
          ++curEvents.deaths[2];
        }

        if(nOnesAlg == nOnesGT && nOnesAlg == 1
           && columnFillAlg[indexIdMapAlg0[c]][0]
                == columnFillGT[indexIdMapGT0[c]][0])
          ++curEvents.continuations[2];

        // To check correctness of the split: check that the feature ids
        // involved in the split matches
        if(nOnesAlg == nOnesGT && nOnesAlg >= 2) {
          bool completeMatch = true;
          bool atLeastOneMatch = false;
          for(int i = 0; i < nOnesAlg; i++) {
            if(columnFillAlg[indexIdMapAlg0[c]][i]
               != columnFillGT[indexIdMapGT0[c]][i])
              completeMatch = false;
            else
              atLeastOneMatch = true;
          }
          if(completeMatch)
            ++curEvents.splits[2];
          if(atLeastOneMatch)
            ++curEvents.continuations[2];
        }

        // Check if a continuation occurs within the split, even though it
        // initself is not correct
        if(nOnesAlg >= 2 && nOnesGT < nOnesAlg) {
          bool atLeastOneMatch = false;
          for(int i = 0; i < nOnesAlg; i++) {
            for(int j = 0; j < nOnesGT; j++) {
              if(columnFillAlg[indexIdMapAlg0[c]][i]
                 == columnFillGT[indexIdMapGT0[c]][j]) {
                atLeastOneMatch = true;
              }
            }
          }
          if(atLeastOneMatch)
            ++curEvents.continuations[2];
        }

        if(nOnesGT >= 2 && nOnesAlg < nOnesGT) {
          bool atLeastOneMatch = false;
          for(int i = 0; i < nOnesGT; i++) {
            for(int j = 0; j < nOnesAlg; j++) {
              if(columnFillAlg[indexIdMapAlg0[c]][j]
                 == columnFillGT[indexIdMapGT0[c]][i]) {
                atLeastOneMatch = true;
              }
            }
          }
          if(atLeastOneMatch)
            ++curEvents.continuations[2];
        }
      }

      // In rows we can check for births and merges
      // {alg, gt, correct}
      for(long long int r = 0; r < dims[1]; r++) {
        int nOnesAlg = rowFillAlg[indexIdMapAlg1[r]].size();
        int nOnesGT = rowFillGT[indexIdMapGT1[r]].size();

        // Count events possible in columns
        if(nOnesAlg == 0)
          ++curEvents.births[0];
        else if(nOnesAlg >= 2) {
          ++curEvents.merges[0];
        }

        if(nOnesGT == 0)
          ++curEvents.births[1];
        else if(nOnesGT >= 2) {
          ++curEvents.merges[1];
        }

        // Check correctness
        if(nOnesAlg == nOnesGT && nOnesAlg == 0) {
          ++curEvents.births[2];
        }

        // To check correctness, check that the feature ids involved in the
        // merge matches
        if(nOnesAlg == nOnesGT && nOnesAlg >= 2) {
          bool match = true;
          for(int i = 0; i < nOnesAlg; i++) {
            if(rowFillAlg[indexIdMapAlg1[r]][i]
               != rowFillGT[indexIdMapGT1[r]][i])
              match = false;
          }
          if(match) {
            ++curEvents.merges[2];
          }
        }
      }

      return 1;
    }

  }; // ComparingSimilarityMatrices class

} // namespace ttk
