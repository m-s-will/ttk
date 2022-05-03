/// \ingroup base
/// \class ttk::BaseMPIClass
/// \author Michael Will <mswill@rhrk.uni-kl.de>
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date April 2022
///
/// \brief Base Class and utilities for MPI implementation.

#pragma once
#include <BaseClass.h>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>

#if TTK_ENABLE_MPI
#include <mpi.h>

namespace ttk {
  COMMON_EXPORTS extern int MPIrank_;

  inline MPI_Datatype getMPIType(const float ttkNotUsed(val)) {
    return MPI_FLOAT;
  };
  inline MPI_Datatype getMPIType(const int ttkNotUsed(val)) {
    return MPI_INT;
  };
  inline MPI_Datatype getMPIType(const unsigned int ttkNotUsed(val)) {
    return MPI_UNSIGNED;
  };
  inline MPI_Datatype getMPIType(const double ttkNotUsed(val)) {
    return MPI_DOUBLE;
  };
  inline MPI_Datatype getMPIType(const long double ttkNotUsed(val)) {
    return MPI_LONG_DOUBLE;
  };
  inline MPI_Datatype getMPIType(const long ttkNotUsed(val)) {
    return MPI_LONG;
  };
  inline MPI_Datatype getMPIType(const unsigned long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG;
  };
  inline MPI_Datatype getMPIType(const long long ttkNotUsed(val)) {
    return MPI_LONG_LONG;
  };
  inline MPI_Datatype getMPIType(const unsigned long long ttkNotUsed(val)) {
    return MPI_UNSIGNED_LONG_LONG;
  };

  inline bool isRunningWithMPI() {
    int flag_i;
    MPI_Initialized(&flag_i);
    return flag_i;
  }

  template <typename IT>
  void getGhostCellInfo(std::vector<std::queue<IT>> &rankQueues,
                        const int *const rankArray,
                        const int rank,
                        const IT nVerts) {
    for(IT i = 0; i < nVerts; i++) {
      if(rank != rankArray[i]) {
        rankQueues[rankArray[i]].push(i);
      }
    }
  }

  template <typename DT, typename IT>
  void sendGhostCellInfo(DT *scalarArray,
                         const int *const rankArray,
                         const IT *const globalIds,
                         const std::unordered_map<IT, IT> gidToLidMap,
                         const int rankToSend,
                         const int thisRank,
                         const int nRanks,
                         const IT nVerts) {
    MPI_Datatype MPI_DT = getMPIType(static_cast<DT>(0));
    MPI_Datatype MPI_IT = getMPIType(static_cast<IT>(0));
    int amountTag = 101;
    int idsTag = 102;
    int valuesTag = 103;
    if(rankToSend == thisRank) {
      std::vector<std::vector<IT>> rankVectors;
      rankVectors.resize(nRanks);
      // aggregate the needed ids
      for(IT i = 0; i < nVerts; i++) {
        if(rankToSend != rankArray[i]) {
          rankVectors[rankArray[i]].push_back(globalIds[i]);
        }
      }

      // send the amount of ids and the needed ids themselves
      for(int r = 0; r < nRanks; r++) {
        if(rankToSend != r) {
          IT nValues = rankVectors[r].size();
          MPI_Send(&nValues, 1, MPI_IT, r, amountTag, MPI_COMM_WORLD);
          if(nValues > 0) {
            MPI_Send(rankVectors[r].data(), nValues, MPI_IT, r, idsTag,
                     MPI_COMM_WORLD);
          }
        }
      }

      // receive the scalar values
      for(int r = 0; r < nRanks; r++) {
        if(rankToSend != r) {
          IT nValues = rankVectors[r].size();
          std::vector<DT> receivedValues(nValues);
          if(nValues > 0) {
            MPI_Recv(receivedValues.data(), nValues, MPI_DT, r, valuesTag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(IT i = 0; i < nValues; i++) {
              DT receivedVal = receivedValues[i];
              IT globalId = rankVectors[r][i];
              IT localId = gidToLidMap.at(globalId);
              scalarArray[localId] = receivedVal;
            }
          }
        }
      }

    } else {
      // receive the amount of ids and the needed ids themselves
      IT nValues;
      MPI_Recv(&nValues, 1, MPI_IT, rankToSend, amountTag, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      if(nValues > 0) {
        std::vector<IT> receivedIds(nValues);
        MPI_Recv(receivedIds.data(), nValues, MPI_IT, rankToSend, idsTag, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        // assemble the scalar values
        std::vector<DT> valuesToSend(nValues);
        for(IT i = 0; i < nValues; i++) {
          IT globalId = receivedIds[i];
          IT localId = gidToLidMap.at(globalId);
          valuesToSend[i] = scalarArray[localId];
        }
        // send the scalar values
        MPI_Send(valuesToSend.data(), nValues, MPI_DT, rankToSend, valuesTag,
                 MPI_COMM_WORLD);
      }
    }
  }

  class BaseMPIClass : public BaseClass {

  public:
    BaseMPIClass();
    virtual ~BaseMPIClass() = default;
  };
} // namespace ttk

#endif