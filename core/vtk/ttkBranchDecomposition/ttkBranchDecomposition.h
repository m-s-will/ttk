/// \ingroup vtk
/// \class ttkBranchDecomposition
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 29.07.2021
///
/// \brief This module assigns to each vertex of a tracking graph a branch id
/// based on a given attribute.
///
/// This module assigns to each vertex of a tracking graph a branch id based on
/// a given attribute. First, all birth nodes are assigned a unique branch id,
/// and then the algorithm iterates over every vertex in order of time and then
/// either inherits the branch id of its largest predecessor (but only if the
/// current vertex is also the largest successor of this predecessor), or the
/// vertex gets assinged a new unique branch id.
///
/// \param Input vtkPolyData or vtkUnstructuredGrid representing a tracking
/// graph. \param Output A shallow copy of the input data object augmented with
/// a BranchId point data array.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// This module requires two scalars fields (time and attribute) that need to be
/// specified via the standard VTK call vtkAlgorithm::SetInputArrayToProcess().
///
/// The time array needs the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires (time))
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the time array)
///
/// The attribute array needs the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires (attribute))
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 or 1 (DYNAMIC: either point or cell data)
/// \param arrayName (DYNAMIC: string identifier of the attribute array)
///
/// \sa ttk::BranchDecomposition
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkBranchDecompositionModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <BranchDecomposition.h>

class TTKBRANCHDECOMPOSITION_EXPORT ttkBranchDecomposition
  : public ttkAlgorithm,
    protected ttk::BranchDecomposition {
public:
  static ttkBranchDecomposition *New();
  vtkTypeMacro(ttkBranchDecomposition, ttkAlgorithm);

protected:
  ttkBranchDecomposition();
  ~ttkBranchDecomposition() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};