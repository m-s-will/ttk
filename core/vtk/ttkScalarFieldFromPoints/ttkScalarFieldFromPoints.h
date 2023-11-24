/// \ingroup vtk
/// \class ttkScalarFieldFromPoints
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-10-11.
///
/// \brief TTK VTK-filter that wraps the ttk::ScalarFieldFromPoints module.
///
/// This VTK filter uses the ttk::ScalarFieldFromPoints module to compute the
/// scalar field from a set of points using a kernel to generate scalar values.
///
/// \param Input vtkPointSet.
/// \param Output vtkImageData.
///
/// \sa ttk::ScalarFieldFromPoints
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkScalarFieldFromPointsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <ScalarFieldFromPoints.h>

class TTKSCALARFIELDFROMPOINTS_EXPORT ttkScalarFieldFromPoints
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ScalarFieldFromPoints // and we inherit from the base
                                         // class
{
private:
  double ImageBounds[6]{0, 1, 0, 1, 0, 1};
  int Resolution[3]{2, 2, 2};
  int Kernel{0};

public:
  vtkSetVector6Macro(ImageBounds, double);
  vtkGetVector6Macro(ImageBounds, double);
  vtkSetVector3Macro(Resolution, int);
  vtkGetVector3Macro(Resolution, int);
  vtkSetMacro(Kernel, int);
  vtkGetMacro(Kernel, int);

  static ttkScalarFieldFromPoints *New();
  vtkTypeMacro(ttkScalarFieldFromPoints, ttkAlgorithm);

protected:
  ttkScalarFieldFromPoints();
  ~ttkScalarFieldFromPoints() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector);
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
