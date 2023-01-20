#pragma once

// VTK Module
#include <ttkScalarFieldCriticalPoints2Module.h>

// VTK Includes
#include <ttkAlgorithm.h>

#include <ScalarFieldCriticalPoints2.h>

class TTKSCALARFIELDCRITICALPOINTS2_EXPORT ttkScalarFieldCriticalPoints2
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ScalarFieldCriticalPoints2 // and we inherit from the base class
{
private:
  std::string OutputArrayName{"AveragedScalarField"};

public:
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  static ttkScalarFieldCriticalPoints2 *New();
  vtkTypeMacro(ttkScalarFieldCriticalPoints2, ttkAlgorithm);

protected:
  ttkScalarFieldCriticalPoints2();
  ~ttkScalarFieldCriticalPoints2() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
