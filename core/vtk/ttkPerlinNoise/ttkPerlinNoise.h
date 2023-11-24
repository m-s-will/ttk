/// \ingroup vtk
/// \class ttkPerlinNoise
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-06-04
///
/// This filter creates a Perlin noise scalar field of chosen spatio-temporal
/// dimension.
///
/// \sa ttk::PerlinNoise
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkPerlinNoiseModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkImageData.h>

// TTK Base Includes
#include <PerlinNoise.h>

class TTKPERLINNOISE_EXPORT ttkPerlinNoise : public ttkAlgorithm,
                                             protected ttk::PerlinNoise {
private:
  int Resolution[3]{0, 0, 0};
  double Scale{0};
  int Frequency{0};
  int nOctaves{0};
  double Persistence{0.0};

  int TimeProp{0};
  double TimeStep{0.0};
  int TimeSeries{0};
  double Interval{0};

public:
  vtkSetVector3Macro(Resolution, int);
  vtkGetVector3Macro(Resolution, int);

  vtkSetMacro(Scale, double);
  vtkGetMacro(Scale, double);

  vtkSetMacro(Frequency, int);
  vtkGetMacro(Frequency, int);

  vtkSetMacro(nOctaves, int);
  vtkGetMacro(nOctaves, int);

  vtkSetMacro(Persistence, double);
  vtkGetMacro(Persistence, double);

  vtkSetMacro(TimeProp, int);
  vtkGetMacro(TimeProp, int);

  vtkSetMacro(TimeStep, double);
  vtkGetMacro(TimeStep, double);

  vtkSetMacro(TimeSeries, int);
  vtkGetMacro(TimeSeries, int);

  vtkSetMacro(Interval, double);
  vtkGetMacro(Interval, double);

  static ttkPerlinNoise *New();
  vtkTypeMacro(ttkPerlinNoise, ttkAlgorithm);

protected:
  ttkPerlinNoise();
  ~ttkPerlinNoise() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;

  int initializeOutput(vtkImageData *img,
                       int extent[6],
                       const int nTuples,
                       const double t);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
