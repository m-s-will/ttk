#include <ttkPerlinNoise.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPerlinNoise);

ttkPerlinNoise::ttkPerlinNoise() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkPerlinNoise::~ttkPerlinNoise() {
}

int ttkPerlinNoise::FillInputPortInformation(int, vtkInformation *) {
  return 0;
}

int ttkPerlinNoise::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    if(this->TimeProp < 2) {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    } else {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    }
  } else
    return 0;
  return 1;
}

int ttkPerlinNoise::RequestInformation(vtkInformation *,
                                       vtkInformationVector **,
                                       vtkInformationVector *outputVector) {
  // For vtkImageData output we have to set extent, spacing and origin already
  // in request information
  if(this->TimeProp < 2) {
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    int dimZ = 0;

    if(this->Resolution[2] == 0) {
      dimZ = 0;
    } else if(this->Resolution[2] > 0) {
      dimZ = this->Resolution[2] - 1;
    } else {
      this->printErr("Perlin noise dimension is off.");
      return 0;
    }

    int extent[6]
      = {0, this->Resolution[0] - 1, 0, this->Resolution[1] - 1, 0, dimZ};
    double spacing[3] = {1.0, 1.0, 1.0};
    double origin[3] = {0.0, 0.0, 0.0};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
    outInfo->Set(vtkDataObject::SPACING(), spacing, 3);
    outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
  }

  return 1;
}

int ttkPerlinNoise::initializeOutput(vtkImageData *img,
                                     int extent[6],
                                     const int nTuples,
                                     const double t) {

  // Set dimensions and extent for image data
  img->SetOrigin(0, 0, 0);
  img->SetExtent(extent);
  const double spacing[3]{
    this->Resolution[0] > 1 ? 1.0 / (this->Resolution[0] - 1) : 0,
    this->Resolution[1] > 1 ? 1.0 / (this->Resolution[1] - 1) : 0,
    this->Resolution[2] > 1 ? 1.0 / (this->Resolution[2] - 1) : 0};
  img->SetSpacing(spacing);

  auto noiseArray = vtkSmartPointer<vtkDoubleArray>::New();
  noiseArray->SetName("Field");
  noiseArray->SetNumberOfComponents(1);
  noiseArray->SetNumberOfTuples(nTuples);
  img->GetPointData()->AddArray(noiseArray);

  if(this->TimeProp > 0) {

    // Create an array to store time-step in
    auto tsArray = vtkSmartPointer<vtkDoubleArray>::New();
    tsArray->SetName("Time");
    tsArray->SetNumberOfComponents(1);
    tsArray->InsertNextValue(t);
    img->GetFieldData()->AddArray(tsArray);
  }

  return 1;
}

int ttkPerlinNoise::RequestData(vtkInformation *,
                                vtkInformationVector **,
                                vtkInformationVector *outputVector) {

  // Set dimensions, extent and number of tuples with values
  int dimZ = 0;
  int nTuples = 0;

  if(this->Resolution[2] == 1) {
    nTuples = this->Resolution[0] * this->Resolution[1];
  } else if(this->Resolution[2] > 1) {
    dimZ = this->Resolution[2] - 1;
    nTuples = this->Resolution[0] * this->Resolution[1] * this->Resolution[2];
  } else {
    this->printErr("Perlin noise dimension is invalid.");
    return 0;
  }
  int extent[6]
    = {0, this->Resolution[0] - 1, 0, this->Resolution[1] - 1, 0, dimZ};

  // Create timer for measuring execution
  ttk::Timer timer;

  // Switch on type of time-depdendency for the Perlin noise
  switch(this->TimeProp) {
    case 0: { // no time
      auto output = vtkImageData::GetData(outputVector);
      if(!output) {
        this->printErr("No output data available.");
        return 0;
      }

      this->printMsg("Creating Perlin noise image ["
                       + std::to_string(this->Resolution[0]) + ", "
                       + std::to_string(this->Resolution[1]) + ", "
                       + std::to_string(this->Resolution[2]) + "]",
                     0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      initializeOutput(output, extent, nTuples, -1);

      if(this->Resolution[2] == 1) {
        // Calculate 2D noise for image
        int dims[2] = {this->Resolution[0], this->Resolution[1]};
        switch(VTK_DOUBLE) {
          vtkTemplateMacro(this->perlin2Daux<VTK_TT>(
            dims, this->nOctaves, this->Scale, this->Frequency,
            this->Persistence,
            static_cast<VTK_TT *>(
              ttkUtils::GetVoidPointer(output->GetPointData()->GetArray(0)))));
        }
      } else {
        // Calculate 3D noise for image
        switch(VTK_DOUBLE) {
          vtkTemplateMacro(this->perlin3Daux<VTK_TT>(
            this->Resolution, this->nOctaves, this->Scale, this->Frequency,
            this->Persistence,
            static_cast<VTK_TT *>(
              ttkUtils::GetVoidPointer(output->GetPointData()->GetArray(0)))));
        }
      }

      this->printMsg("Creating Perlin noise image ["
                       + std::to_string(this->Resolution[0]) + ", "
                       + std::to_string(this->Resolution[1]) + ", "
                       + std::to_string(this->Resolution[2]) + "]",
                     1, timer.getElapsedTime(), this->threadNumber_);
      timer.reStart();

      break;
    }
    case 1: { // singular time-step
      auto output = vtkImageData::GetData(outputVector);
      if(!output) {
        this->printErr("No output data available.");
        return 0;
      }

      this->printMsg("Creating Perlin noise image ["
                       + std::to_string(this->Resolution[0]) + ", "
                       + std::to_string(this->Resolution[1]) + ", "
                       + std::to_string(this->Resolution[2])
                       + "], t = " + std::to_string(this->TimeStep),
                     0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      initializeOutput(output, extent, nTuples, this->TimeStep);

      // Check perlin dimension
      if(this->Resolution[2] == 1) {
        // Get 2D+T noise for chosen time-step
        int dims[2] = {this->Resolution[0], this->Resolution[1]};
        switch(VTK_DOUBLE) {
          vtkTemplateMacro(this->perlin2DTaux<VTK_TT>(
            dims, this->TimeStep, this->nOctaves, this->Scale, this->Frequency,
            this->Persistence,
            static_cast<VTK_TT *>(
              ttkUtils::GetVoidPointer(output->GetPointData()->GetArray(0)))));
        }
      } else {
        // Get 3D+T noise for chosen time-step
        switch(VTK_DOUBLE) {
          vtkTemplateMacro(this->perlin3DTaux<VTK_TT>(
            this->Resolution, this->TimeStep, this->nOctaves, this->Scale,
            this->Frequency, this->Persistence,
            static_cast<VTK_TT *>(
              ttkUtils::GetVoidPointer(output->GetPointData()->GetArray(0)))));
        }
      }

      this->printMsg("Creating Perlin noise image ["
                       + std::to_string(this->Resolution[0]) + ", "
                       + std::to_string(this->Resolution[1]) + ", "
                       + std::to_string(this->Resolution[2])
                       + "], t = " + std::to_string(this->TimeStep),
                     1, timer.getElapsedTime(), this->threadNumber_);
      timer.reStart();

      break;
    }
    case 2: { // Outputs a time-series
      auto output = vtkMultiBlockDataSet::GetData(outputVector);
      if(!output) {
        this->printErr("No output data available.");
        return 0;
      }

      // Loop through all time-steps
      this->printMsg(
        "Creating time series of " + std::to_string(this->TimeSeries)
          + " timesteps, interval = " + std::to_string(this->Interval)
          + ", of Perlin noise images [" + std::to_string(this->Resolution[0])
          + ", " + std::to_string(this->Resolution[1]) + ", "
          + std::to_string(this->Resolution[2]) + "]",
        0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      for(int t = 0; t < this->TimeSeries; t++) {
        // Calculate actual time with Interval
        double time = this->Interval * t;

        // Create a VTK image for current time-step using smart pointers
        auto image = vtkSmartPointer<vtkImageData>::New();

        initializeOutput(image, extent, nTuples, time);

        // Check perlin dimension
        if(this->Resolution[2] == 1) {
          int dims[2] = {this->Resolution[0], this->Resolution[1]};
          // Execute perlin for time-step
          switch(VTK_DOUBLE) {
            vtkTemplateMacro(this->perlin2DTaux<VTK_TT>(
              dims, time, this->nOctaves, this->Scale, this->Frequency,
              this->Persistence,
              static_cast<VTK_TT *>(
                ttkUtils::GetVoidPointer(image->GetPointData()->GetArray(0)))));
          }
        } else {
          // Execute perlin for time-step
          switch(VTK_DOUBLE) {
            vtkTemplateMacro(this->perlin3DTaux<VTK_TT>(
              this->Resolution, time, this->nOctaves, this->Scale,
              this->Frequency, this->Persistence,
              static_cast<VTK_TT *>(
                ttkUtils::GetVoidPointer(image->GetPointData()->GetArray(0)))));
          }
        }

        // Set image to a block in the output dataset
        size_t nBlocks = output->GetNumberOfBlocks();
        output->SetBlock(nBlocks, image);
      }
      this->printMsg(
        "Creating time series of " + std::to_string(this->TimeSeries)
          + " timesteps, interval = " + std::to_string(this->Interval)
          + ", of Perlin noise images [" + std::to_string(this->Resolution[0])
          + ", " + std::to_string(this->Resolution[1]) + ", "
          + std::to_string(this->Resolution[2]) + "]",
        1, timer.getElapsedTime(), this->threadNumber_);
      timer.reStart();
      break;
    }
    default: {
      this->printErr("No time option selected.");
      break;
    }
  }

  // return success
  return 1;
}
