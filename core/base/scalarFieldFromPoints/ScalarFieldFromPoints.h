/// \ingroup base
/// \class ttk::ScalarFieldFromPoints
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-10-11.
///
/// This module defines the %ScalarFieldFromPoints class that computes a 2D or
/// 3D the scalar field by using kernels on input points.
///
/// \b Related \b publication: \n
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <math.h>

namespace ttk {

  class ScalarFieldFromPoints : virtual public Debug {

  public:
    typedef double (*KERNEL)(const double &, const double &, const double &);

    static double
      Gaussian(const double &u, const double &bandwidth, const double &weight) {
      return weight * exp(-0.5 * (u / bandwidth));
    };

    ScalarFieldFromPoints() {
      this->setDebugMsgPrefix("ScalarFieldFromPoints");
    };
    ~ScalarFieldFromPoints(){};

    template <KERNEL k, typename DT>
    int computeScalarField3D(double *maxValData,
                             double *addValData,
                             int *maxIdData,
                             const DT *pointCoordiantes,
                             const int *pointIds,
                             const double *weights,
                             const double *constants,
                             const double *bounds,
                             const double *spacing,
                             const int *dims,
                             const size_t &nPoints,
                             const size_t &nPixels) const {
      ttk::Timer timer;

      this->printMsg("Computing Scalar Field 3D",
                     0, // progress form 0-1
                     0, // elapsed time so far
                     this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const double dx = spacing[0];
      const double dy = spacing[1];
      const double dz = spacing[2];
      const double dx2 = dx / 2.0;
      const double dy2 = dy / 2.0;
      const double dz2 = dz / 2.0;
      const int width = dims[0];
      const int height = dims[1];
      const int depth = dims[2];

      const double xBound = bounds[0] - dx2;
      const double yBound = bounds[2] - dy2;
      const double zBound = bounds[4] - dz2;

      // create locks
      omp_lock_t lock[nPixels];

      // clear data and init locks
      for(size_t i = 0; i < nPixels; i++) {
        maxValData[i] = 0.0;
        addValData[i] = 0.0;
        maxIdData[i] = -1;
        omp_init_lock(&(lock[i]));
      }

// compute scalar field
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(size_t i = 0; i < nPoints; i++) {
        const double xP = pointCoordiantes[i * 3 + 0];
        const double yP = pointCoordiantes[i * 3 + 1];
        const double zP = pointCoordiantes[i * 3 + 2];

        const int xi
          = std::min(width, std::max(0, (int)floor((xP - xBound) / dx)));
        const int yi
          = std::min(height, std::max(0, (int)floor((yP - yBound) / dy)));
        const int zi
          = std::min(depth, std::max(0, (int)floor((zP - zBound) / dz)));

        const int kdx = floor(3 * sqrt(constants[i]) / dx + 0.5);
        const int kdy = floor(3 * sqrt(constants[i]) / dy + 0.5);
        const int kdz = floor(3 * sqrt(constants[i]) / dz + 0.5);

        const int x0 = std::max(0, std::min(width - 1, xi - kdx));
        const int x1 = std::max(0, std::min(width - 1, xi + kdx));
        const int y0 = std::max(0, std::min(height - 1, yi - kdy));
        const int y1 = std::max(0, std::min(height - 1, yi + kdy));
        const int z0 = std::max(0, std::min(depth - 1, zi - kdz));
        const int z1 = std::max(0, std::min(depth - 1, zi + kdz));

        // for all points in the bandwidth interval, calculate scalar value
        for(int x = x0; x <= x1; x++) {
          for(int y = y0; y <= y1; y++) {
            for(int z = z0; z <= z1; z++) {
              double xxx = (x - xi) * dx;
              double yyy = (y - yi) * dy;
              double zzz = (z - zi) * dz;
              const double u = (xxx * xxx + yyy * yyy + zzz * zzz);
              const double ku = k(u, constants[i], weights[i]);

              int pixelIndex = z * width * height + y * width + x;
              omp_set_lock(&(lock[pixelIndex]));
              // Max mixture
              if(ku > maxValData[pixelIndex]) {
                maxValData[pixelIndex] = ku;
                maxIdData[pixelIndex] = pointIds[i];
              }

              // Additive mixture
              addValData[pixelIndex] += ku;
              omp_unset_lock(&(lock[pixelIndex]));
            }
          }
        }
      }

      // destroy all locks
      for(int i = 0, j = nPixels; i < j; i++) {
        omp_destroy_lock(&(lock[i]));
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing Scalar Field 3D",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

  }; // ScalarFieldFromPoints class

} // namespace ttk
