/// \ingroup base
/// \class ttk::InputPointAdvection
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2022-06-07.
///
///
/// \b Related \b publication: \n
/// '???'
///
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <PerlinNoise.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The InputPointAdvection class provides methods to advect points in a vector
   * field.
   */
  class InputPointAdvection : virtual public Debug {

  public:
    InputPointAdvection();

    // Struct Point used for integration along the vector field
    struct Point {
      int pointId{-1};
      int timestep{-1};
      int birth{-1};
      int death{-1};
      double x{0.0};
      double y{0.0};
      double z{0.0};
      double v[3]{0.0, 0.0, 0.0};
      double amplitude{0.0};
      double variance{0.0};
      double rate{0.0};
      bool outsideDomain{false};

      Point() {
      }

      Point(const double &px, const double &py, const double &pz) {
        x = px;
        y = py;
        z = pz;
      }

      Point(const Point &p) {
        pointId = p.pointId;
        timestep = p.timestep;
        birth = p.birth;
        death = p.death;
        x = p.x;
        y = p.y;
        z = p.z;
        v[0] = p.v[0];
        v[1] = p.v[1];
        v[2] = p.v[2];
        amplitude = p.amplitude;
        variance = p.variance;
        rate = p.rate;
        outsideDomain = p.outsideDomain;
      }

      Point &operator=(const Point &p) {
        pointId = p.pointId;
        timestep = p.timestep;
        birth = p.birth;
        death = p.death;
        x = p.x;
        y = p.y;
        z = p.z;
        v[0] = p.v[0];
        v[1] = p.v[1];
        v[2] = p.v[2];
        amplitude = p.amplitude;
        variance = p.variance;
        rate = p.rate;
        outsideDomain = p.outsideDomain;
        return *this;
      }

      Point operator+(const Point &a) const {
        Point p;
        p.pointId = pointId;
        p.timestep = timestep;
        p.birth = birth;
        p.death = death;
        p.amplitude = amplitude;
        p.variance = variance;
        p.rate = rate;
        p.outsideDomain = outsideDomain;
        p.x = x + a.x;
        p.y = y + a.y;
        p.z = z + a.z;

        return p;
      }

      Point operator*(double k) {
        Point p;
        p.pointId = pointId;
        p.timestep = timestep;
        p.birth = birth;
        p.death = death;
        p.amplitude = amplitude;
        p.variance = variance;
        p.rate = rate;
        p.outsideDomain = outsideDomain;
        p.x = k * x;
        p.y = k * y;
        p.z = k * z;
        return p;
      }

      void setVelocity(const double vel[3]) {
        v[0] = vel[0];
        v[1] = vel[1];
        v[2] = vel[2];
      }
    };

    // Class used to determine which vector fields exist and can be used
    enum class VectorField {
      PerlinPerturbed,
      PerlinGradient,
      PosDiagonal,
      PosX,
      Sink,
      Saddle
    };

    void setStepLength(const double stepLength) {
      h_ = stepLength;
    }

    void setPerlinScaleFactor(const double psf) {
      psf_[0] = psf;
      psf_[1] = psf;
      psf_[2] = psf;
    }

    void setVectorField(const InputPointAdvection::VectorField &vf) {
      vf_ = vf;
    }

    Point sampleVectorField(const Point &p, double t) {
      Point dv(p);
      switch(vf_) {
        case InputPointAdvection::VectorField::PerlinPerturbed: {
          // Use perturbation in the 3D domain for the vector field
          pn.perlin4D<double>(
            p.x / psf_[0], p.y / psf_[1], p.z / psf_[2], t, dv.x);
          pn.perlin4D<double>(
            p.z / psf_[2], p.x / psf_[0], p.y / psf_[1], t, dv.y);
          pn.perlin4D<double>(
            p.y / psf_[1], p.z / psf_[2], p.x / psf_[0], t, dv.z);
          break;
        }
        case InputPointAdvection::VectorField::PerlinGradient: {
          // Calculate gradient using finite differences
          double s1, s2 = 0.0;
          pn.perlin4D<double>(
            (p.x - h_) / psf_[0], p.y / psf_[1], p.z / psf_[2], t, s1);
          pn.perlin4D<double>(
            (p.x + h_) / psf_[0], p.y / psf_[1], p.z / psf_[2], t, s2);
          dv.x = (s2 - s1) / (2 * h_);

          pn.perlin4D<double>(
            p.x / psf_[0], (p.y - h_) / psf_[1], p.z / psf_[2], t, s1);
          pn.perlin4D<double>(
            p.x / psf_[0], (p.y + h_) / psf_[1], p.z / psf_[2], t, s2);
          dv.y = (s2 - s1) / (2 * h_);

          pn.perlin4D<double>(
            p.x / psf_[0], p.y / psf_[1], (p.z - h_) / psf_[2], t, s1);
          pn.perlin4D<double>(
            p.x / psf_[0], p.y / psf_[1], (p.z + h_) / psf_[2], t, s2);
          dv.z = (s2 - s1) / (2 * h_);
          break;
        }
        case InputPointAdvection::VectorField::PosDiagonal: {
          // Go (1, 1, 1) along positive diagonal
          dv.x = 1.0;
          dv.y = 1.0;
          dv.z = 1.0;
          break;
        }
        case InputPointAdvection::VectorField::PosX: {
          // Go (1, 0, 0)
          dv.x = 1.0;
          dv.y = 0.0;
          dv.z = 0.0;
          break;
        }
        case InputPointAdvection::VectorField::Sink: {
          dv.x = -dv.x;
          dv.y = -dv.y;
          dv.z = -dv.z;
          break;
        }
        case InputPointAdvection::VectorField::Saddle: {
          dv.x = dv.x;
          dv.y = -dv.y;
          dv.z = 0;
          break;
        }
      }

      return dv;
    }

    int RK4(Point &prevP, Point &newP, double time) {
      // RK4 integration
      Point q1 = sampleVectorField(prevP, time) * h_;
      Point q2 = sampleVectorField(prevP + (q1 * 0.5), time) * h_;
      Point q3 = sampleVectorField(prevP + (q2 * 0.5), time) * h_;
      Point q4 = sampleVectorField(prevP + q3, time) * h_;

      Point vel = (q1 + q2 * 2 + q3 * 2 + q4) * (1.0 / 6);
      newP = prevP + vel;

      // Set velocity of previous point
      double v[3] = {vel.x, vel.y, vel.z};
      prevP.setVelocity(v);

      return 1;
    }

    int setVariables(const double stepLength,
                     const double psf,
                     const VectorField &vf) {
      // Set class variables
      setStepLength(stepLength);
      setPerlinScaleFactor(psf);
      setVectorField(vf);

      return 1;
    }

    int advect(std::vector<Point> &points,
               const std::vector<int> &ids,
               const int timestep,
               const double timeInterval) {
      // current time
      double time = timestep * timeInterval;
      // for all points to be advected
      for(unsigned int i = 0; i < ids.size(); i++) {
        Point newP;
        auto &curP = points[ids[i]];

        // Integrate using RK4
        RK4(curP, newP, time);

        // Change time-step
        newP.timestep = timestep + 1;
        newP.pointId = curP.pointId;
        curP = newP;
      }
      return 1;
    }

  protected:
    double h_{};
    double psf_[3]{};
    VectorField vf_{};
    PerlinNoise pn;
  }; // InputPointAdvection class

} // namespace ttk
