/// \ingroup base
/// \class ttk::PerlinNoise
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-06-04.
///
/// This module defines the %PerlinNoise class that computes a perlin noise
/// scalar field for the chosen dimensions.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {
  class PerlinNoise : virtual public Debug {

  public:
    PerlinNoise();
    /* HELPER FUNCTIONS */
    template <class dataType>
    dataType lerp(dataType t, dataType a, dataType b) const;

    template <class dataType>
    dataType fade(dataType t) const;

    /* DOT PRODUCT FUNCTIONS FOR 2,3,4 DIM */
    template <class dataType>
    dataType dot2D(int g[2], dataType o[2]) const;

    template <class dataType>
    dataType dot3D(int g[3], dataType o[3]) const;

    template <class dataType>
    dataType dot4D(int g[4], dataType o[4]) const;

    /* PERLIN FUNCTIONS */
    template <class dataType>
    int perlin2D(dataType x, dataType y, dataType &noise);

    template <class dataType>
    int perlin3D(dataType x, dataType y, dataType z, dataType &noise);

    template <class dataType>
    int
      perlin4D(dataType x, dataType y, dataType z, dataType t, dataType &noise);

    /* AUX FUNCTIONS TO ACCES PERLIN FUNCTIONS EASILY */
    template <class dataType>
    int perlin2Daux(int dims[2],
                    int nOctaves,
                    double scale,
                    int frequency,
                    double persistence,
                    dataType *outputData);

    template <class dataType>
    int perlin2DTaux(int dims[2],
                     dataType timeStep,
                     int nOctaves,
                     double scale,
                     int frequency,
                     double persistence,
                     dataType *outputData);

    template <class dataType>
    int perlin3Daux(int dims[3],
                    int nOctaves,
                    double scale,
                    int frequency,
                    double persistence,
                    dataType *outputData);

    template <class dataType>
    int perlin3DTaux(int dims[3],
                     dataType timeStep,
                     int nOctaves,
                     double scale,
                     int frequency,
                     double persistence,
                     dataType *outputData);

  private:
    std::vector<int> pT;

    // gradient vectors
    int grad2b[4][2] = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
    int grad3[12][3] = {{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
                        {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
                        {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}};
    int grad4[32][4] = {
      {0, 1, 1, 1},  {0, 1, 1, -1},  {0, 1, -1, 1},  {0, 1, -1, -1},
      {0, -1, 1, 1}, {0, -1, 1, -1}, {0, -1, -1, 1}, {0, -1, -1, -1},
      {1, 0, 1, 1},  {1, 0, 1, -1},  {1, 0, -1, 1},  {1, 0, -1, -1},
      {-1, 0, 1, 1}, {-1, 0, 1, -1}, {-1, 0, -1, 1}, {-1, 0, -1, -1},
      {1, 1, 0, 1},  {1, 1, 0, -1},  {1, -1, 0, 1},  {1, -1, 0, -1},
      {-1, 1, 0, 1}, {-1, 1, 0, -1}, {-1, -1, 0, 1}, {-1, -1, 0, -1},
      {1, 1, 1, 0},  {1, 1, -1, 0},  {1, -1, 1, 0},  {1, -1, -1, 0},
      {-1, 1, 1, 0}, {-1, 1, -1, 0}, {-1, -1, 1, 0}, {-1, -1, -1, 0},
    };

  }; // PerlinNoise class
} // namespace ttk

template <class dataType>
dataType ttk::PerlinNoise::lerp(dataType t, dataType a, dataType b) const {
  return (1 - t) * a + t * b;
}

template <class dataType>
dataType ttk::PerlinNoise::fade(dataType t) const {
  return t * t * t * (t * (t * 6 - 15) + 10);
}

template <class dataType>
dataType ttk::PerlinNoise::dot2D(int g[2], dataType o[2]) const {
  return (g[0] * o[0]) + (g[1] * o[1]);
}

template <class dataType>
dataType ttk::PerlinNoise::dot3D(int g[3], dataType o[3]) const {
  return (g[0] * o[0]) + (g[1] * o[1]) + (g[2] * o[2]);
}

template <class dataType>
dataType ttk::PerlinNoise::dot4D(int g[4], dataType o[4]) const {
  return (g[0] * o[0]) + (g[1] * o[1]) + (g[2] * o[2]) + (g[3] * o[3]);
}

template <class dataType>
int ttk::PerlinNoise::perlin2D(dataType x, dataType y, dataType &noise) {
  // Find unit cube that contains point (from 256 unit cubes)
  int X = (int)floor(x) & 255; // like modulus
  int Y = (int)floor(y) & 255;

  // Relative coords in current unit cube
  dataType xRel = x - dataType(floor(x));
  dataType yRel = y - dataType(floor(y));

  dataType u = fade<dataType>(xRel);
  dataType v = fade<dataType>(yRel);

  // 01___11
  // |     |
  // 00___10

  // Calculate hash indices to get gradients (& 3 because we have four different
  // possible gradients)
  int id00 = pT[pT[X] + Y] & 3;
  int id10 = pT[pT[X + 1] + Y] & 3;
  int id01 = pT[pT[X] + Y + 1] & 3;
  int id11 = pT[pT[X + 1] + Y + 1] & 3;

  // Get offset vectors from vertex points to current point
  dataType off00[2] = {xRel, yRel};
  dataType off10[2] = {dataType(xRel - 1), yRel};
  dataType off01[2]{xRel, dataType(yRel - 1)};
  dataType off11[2] = {dataType(xRel - 1), dataType(yRel - 1)};

  // Calculate dot products
  dataType dot00 = dot2D<dataType>(grad2b[id00], off00);
  dataType dot10 = dot2D<dataType>(grad2b[id10], off10);
  dataType dot01 = dot2D<dataType>(grad2b[id01], off01);
  dataType dot11 = dot2D<dataType>(grad2b[id11], off11);

  // Linear interpolation along x-axis
  dataType resy0 = lerp<dataType>(u, dot00, dot10);
  dataType resy1 = lerp<dataType>(u, dot01, dot11);

  // Linear interpolation along y-axis
  noise = lerp<dataType>(v, resy0, resy1);

  return 1;
}

template <class dataType>
int ttk::PerlinNoise::perlin3D(dataType x,
                               dataType y,
                               dataType z,
                               dataType &noise) {
  // Find unit cube that contains point (from 256 unit cubes)
  int X = (int)floor(x) & 255; // like modulus
  int Y = (int)floor(y) & 255;
  int Z = (int)floor(z) & 255;

  // Relative coords in current unit cube
  dataType xRel = x - floor(x);
  dataType yRel = y - floor(y);
  dataType zRel = z - floor(z);

  dataType u = fade(xRel);
  dataType v = fade(yRel);
  dataType w = fade(zRel);

  // Calculate hash indices to get gradients (% 12 because we have listed 12
  // different possible gradients)
  int id000 = pT[pT[pT[X] + Y] + Z] % 12;
  int id100 = pT[pT[pT[X + 1] + Y] + Z] % 12;
  int id010 = pT[pT[pT[X] + Y + 1] + Z] % 12;
  int id110 = pT[pT[pT[X + 1] + Y + 1] + Z] % 12;
  int id001 = pT[pT[pT[X] + Y] + Z + 1] % 12;
  int id101 = pT[pT[pT[X + 1] + Y] + Z + 1] % 12;
  int id011 = pT[pT[pT[X] + Y + 1] + Z + 1] % 12;
  int id111 = pT[pT[pT[X + 1] + Y + 1] + Z + 1] % 12;

  // Get offset vectors from vertex points to current point
  dataType xRelMin = xRel - 1; // Set as dataType to avoid narrowing conversions
  dataType yRelMin = yRel - 1; // Set as dataType to avoid narrowing conversions
  dataType zRelMin = zRel - 1; // Set as dataType to avoid narrowing conversions
  dataType off000[3] = {xRel, yRel, zRel};
  dataType off100[3] = {xRelMin, yRel, zRel};
  dataType off010[3] = {xRel, yRelMin, zRel};
  dataType off110[3] = {xRelMin, yRelMin, zRel};
  dataType off001[3] = {xRel, yRel, zRelMin};
  dataType off101[3] = {xRelMin, yRel, zRelMin};
  dataType off011[3] = {xRel, yRelMin, zRelMin};
  dataType off111[3] = {xRelMin, yRelMin, zRelMin};

  // Calculate dot products
  dataType dot000 = dot3D(grad3[id000], off000);
  dataType dot100 = dot3D(grad3[id100], off100);
  dataType dot010 = dot3D(grad3[id010], off010);
  dataType dot110 = dot3D(grad3[id110], off110);
  dataType dot001 = dot3D(grad3[id001], off001);
  dataType dot101 = dot3D(grad3[id101], off101);
  dataType dot011 = dot3D(grad3[id011], off011);
  dataType dot111 = dot3D(grad3[id111], off111);

  // Linear interpolation along x-axis
  dataType resy0z0 = lerp(u, dot000, dot100);
  dataType resy1z0 = lerp(u, dot010, dot110);
  dataType resy0z1 = lerp(u, dot001, dot101);
  dataType resy1z1 = lerp(u, dot011, dot111);

  // Linear interpolation along y-axis
  dataType resz0 = lerp(v, resy0z0, resy1z0);
  dataType resz1 = lerp(v, resy0z1, resy1z1);

  // Linear interpolation along z-axis
  noise = lerp(w, resz0, resz1);

  return 1;
}

template <class dataType>
int ttk::PerlinNoise::perlin4D(
  dataType x, dataType y, dataType z, dataType t, dataType &noise) {
  // Find unit cube that contains point (from 256 unit cubes)
  int X = (int)floor(x) & 255; // like modulus
  int Y = (int)floor(y) & 255;
  int Z = (int)floor(z) & 255;
  int T = (int)floor(t) & 255;

  // Relative coords in current unit cube
  dataType xRel = x - floor(x);
  dataType yRel = y - floor(y);
  dataType zRel = z - floor(z);
  dataType tRel = t - floor(t);

  dataType u = fade(xRel);
  dataType v = fade(yRel);
  dataType w = fade(zRel);
  dataType r = fade(tRel);

  // Calculate hash indices to get gradients (% 32 because we have listed 32
  // different possible gradients)
  int id0000 = pT[pT[pT[pT[X] + Y] + Z] + T] % 32;
  int id1000 = pT[pT[pT[pT[X + 1] + Y] + Z] + T] % 32;
  int id0100 = pT[pT[pT[pT[X] + Y + 1] + Z] + T] % 32;
  int id1100 = pT[pT[pT[pT[X + 1] + Y + 1] + Z] + T] % 32;
  int id0010 = pT[pT[pT[pT[X] + Y] + Z + 1] + T] % 32;
  int id1010 = pT[pT[pT[pT[X + 1] + Y] + Z + 1] + T] % 32;
  int id0110 = pT[pT[pT[pT[X] + Y + 1] + Z + 1] + T] % 32;
  int id1110 = pT[pT[pT[pT[X + 1] + Y + 1] + Z + 1] + T] % 32;
  int id0001 = pT[pT[pT[pT[X] + Y] + Z] + T + 1] % 32;
  int id1001 = pT[pT[pT[pT[X + 1] + Y] + Z] + T + 1] % 32;
  int id0101 = pT[pT[pT[pT[X] + Y + 1] + Z] + T + 1] % 32;
  int id1101 = pT[pT[pT[pT[X + 1] + Y + 1] + Z] + T + 1] % 32;
  int id0011 = pT[pT[pT[pT[X] + Y] + Z + 1] + T + 1] % 32;
  int id1011 = pT[pT[pT[pT[X + 1] + Y] + Z + 1] + T + 1] % 32;
  int id0111 = pT[pT[pT[pT[X] + Y + 1] + Z + 1] + T + 1] % 32;
  int id1111 = pT[pT[pT[pT[X + 1] + Y + 1] + Z + 1] + T + 1] % 32;

  // Get offset vectors from vertex points to current point
  dataType xRelMin = xRel - 1; // Set as dataType to avoid narrowing conversions
  dataType yRelMin = yRel - 1; // Set as dataType to avoid narrowing conversions
  dataType zRelMin = zRel - 1; // Set as dataType to avoid narrowing conversions
  dataType tRelMin = tRel - 1; // Set as dataType to avoid narrowing conversions
  dataType off0000[4] = {xRel, yRel, zRel, tRel};
  dataType off1000[4] = {xRelMin, yRel, zRel, tRel};
  dataType off0100[4] = {xRel, yRelMin, zRel, tRel};
  dataType off1100[4] = {xRelMin, yRelMin, zRel, tRel};
  dataType off0010[4] = {xRel, yRel, zRelMin, tRel};
  dataType off1010[4] = {xRelMin, yRel, zRelMin, tRel};
  dataType off0110[4] = {xRel, yRelMin, zRelMin, tRel};
  dataType off1110[4] = {xRelMin, yRelMin, zRelMin, tRel};
  dataType off0001[4] = {xRel, yRel, zRel, tRelMin};
  dataType off1001[4] = {xRelMin, yRel, zRel, tRelMin};
  dataType off0101[4] = {xRel, yRelMin, zRel, tRelMin};
  dataType off1101[4] = {xRelMin, yRelMin, zRel, tRelMin};
  dataType off0011[4] = {xRel, yRel, zRelMin, tRelMin};
  dataType off1011[4] = {xRelMin, yRel, zRelMin, tRelMin};
  dataType off0111[4] = {xRel, yRelMin, zRelMin, tRelMin};
  dataType off1111[4] = {xRelMin, yRelMin, zRelMin, tRelMin};

  // Calculate dot products
  dataType dot0000 = dot4D(grad4[id0000], off0000);
  dataType dot1000 = dot4D(grad4[id1000], off1000);
  dataType dot0100 = dot4D(grad4[id0100], off0100);
  dataType dot1100 = dot4D(grad4[id1100], off1100);
  dataType dot0010 = dot4D(grad4[id0010], off0010);
  dataType dot1010 = dot4D(grad4[id1010], off1010);
  dataType dot0110 = dot4D(grad4[id0110], off0110);
  dataType dot1110 = dot4D(grad4[id1110], off1110);
  dataType dot0001 = dot4D(grad4[id0001], off0001);
  dataType dot1001 = dot4D(grad4[id1001], off1001);
  dataType dot0101 = dot4D(grad4[id0101], off0101);
  dataType dot1101 = dot4D(grad4[id1101], off1101);
  dataType dot0011 = dot4D(grad4[id0011], off0011);
  dataType dot1011 = dot4D(grad4[id1011], off1011);
  dataType dot0111 = dot4D(grad4[id0111], off0111);
  dataType dot1111 = dot4D(grad4[id1111], off1111);

  // Linear interpolation along x-axis
  dataType resy0z0t0 = lerp(u, dot0000, dot1000);
  dataType resy1z0t0 = lerp(u, dot0100, dot1100);
  dataType resy0z1t0 = lerp(u, dot0010, dot1010);
  dataType resy1z1t0 = lerp(u, dot0110, dot1110);
  dataType resy0z0t1 = lerp(u, dot0001, dot1001);
  dataType resy1z0t1 = lerp(u, dot0101, dot1101);
  dataType resy0z1t1 = lerp(u, dot0011, dot1011);
  dataType resy1z1t1 = lerp(u, dot0111, dot1111);

  // Linear interpolation along y-axis
  dataType resz0t0 = lerp(v, resy0z0t0, resy1z0t0);
  dataType resz1t0 = lerp(v, resy0z1t0, resy1z1t0);
  dataType resz0t1 = lerp(v, resy0z0t1, resy1z0t1);
  dataType resz1t1 = lerp(v, resy0z1t1, resy1z1t1);

  // Linear interpolation along z-axis
  dataType resw0 = lerp(w, resz0t0, resz1t0);
  dataType resw1 = lerp(w, resz0t1, resz1t1);

  // Linear interpolation along t-axis
  noise = lerp(r, resw0, resw1);

  return 1;
}

template <class dataType>
int ttk::PerlinNoise::perlin2Daux(int dims[2],
                                  int nOctaves,
                                  double scale,
                                  int frequency,
                                  double persistence,
                                  dataType *outputData) {
  // Get dimensions and calculate divider used to make coordinates between grid
  // points
  int dimX = dims[0];
  int dimY = dims[1];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int y = 0; y < dimY; y++) {
    for(int x = 0; x < dimX; x++) {
      // Acccumulate noise values over the number of octaves
      dataType accNoise = 0;
      dataType xD = ((dataType)x) / scale;
      dataType yD = ((dataType)y) / scale;

      for(int o = 0; o < nOctaves; o++) {
        dataType noise;
        this->perlin2D<dataType>(
          std::pow(frequency, o) * xD, std::pow(frequency, o) * yD, noise);
        accNoise += std::pow(persistence, o) * noise;
      }

      // Calculate index and set noise to it in the array
      int idx = x + dimX * y;
      outputData[idx] = accNoise;
    }
  }
  return 1;
}

template <class dataType>
int ttk::PerlinNoise::perlin2DTaux(int dims[2],
                                   dataType timeStep,
                                   int nOctaves,
                                   double scale,
                                   int frequency,
                                   double persistence,
                                   dataType *outputData) {
  // Get dimensions and calculate divider used to make coordinates between grid
  // points
  int dimX = dims[0];
  int dimY = dims[1];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int y = 0; y < dimY; y++) {
    for(int x = 0; x < dimX; x++) {
      // Acccumulate noise values over the number of octaves
      dataType accNoise = 0;
      dataType xD = ((dataType)x) / scale;
      dataType yD = ((dataType)y) / scale;

      for(int o = 0; o < nOctaves; o++) {
        dataType noise;
        this->perlin3D<dataType>(std::pow(frequency, o) * xD,
                                 std::pow(frequency, o) * yD, timeStep, noise);
        accNoise += std::pow(persistence, o) * noise;
      }

      // Calculate index and set noise to it in the array
      int idx = x + dimX * y;
      outputData[idx] = accNoise;
    }
  }
  return 1;
}

template <class dataType>
int ttk::PerlinNoise::perlin3Daux(int dims[3],
                                  int nOctaves,
                                  double scale,
                                  int frequency,
                                  double persistence,
                                  dataType *outputData) {
  // Get dimensions and calculate divider used to make coordinates between grid
  // points
  int dimX = dims[0];
  int dimY = dims[1];
  int dimZ = dims[2];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int z = 0; z < dimZ; z++) {
    for(int y = 0; y < dimY; y++) {
      for(int x = 0; x < dimX; x++) {
        // Acccumulate noise values over the number of octaves
        dataType accNoise = 0;
        dataType xD = ((dataType)x) / scale;
        dataType yD = ((dataType)y) / scale;
        dataType zD = ((dataType)z) / scale;

        for(int o = 0; o < nOctaves; o++) {
          dataType noise;
          this->perlin3D<dataType>(std::pow(frequency, o) * xD,
                                   std::pow(frequency, o) * yD,
                                   std::pow(frequency, o) * zD, noise);
          accNoise += std::pow(persistence, o) * noise;
        }

        // Calculate index and set noise to it in the array
        int idx = x + dimX * y + dimX * dimY * z;
        outputData[idx] = accNoise;
      }
    }
  }
  return 1;
}

template <class dataType>
int ttk::PerlinNoise::perlin3DTaux(int dims[3],
                                   dataType timeStep,
                                   int nOctaves,
                                   double scale,
                                   int frequency,
                                   double persistence,
                                   dataType *outputData) {
  // Get dimensions and calculate divider used to make coordinates between grid
  // points
  int dimX = dims[0];
  int dimY = dims[1];
  int dimZ = dims[2];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int z = 0; z < dimZ; z++) {
    for(int y = 0; y < dimY; y++) {
      for(int x = 0; x < dimX; x++) {
        // Acccumulate noise values over the number of octaves
        dataType accNoise = 0;
        dataType xD = ((dataType)x) / scale;
        dataType yD = ((dataType)y) / scale;
        dataType zD = ((dataType)z) / scale;

        for(int o = 0; o < nOctaves; o++) {
          dataType noise;
          this->perlin4D<dataType>(
            std::pow(frequency, o) * xD, std::pow(frequency, o) * yD,
            std::pow(frequency, o) * zD, timeStep, noise);
          accNoise += std::pow(persistence, o) * noise;
        }

        // Calculate index and set noise to it in the array
        int idx = x + dimX * y + dimX * dimY * z;
        outputData[idx] = accNoise;
      }
    }
  }
  return 1;
}