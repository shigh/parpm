
#pragma once

#include "types.hpp"

struct Grid {

  // Number of grid points
  int Nx;
  int Ny;
  int Nz;

  // Starting grid points
  int Nx0;
  int Ny0;
  int Nz0;

  // Length of domain in this grid
  // Domain in this grid is [L0,L0+L]
  FLOAT Lx;
  FLOAT Ly;
  FLOAT Lz;

  FLOAT Lx0;
  FLOAT Ly0;
  FLOAT Lz0;

  // MPI info
  int rank;
  int size;
  int[3][3][3] neighbors;

};
