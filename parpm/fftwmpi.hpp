
#pragma once

#include <fftw3.h>
#include <fftw3-mpi.h>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <vector>

void setup_fftw_mpi();

void cleanup_fftw_mpi();

class FFTWPoisson3DMPI {

private:

  // fftw structures
  fftw_complex *out;
  double *in;
  fftw_plan p, pi;

  // Global problem dimensions
  ptrdiff_t N0, N1, N2;
  ptrdiff_t local_n0, local_0_start, alloc_local;
  double Lz, Ly, Lx;

  // Local problem dimensions
  int nx, ny, nz, nx_pad, z0;
  int nzc, nyc, nxc;

  // local k values
  std::vector<double> kz2_vals, ky2_vals, kx2_vals;

public:

  FFTWPoisson3DMPI(ptrdiff_t N0, double Lz,
                   ptrdiff_t N1, double Ly,
                   ptrdiff_t N2, double Lx);

  void solve(double *x);

  int get_nx();
  int get_ny();
  int get_nz();
  int get_z0();

  ~FFTWPoisson3DMPI();

};
