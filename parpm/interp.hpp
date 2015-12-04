
#pragma once

#include <math.h>
#include <omp.h>
#include <vector>
#include <cassert>

void interp_cic_par(const int nz, const int ny, const int nx, const double *vals,
                    const int N, const double *z, const double dz,
                    const double *y, const double dy,
                    const double *x, const double dx, double *c);

void weight_cic_par(const int nz, const int ny, const int nx, double *grid,
                    const int N, const double *z, const double dz, const double *y, const double dy,
                    const double *x, const double dx, const double *q);
