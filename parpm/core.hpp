
#pragma once
#define _USE_MATH_DEFINES
#include <omp.h>
#include <math.h>
#include <stdlib.h>

void move_par(const int N, const double dt,
              double *zp, const double *vz,
              double *yp, const double *vy,
              double *xp, const double *vx);

void accel_par(const int N, const double dt, const double *qm,
               const double *Ez, double *vz,
               const double *Ey, double *vy,
               const double *Ex, double *vx);

void scale_array(const int N, double *x, double sx);

void scale_array_3_copy(const int N,
                        const double *z, const double sz, double *zc,
                        const double *y, const double sy, double *yc,
                        const double *x, const double sx, double *xc);
