
#include "par_core.hpp"

void move_par(const int N, const double dt,
              double *zp, const double *vz,
              double *yp, const double *vy,
              double *xp, const double *vx) {

  int i;
#pragma omp parallel for
  for(i=0; i<N; i++) {
    zp[i] = zp[i] + dt*vz[i];
    yp[i] = yp[i] + dt*vy[i];
    xp[i] = xp[i] + dt*vx[i];
  }

}


void accel_par(const int N, const double dt, const double *qm,
               const double *Ez, double *vz,
               const double *Ey, double *vy,
               const double *Ex, double *vx) {

  int i;
#pragma omp parallel for
  for(i=0; i<N; i++) {
    double dtqmi = dt*qm[i];
    vz[i] = vz[i] + dtqmi*Ez[i];
    vy[i] = vy[i] + dtqmi*Ey[i];
    vx[i] = vx[i] + dtqmi*Ex[i];
	
  }
	
}


void scale_array(const int N, double *x, double sx) {

  int i;
#pragma omp parallel for
  for(i=0; i<N; i++)
    x[i] = x[i]*sx;

}


void scale_array_3_copy(const int N,
                        const double *z, const double sz, double *zc,
                        const double *y, const double sy, double *yc,
                        const double *x, const double sx, double *xc) {

  int i;
#pragma omp parallel for
  for(i=0; i<N; i++) {
    zc[i] = z[i]*sz;
    yc[i] = y[i]*sy;
    xc[i] = x[i]*sx;
  }

}
