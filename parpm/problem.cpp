
#include "problem.hpp"


void build_problem(double *x,
                   int zstart, int nz, double dz,
                   int ystart, int ny, double dy,
                   int xstart, int nx, double dx,
                   double k) {

  const double kdx = k*dx;
  const double kdy = k*dy;
  const double kdz = k*dz;

  for(int m=0; m<nz; m++)
    for(int j=0; j<ny; j++)
      for(int i=0; i<nx; i++)
        x[(m*ny+j)*nx+i] = sin(kdx*(xstart+i))*sin(kdy*(ystart+j))*sin(kdz*(zstart+m));

}

void build_solution(double *x,
                    int zstart, int nz, double dz,
                    int ystart, int ny, double dy,
                    int xstart, int nx, double dx,
                    double k) {

  const double kdx = k*dx;
  const double kdy = k*dy;
  const double kdz = k*dz;
  const double s   = 1./(3*k*k);

  for(int m=0; m<nz; m++)
    for(int j=0; j<ny; j++)
      for(int i=0; i<nx; i++)
        x[(m*ny+j)*nx+i] = s*sin(kdx*(xstart+i))*sin(kdy*(ystart+j))*sin(kdz*(zstart+m));

}
