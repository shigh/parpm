
#include <vector>
#include <math.h>

// 3D test problem setup
void build_problem(double *x,
                   int zstart, int nz, double dz,
                   int ystart, int ny, double dy,
                   int xstart, int nx, double dx,
                   double k);

void build_solution(double *x,
                    int zstart, int nz, double dz,
                    int ystart, int ny, double dy,
                    int xstart, int nx, double dx,
                    double k);
