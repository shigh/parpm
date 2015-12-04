
#include <vector>
#include <cassert>

#include "types.hpp"
#include "fftwmpi.hpp"
#include "particle.hpp"

int pic_fft(int argc, char* argv[]) {

  // Temp test params
  const int nt  = 1;
  const int ppc = 1;

  // Global problem setup
  const FLOAT Lz = 2*M_PI;
  const FLOAT Ly = 2*M_PI;
  const FLOAT Lx = 2*M_PI;

  const ptrdiff_t N0 = 100;
  const ptrdiff_t N1 = 100;
  const ptrdiff_t N2 = 100;

  const FLOAT dz = Lz/N0;
  const FLOAT dy = Ly/N1;
  const FLOAT dx = Lx/N2;
  const FLOAT V  = dz*dy*dx;
  const FLOAT Vi = dz*dy*dx;

  // Build solver and get local params
  FFTWPoisson3DMPI solver(N0, Lz, N1, Ly, N2, Lx);
  const int nz = solver.get_nz();
  const int nx = solver.get_nx();
  const int ny = solver.get_ny();
  const int z0 = solver.get_z0();
  // Local domain length
  const FLOAT Lxl = nx*dx;
  const FLOAT Lyl = ny*dy;
  const FLOAT Lzl = nz*dz;

  // Field arrays
  std::vector<FLOAT> phi(nx*ny*nz, 0);
  std::vector<FLOAT> Ex(nx*ny*nz, 0);
  std::vector<FLOAT> Ey(nx*ny*nz, 0);
  std::vector<FLOAT> Ez(nx*ny*nz, 0);

  // Get MPI params
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int right = (rank+1)%size;
  const int left  = rank>0 ? rank-1 : size-1;
  assert((left>=0) && (right<size));

  // Communication buffers
  std::vector<FLOAT> send_right(nx*ny, 0);
  std::vector<FLOAT> send_left(nx*ny, 0);
  std::vector<FLOAT> recv_right(nx*ny, 0);
  std::vector<FLOAT> recv_left(nx*ny, 0);
  std::vector<FLOAT> send_parts_right;
  std::vector<FLOAT> recv_parts_right;
  std::vector<FLOAT> send_parts_left;
  std::vector<FLOAT> recv_parts_left;

  // Construct init particles
  ParticleList particles(ppc*ppc*ppc*nx*ny*nz);
  {
    const FLOAT pdx = dx/(ppc+1);
    const FLOAT pdy = dy/(ppc+1);
    const FLOAT pdz = dz/(ppc+1);
    FLOAT x0, y0, z0;
    for(int iz=0; iz<nz; ++iz)
      for(int iy=0; iy<ny; ++iy)
        for(int ix=0; ix<nx; ++ix) {
          x0 = ix*dx;
          y0 = iy*dy;
          z0 = iz*dz;
          for(int pz=0; pz<ppc; ++pz)
            for(int py=0; py<ppc; ++py)
              for(int px=0; px<ppc; ++px) {
                particles.xp_.push_back(x0+px*pdx);
                particles.yp_.push_back(y0+py*pdy);
                particles.zp_.push_back(z0+pz*pdz);

                // TODO: Init particle velocities
                particles.vx_.push_back(1.0);
                particles.vy_.push_back(1.0);
                particles.vz_.push_back(1.0);

                particles.q_.push_back(1.0);
              }

        }

  }

  // Main time stepping loop
  for(int it=0; it<nt; ++it) {

    // Weight particles to mesh

    // Solve for phi
    solver.solve(&phi[0]);

    // Solve for Exyz

    // Weight fields to particles

    // Accel and move

    // Communicate particles

    // Diagnostics

  }

  
}

int main(int argc, char* argv[]) {

  MPI_Init(&argc, &argv);
  fftw_mpi_init();
  
  pic_fft(argc, argv);

  MPI_Finalize();  

}
