
#include <vector>
#include <algorithm>
#include <cassert>
#include <fstream>
#include "omp.h"
#include "mpi.h"

#include "types.hpp"
#include "fftwmpi.hpp"
#include "particle.hpp"
#include "interp.hpp"
#include "core.hpp"

int pic_fft(int argc, char* argv[]) {

  // Temp test params

  assert(argc==4);
  const int nt       = atoi(argv[1]);
  const ptrdiff_t N0 = atoi(argv[2]);
  const ptrdiff_t N1 = N0;
  const ptrdiff_t N2 = N0;
  const int ppc      = atoi(argv[3]);

  int _n_threads;
#pragma omp parallel
  {
    _n_threads  = omp_get_num_threads();
  }
  const int n_threads = _n_threads;
  const FLOAT dt = 0.1;

  // Global problem setup
  const FLOAT Lz = 2*M_PI;
  const FLOAT Ly = 2*M_PI;
  const FLOAT Lx = 2*M_PI;

  const FLOAT dz = Lz/N0;
  const FLOAT dy = Ly/N1;
  const FLOAT dx = Lx/N2;
  const FLOAT V  = dz*dy*dx;
  const FLOAT Vi = 1.0/dz*dy*dx;

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
  const FLOAT Lz0 = z0*dz;

  // Field arrays
  std::vector<FLOAT> phi(nx*ny*nz, 0);
  std::vector<FLOAT> Ex(nx*ny*nz, 0);
  std::vector<FLOAT> Ey(nx*ny*nz, 0);
  std::vector<FLOAT> Ez(nx*ny*nz, 0);


  // Communication buffers
  std::vector<FLOAT> send_parts_right;
  std::vector<FLOAT> recv_parts_right;
  std::vector<FLOAT> send_parts_left;
  std::vector<FLOAT> recv_parts_left;

  // Construct particles
  ParticleList particles(ppc*ppc*ppc*nx*ny*nz);
  {
    const FLOAT pdx = dx/(ppc+1);
    const FLOAT pdy = dy/(ppc+1);
    const FLOAT pdz = dz/(ppc+1);
    FLOAT xl, yl, zl, v;
    for(int iz=0; iz<nz; ++iz)
      for(int iy=0; iy<ny; ++iy)
        for(int ix=0; ix<nx; ++ix) {
          xl = ix*dx;
          yl = iy*dy;
          zl = iz*dz;
          for(int pz=0; pz<ppc; ++pz)
            for(int py=0; py<ppc; ++py)
              for(int px=0; px<ppc; ++px) {
                particles.xp_.push_back(xl+px*pdx);
                particles.yp_.push_back(yl+py*pdy);
                particles.zp_.push_back(zl+pz*pdz);

                v = (((FLOAT)rand())/RAND_MAX-.5)*2.0;
                particles.vx_.push_back(v);
                v = (((FLOAT)rand())/RAND_MAX-.5)*2.0;
                particles.vy_.push_back(v);
                v = (((FLOAT)rand())/RAND_MAX-.5)*2.0;
                particles.vz_.push_back(v);
                particles.q_.push_back(1.0);
              }

        }

  }

  int ind;
  std::vector<FLOAT> Exp, Eyp, Ezp;
  std::vector<FLOAT> phi_send_r(ny*nx, 0);
  std::vector<FLOAT> phi_send_l(ny*nx, 0);
  std::vector<FLOAT> phi_recv_r(ny*nx, 0);
  std::vector<FLOAT> phi_recv_l(ny*nx, 0);
  
  // Diagnostic vectors
  std::vector<int> vn_send_r(nt,0), vn_send_l(nt,0);
  std::vector<int> vn_recv_r(nt,0), vn_recv_l(nt,0);
  std::vector<FLOAT> t_field_solve(nt,0), t_weight(nt,0), t_interp(nt,0);
  std::vector<FLOAT> t_calc_E(nt,0), t_comm(nt,0), t_loop(nt,0);
  std::vector<FLOAT> t_accel(nt,0), t_move(nt,0);
  std::vector<int> n_parts_local(nt,0);
  FLOAT t_start, t_end;
  FLOAT t_loop_start, t_loop_end;

  // Main time stepping loop
  for(int it=0; it<nt; ++it) {

    t_loop_start = MPI_Wtime();
    // Weight particles to mesh
    t_start = MPI_Wtime();
    weight_cic_par(nz, ny, nx, &phi[0], particles.size(),
                   &particles.zp_[0], dz,
                   &particles.yp_[0], dy,
                   &particles.xp_[0], dx,
                   &particles.q_[0]);
    t_end = MPI_Wtime();
    t_weight.at(it) = t_end-t_start;

    // Solve for phi
    t_start = MPI_Wtime();
    for(int i=0; i<phi.size(); ++i)
      phi.at(i) *= Vi;

    solver.solve(&phi[0]);
    t_end = MPI_Wtime();
    t_field_solve.at(it) = t_end-t_start;


    // Calc Exyz at grid points
    t_start = MPI_Wtime();

    for(int iz=1; iz<nz-1; ++iz)
      for(int iy=1; iy<ny-1; ++iy)
        for(int ix=1; ix<nx-1; ++ix) {
          ind = iz*ny*nx+iy*nx+ix;
          Ex.at(ind) = (phi.at(ind+1)-phi.at(ind-1))/(2*dx);
          Ey.at(ind) = (phi.at(ind+nx)-phi.at(ind-nx))/(2*dy);
          Ez.at(ind) = (phi.at(ind+nx*ny)-phi.at(ind-nx*ny))/(2*dz);
        }

    for(int iy=1; iy<ny-1; ++iy)
      for(int ix=1; ix<nx-1; ++ix) {
        ind = (nz-2)*ny*nx;
        Ez.at(ind+iy*nx+ix) = (phi.at(ind+nx*ny)-phi.at(ind-nx*ny))/(2*dz);
        ind = nx*ny;
        Ez.at(ind+iy*nx+ix) = (phi.at(ind+nx*ny)-phi.at(ind-nx*ny))/(2*dz);
      }
    t_end = MPI_Wtime();
    t_calc_E.at(it) = t_end-t_start;

    // Weight Exyz to particles
    t_start = MPI_Wtime();
    Exp.resize(particles.size());
    Eyp.resize(particles.size());
    Ezp.resize(particles.size());

    interp_cic_par(nz, ny, nx, &Ex[0],
                   particles.size(),
                   &particles.zp_[0], dz,
                   &particles.yp_[0], dy,
                   &particles.xp_[0], dx,
                   &Exp[0]);
    interp_cic_par(nz, ny, nx, &Ey[0],
                   particles.size(),
                   &particles.zp_[0], dz,
                   &particles.yp_[0], dy,
                   &particles.xp_[0], dx,
                   &Eyp[0]);
    interp_cic_par(nz, ny, nx, &Ez[0],
                   particles.size(),
                   &particles.zp_[0], dz,
                   &particles.yp_[0], dy,
                   &particles.xp_[0], dx,
                   &Ezp[0]);
    t_end = MPI_Wtime();
    t_interp.at(it) = t_end-t_start;

    // Accel and move
    t_start = MPI_Wtime();
    accel_par(particles.size(), dt, &particles.q_[0],
              &Ezp[0], &particles.vz_[0],
              &Eyp[0], &particles.vy_[0],
              &Exp[0], &particles.vx_[0]);
    t_end = MPI_Wtime();
    t_accel.at(it) = t_end-t_start;


    t_start = MPI_Wtime();
    move_par(particles.size(), dt,
             &particles.zp_[0], &particles.vz_[0],
             &particles.yp_[0], &particles.vy_[0],
             &particles.xp_[0], &particles.vx_[0]);
    t_end = MPI_Wtime();
    t_move.at(it) = t_end-t_start;

    // Find particles to send
    t_start = MPI_Wtime();

    for(int ipart=0; ipart<particles.size(); ++ipart) {
      // Normalize particle positions
      if(particles.xp_.at(ipart)<0)
        particles.xp_.at(ipart) += Lxl;
      if(particles.yp_.at(ipart)<0)
        particles.yp_.at(ipart) += Lyl;
      if(particles.zp_.at(ipart)<0)
        particles.zp_.at(ipart) += Lzl;
      if(particles.xp_.at(ipart)>=Lxl)
        particles.xp_.at(ipart) -= Lxl;
      if(particles.yp_.at(ipart)>=Lyl)
        particles.yp_.at(ipart) -= Lyl;
      if(particles.zp_.at(ipart)>=Lzl)
        particles.zp_.at(ipart) -= Lzl;

      assert(particles.xp_.at(ipart)>=0 && particles.xp_.at(ipart)<Lxl);
      assert(particles.yp_.at(ipart)>=0 && particles.yp_.at(ipart)<Lyl);
      assert(particles.zp_.at(ipart)>=0 && particles.zp_.at(ipart)<Lzl);

    }

    t_end = MPI_Wtime();
    t_comm.at(it) = t_end-t_start;


    t_loop_end = MPI_Wtime();
    t_loop.at(it) = t_loop_end-t_loop_start;

    std::cout << it+1 << ' ' << t_loop.at(it) << std::endl;;

  }

  // Save diagnostic results
  if(true) {
    std::ofstream odiag;
    odiag.open("./diag.dat");
    // Write header
    odiag << "rank" << ' ' << "iter" << ' '
          << "size" << ' ' << "n_threads" <<  ' '
          << "ppc" << ' ' << "nt" << ' '
          << "t_loop" << ' '
          << "t_field_solve" << ' '
          << "t_weight" << ' '
          << "t_interp" << ' '
          << "t_calc_E" << ' '
          << "t_comm" << ' '
          << "t_accel" << ' '
          << "t_move" << ' '
          << "n_local" << ' '
          << "n_sent" << std::endl;

    for(int i=0; i<nt; ++i) {
      odiag << 0 << ' '
            << i+1 << ' '<< 1 << ' '
            << n_threads << ' '
            << ppc << ' ' << nt << ' '
            << t_loop.at(i) << ' '
            << t_field_solve.at(i) << ' '
            << t_weight.at(i) << ' '
            << t_interp.at(i) << ' '
            << t_calc_E.at(i) << ' '
            << t_comm.at(i) << ' '
            << t_accel.at(i) << ' '
            << t_move.at(i) << ' '
            << particles.size() << ' '
            << 0 << std::endl;
    }

    odiag.close();
  }

}

int main(int argc, char* argv[]) {

  MPI_Init(&argc, &argv);
  fftw_mpi_init();
  
  pic_fft(argc, argv);

  MPI_Finalize();  

}

