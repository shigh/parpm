
#include <vector>
#include <algorithm>
#include <cassert>
#include <fstream>
#include "mpi.h"
#include "omp.h"

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

  // Get MPI params
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int right = (rank+1)%size;
  const int left  = rank>0 ? rank-1 : size-1;
  assert((left>=0) && (right<size));

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

  // Get global initial number of particles
  int init_n_parts;
  {
    int n_parts;
    n_parts = particles.size();
    MPI_Reduce(&n_parts, &init_n_parts, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
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
  MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Request comm_phi[2];
    for(int i=0; i<ny*nx; ++i) {
      phi_send_r.at(i) = phi.at(ny*nx*(nz-1)+i);
      phi_send_l.at(i) = phi.at(i);
    }
    MPI_Isend(&phi_send_r[0], ny*nx, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &comm_phi[0]);
    MPI_Isend(&phi_send_l[0], ny*nx, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &comm_phi[1]);

    for(int iz=1; iz<nz-1; ++iz)
      for(int iy=1; iy<ny-1; ++iy)
        for(int ix=1; ix<nx-1; ++ix) {
          ind = iz*ny*nx+iy*nx+ix;
          Ex.at(ind) = (phi.at(ind+1)-phi.at(ind-1))/(2*dx);
          Ey.at(ind) = (phi.at(ind+nx)-phi.at(ind-nx))/(2*dy);
          Ez.at(ind) = (phi.at(ind+nx*ny)-phi.at(ind-nx*ny))/(2*dz);
        }

    MPI_Recv(&phi_recv_l[0], ny*nx, MPI_DOUBLE, left, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&phi_recv_r[0], ny*nx, MPI_DOUBLE, right, 1,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(int iy=1; iy<ny-1; ++iy)
      for(int ix=1; ix<nx-1; ++ix) {
        ind = (nz-2)*ny*nx;
        Ez.at(ind+iy*nx+ix) = (phi.at(ind+nx*ny)-phi.at(ind-nx*ny))/(2*dz);
        ind = nx*ny;
        Ez.at(ind+iy*nx+ix) = (phi.at(ind+nx*ny)-phi.at(ind-nx*ny))/(2*dz);
      }
    MPI_Waitall(2, comm_phi, MPI_STATUS_IGNORE);
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
    std::vector<int> to_send_right, to_send_left;
    for(int ipart=0; ipart<particles.size(); ++ipart)
      if(particles.zp_.at(ipart)>=Lzl)
        to_send_right.push_back(ipart);
      else if(particles.zp_.at(ipart)<0)
        to_send_left.push_back(ipart);

    // Communicate particles
    send_parts_right.resize(std::max<int>(1, to_send_right.size()*partFields));
    send_parts_left.resize(std::max<int>(1, to_send_left.size()*partFields));
    int ipart;
    for(int i=0; i<to_send_right.size(); ++i) {
      ipart = to_send_right.at(i);
      send_parts_right.at(i*partFields+0) = particles.xp_.at(ipart);
      send_parts_right.at(i*partFields+1) = particles.yp_.at(ipart);
      if(rank==size-1)
        // Periodic wrap around
        send_parts_right.at(i*partFields+2) = particles.zp_.at(ipart)-Lzl;
      else
        send_parts_right.at(i*partFields+2) = particles.zp_.at(ipart)+Lz0;
      send_parts_right.at(i*partFields+3) = particles.vx_.at(ipart);
      send_parts_right.at(i*partFields+4) = particles.vy_.at(ipart);
      send_parts_right.at(i*partFields+5) = particles.vz_.at(ipart);
      send_parts_right.at(i*partFields+6) = particles.q_.at(ipart);
    }
    for(int i=0; i<to_send_left.size(); ++i) {
      ipart = to_send_left.at(i);
      send_parts_left.at(i*partFields+0) = particles.xp_.at(ipart);
      send_parts_left.at(i*partFields+1) = particles.yp_.at(ipart);
      if(rank==0)
        send_parts_left.at(i*partFields+2) = particles.zp_.at(ipart)+Lz;
      else
        send_parts_left.at(i*partFields+2) = particles.zp_.at(ipart)+Lz0;
      send_parts_left.at(i*partFields+3) = particles.vx_.at(ipart);
      send_parts_left.at(i*partFields+4) = particles.vy_.at(ipart);
      send_parts_left.at(i*partFields+5) = particles.vz_.at(ipart);
      send_parts_left.at(i*partFields+6) = particles.q_.at(ipart);
    }

    int n_send_r, n_send_l;
    int n_recv_r, n_recv_l;
    MPI_Request send_size[2], send_parts[2];

    n_send_r = to_send_right.size();
    n_send_l = to_send_left.size();    
    MPI_Isend(&n_send_r, 1, MPI_INT, right, 0, MPI_COMM_WORLD, &send_size[0]);
    MPI_Isend(&n_send_l, 1, MPI_INT, left, 0, MPI_COMM_WORLD, &send_size[1]);

    MPI_Isend(&send_parts_right[0], n_send_r*partFields, MPI_DOUBLE, right,
              1, MPI_COMM_WORLD, &send_parts[0]);
    MPI_Isend(&send_parts_left[0], n_send_l*partFields, MPI_DOUBLE, left,
              1, MPI_COMM_WORLD, &send_parts[1]);

    MPI_Recv(&n_recv_r, 1, MPI_INT, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&n_recv_l, 1, MPI_INT, left,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    recv_parts_right.resize(std::max<int>(1, n_recv_r*partFields));
    recv_parts_left.resize(std::max<int>(1, n_recv_l*partFields));
    MPI_Recv(&recv_parts_right[0], n_recv_r*partFields, MPI_DOUBLE, right,
             1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&recv_parts_left[0], n_recv_l*partFields, MPI_DOUBLE, left,
             1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Waitall(2, send_size,  MPI_STATUS_IGNORE);
    MPI_Waitall(2, send_parts, MPI_STATUS_IGNORE);

    // Update data structures
    std::vector<int> to_remove = to_send_left;
    for(int i=0; i<to_send_right.size(); ++i)
      to_remove.push_back(to_send_right.at(i));
    particles.RemoveParticles(to_remove);
    
    for(int ipart=0; ipart<particles.size(); ++ipart)
      assert(particles.zp_.at(ipart)>=0 && particles.zp_.at(ipart)<Lzl);

    for(int ipart=0; ipart<n_recv_r; ++ipart) {
      particles.xp_.push_back(recv_parts_right.at(ipart*partFields+0));
      particles.yp_.push_back(recv_parts_right.at(ipart*partFields+1));
      particles.zp_.push_back(recv_parts_right.at(ipart*partFields+2)-Lz0);
      particles.vx_.push_back(recv_parts_right.at(ipart*partFields+3));
      particles.vy_.push_back(recv_parts_right.at(ipart*partFields+4));
      particles.vz_.push_back(recv_parts_right.at(ipart*partFields+5));
      particles.q_.push_back(recv_parts_right.at(ipart*partFields+6));
    }
    for(int ipart=0; ipart<n_recv_l; ++ipart) {
      particles.xp_.push_back(recv_parts_left.at(ipart*partFields+0));
      particles.yp_.push_back(recv_parts_left.at(ipart*partFields+1));
      particles.zp_.push_back(recv_parts_left.at(ipart*partFields+2)-Lz0);
      particles.vx_.push_back(recv_parts_left.at(ipart*partFields+3));
      particles.vy_.push_back(recv_parts_left.at(ipart*partFields+4));
      particles.vz_.push_back(recv_parts_left.at(ipart*partFields+5));
      particles.q_.push_back(recv_parts_left.at(ipart*partFields+6));
    }

    for(int ipart=0; ipart<particles.size(); ++ipart) {
      // Normalize particle positions
      if(particles.xp_.at(ipart)<0)
        particles.xp_.at(ipart) += Lxl;
      if(particles.yp_.at(ipart)<0)
        particles.yp_.at(ipart) += Lyl;
      if(particles.xp_.at(ipart)>=Lxl)
        particles.xp_.at(ipart) -= Lxl;
      if(particles.yp_.at(ipart)>=Lyl)
        particles.yp_.at(ipart) -= Lyl;

      assert(particles.xp_.at(ipart)>=0 && particles.xp_.at(ipart)<Lxl);
      assert(particles.yp_.at(ipart)>=0 && particles.yp_.at(ipart)<Lyl);
      assert(particles.zp_.at(ipart)>=0 && particles.zp_.at(ipart)<Lzl);

    }

    t_end = MPI_Wtime();
    t_comm.at(it) = t_end-t_start;


    // Diagnostics
    {
      // Check total particle counts
      int part_sum, n_parts;
      n_parts = particles.size();
      n_parts_local.at(it) = n_parts;
      MPI_Reduce(&n_parts, &part_sum, 1, MPI_INT, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      if(rank==0)
        assert(part_sum==init_n_parts);

      // Save comm counts
      vn_send_r.at(it) = n_send_r;
      vn_send_l.at(it) = n_send_l;
      vn_recv_r.at(it) = n_recv_r;
      vn_recv_l.at(it) = n_recv_l;

    }

    t_loop_end = MPI_Wtime();
    t_loop.at(it) = t_loop_end-t_loop_start;

    if(rank==0)
      std::cout << it+1 << ' ' << t_loop.at(it) << std::endl;;

  }

  // Save diagnostic results
  if(rank==0) {
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

    for(int irank=0; irank<size; ++irank) {

      if(irank>0) {
        MPI_Recv(&t_loop[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_field_solve[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_weight[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_interp[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_calc_E[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_comm[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_accel[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_move[0], nt, MPI_DOUBLE, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&n_parts_local[0], nt, MPI_INT, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&vn_send_l[0], nt, MPI_INT, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&vn_send_r[0], nt, MPI_INT, irank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      for(int i=0; i<nt; ++i) {
        odiag << irank << ' '
              << i+1 << ' '<< size << ' '
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
              << n_parts_local.at(i) << ' '
              << vn_send_r.at(i)+vn_send_l.at(i) << std::endl;
      }

    }
    odiag.close();
  }
  else {
    MPI_Send(&t_loop[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_field_solve[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_weight[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_interp[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_calc_E[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_comm[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_accel[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&t_move[0], nt, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&n_parts_local[0], nt, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&vn_send_l[0], nt, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&vn_send_r[0], nt, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

}

int main(int argc, char* argv[]) {

  MPI_Init(&argc, &argv);
  fftw_mpi_init();
  
  pic_fft(argc, argv);

  MPI_Finalize();  

}

