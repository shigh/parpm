
#include "particle.hpp"

ParticleList::ParticleList(int n_resv) {

  xp_.reserve(n_resv);
  yp_.reserve(n_resv);
  zp_.reserve(n_resv);
  vx_.reserve(n_resv);
  vy_.reserve(n_resv);
  vz_.reserve(n_resv);
  q_.reserve(n_resv);

}

void ParticleList::RemoveParticles(std::vector<int> particles) {

  // Remove particles in reverse order to avoid index problems
  std::sort(particles.begin(), particles.end());
  int ipart;
  int last = xp_.size()-1;  
  for(int i=particles.size()-1; i>=0; --i) {
    ipart = particles.at(i);
    assert(ipart>=0);
    assert(ipart<size());
    assert(last>=0);
    if(ipart!=last) {
      xp_.at(ipart) = xp_.at(last);
      yp_.at(ipart) = yp_.at(last);
      zp_.at(ipart) = zp_.at(last);
      vx_.at(ipart) = vx_.at(last);
      vy_.at(ipart) = vy_.at(last);
      vz_.at(ipart) = vz_.at(last);
      q_.at(ipart)  = q_.at(last);
    }
    --last;
  }

  xp_.resize(last+1);
  yp_.resize(last+1);
  zp_.resize(last+1);
  vx_.resize(last+1);
  vy_.resize(last+1);
  vz_.resize(last+1);
  q_.resize(last+1);

}
