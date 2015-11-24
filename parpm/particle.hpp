
#pragma once

#include "types.hpp"

struct Particle {

  FLOAT x;
  FLOAT y;
  FLOAT z;

  FLOAT vx;
  FLOAT vy;
  FLOAT vz;

  FLOAT q;
  
};


class ParticleList {

private:

  std::vector<FLOAT> xp_;
  std::vector<FLOAT> yp_;
  std::vector<FLOAT> zp_;

  std::vector<FLOAT> vx_;
  std::vector<FLOAT> vy_;
  std::vector<FLOAT> vz_;

  std::vector<FLOAT> q_;

  int n_parts_;

public:

  Particle get_particle(int i) const;
  void AddParticle(Particle particle);
  void AddParticle(Particle* particles, int n);
  void RemoveParticle(int i);
  void RemoveParticle(std::vector<int> particles);

};
