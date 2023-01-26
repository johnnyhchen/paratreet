#ifndef PARATREET_FOFDATA_H_
#define PARATREET_FOFDATA_H_

#include "common.h"
#include <vector>
#include <queue>
#include "Particle.h"
#include "OrientedBox.h"
#include "MultipoleMoments.h"

struct pqSmoothNode {
  Real fKey = 0.;// distance^2 -> place in priority queue
  const Particle* pPtr = nullptr;

  inline bool operator<(const pqSmoothNode& n) const {
    return fKey < n.fKey;
  }

  void pup(PUP::er& p) {
    p|fKey;
  }
};

struct FoFData {
  Vector3D<Real> moment;
  Real sum_mass;
  Vector3D<Real> centroid; // too slow to compute this on the fly
  OrientedBox<Real> box;
  int count;

  FoFData() :
  moment(Vector3D<Real> (0,0,0)), sum_mass(0), count(0) {}

  /// Construct centroid from particles.
  FoFData(const Particle* particles, int n_particles, int depth) : FoFData() {
    for (int i = 0; i < n_particles; i++) {
      moment += particles[i].mass * particles[i].position;
      sum_mass += particles[i].mass;
      box.grow(particles[i].position);
    }
    centroid = moment / sum_mass;
    count = n_particles;
  }

  const FoFData& operator+=(const FoFData& cd) { // needed for upward traversal
    moment += cd.moment;
    sum_mass += cd.sum_mass;
    centroid = moment / sum_mass;
    box.grow(cd.box);
    count += cd.count;
    return *this;
  }

  FoFData& operator=(const FoFData&) = default;

  void pup(PUP::er& p) {
    p|moment;
    p|sum_mass;
    p|centroid;
    p|box;
    p|count;
  }

};

#endif // PARATREET_SEARCH_H_
