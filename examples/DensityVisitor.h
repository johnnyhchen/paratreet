#ifndef PARATREET_DENSITYVISITOR_H_
#define PARATREET_DENSITYVISITOR_H_

#include "paratreet.decl.h"
#include "common.h"
#include "Space.h"
#include <cmath>
#include <vector>
#include <queue>

struct DensityVisitor {
public:
  static constexpr const bool CallSelfLeaf = true;

  void pup(PUP::er& p) {}

// in leaf check for not same particle plz
private:
  static constexpr const int k = 32;

public:
  static bool open(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    // Check if any of the target balls intersect the source volume
    for (int i = 0; i < target.n_particles; i++) {
      if (target.data.pps[i].neighbors.size() < k) return true;
      if(Space::intersect(source.data.box, target.particles()[i].position, target.data.pps[i].neighbors[0].fKey))
        return true;
    }
    return false;
  }

  static void node(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {}

  static void leaf(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    for (int i = 0; i < target.n_particles; i++) {
      auto& Q = target.data.pps[i].neighbors;
      for (int j = 0; j < source.n_particles; j++) {
        const auto& sp = source.particles()[j]; //source particle
        Vector3D<Real> dr = target.particles()[i].position - sp.position;
        auto dsq = dr.lengthSquared();
        // Remove the most distant neighbor if this one is closer and the list is full
        if (Q.size() == k) {
          if (dsq < Q[0].fKey) {
            std::pop_heap(&(Q[0]) + 0, &(Q)[0] + k);
            Q.resize(k-1);
          }
        }
        // Add the particle to the neighbor list if it isnt filled up
        if (Q.size() < k) {
          pqSmoothNode pqNew;
          pqNew.pPtr = &sp;
          pqNew.fKey = dsq;
          Q.push_back(pqNew);
          std::push_heap(&(Q)[0] + 0, &(Q)[0] + Q.size());
        }
      }
    }
  }
};

#endif // PARATREET_DENSITYVISITOR_H_
