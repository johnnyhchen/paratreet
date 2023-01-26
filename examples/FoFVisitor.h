#ifndef PARATREET_FOFVISITOR_H_ // does this macro need to be defined in some configeration file as well? J.C.
#define PARATREET_FOFVISITOR_H_

#include "paratreet.decl.h"
#include "common.h"
#include "Space.h"
#include <cmath>
#include <vector>
#include <queue>

extern Real max_timestep;

struct FoFVisitor {
public:
  // static constexpr const bool CallSelfLeaf = true;
  void pup(PUP::er& p) {}

public:
  FoFVisitor() { // Had a vector param here initially but don't think we need it like in GravityVisitor.h
    // TODO: Write constructor?
  }

  // TODO: update with "FoFData" if we decide we need a more parred down version
  bool open(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    Real r_bucket = target.data.size_sm + 0.2; // TODO: verify instead of max_rad, have 0.2 for FoF red sphere len
    if (!Space::intersect(source.data.box, target.data.box.center(), r_bucket*r_bucket))
      return false;

    // Check if any of the target balls intersect the source volume
    for (int i = 0; i < target.n_particles; i++) {
      Real ballSq = 0.2 * 0.2; // TODO: make global variable for linking length
      if(Space::intersect(source.data.box, target.particles()[i].position, ballSq))
        return true;
    }
    return false;
  }

  void node(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {

  }

  void leaf(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    for (int i = 0; i < target.n_particles; i++) {
      for (int j = 0; j < source.n_particles; j++) {
        auto& sp = source.particles()[j];
        auto& tp = target.particles()[i];
        Real distance = (tp.position - sp.position).length();
        Real linkingLength = 0.2;
        // if (source and target ptcl linking length spheres intersect, and order not equal)
        // using less than for order comparison to avoid adding same pair twice
        if (distance < linkingLength && sp.order < tp.order) {
          // TODO: Call UnionFind here and determine if particles in this leaf are part of other particle's halo
          
        }
      }
    }
  }
};

#endif // PARATREET_FOFVISITOR_H_
