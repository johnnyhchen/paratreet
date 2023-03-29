#ifndef PARATREET_FOFVISITOR_H_ // does this macro need to be defined in some configeration file as well? J.C.
#define PARATREET_FOFVISITOR_H_

#include "paratreet.decl.h"
#include "common.h"
#include "Space.h"
#include "CentroidData.h"
#include <cmath>
#include <vector>
#include <queue>
#include "unionFindLib.h"
#include "Partition.h"
// #include "FoF.h"

extern CProxy_UnionFindLib libProxy;
extern CProxy_Partition<CentroidData> partitionProxy;
extern Real LINKING_LENGTH;

struct FoFVisitor {
public:
  static constexpr const bool CallSelfLeaf = true;
  void pup(PUP::er& p) {}

public:
  FoFVisitor() {
  }

  bool open(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    Real r_bucket = target.data.size_sm + LINKING_LENGTH;
    if (!Space::intersect(source.data.box, target.data.box.center(), r_bucket*r_bucket))
      return false;

    // Check if any of the target balls intersect the source volume
    for (int i = 0; i < target.n_particles; i++) {
      Real ballSq = LINKING_LENGTH * LINKING_LENGTH;
      if(Space::intersect(source.data.box, target.particles()[i].position, ballSq))
        return true;
    }
    return false;
  }

  void node(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {}

  void leaf(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    for (int i = 0; i < target.n_particles; i++) {
      for (int j = 0; j < source.n_particles; j++) {
        const Particle& sp = source.particles()[j];
        const Particle& tp = target.particles()[i];
        Real distance = (tp.position - sp.position).length();
        // union two particles if source and target particles linking length spheres intersect. avoid union of same pair twice by comapring particle orders (ID)
        if (distance < LINKING_LENGTH && sp.order < tp.order) {
          CkPrintf("union_requst on vid1=%ld, vid2=%ld\n", sp.order, tp.order); // TODO: remove debugging printf
          libProxy[tp.partition_idx].ckLocal()->union_request(sp.vertex_id, tp.vertex_id);
        }
      }
    }
  }
};

#endif // PARATREET_FOFVISITOR_H_
