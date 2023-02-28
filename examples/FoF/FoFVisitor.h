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
extern int leafNum;
extern int openNum;

struct FoFVisitor {
private:
// TODO: make this variable command line arg (so can try options when running)
  static constexpr const Real linkingLength = 0.00417; // 0.2; // can't do static constexpr const Real like in CollisionVisitor.h

public:
  static constexpr const bool CallSelfLeaf = true;
  void pup(PUP::er& p) {}

public:
  FoFVisitor() { // Had a vector param here initially but don't think we need it like in GravityVisitor.h
    // TODO: Write constructor?
  }

  // TODO: update with "FoFData" if we decide we need a more parred down version
  bool open(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    CkPrintf("open %d\n", openNum+=1);
    Real r_bucket = target.data.size_sm + linkingLength; // TODO: verify instead of max_rad, have 0.2 for FoF red sphere len
    if (!Space::intersect(source.data.box, target.data.box.center(), r_bucket*r_bucket))
      return false;

    // Check if any of the target balls intersect the source volume
    for (int i = 0; i < target.n_particles; i++) {
      Real ballSq = linkingLength * linkingLength;
      if(Space::intersect(source.data.box, target.particles()[i].position, ballSq))
        return true;
    }
    return false;
  }

  void node(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {

  }

  void leaf(const SpatialNode<CentroidData>& source, SpatialNode<CentroidData>& target) {
    CkPrintf("leaf %d\n", leafNum+=1);
    for (int i = 0; i < target.n_particles; i++) {
      for (int j = 0; j < source.n_particles; j++) {
        const Particle& sp = source.particles()[j];
        const Particle& tp = target.particles()[i];
        Real distance = (tp.position - sp.position).length();
        // if (source and target ptcl linking length spheres intersect, and order not equal)
        // using less than for order comparison to avoid adding same pair twice
        if (distance < linkingLength && sp.order < tp.order) {
          // TODO: Call UnionFind here and determine if particles in this leaf are part of other particle's halo
          // Problem: FoFVisitor is not of type partition so does not have thisIndex
          // Solution?: Call partition proxy? nope, this is a broadcast (sends to all elements) need to send to ONE element with the particles
          // possible solution: store partition index in SpatialNode struct so it knows what partition it sits on (requires changing init for SpatialNodes or Partitions?)
          
          partitionProxy[tp.partition_idx].unionRequest(sp.order, tp.order);
        }
      }
    }
  }
};

#endif // PARATREET_FOFVISITOR_H_
