#include "Main.h"
#include "Paratreet.h"
#include <iostream>
#include "../unionfind/unionFindLib.h"
#include "graph.decl.h"
#include "graph-io.h"
#include "FoFVisitor.h"

extern bool verify;
extern Real max_timestep;

/*readonly*/ CProxy_UnionFindLib libProxy;
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int NUM_VERTICES;
/*readonly*/ int NUM_EDGES;
/*readonly*/ int NUM_TREEPIECES;
/*readonly*/ long int lastChareBegin;

  using namespace paratreet;

  // following Collision.C example

  // Do we need to override main() as well?

  void ExMain::preTraversalFn(ProxyPack<FoFData>& proxy_pack) {
    // tells the Driver to load the Cache Manager with a starter pack of data, specified in Configuration.cache_share_depth
    proxy_pack.driver.loadCache(CkCallbackResumeThread()); // TODO: need to make config file and specify Configuration.cache_share_depth
  }

  void ExMain::traversalFn(BoundingBox& universe, ProxyPack<FoFData>& proxy_pack, int iter) {
    // start a traversal
    proxy_pack.partition.template startDown<FoFVisitor>(FoFVisitor()); // TODO: What is theta? Also add params to FoFVisitor
  }

  void ExMain::postIterationFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) {
    // output results from UnionFind for halos somehow
    // what data format should we output?
  }

  // TODO: Do we need this?
  Real ExMain::getTimestep(BoundingBox& universe, Real max_velocity) {
    Real universe_box_len = universe.box.greater_corner.x - universe.box.lesser_corner.x;
    Real temp = universe_box_len / max_velocity / std::cbrt(universe.n_particles);
    return std::min(temp, max_timestep);
  }

  /*
  void run(CkCallback cb) {
    // not sure where to set config
    auto& config = getConfiguration();
    config.num_iterations = 1;  // analyzing a snapshot for FoF, no timesteps needed
  }
  */