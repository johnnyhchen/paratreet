#include <iostream>
#include "Paratreet.h"
#include "CentroidData.h"
#include "unionFindLib.h"
#include "FoFVisitor.h"
#include "FoF.decl.h"

/* readonly */ CProxy_UnionFindLib libProxy;
/* readonly */ CProxy_Partition<CentroidData> partitionProxy;
/* readonly */ Real max_timestep;

using namespace paratreet;

static void initialize() {
  BoundingBox::registerReducer();
}

class FoF : public paratreet::Main<CentroidData> { 
  public:
  void main(CkArgMsg* m) override {    // Initialize readonly variables
    max_timestep = 1e-5;

    
    // Process command line arguments
    int c;
    std::string input_str;

    while ((c = getopt(m->argc, m->argv, "mec:j:")) != -1) {
      switch (c) {
        case 'j':
          max_timestep = atof(optarg);
          break;
      }
    }
    //main::initializeDriver() will be run after main exits. After that main::run() is ran. See Paratreet.C::MainChare class
  }


  void setDefaults(void) {
    // TODO: adjust values as needed
    conf.min_n_subtrees = CkNumPes() * 8; // default from ChaNGa
    conf.min_n_partitions = CkNumPes() * 8;
    conf.max_particles_per_leaf = 12; // default from ChaNGa
    conf.decomp_type = paratreet::DecompType::eBinaryOct;
    conf.tree_type = paratreet::TreeType::eBinaryOct;
    conf.num_iterations = 1; // for FoF, only 1 iteration
    conf.num_share_nodes = 0; // 3;
    conf.cache_share_depth = 3;
    conf.pool_elem_size;
    conf.flush_period = 0;
    conf.flush_max_avg_ratio = 10.;
    conf.lb_period = 5;
    conf.request_pause_interval = 20;
    conf.iter_pause_interval = 1000;
  }

  // -------------------
  // Traversal functions
  // -------------------
  void preTraversalFn(ProxyPack<CentroidData>& proxy_pack) override {
    // tells the Driver to load the Cache Manager with a starter pack of data, specified in Configuration.cache_share_depth
    proxy_pack.driver.loadCache(CkCallbackResumeThread());
    
    // store proxies for use in FoFVisitor
    libProxy = proxy_pack.libProxy;
    partitionProxy = proxy_pack.partition;
  }

  void traversalFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) override {
    proxy_pack.partition.template startDown<FoFVisitor>(FoFVisitor());
  }

  //template <typename Data>
  void postIterationFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) override {
    // output results from UnionFind for halos somehow
    // what data format should we output?

    // CkCallback cb(CkIndex_Main::done(), mainProxy);
    partitionProxy.ckGetArrayID();
    libProxy.find_components(CkCallbackResumeThread());
    // print libProxy results: to access results, iterate through libVertices (externed above, remove if we don't access in this file)
    // TODO: figure out output format: what format to print to and what API (see Writer.h)
    partitionProxy.getConnectedComponents(CkCallbackResumeThread());
    paratreet::outputParticleAccelerations(universe, partitionProxy);
  }

  // TODO: Figure out what prune components means
  /*
  void doneFindComponents() {
      CkPrintf("[Main] Components identified, prune unecessary ones now\n");
      CkPrintf("[Main] Components detection time: %f\n", CkWallTimer()- start_time);
      // callback for library to report to after pruning
      CkCallback cb(CkIndex_TreePiece::requestVertices(), partitions); // TODO: define requestVertices() function in Partitions.h
      libProxy.prune_components(1, cb);  // can imagine will prune_components that i.e. have only 1 element. also reports # of components in parent tree (see UnionFindLib.C)
      // requestVertices() called when prune_components done
  }

  void donePrinting() {
    // print output here
    partitionProxy.getConnectedComponents();
  }
  */
  Real getTimestep(BoundingBox& universe, Real max_velocity) {
    Real universe_box_len = universe.box.greater_corner.x - universe.box.lesser_corner.x;
    Real temp = universe_box_len / max_velocity / std::cbrt(universe.n_particles);
    return std::min(temp, max_timestep);
  }


  void run() {
    driver.run(CkCallbackResumeThread());
    CkExit();
  }
};

PARATREET_REGISTER_MAIN(FoF);
#include "templates.h"
#include "FoF.def.h"
