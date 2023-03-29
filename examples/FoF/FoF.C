#include <iostream>
#include "Paratreet.h"
#include "CentroidData.h"
#include "unionFindLib.h"
#include "FoFVisitor.h"
#include "FoF.decl.h"

/* readonly */ bool outputFileConfigured;
/* readonly */ CProxy_UnionFindLib libProxy;
/* readonly */ CProxy_Partition<CentroidData> partitionProxy;
/* readonly */ Real max_timestep;
/* readonly */ int peanoKey;
/* readonly */ Real LINKING_LENGTH;

using namespace paratreet;

static void initialize() {
  BoundingBox::registerReducer();
}

class FoF : public paratreet::Main<CentroidData> { 
  void main(CkArgMsg* m) override {
    // Initialize readonly variables
    outputFileConfigured = !conf.output_file.empty();
    peanoKey = 3;
    max_timestep = 1e-5;
    LINKING_LENGTH = conf.linking_length;

    // Process command line arguments
    int c;
    std::string input_str;

    while ((c = getopt(m->argc, m->argv, "mec:j:")) != -1) {
      switch (c) {
        case 'm':
          peanoKey = 0; // morton
          break;
        case 'j':
          max_timestep = atof(optarg);
          break;
        
        default:
          CkPrintf("Usage: %s\n", m->argv[0]);
          CkPrintf("\t-f [input file]\n");
          CkPrintf("\t-n [number of treepieces]\n");
          CkPrintf("\t-p [maximum number of particles per treepiece]\n");
          CkPrintf("\t-l [maximum number of particles per leaf]\n");
          CkPrintf("\t-d [decomposition type: oct, sfc, kd]\n");
          CkPrintf("\t-t [tree type: oct, bin, kd]\n");
          CkPrintf("\t-i [number of iterations]\n");
          CkPrintf("\t-s [number of shared tree levels]\n");
          CkPrintf("\t-u [flush period]\n");
          CkPrintf("\t-r [flush threshold for Subtree max_average ratio]\n");
          CkPrintf("\t-b [load balancing period]\n");
          CkPrintf("\t-v [filename prefix]\n");
          CkPrintf("\t-j [max timestep]\n");
          CkPrintf("\t-ll [linking length]\n");
          CkExit();
      }
    }
    delete m;

    // Print configuration
    CkPrintf("\n[PARATREET]\n");
    if (conf.input_file.empty()) CkAbort("Input file unspecified");
    CkPrintf("Input file: %s\n", conf.input_file.c_str());
    CkPrintf("Decomposition type: %s\n", paratreet::asString(conf.decomp_type).c_str());
    CkPrintf("Tree type: %s\n", paratreet::asString(conf.tree_type).c_str());
    CkPrintf("Minimum number of subtrees: %d\n", conf.min_n_subtrees);
    CkPrintf("Minimum number of partitions: %d\n", conf.min_n_partitions);
    CkPrintf("Maximum number of particles per leaf: %d\n", conf.max_particles_per_leaf);
    CkPrintf("Linking length for friends-of-friends: %f\n", conf.linking_length);
    
    //main::initializeDriver() will be run after main exits. After that main::run() is ran. See Paratreet.C::MainChare class
  }


  void setDefaults(void) {
    conf.min_n_subtrees = CkNumPes() * 8; // default from ChaNGa
    conf.min_n_partitions = CkNumPes() * 8;
    conf.max_particles_per_leaf = 12; // default from ChaNGa
    conf.decomp_type = paratreet::DecompType::eBinaryOct;
    conf.tree_type = paratreet::TreeType::eBinaryOct;
    conf.num_iterations = 1;
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
    // The size of the starter pack of data loaded by the cache manager is specified in Configuration.cache_share_depth
    proxy_pack.driver.loadCache(CkCallbackResumeThread());
    
    libProxy = proxy_pack.libProxy;
    partitionProxy = proxy_pack.partition;

  }

  void traversalFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) override {
    proxy_pack.partition.template startDown<FoFVisitor>(FoFVisitor());
  }

  void postIterationFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) override {
    CkPrintf("[Main] Inverted trees constructed for unionFindLib. Performing components detection\n");
    int startTime = CkWallTimer();
    libProxy.find_components(CkCallbackResumeThread());

    CkPrintf("[Main] Components identified, prune unecessary ones now\n");
    CkPrintf("[Main] Components detection time: %f\n", CkWallTimer()- startTime);
    int minVerticesPerComponent = 1; // strictly greater than
    libProxy.prune_components(minVerticesPerComponent, CkCallbackResumeThread());

    partitionProxy.getConnectedComponents(CkCallbackResumeThread());
    
    //if (iter == 0 && outputFileConfigured) {
    paratreet::outputParticleAccelerations(universe, partitionProxy);
    //}
  }

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
