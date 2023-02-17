#include "Main.h"
#include "Paratreet.h"
#include <iostream>
#include "../unionfind/unionFindLib.h"
#include "FoFVisitor.h"

extern bool verify;
extern Real max_timestep;
extern unionFindVertex *libVertices;

// IS it ok for us to declare a global variable here, and then define it in another file that extern's it?
/*readonly*/ CProxy_UnionFindLib libProxy;
/*readonly*/ CProxy_Partition<CentroidData> partitionProxy;

  using namespace paratreet;

  template <typename Data>
  class FoF : public ExMain {
    public:
    BoundingBox universe;
    int n_partitions;
    // additional fields
    int n_subtrees;
    std::vector<int> partition_locations;
    CProxy_Partition<Data> partitions;
    CProxy_CacheManager<Data> cache_manager;
    CProxy_Subtree<Data> subtrees;
    CProxy_Resumer<Data> resumer;
    CProxy_TreeCanopy<Data> canopy;  // called calculator in Driver.h code
    double start_time;  // for benchmarking

    void main(CkArgMsg* msg) override {
      /*
      // I think we might not need to set the config since we are not using the iteration code in Driver.h
      // set config: todo check driver code, there's a section about dropping const-ness somewhere which we need to do here
      auto& config = getConfiguration();
      config.num_iterations = 1;  // analyzing a snapshot for FoF, no timesteps needed
      */

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
      proxy_pack.driver.loadCache(CkCallbackResumeThread()); // TODO: need to make config file and specify Configuration.cache_share_depth
    }

    void traversalFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) override {
      proxy_pack.partition.template startDown<FoFVisitor>(FoFVisitor());
    }


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

    // TODO: Do we need this?
    Real getTimestep(BoundingBox& universe, Real max_velocity) {
      Real universe_box_len = universe.box.greater_corner.x - universe.box.lesser_corner.x;
      Real temp = universe_box_len / max_velocity / std::cbrt(universe.n_particles);
      return std::min(temp, max_timestep);
    }*/
  };