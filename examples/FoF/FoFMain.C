#include "Main.h"
#include "Paratreet.h"
#include <iostream>
#include "../unionfind/unionFindLib.h"
#include "graph.decl.h"
#include "graph-io.h"
#include "FoFVisitor.h"

extern bool verify;
extern Real max_timestep;
extern unionFindVertex *libVertices;

/*readonly*/ CProxy_Main mainProxy; // TODO: how to populate this field so Partitions::initializeVertices() can access this?
/*readonly*/ CProxy_UnionFindLib libProxy;
/*readonly*/ CProxy_Partition<CentroidData> partitionProxy;
/*readonly*/ UnionFindLib *libPtr;
/*readonly*/ int NUM_VERTICES;
/*readonly*/ int NUM_EDGES;
/*readonly*/ int NUM_TREEPIECES;
/*readonly*/ long int lastChareBegin;

  using namespace paratreet;

  /*
  class FoFMain : public ExMain {
    
  };*/

  // Question: since readonly variables are supposed to not be modifed after we exit main, and we are depending on Paratreet main() as our main(). how can we define UnionFindProxy as a global variable?
  // - are we using Paratreet main()? If so, what is the Main() in the example programs?

  void ExMain::preTraversalFn(ProxyPack<CentroidData>& proxy_pack) {
    // tells the Driver to load the Cache Manager with a starter pack of data, specified in Configuration.cache_share_depth
    proxy_pack.driver.loadCache(CkCallbackResumeThread()); // TODO: need to make config file and specify Configuration.cache_share_depth
  }

  void ExMain::traversalFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) {
    // start a traversal
    NUM_VERTICES = universe.n_particles; // TODO: verify this works
    
    // putting below init section in traversalFn because we need the number of particles (vertices) to do initializeLibVertices

    // mainProxy = this->thisProxy; // TODO: does this work?
    partitionProxy = proxy_pack.partition;
    // TODO: figure out if the partition proxy and subtree proxy are bound arrays (so we can use the thisIndex trick)
    // also find out how to get `thisIndex` variable in scope

    // Initialize UnionFindLib
    libProxy = UnionFindLib::unionFindInit(partitionProxy, proxy_pack.subtree[0].ckLocal()->n_partitions); // can also use proxy_pack.partition[0].ckLocal()->n_partitions;

    // Add vertices
    //libProxy[0].register_phase_one_cb(CkCallbackResumeThread()); // TODO: do we need this since we do nothing in the callback?

    // how did Graph.C define a initializeLibVertices() function in TreePiece class, but was able to call the function on the tpProxy?
    // It think charm++ will autogenerate the member function for a proxy in the FoFMain.def.h file once we compile it...?
    partitionProxy.initializeLibVertices(CkCallbackResumeThread());

    // how to get array of all particles? do this in main?
    // what passes the traversalFn the universe field? no results so far...
    // (*libPtr).initialize_vertices( , NUM_VERTICES);// TODO how to get other field (all particles in system)?
    // libProxy.initialize_vertices(proxy_pack.cache, NUM_VERTICES);

    proxy_pack.partition.template startDown<FoFVisitor>(FoFVisitor());
  }

  void ExMain::postIterationFn(BoundingBox& universe, ProxyPack<CentroidData>& proxy_pack, int iter) {
    // output results from UnionFind for halos somehow
    // what data format should we output?

    // CkCallback cb(CkIndex_Main::done(), mainProxy);
    partitionProxy.ckGetArrayID();
    libProxy.find_components(CkCallbackResumeThread());
    // print libProxy results: to access results, iterate through libVertices (externed above, remove if we don't access in this file)
    // TODO: figure out output format: what format to print to and what API (see Writer.h)
    partitionProxy.getConnectedComponents();
  }

  /*
  void ExMain::done() {
    // does nothing
  }
  */


  void doneFindComponents() {
      CkPrintf("[Main] Components identified, prune unecessary ones now\n");
      CkPrintf("[Main] Components detection time: %f\n", CkWallTimer()-startTime);
      // callback for library to report to after pruning
      CkCallback cb(CkIndex_TreePiece::requestVertices(), tpProxy); // TODO: define requestVertices() function in Partitions.h
      libProxy.prune_components(1, cb);  // can imagine will prune_components that i.e. have only 1 element. also reports # of components in parent tree (see UnionFindLib.C)
      // requestVertices() called when prune_components done
  }

  void donePrinting() {
    // print output here
    partitionProxy.getConnectedComponents();
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