#include "Main.h"
#include "Paratreet.h"
#include <iostream>
#include "../unionfind/unionFindLib.h"
#include "FoFVisitor.h"

extern bool verify;
extern Real max_timestep;
extern unionFindVertex *libVertices;

/*readonly*/ CProxy_UnionFindLib libProxy;
/*readonly*/ CProxy_Partition<CentroidData> partitionProxy;
/*readonly*/ UnionFindLib *libPtr;
/*readonly*/ int NUM_VERTICES;
/*readonly*/ int NUM_EDGES;
/*readonly*/ int NUM_TREEPIECES;
/*readonly*/ long int lastChareBegin;

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

    void initializeDriver(const CkCallback& cb) override {
      initialize(cb);
    }

    void initialize(const CkCallback& cb) {
        // Create readers
        n_readers = CkNumPes();
        readers = CProxy_Reader::ckNew();
        treespec = CProxy_TreeSpec::ckNew();
        thread_state_holder = CProxy_ThreadStateHolder::ckNew();

        // Create library chares
        // TODO: do we need to set these as global variables somewhere?
        canopy = CProxy_TreeCanopy<Data>::ckNew();
        canopy.doneInserting();
        cache_manager= CProxy_CacheManager<Data>::ckNew();
        resumer = CProxy_Resumer<Data>::ckNew();
        /*
        // commenting out for no driver code
        CProxy_Driver<Data> driver = CProxy_Driver<Data>::ckNew(cache, resumer, canopy, CkMyPe());
        // Call the driver initialization routine (performs decomposition)
        auto& cfg = const_cast<paratreet::Configuration&>(paratreet::getConfiguration());
        driver.init(cb, CkReference<Configuration>(cfg));
        
        return driver;
        */

        auto& cfg = const_cast<paratreet::Configuration&>(paratreet::getConfiguration());
        init(cb, CkReference<Configuration>(cfg)); // once init() is done it will call run() as its callback
    }


    // Performs initial decomposition
    void init(const CkCallback& cb, const paratreet::Configuration& cfg) {
      // Ensure all treespecs have been created
      CkPrintf("* Validating tree specifications.\n");
      treespec.receiveConfiguration(
        CkCallbackResumeThread(),
        const_cast<paratreet::Configuration*>(&cfg)
      );
      // Then, initialize the cache managers
      CkPrintf("* Initializing cache managers.\n");
      cache_manager.initialize(CkCallbackResumeThread());
      // Useful particle keys
      CkPrintf("* Initialization\n");
      decompose(0);
      cb.send();
    }


    // Performs decomposition by distributing particles among Subtrees,
    // by either loading particle information from input file or re-computing
    // the universal bounding box
    void decompose(int iter) {
      auto& config = paratreet::getConfiguration();
      double decomp_time = CkWallTimer();
      // Build universe
      double start_time = CkWallTimer();
      CkReductionMsg* result;
      readers.load(config.input_file, CkCallbackResumeThread((void*&)result));
      CkPrintf("Loading Tipsy data and building universe: %.3lf ms\n",
          (CkWallTimer() - start_time) * 1000);
      
      if(config.origin_of("dSoft") != paratreet::FieldOrigin::Unknown) {
          CkPrintf("Setting softening to %f \n", config.dSoft);
          // Softening is specified: set it for all particles.
          readers.setSoft(config.dSoft, CkCallbackResumeThread());
      }
      universe = *((BoundingBox*)result->getData());
      delete result;
      remakeUniverse();
      if (config.min_n_subtrees < CkNumPes() || config.min_n_partitions < CkNumPes()) {
        CkPrintf("WARNING: Consider increasing min_n_subtrees and min_n_partitions to at least #pes\n");
      }
      // Assign keys and sort particles locally
      start_time = CkWallTimer();
      readers.assignKeys(universe, CkCallbackResumeThread());
      CkPrintf("Assigning keys and sorting particles: %.3lf ms\n",
        (CkWallTimer() - start_time) * 1000);

      bool matching_decomps = config.decomp_type == paratreet::subtreeDecompForTree(config.tree_type);
      // Set up splitters for decomposition
      start_time = CkWallTimer();
      n_partitions = treespec.ckLocalBranch()->getPartitionDecomposition()->findSplitters(universe, readers, config.min_n_partitions);
      partition_locations.resize(n_partitions);
      treespec.receiveDecomposition(CkCallbackResumeThread(),
          CkPointer<Decomposition>(treespec.ckLocalBranch()->getPartitionDecomposition()), false);
      CkPrintf("Setting up splitters for particle decompositions: %.3lf ms\n",
          (CkWallTimer() - start_time) * 1000);

      // Create Partitions
      CkArrayOptions partition_opts(n_partitions);
      treespec.ckLocalBranch()->getPartitionDecomposition()->setArrayOpts(partition_opts, {}, false);
      partitions = CProxy_Partition<Data>::ckNew(
        n_partitions, cache_manager, resumer, canopy,
        this->thisProxy, matching_decomps, partition_opts
        );
      CkPrintf("Created %d Partitions: %.3lf ms\n", n_partitions,
          (CkWallTimer() - start_time) * 1000);

      start_time = CkWallTimer();
      readers.assignPartitions(n_partitions, partitions);
      CkStartQD(CkCallbackResumeThread());
      CkPrintf("Assigning particles to Partitions: %.3lf ms\n",
          (CkWallTimer() - start_time) * 1000);

      start_time = CkWallTimer();
      if (matching_decomps) {
        n_subtrees = n_partitions;
        CkPrintf("Using same decomposition for subtrees and partitions\n");
        treespec.receiveDecomposition(CkCallbackResumeThread(),
          CkPointer<Decomposition>(treespec.ckLocalBranch()->getPartitionDecomposition()), true);
      }
      else {
        n_subtrees = treespec.ckLocalBranch()->getSubtreeDecomposition()->findSplitters(universe, readers, config.min_n_subtrees);
        treespec.receiveDecomposition(CkCallbackResumeThread(),
          CkPointer<Decomposition>(treespec.ckLocalBranch()->getSubtreeDecomposition()), true);
        CkPrintf("Setting up splitters for subtree decompositions: %.3lf ms\n",
            (CkWallTimer() - start_time) * 1000);
      }

      // Create Subtrees
      start_time = CkWallTimer();
      CkArrayOptions subtree_opts(n_subtrees);
      if (matching_decomps) subtree_opts.bindTo(partitions);
      treespec.ckLocalBranch()->getSubtreeDecomposition()->setArrayOpts(subtree_opts, partition_locations, !matching_decomps);
      subtrees = CProxy_Subtree<Data>::ckNew(
        CkCallbackResumeThread(),
        universe.n_particles, n_subtrees, n_partitions,
        canopy, resumer,
        cache_manager, this->thisProxy, matching_decomps, subtree_opts
        );
      CkPrintf("Created %d Subtrees: %.3lf ms\n", n_subtrees,
          (CkWallTimer() - start_time) * 1000);

      start_time = CkWallTimer();
      readers.flush(n_subtrees, subtrees);
      CkStartQD(CkCallbackResumeThread());
      CkPrintf("Flushing particles to Subtrees: %.3lf ms\n",
          (CkWallTimer() - start_time) * 1000);
      CkPrintf("**Total Decomposition time: %.3lf ms\n",
          (CkWallTimer() - decomp_time) * 1000);

      // Initialize unionFindLib

      NUM_VERTICES = universe.n_particles;
      partitionProxy = partitions; // TODO: delete partitions local variable, just assign to partitionProxy global var

      libProxy = UnionFindLib::unionFindInit(partitionProxy, subtrees[0].ckLocal()->n_partitions); // can also use proxy_pack.partition[0].ckLocal()->n_partitions;

      // Add vertices
      //libProxy[0].register_phase_one_cb(CkCallbackResumeThread()); // TODO: do we need this since we do nothing in the callback?
      partitionProxy.initializeLibVertices(CkCallbackResumeThread());
      CkPrintf("Initialized %d vertices in UnionFindLib\n", NUM_VERTICES);
    }

    void remakeUniverse() {
      Vector3D<Real> bsize = universe.box.size();
      Real max = (bsize.x > bsize.y) ? bsize.x : bsize.y;
      max = (max > bsize.z) ? max : bsize.z;
      Vector3D<Real> bcenter = universe.box.center();
      // The magic number below is approximately 2^(-19)
      const Real fEps = 1.0 + 1.91e-6;  // slop to ensure keys fall between 0 and 1.
      bsize = Vector3D<Real>(fEps*0.5*max);
      universe.box = OrientedBox<Real>(bcenter-bsize, bcenter+bsize);
      thread_state_holder.setUniverse(universe);

      std::cout << "Universal bounding box: " << universe << " with volume "
        << universe.box.volume() << std::endl;
    }

    void run() override {
      // TODO: where are the fields n_partitions, partition_proxy, and n_particles (and universe with particle ids) populated? Need to init UnionFind AFTER that
      // - check decompose() then run()
      // n_partitions and partitions set in decompose()
      // n_particles set in


      auto& config = paratreet::getConfiguration(); // The & operator ensures we store a reference to the object and not a copy
      double total_time = 0;

      // below is iteration code orinally in for loop of Driver.h
      int iter = 0; // TODO: figure out what logic below is doing with iter variable (always 0 for our case) and see if we can remove those if/ternary statements
      
      double iter_start_time = CkWallTimer();
      // Start tree build in Subtrees
      start_time = CkWallTimer();
      CkCallback timeCb (CkReductionTarget(Driver<Data>, reportTime), this->thisProxy);
      subtrees.buildTree(partitions, timeCb);
      CkWaitQD();
      CkPrintf("Tree build and sending leaves: %.3lf ms\n", (CkWallTimer() - start_time) * 1000);

      // Meta data collections, first for max velo
      CkReductionMsg * msg, *msg2;
      subtrees.collectMetaData(CkCallbackResumeThread((void *&) msg));
      // Parse Subtree reduction message
      int numRedn = 0, numRedn2 = 0;
      CkReduction::tupleElement* res = nullptr, *res2 = nullptr;
      msg->toTuple(&res, &numRedn);
      Real max_velocity = *(Real*)(res[0].data); // avoid max_velocity = 0.0
      Real timestep_size = paratreet::getTimestep(universe, max_velocity);

      ProxyPack<Data> proxy_pack (this->thisProxy, subtrees, partitions, cache_manager);

      // Prefetch into cache
      start_time = CkWallTimer();
      // use exactly one of these three commands to load the software cache
      paratreet::preTraversalFn(proxy_pack);
      CkWaitQD();
      CkPrintf("TreeCanopy cache loading: %.3lf ms\n",
          (CkWallTimer() - start_time) * 1000);

      // Perform traversals
      start_time = CkWallTimer();
      paratreet::traversalFn(universe, proxy_pack, iter);
      CkWaitQD();
      CkPrintf("Tree traversal: %.3lf ms\n", (CkWallTimer() - start_time) * 1000);

      start_time = CkWallTimer();

      // Move the particles in Partitions
      partitions.kick(timestep_size, CkCallbackResumeThread());  // TODO: determine if I should remove this step

      // Now track PE imbalance for memory reasons
      thread_state_holder.collectMetaData(CkCallbackResumeThread((void *&) msg2));
      msg2->toTuple(&res2, &numRedn2);
      int numParticleCopies = *(int*)(res2[2].data);
      int numParticleShares = *(int*)(res2[3].data);
      int maxPESize = *(int*)(res2[0].data);
      int sumPESize = *(int*)(res2[1].data);
      float avgPESize = (float) universe.n_particles / (float) CkNumPes();
      float ratio = (float) maxPESize / avgPESize;
      bool complete_rebuild = false;  // only doing one iteration so no need to do complete rebuild
      CkPrintf("[Meta] n_subtree = %d; timestep_size = %f; numPSParticleCopies = %d; numPSParticleShares = %d; sumPESize = %d; maxPESize = %d, avgPESize = %f; ratio = %f; maxVelocity = %f; rebuild = %s\n", n_subtrees, timestep_size, numParticleCopies, numParticleShares, sumPESize, maxPESize, avgPESize, ratio, max_velocity, (complete_rebuild? "yes" : "no"));
      //End Subtree reduction message parsing

      paratreet::postIterationFn(universe, proxy_pack, iter);

      CkReductionMsg* result;
      partitions.perturb(timestep_size, CkCallbackResumeThread((void *&)result));
      universe = *((BoundingBox*)result->getData());
      delete result;
      remakeUniverse();
      partitions.rebuild(universe, subtrees, complete_rebuild); // 0.1s for example
      CkWaitQD();
      CkPrintf("Perturbations: %.3lf ms\n", (CkWallTimer() - start_time) * 1000);
      if (!complete_rebuild && config.lb_period > 0 && iter % config.lb_period == config.lb_period - 1){
        start_time = CkWallTimer();
        //subtrees.pauseForLB(); // move them later
        partitions.pauseForLB();
        CkWaitQD();
        CkPrintf("Load balancing: %.3lf ms\n", (CkWallTimer() - start_time) * 1000);
      }
      // Destroy subtrees and perform decomposition from scratch
      resumer.reset();
      if (complete_rebuild) {
        treespec.reset();
        subtrees.destroy();
        partitions.destroy();
        decompose(iter+1);
      } else {
        partitions.reset();
        subtrees.reset();
      }

      // Clear cache and other storages used in this iteration
      cache_manager.destroy(true);
      CkCallback statsCb (CkReductionTarget(Driver<Data>, countInts), this->thisProxy);
      thread_state_holder.collectAndResetStats(statsCb);
      storage.clear();
      storage_sorted = false;
      CkWaitQD();
      double iter_time = CkWallTimer() - iter_start_time;
      total_time += iter_time;
      CkPrintf("Iteration %d time: %.3lf ms\n", iter, iter_time * 1000);
      if (iter == config.num_iterations-1) {
        CkPrintf("Average iteration time: %.3lf ms\n", total_time / config.num_iterations * 1000);
      }

      cb.send();
    }

    // -------------------
    // Auxiliary functions
    // -------------------
    void reportTime(){
      CkPrintf("Tree build: %.3lf ms\n", (CkWallTimer() - start_time) * 1000);
    }

    void countInts(unsigned long long* intrn_counts) {
      CkPrintf("%llu node-particle interactions, %llu bucket-particle interactions %llu node opens, %llu node closes\n", intrn_counts[0], intrn_counts[1], intrn_counts[2], intrn_counts[3]);
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
    }
  };