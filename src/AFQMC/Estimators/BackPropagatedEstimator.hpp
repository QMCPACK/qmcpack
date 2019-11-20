#ifndef QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_HPP
#define QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_HPP

#include "AFQMC/config.h"
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include "io/hdf_multi.h"
#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/NewTimer.h"
#include "Utilities/Timer.h"

#include "AFQMC/Estimators/FullObsHandler.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Propagators/Propagator.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

namespace qmcplusplus
{
namespace afqmc
{

/*
 * Top class for back propagated estimators. 
 * An instance of this class will manage a set of back propagated observables.
 * The main task of this class is to coordinate the generation of left-handed
 * states during the back propagation algorithm. The calculation and accumulation
 * of actual observables is handled by an object of the variant FullObsHandler.
 * BackPropagatedEstimator provides a list of back propagated references to 
 * FullObsHandler, whose job it is to calculate, accumulate and print quantities.
 */
class BackPropagatedEstimator: public EstimatorBase
{
  // allocators
  using Allocator = device_allocator<ComplexType>;

  // type defs
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;

  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;
  using CVector = boost::multi::array<ComplexType,1,Allocator>;
  using CMatrix = boost::multi::array<ComplexType,2,Allocator>;
  using stdCVector_ref = boost::multi::array_ref<ComplexType,1>;
  using stdCVector = boost::multi::array<ComplexType,1>;
  using stdCMatrix = boost::multi::array<ComplexType,2>;
  using stdCTensor = boost::multi::array<ComplexType,3>;
  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;

  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string name, xmlNodePtr cur, WALKER_TYPES wlk, WalkerSet& wset, 
        Wavefunction& wfn, Propagator& prop, bool impsamp_=true) :
                                      EstimatorBase(info),TG(tg_), walker_type(wlk),
                                      writer(false), 
                                      Refs({0,0,0},shared_allocator<ComplexType>{TG.TG_local()}),
                                      observ0(TG,info,name,cur,wlk,wfn), wfn0(wfn), prop0(prop),
                                      max_nback_prop(10),
                                      nStabalize(10), block_size(1), path_restoration(false),
                                      importanceSampling(impsamp_),first(true)
  {
    int nave(1);
    if(cur != NULL) {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(max_nback_prop, "nsteps", "int");
      m_param.add(nave, "naverages", "int");
      m_param.add(restore_paths, "path_restoration", "std::string");
      m_param.add(block_size, "block_size", "int");
      m_param.add(nblocks_skip, "nskip", "int");
      m_param.put(cur);
      if(restore_paths == "true") {
        path_restoration = true;
      } else {
        path_restoration = false;
      }
    }

    if(nave <= 0)
      APP_ABORT("naverages <= 0 is not allowed.\n");

    nback_prop_steps.reserve(nave);
    for(int i=1; i<nave; i++)
      nback_prop_steps.push_back(i*max_nback_prop/nave);
    nback_prop_steps.push_back(max_nback_prop);

    // sort the requested number of steps
    std::sort(nback_prop_steps.begin(),nback_prop_steps.end());

    if(max_nback_prop <= 0)
      APP_ABORT("max_nback_prop <= 0 is not allowed.\n");

    int ncv(prop0.global_number_of_cholesky_vectors());
    int nref(wfn0.number_of_references_for_back_propagation());
    wset.resize_bp(max_nback_prop,ncv,nref);
    wset.setBPPos(0);
    // set SMN in case BP begins right away
    if(nblocks_skip==0)
      for(auto it=wset.begin(); it<wset.end(); ++it)
        it->setSlaterMatrixN(); 

    if(TG.getNCoresPerTG() > 1)
      APP_ABORT("ncores > 1 is broken with back propagation. Fix this.");
    writer = (TG.getGlobalRank()==0);
  }
  
  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    accumulated_in_last_block=false;
    int bp_step = wset.getBPPos();
    if(bp_step <=0)
      APP_ABORT(" Error: Found bp_step <=0 in BackPropagate::accumulate_block. \n"); 
    if(bp_step > max_nback_prop)
      APP_ABORT(" Error: max_nback_prop in back propagation estimator must be conmensuate with nStep*nSubStep.\n");
    if(max_nback_prop > wset.NumBackProp())
      APP_ABORT(" Error: max_nback_prop > wset.NumBackProp() \n");

    // check if measurement is needed 
    int iav(-1);
    for(int i=0; i<nback_prop_steps.size(); i++) {
      if(bp_step == nback_prop_steps[i]) {
        iav = i; 
        break;
      }    
    }
    if(iav < 0) return;

    using std::fill_n;
    // 0. skip if requested  
    if(iblock < nblocks_skip) {
      if(bp_step == max_nback_prop) {
        if( iblock+1 == nblocks_skip )
          for(auto it=wset.begin(); it<wset.end(); ++it)
            it->setSlaterMatrixN(); 
        iblock++;
        wset.setBPPos(0);
      }  
      return;
    }

    AFQMCTimers[back_propagate_timer]->start();
    int nrefs = wfn0.number_of_references_for_back_propagation();
    int nrow(NMO*((walker_type==NONCOLLINEAR)?2:1));
    int ncol(NAEA+((walker_type==CLOSED)?0:NAEB));
    int nx((walker_type==COLLINEAR)?2:1);

    // 1. check structures
    if(Refs.size(0) != wset.size() || Refs.size(1) != nrefs || Refs.size(2) !=nrow*ncol) 
      Refs = std::move(mpi3CTensor({wset.size(),nrefs,nrow*ncol},Refs.get_allocator()));
    if(detR.size(0) != wset.size() || detR.size(1) != nx*nrefs)
      detR.reextent({wset.size(),nrefs*nx}); 

    int n0,n1;
    std::tie(n0,n1) = FairDivideBoundary(TG.getLocalTGRank(),int(Refs.size(2)),TG.getNCoresPerTG());
    boost::multi::array_ref<ComplexType,3> Refs_(to_address(Refs.origin()),Refs.extensions());

    // 2. setup back propagated references
    wfn0.getReferencesForBackPropagation(Refs_[0]);
    for(int iw=1; iw<wset.size(); ++iw)
      for(int ref=0; ref<nrefs; ++ref)
       copy_n(Refs_[0][ref].origin()+n0,n1-n0,Refs_[iw][ref].origin()+n0);
    TG.TG_local().barrier();

    //3. propagate backwards the references
    prop0.BackPropagate(bp_step,nStabalize,wset,Refs_,detR);

    //4. calculate properties 
    // adjust weights here is path restoration
    stdCVector wgt(iextensions<1u>{wset.size()});
    wset.getProperty(WEIGHT,wgt);
    if(path_restoration) {
      auto&& factors(wset.getWeightFactors());
      for(int k=0; k<bp_step; k++)  
        for(int i=0; i<wgt.size(); i++) 
          wgt[i] *= factors[k][i];
    } else if(!importanceSampling) {
      stdCVector phase(iextensions<1u>{wset.size()});
      wset.getProperty(PHASE,phase);
      for(int i=0; i<wgt.size(); i++) wgt[i] *= phase[i];
    }
    observ0.accumulate(iav, wset, Refs_, wgt, detR, importanceSampling);

    if(bp_step == max_nback_prop) {
      // 5. setup for next block 
      for(auto it=wset.begin(); it<wset.end(); ++it)
        it->setSlaterMatrixN(); 
      wset.setBPPos(0);

      // 6. increase block counter
      iblock++;
      accumulated_in_last_block=true;
    }
    AFQMCTimers[back_propagate_timer]->stop();
  }

  void tags(std::ofstream& out) {
    if(writer) 
      out<<"BP_timer ";
  }

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset)
  {
    // I doubt we will ever collect a billion blocks of data.
    if(writer) {
      out<<std::setprecision(5)
                    <<AFQMCTimers[back_propagate_timer]->get_total() <<" ";
      AFQMCTimers[back_propagate_timer]->reset();  
    }
    if(accumulated_in_last_block) {
      if(writer && first) { 
        first=false;
        int nave(nback_prop_steps.size());
        int nref(detR.size(1));
        if(write_metadata) {
          dump.push("Observables");
          dump.push("BackPropagated");
          dump.push("Metadata");
          dump.write(nback_prop_steps, "BackPropSteps");
          dump.write(nave, "NumAverages");
          dump.write(nref, "NumReferences");
          dump.pop();
          dump.pop();
          dump.pop();
          write_metadata = false;
        }
      }
      if(writer) { 
        dump.push("Observables");
        dump.push("BackPropagated");
      }
      observ0.print(iblock,dump);
      if(writer) { 
        dump.pop();
        dump.pop();
      }  
    }
  }

  private:

  TaskGroup_& TG;

  WALKER_TYPES walker_type;

  bool writer;
  bool accumulated_in_last_block;

  mpi3CTensor Refs;

  FullObsHandler observ0; 

  Wavefunction& wfn0;

  Propagator& prop0;

  int max_nback_prop;
  std::vector<int> nback_prop_steps;  

  boost::multi::array<ComplexType,2> detR;

  RealType weight, weight_sub;
  RealType targetW = 1;
  int iblock = 0;
  int nblocks_skip = 0;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one = ComplexType(1.0, 0.0);

  // Frequency of reorthogonalisation.
  int nStabalize;
  // Block size over which RDM will be averaged.
  int block_size;
  // Whether to restore cosine projection and real local energy apprximation for weights
  // along back propagation path.
  bool path_restoration, importanceSampling;

  int first;

  bool write_metadata = true;

};
}
}

#endif
