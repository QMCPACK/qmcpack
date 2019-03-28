#ifndef QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_H

#include <Message/MPIObjectBase.h>
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

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Propagators/Propagator.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class BackPropagatedEstimator: public EstimatorBase
{

  // allocators
  using Allocator = device_allocator<ComplexType>;

  // type defs
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;

  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;
  using CVector = boost::multi::array<ComplexType,1,Allocator>;
  using stdCVector = boost::multi::array<ComplexType,1>;
  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string title, xmlNodePtr cur, WALKER_TYPES wlk, WalkerSet& wset, 
        Wavefunction& wfn, Propagator& prop, bool impsamp_=true) :
                                            EstimatorBase(info),TG(tg_), wfn0(wfn), prop0(prop),
                                            greens_function(false),nback_prop(10),
                                            nStabalize(1), path_restoration(false), block_size(1),
                                            writer(false), importanceSampling(impsamp_)
  {

    if(cur != NULL) {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(nback_prop, "nsteps", "int");
      m_param.add(restore_paths, "path_restoration", "std::string");
      m_param.add(block_size, "block_size", "int");
      m_param.put(cur);
      if(restore_paths == "true") {
        path_restoration = true;
      } else {
        path_restoration = false;
      }
    }

    if(nback_prop <= 0)
      APP_ABORT("nback_prop == 0 is not allowed.\n");

    int ncv(prop0.global_number_of_cholesky_vectors());
    int nref(wfn0.number_of_references_for_back_propagation());
    wset.resize_bp(nback_prop,ncv,nref);
    // set SMN in case BP begins right away
    for(auto it=wset.begin(); it<wset.end(); ++it)
      it->setSlaterMatrixN(); 

    ncores_per_TG = TG.getNCoresPerTG();
    if(ncores_per_TG > 1)
      APP_ABORT("ncores > 1 is broken with back propagation. Fix this.");
    core_rank = TG.getLocalTGRank();
    writer = (TG.getGlobalRank()==0);
    if(wlk == CLOSED) {
      dm_size = NMO*NMO;
      dm_dims = {NMO,NMO};
    } else if(wlk == COLLINEAR) {
      dm_size = 2*NMO*NMO;
      dm_dims = {2*NMO,NMO};
    } else if(wlk == NONCOLLINEAR) {
      dm_size = 4*NMO*NMO;
      dm_dims = {2*NMO,2*NMO};
    }
    if(DMBuffer.size() < dm_size) {
      DMBuffer.reextent(iextensions<1u>{dm_size});
    }
    if(DMAverage.size() < dm_size) {
      DMAverage.reextent(iextensions<1u>{dm_size});
    }
    using std::fill_n;
    fill_n(DMBuffer.origin(), DMBuffer.num_elements(), ComplexType(0.0,0.0));
    fill_n(DMAverage.origin(), DMAverage.num_elements(), ComplexType(0.0,0.0));
    denom.reextent({1});
    denom_average.reextent({1});
    denom_average[0] = 0;
  }

  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
/*
    // check to see whether we should be accumulating estimates.
    bool back_propagate = wset[0].isBMatrixBufferFull();
    if(back_propagate) {
      CMatrix_ref BackPropDM(DMBuffer.origin(), {dm_dims.first,dm_dims.second});
      // Computes GBP(i,j)
      using std::fill_n;
      fill_n(denom.origin(),1,ComplexType(0.0,0.0));
      fill_n(DMBuffer.origin(), DMBuffer.num_elements(), ComplexType(0.0,0.0));
      wfn0.WalkerAveragedDensityMatrix(wset, BackPropDM, denom, path_restoration, !importanceSampling, back_propagate); 
      for(int iw = 0; iw < wset.size(); iw++) {
        // Resets B matrix buffer to identity, copies current wavefunction and resets weight
        // factor.
        if( iw%TG.TG_local().size() != TG.TG_local().rank() ) continue;
        wset[iw].resetForBackPropagation();
      }
      iblock++;
    }
*/
  }

  void tags(std::ofstream& out) {}

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset)
  {
    // I doubt we will ever collect a billion blocks of data.
    int n_zero = 9;
    if(writer) {
      if(write_metadata) {
        dump.push("Metadata");
        dump.write(nback_prop, "NumBackProp");
        dump.pop();
        write_metadata = false;
      }
      bool write = false; //wset[0].isBMatrixBufferFull();
      if(write) {
//        for(int i = 0; i < DMBuffer.size(); i++)
//          DMAverage[i] += DMBuffer[i];
// MAM: make a wrapper for this type of operation
// e.g. auto reference_or_copy<stdCVector>(DMBuffer);
#ifdef QMC_CUDA
        stdCVector buff(DMBuffer);
#else
        CVector& buff(DMBuffer);
#endif
        ma::axpy(ComplexType(1.0),buff,DMAverage);
        denom_average[0] += ComplexType(*denom.origin());
        if(iblock%block_size == 0) {
          for(int i = 0; i < DMAverage.size(); i++)
            DMAverage[i] /= block_size;
          denom_average[0] /= block_size;
          dump.push("BackPropagated");
          std::string padded_iblock = std::string(n_zero-std::to_string(iblock).length(),'0')+std::to_string(iblock);
          dump.write(DMAverage, "one_rdm_"+padded_iblock);
          dump.write(denom_average, "one_rdm_denom_"+padded_iblock);
          dump.pop();
          using std::fill;  
          fill(DMAverage.begin(), DMAverage.end(), ComplexType(0.0,0.0));
          fill(denom_average.begin(), denom_average.end(), ComplexType(0.0,0.0));
        }
      }
    }
  }

  private:

  TaskGroup_& TG;

  bool writer;

  int nback_prop;

  Wavefunction& wfn0;

  Propagator& prop0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  CVector DMBuffer;
  stdCVector DMAverage;

  RealType weight, weight_sub;
  RealType targetW = 1;
  int core_rank;
  int ncores_per_TG;
  int iblock = 0;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one = ComplexType(1.0, 0.0);

  // Print the full mixed estimator for the one-particle reduced density matrix.
  bool greens_function;
  // Frequency of reorthogonalisation.
  int nStabalize;
  // Block size over which RDM will be averaged.
  int block_size;
  // Whether to restore cosine projection and real local energy apprximation for weights
  // along back propagation path.
  bool path_restoration, importanceSampling;
  std::vector<ComplexType> weights;
  int dm_size;
  std::pair<int,int> dm_dims;
  CVector denom;
  stdCVector denom_average;
  bool write_metadata = true;

};
}
}

#endif
