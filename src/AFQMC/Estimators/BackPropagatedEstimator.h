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

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class BackPropagatedEstimator: public EstimatorBase
{

  using CMatrix_ref = boost::multi::array_ref<ComplexType,2>;
  using CVector = boost::multi::array<ComplexType,1>;
  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string title, xmlNodePtr cur, WALKER_TYPES wlk, Wavefunction& wfn, bool impsamp_=true) :
                                            EstimatorBase(info),TG(tg_), wfn0(wfn),
                                            greens_function(false),
                                            nStabalize(1), path_restoration(false), block_size(1),
                                            writer(false), importanceSampling(impsamp_)
  {

    if(cur != NULL) {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(restore_paths, "path_restoration", "std::string");
      m_param.add(block_size, "block_size", "int");
      m_param.put(cur);
      if(restore_paths == "true") {
        path_restoration = true;
      } else {
        path_restoration = false;
      }
    }

    ncores_per_TG = TG.getNCoresPerTG();
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
    std::fill(DMBuffer.begin(), DMBuffer.end(), ComplexType(0.0,0.0));
    std::fill(DMAverage.begin(), DMAverage.end(), ComplexType(0.0,0.0));
    denom.reextent({1});
    denom_average.reextent({1});
    denom_average[0] = 0;
  }

  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    // check to see whether we should be accumulating estimates.
    bool back_propagate = wset[0].isBMatrixBufferFull();
    if(back_propagate) {
      CMatrix_ref BackPropDM(DMBuffer.data(), {dm_dims.first,dm_dims.second});
      // Computes GBP(i,j)
      denom[0] = ComplexType(0.0,0.0);
      std::fill(DMBuffer.begin(), DMBuffer.end(), ComplexType(0.0,0.0));
      wfn0.WalkerAveragedDensityMatrix(wset, BackPropDM, denom, path_restoration, !importanceSampling, back_propagate);
      for(int iw = 0; iw < wset.size(); iw++) {
        // Resets B matrix buffer to identity, copies current wavefunction and resets weight
        // factor.
        if( iw%TG.TG_local().size() != TG.TG_local().rank() ) continue;
        wset[iw].resetForBackPropagation();
      }
      write_back_prop = true;
      iblock++;
    }
  }

  void tags(std::ofstream& out) {}

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset)
  {
    if(writer) {
      if(write_back_prop) {
        for(int i = 0; i < DMBuffer.size(); i++)
          DMAverage[i] += DMBuffer[i];
        denom_average[0] += denom[0];
        if(iblock%block_size == 0) {
          for(int i = 0; i < DMAverage.size(); i++)
            DMAverage[i] /= block_size;
          denom_average[0] /= block_size;
          dump.push("BackPropagated");
          dump.write(DMAverage, "one_rdm_"+std::to_string(iblock));
          dump.write(denom_average, "one_rdm_denom_"+std::to_string(iblock));
          dump.pop();
          std::fill(DMAverage.begin(), DMAverage.end(), ComplexType(0.0,0.0));
          std::fill(denom_average.begin(), denom_average.end(), ComplexType(0.0,0.0));
        }
        write_back_prop = false;
      }
    }
  }

  private:

  TaskGroup_& TG;

  bool writer;

  Wavefunction& wfn0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  CVector DMBuffer, DMAverage;

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
  bool write_back_prop = false;
  std::pair<int,int> dm_dims;
  CVector denom, denom_average;

};
}
}

#endif
