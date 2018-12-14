#ifndef QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_BACKPROPAGATEDESTIMATOR_H

#include <Message/MPIObjectBase.h>
#include "AFQMC/config.h"
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

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

  using CMatrix_ref = boost::multi_array_ref<ComplexType,2>;
  using CVector = boost::multi_array<ComplexType,1>;
  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string title, xmlNodePtr cur, WALKER_TYPES wlk, Wavefunction& wfn, bool impsamp_=true) :
                                            EstimatorBase(info),TG(tg_), wfn0(wfn),
                                            greens_function(false),
                                            nStabalize(1), path_restoration(false),
                                            writer(false), importanceSampling(impsamp_)
  {

    if(cur != NULL) {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(restore_paths, "path_restoration", "std::string");
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
      DMBuffer.resize(extents[dm_size]);
    }
    std::fill(DMBuffer.begin(), DMBuffer.end(), ComplexType(0.0,0.0));
    denom.resize(extents[1]);
  }

  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    bool transpose = false, compact = false;
    // check to see whether we should be accumulating estimates.
    bool back_propagate = wset[0].isBMatrixBufferFull();
    if(back_propagate) {
      CMatrix_ref BackPropDM(DMBuffer.data(), extents[dm_dims.first][dm_dims.second]);
      // Computes GBP(i,j)
      denom[0] = ComplexType(0.0,0.0);
      wfn0.BackPropagatedDensityMatrix(wset, BackPropDM, denom, path_restoration, !importanceSampling);
      for(int iw = 0; iw < wset.size(); iw++) {
        // Resets B matrix buffer to identity, copies current wavefunction and resets weight
        // factor.
        if( iw%TG.TG_local().size() != TG.TG_local().rank() ) continue;
        wset[iw].resetForBackPropagation();
      }
    }
  }

  //template <typename T>
  //ComplexType back_propagate_wavefunction(const boost::multi_array<ComplexType,2>& trialSM, SMType& psiBP, T& walker, int nback_prop)
  //{
    //ComplexType detR = one;
    //walker.decrementBMatrix();
    //SMType B = walker.BMatrix();
    //ma::product(ma::H(B), trialSM, T1);
    //for (int i = 0; i < nback_prop-1; i++) {
      //walker.decrementBMatrix();
      //SMType B = walker.BMatrix();
      //ma::product(ma::H(B), T1, std::forward<SMType>(psiBP));
      //T1 = psiBP;
      //if ((i != 0) && (i%nStabalize == 0)) {
        ////detR *= DenseMatrixOperators::GeneralizedGramSchmidt(T1.data(),NAEA,NMO,NAEA);
      //}
    //}
    //return detR;
  //}

  void tags(std::ofstream& out)
  {
    if(TG.getGlobalRank() == 0) {
      for(int i = 0; i < dm_dims.first; i++) {
        for(int j = 0; j < dm_dims.second; j++) {
          out << "ReG_" << i << "_"  << j << " " << "ImG_" << i << "_" << j << " ";
        }
      }
      if(!importanceSampling) out << "ReG_denom ImG_denom ";
    }
  }

  void print(std::ofstream& out, WalkerSet& wset)
  {
    if(writer) {
      if(importanceSampling) {
        for(int i = 0; i < DMBuffer.size(); i++) {
          out << std::setprecision(16) << " " << DMBuffer[i].real() << " ";
        }
      } else {
        for(int i = 0; i < DMBuffer.size(); i++) {
          //RealType re_num = DMBuffer[i].real()*denom[0].real() + DMBuffer[i].imag()*denom[0].imag();
          out << std::setprecision(16) << " " << DMBuffer[i].real() << " " << DMBuffer[i].imag() << " ";
        }
        out << denom[0].real() << " " << denom[0].imag() << " ";
      }
    }
    // Zero our estimator array.
    std::fill(DMBuffer.begin(), DMBuffer.end(), ComplexType(0.0,0.0));
    denom[0] = ComplexType(0.0,0.0);
  }

  private:

  //template <typename T>
  //ComplexType reweighting_factor(T walker, int nback_prop)
  //{
    //ComplexType factor = one;
    //// Get final entry in back propagation buffers = first entry for back propagation.
    //for (int i = 0; i < nback_prop; i++) {
      //// undo cosine projection.
      //walker.decrementBMatrix();
      //factor *= walker.weightFactor();
      //if (fully_restore_weights) {
        //// restore imaginary part of local energy.
        //factor /= walker.cosineFactor();
      //}
    //}
    //return factor;
  //}

  TaskGroup_& TG;

  bool writer;

  Wavefunction& wfn0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  CVector DMBuffer;

  RealType weight, weight_sub;
  RealType targetW = 1;
  int core_rank;
  int ncores_per_TG;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one = ComplexType(1.0, 0.0);

  // Print the full mixed estimator for the one-particle reduced density matrix.
  bool greens_function;
  // Frequency of reorthogonalisation.
  int nStabalize;
  // Whether to restore cosine projection and real local energy apprximation for weights
  // along back propagation path.
  bool path_restoration, importanceSampling;
  std::vector<ComplexType> weights;
  int dm_size;
  std::pair<int,int> dm_dims;
  CVector denom;

};
}
}

#endif
