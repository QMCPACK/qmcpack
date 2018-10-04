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
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class BackPropagatedEstimator: public EstimatorBase
{

  using SMType = boost::multi_array_ref<ComplexType,2>;
  using const_SMType = boost::multi_array_ref<const ComplexType,2>;
  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string title, xmlNodePtr cur, WALKER_TYPES wlk, Wavefunction& wfn) :
                                            EstimatorBase(info),TG(tg_), wfn0(wfn),
                                            greens_function(false),
                                            nStabalize(1), fully_restore_weights(false),
                                            writer(false), partially_restore_weights(false),
                                            SDetOp(((wlk!=2)?(NMO):(2*NMO)),((wlk!=2)?(NAEA):(NAEA+NAEB)))
  {

    if(cur != NULL) {
      ParameterSet m_param;
      std::string weight_restoration;
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(weight_restoration, "restore_weights", "std::string");
      m_param.put(cur);
      if (weight_restoration == "partial") {
        partially_restore_weights = true;
      } else if (weight_restoration == "full") {
        fully_restore_weights = true;
      } else {
        // Other choices?
      }
    }

    ncores_per_TG = TG.getNCoresPerTG(); 
    core_rank = TG.getLocalTGRank();
    writer = (TG.getGlobalRank()==0);

    // Back propagated walker Slater Matrix.
    int NEL = (wlk==CLOSED)?(NAEA):(NAEA+NAEB);
    int NMO2 = (wlk==CLOSED)?(NMO):(2*NMO);
    backPropSM.resize(extents[NMO][NEL]);
    walkerDM.resize(extents[NMO2][NMO]);
    T1.resize(extents[NMO2][NAEA]);
    // (weight, GF)
    data.resize(1+NMO2*NMO);
    std::fill(data.begin(), data.end(), zero);
  }

  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
    size_t nwalk = wset.size();
    int nback_prop = wset.getNBackProp();
    std::fill(data.begin(), data.end(), 0.0);
    const std::vector<PsiT_Matrix> trialSM2 = wfn0.getOrbMat();
    const std::vector<PsiT_Matrix> trialSM;
    auto walker_type = wset.getWalkerType();
    for (int i = 0; i < nwalk; i++) {
      if (i%ncores_per_TG == core_rank) {
        // Initial expansion coefficients for back propagation. These may be modified during
        // the QR decomposition (stabalisation) step.
        ComplexType weight = wi.weight();
        if (std::abs(weight) <= 1e-6 || std::isnan(weight.real())
            || !wi.isBMatrixBufferFull()) continue;
        // 1. Back propagate each walker.
        // Initial state for back propagation.
        SMType BPAlpha = SMType(backPropSM.data(), extents[NMO][NAEA]);
        ComplexType detR = back_propagate_wavefunction(trialSM, BPAlpha, wi, nback_prop);
        // 2. Construct Green's function / 1RDM.
        // Note we use the Slater Matrix from timestep N (before back propagation).
        SMType DMAlpha = SMType(walkerDM.data(), extents[NMO][NMO]);
        SDetOp.MixedDensityMatrix_noHerm(BPAlpha, wi.SlaterMatrixN(Alpha), std::forward<SMType>(DMAlpha));
        if (walker_type == COLLINEAR) {
          // Do the same for down spin block.
          SMType BPBeta = SMType(backPropSM.data()+NMO*NAEA, extents[NMO][NAEB]);
          detR *= back_propagate_wavefunction(trialSM, BPBeta, wi, nback_prop);
          SMType DMBeta = SMType(walkerDM.data()+NMO*NMO, extents[NMO][NMO]);
          SDetOp.MixedDensityMatrix_noHerm(BPBeta, wi.SlaterMatrixN(Alpha), std::forward<SMType>(DMBeta));
        }
        // 3. Update averages.
        // 3.a Get walker weight at timestep M+N
        if (fully_restore_weights || partially_restore_weights) {
          weight *= reweighting_factor(wi, nback_prop);
        }
        if (!std::isfinite(weight.real())) continue;
        data[0] += weight;
        // 3.b Update Green's function estimator.
        // Update all the elements of the Green's function estimator.
        for (int i = 0; i < walkerDM.shape()[0]; i++) {
          for (int j = 0; j < walkerDM.shape()[1]; j++) {
            data[i+j*walkerDM.shape()[1]+1] += weight * walkerDM[i][j];
          }
        }
        // 4. Copy current Slater Matrix (timestep M+N) to be the Nth Wavefunction for the
        // next round of back propagation.
        // Warning: This assumes that estimators and back propagation are linked. We would
        // need to move this elsewhere if this link is broken.
        wi.setSlaterMatrixN();
      } // walker loop
    }
    // Crashing code
    TG.Global().all_reduce_in_place_n(data.begin(),data.size(),std::plus<>());
  }

  template <typename T>
  ComplexType back_propagate_wavefunction(const boost::multi_array<ComplexType,2>& trialSM, SMType& psiBP, T& walker, int nback_prop)
  {
    ComplexType detR = one;
    walker.decrementBMatrix();
    SMType B = walker.BMatrix();
    ma::product(ma::H(B), trialSM, T1);
    for (int i = 0; i < nback_prop-1; i++) {
      SMType B = walker.BMatrix();
      ma::product(ma::H(B), T1, std::forward<SMType>(psiBP));
      T1 = psiBP;
      if ((i != 0) && (i%nStabalize == 0)) {
        //detR *= DenseMatrixOperators::GeneralizedGramSchmidt(T1.data(),NAEA,NMO,NAEA);
      }
      walker.decrementBMatrix();
    }
    return detR;
  }

  void tags(std::ofstream& out)
  {
    if (TG.getGlobalRank() == 0) {
      for (int i = 0; i < walkerDM.shape()[0]; i++) {
        for (int j = 0; j < walkerDM.shape()[1]; j++) {
          if (i < NMO) {
            out << "Gup_" << i << "_"  << j << " ";
          } else {
              out << "Gdn_" << i-NMO << "_"  << j << " ";
          }
        }
      }
    }
  }

  void print(std::ofstream& out, WalkerSet& wset)
  {
    if (writer) {
      for (int i = 1; i < data.size(); i++) {
        out << std::setprecision(16) << " " << (data[i]/data[0]).real() << " ";
      }
    }
    // Zero our estimator array.
    std::fill(data.begin(), data.end(), 0);
  }

  private:

  template <typename T>
  ComplexType reweighting_factor(T walker, int nback_prop)
  {
    ComplexType factor = one;
    // Get final entry in back propagation buffers = first entry for back propagation.
    for (int i = 0; i < nback_prop; i++) {
      // undo cosine projection.
      walker.decrementBMatrix();
      factor *= walker.weightFactor();
      if (fully_restore_weights) {
        // restore imaginary part of local energy.
        factor /= walker.cosineFactor();
      }
    }
    return factor;
  }

  TaskGroup_& TG;

  bool writer;

  Wavefunction& wfn0;
  SlaterDetOperations<ComplexType> SDetOp;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  std::vector<ComplexType> data;

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
  bool fully_restore_weights, partially_restore_weights;

  SlaterDetOperations SDet;
  // Temporary storage for a given walker's density matrix.
  boost::multi_array<ComplexType,2> backPropSM, T1, walkerDM;

};
}
}

#endif
