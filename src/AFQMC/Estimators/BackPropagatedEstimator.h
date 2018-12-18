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
#include "AFQMC/Wavefunctions/SlaterDetOperations.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class BackPropagatedEstimator: public EstimatorBase
{

  public:

  BackPropagatedEstimator(afqmc::TaskGroup_& tg_, AFQMCInfo& info,
        std::string title, xmlNodePtr cur, Wavefunction& wfn) : 
                                            EstimatorBase(info),TG(tg_), wfn0(wfn),
                                            greens_function(false),
                                            nStabalize(1), fully_restore_weights(false),
                                            writer(false),partially_restore_weights(false), 
                                            SDet(info) 
  {

    if(cur != NULL) {
      ParameterSet m_param;
      std::string gf;
      std::string weight_restoration;
      m_param.add(gf, "greens_function", "std::string");
      m_param.add(nStabalize, "ortho", "int");
      m_param.add(weight_restoration, "restore_weights", "std::string");
      m_param.put(cur);
      if (gf == "true") {
        greens_function = true;
      }
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

    // Manually construct indexing for some predetermined options.
    if (greens_function) {
      for (int i = 0; i < 2*NMO; i++) {
        for (int j = 0; j < NMO; j++) {
          int ix = twoDimensionalIndexToOneDimension(i, j, NMO);
          elements.push_back(ix);
        }
      }
    }
    // (weight, GF)
    data.resize(1+elements.size());
    tmp_data.resize(data.size());
    std::fill(data.begin(), data.end(), 0.0);
    walker_density_matrix.resize(2*NMO, NMO);

    // Store for back propagated walker's Slater Matrix and expansion coefficients.
APP_ABORT("Finish \n");
/*
    BPCoeff = wfn0.getCiCoeff();
    int nterms = BPCoeff.size();
    closed_shell = wfn0.isClosedShell();
    orbsize = wfn0.getOrbSize();
    if (nterms == 1) {
      // Single-determinant trial wavefunction.
      BPSlaterMat.resize(orbsize,NAEA);
    } else {
      // Multi-determinant trial wavefunction.
      BPSlaterMat.resize(nterms*orbsize,NAEA);
    }
    // Temporary storage
    T1.resize(NMO,NAEA);
    T2.resize(NMO,NAEA);
*/
  } 

  ~BackPropagatedEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) {}

  void accumulate_block(WalkerSet& wset)
  {
APP_ABORT("Finish \n");
/*
    int nW0 = wset->numWalkers(true);
    ComplexType weight, oa, ob, w, ooa, oob, tmp, factor, detR;
    ComplexType *SM, *backPropSM, *trialSM;
    // Initial state for back propagation.
    trialSM = wfn0.getOrbMat();
    // Store for Back-propagated wavefunction.
    backPropSM = BPSlaterMat.data();
    for (int i=0; i<nW0; i++) {
      if (i%ncores_per_TG == core_rank) {
        // Initial expansion coefficients for back propagation. These may be modified during
        // the QR decomposition (stabalisation) step.
        BPCoeff = wfn0.getCiCoeff();
        if (!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6
            || std::isnan(wset->getWeight(i).real()) || !wset->isPMBufferFull(i)) continue;
        // 1. Back propagate each walker.
        if (BPCoeff.size() == 1) {
          back_propagate_wavefunction(trialSM, backPropSM, wset->getNBackProp(), wset, i);
          // 2. Construct Green's function / 1RDM.
          // Note we use the Slater Matrix from timestep N (before back propagation).
          SDet.green_function(backPropSM, wset->getSMN(i), tmp, walker_density_matrix, true, closed_shell);
        } else {
          for (int ic = 0; ic < BPCoeff.size(); ic++) {
            // Construct back propagated wavefunction for each element of the trial wavefunction.
            detR = back_propagate_wavefunction(trialSM+ic*orbsize, backPropSM+ic*orbsize, wset->getNBackProp(), wset, i);
            // Include any additional factor from re-orthogonalisation.
            BPCoeff[ic] *= detR;
          }
          // 2. Construct current walker's Green's function / 1RDM.
          // Note we use the Slater Matrix from timestep N (before back propagation).
          SDet.green_function(backPropSM, wset->getSMN(i), tmp, walker_density_matrix, BPCoeff,
                              orbsize, closed_shell);
        }
        // 3. Evaluate full Green's function.
        // 3.a Get walker weight at timestep M+N
        SM = wset->getWalker(i, w, tmp, ooa, oob);  // "impsampl"
        wfn0.Overlap(wset,ovlp);
        weight = w * ((oa*ob)/(ooa*oob));
        if (fully_restore_weights || partially_restore_weights) {
          factor = reweighting_factor(wset, i);
          weight *= factor;
        }
        if (!std::isfinite(weight.real())) continue;
        data[0] += weight;
        // 3.b Update Green's function estimator.
        SPComplexMatrix::iterator itG = walker_density_matrix.begin();
        std::vector<ComplexType>::iterator itd = data.begin() + 1;
        // Update all the elements of the Green's function estimator.
        for (std::vector<int>::iterator it = elements.begin(); it != elements.end(); ++it) {
          *itd += weight * (*(itG + *it));
          ++itd;
        }
        // 4. Copy current Slater Matrix (timestep M+N) to be the Nth Wavefunction for the
        // next round of back propagation.
        // Warning: This assumes that estimators and back propagation are linked. We would
        // need to move this elsewhere if this link is broken.
        wset->copy_slater_matrix_to_historic_slater_matrix(i);
      } // walker loop
    }
    TG.Global().reduce_in_place_n(data.begin(),data.size(),std::plus<>());

//    std::fill(tmp_data.begin(), tmp_data.end(), 0.0);
//    MPI_Reduce(data.data(), tmp_data.data(), data.size(), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
//    std::swap(data, tmp_data);
*/
  }

  ComplexType back_propagate_wavefunction(ComplexType* trial, ComplexType* BP, int nback_prop, WalkerSet& wset, int n)
  {
APP_ABORT("Finish \n");
/*
    // Back propagate up spin block.
    // Avoid messy repeated copy of trial to T1 due to mismatch in pointer returned from
    // get* routines and full matrix.
    T1ptr = trial;
    ComplexType* PM;
    ComplexType detR = one;
    wset->decrementPM(n);
    for (int i = 0; i < nback_prop-1; i++) {
      PM = wset->getPM(n);
      DenseMatrixOperators::product_AhB(NMO, NAEA, NMO, one, PM, NMO, T1ptr, NAEA, zero, T2.data(), NAEA);
      T1 = T2;
      if ((i != 0) && (i%nStabalize == 0)) {
        detR *= DenseMatrixOperators::GeneralizedGramSchmidt(T1.data(),NAEA,NMO,NAEA);
      }
      T1ptr = T1.data();
      wset->decrementPM(n);
    }
    PM = wset->getPM(n);
    DenseMatrixOperators::product_AhB(NMO, NAEA, NMO, one, PM, NMO, T1ptr, NAEA, zero, BP, NAEA);
    if (closed_shell) {
      detR *= detR;
    } else {
      // Back propagate spin down block.
      int disp = ((spinRestricted)?0:NMO*NMO);
      T1ptr = trial + NMO*NAEA;
      wset->decrementPM(n);
      for (int i = 0; i < nback_prop-1; i++) {
        PM = wset->getPM(n);
        DenseMatrixOperators::product_AhB(NMO, NAEA, NMO, one, PM+disp, NMO, T1ptr, NAEA, zero, T2.data(), NAEA);
        T1 = T2;
        if ((i != 0) && (i%nStabalize == 0)) {
          detR *= DenseMatrixOperators::GeneralizedGramSchmidt(T1.data(),NAEA,NMO,NAEB);
        }
        T1ptr = T1.data();
        wset->decrementPM(n);
      }
      PM = wset->getPM(n);
      DenseMatrixOperators::product_AhB(NMO, NAEA, NMO, one, PM+disp, NMO, T1ptr, NAEA, zero, BP+NMO*NAEA, NAEA);
    }
    return detR;
*/
    return ComplexType(1.0);
  }

  void tags(std::ofstream& out)
  {
    if (TG.getGlobalRank() == 0) {
        std::vector<int> ij;
        ij.resize(2);
        for (std::vector<int>::iterator it=elements.begin(); it!=elements.end(); ++it) {
          oneDimensionalIndexToTwoDimensions(*it, NMO, ij);
          if (*it < NMO*NMO) {
            out << "Gup_" << ij[0] << "_"  << ij[1] << " ";
          } else {
            out << "Gdown_" << ij[0] << "_"  << ij[1] << " ";
          }
        }
    }
  }

  void print(std::ofstream& out, WalkerSet& wset)
  {
    if (writer) {
      for (int i = 1; i < data.size(); i++) {
        out << std::setprecision(16) << " " << (data[i].real()/data[0].real()) << " ";
      }
    }
    // Zero our estimator array.
    std::fill(data.begin(), data.end(), 0);
  }

  private:

  inline void oneDimensionalIndexToTwoDimensions(int ix, int ncols, std::vector<int>& coordinate)
  {
    coordinate[0] = ix / ncols;
    coordinate[1] = ix % ncols;
  }
  inline int twoDimensionalIndexToOneDimension(int x, int y, int ncols)
  {
    return ncols*x + y;
  }
  ComplexType reweighting_factor(WalkerSet& wset, int n)
  {
    ComplexType factor = one;
    // Get final entry in back propagation buffers = first entry for back propagation.
    for (int i = 0; i < wset.getNBackProp(); i++) {
      // undo cosine projection.
      wset.decrementPM(n);
      factor *= wset.getWeightFactor(n);
      if (fully_restore_weights) {
        // restore imaginary part of local energy.
        factor /= wset.getCosineFactor(n);
      }
    }
    return factor;
  }

  TaskGroup_& TG;

  bool writer;

  Wwavefunction& wfn0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  std::vector<ComplexType> data;
  // for MPI communication.
  std::vector<ComplexType> tmp_data;
  ComplexType gdenom = ComplexType(0.0, 0.0);

  RealType weight, weight_sub;
  RealType targetW = 1;
  int core_rank;
  int ncores_per_TG;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one = ComplexType(1.0, 0.0);

  int offset, orbsize;
  bool closed_shell;

  // Print the full mixed estimator for the one-particle reduced density matrix.
  bool greens_function;
  // Frequency of reorthogonalisation.
  int nStabalize;
  // Whether to restore cosine projection and real local energy apprximation for weights
  // along back propagation path.
  bool fully_restore_weights, partially_restore_weights;
  // Note that the user will input a two dimensional list of indices (tuples) {(i,j)}
  // which we then convirt to the one dimensional array `elements'.
  std::vector<int> elements;

  SlaterDetOperations SDet;
  // Temporary storage for a given walker's density matrix.
  SPComplexMatrix walker_density_matrix;
  ComplexMatrix BPSlaterMat, T1, T2;
  std::vector<ComplexType> BPCoeff;
  ComplexType* T1ptr;

};
}
}

#endif
