//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantBatched.h
 * @brief Declaration of DiracDeterminantBatched a specialization of 
 *        for batched walkers of DiracDeterminant with with a 
 *        S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_BATCHED_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_BATCHED_H
#include <typeinfo>
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/SPOSetBatched.h"
#include "QMCWaveFunctions/Fermion/determinant_update.h"
#include "Numerics/CUDA/cuda_inverse.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantEval.h"

namespace qmcplusplus
{

class DiracDeterminantBatched : public DiracDeterminantBase,
				public DiracDeterminantEval<Batching::BATCHED>
{
public:
  using SSTA = SPOSetTypeAliases;
  typedef SSTA::IndexVector_t IndexVector_t;
  typedef SSTA::ValueVector_t ValueVector_t;
  typedef SSTA::ValueMatrix_t ValueMatrix_t;
  typedef SSTA::GradVector_t  GradVector_t;
  typedef SSTA::GradMatrix_t  GradMatrix_t;
  typedef ParticleSet::Walker_t     Walker_t;

  using CudaValueType =  QMCT::CudaValueType;
  
  using SPOSetPtr = SPOSetBatched*;
  SPOSetPtr Phi; //Out local Phi_

  SPOSetPtr getPhi() { return Phi; }
  
  DiracDeterminantBatched(SPOSet* const &spos, int first=0);
  DiracDeterminantBatched(const DiracDeterminantBatched& s) = delete;


public:

  void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  DiracDeterminantBatched* makeCopy(SPOSet* spo) const;

    ///optimizations  are disabled
  virtual inline void checkInVariables(opt_variables_type& active)
  {
    Phi->checkInVariables(active);
    Phi->checkInVariables(myVars);
  }

  virtual inline void checkOutVariables(const opt_variables_type& active)
  {
    Phi->checkOutVariables(active);
    myVars.clear();
    myVars.insertFrom(Phi->myVars);
    myVars.getIndex(active);
  }

  virtual void resetParameters(const opt_variables_type& active)
  {
    Phi->resetParameters(active);
    for(int i=0; i<myVars.size(); ++i)
    {
      int ii=myVars.Index[i];
      if(ii>=0)
        myVars[i]= active[ii];
    }
  }

  virtual void resetTargetParticleSet(ParticleSet& P)
  {
    Phi->resetTargetParticleSet(P);
    targetPtcl = &P;
  }

  void update (std::vector<Walker_t*> &walkers, int iat);
  void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList);

  void recompute (MCWalkerConfiguration &W, bool firstTime);

  void addLog (MCWalkerConfiguration &W, std::vector<RealType> &logPsi);

  void addGradient(MCWalkerConfiguration &W, int iat,
                   std::vector<GradType> &grad);

  void calcGradient(MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);

  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios);

  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad);

  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
              std::vector<ValueType> &lapl);
  void calcRatio (MCWalkerConfiguration &W, int iat,
                  std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
                  std::vector<ValueType> &lapl);
  void addRatio (MCWalkerConfiguration &W, int iat,
                 std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
                 std::vector<ValueType> &lapl);

  void ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
              std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
              std::vector<GradType>  &grad, std::vector<ValueType> &lapl);

  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);

  void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);

  void NLratios_CPU (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                     std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);





};
}
#endif // QMCPLUSPLUS_DIRAC_DETERMINANT_BATCHED_H
