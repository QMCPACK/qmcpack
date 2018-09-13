//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   initially refactored from DiracDeterminantBase.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantSingle.h
 * @brief Declaration of DiracDeterminantSincle specialization of 
 *        DiracDeterminant for single walkers
 *        S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_SINGLE_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_SINGLE_H

#include <typeinfo>
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantEval.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSetBatched.h"
#include "QMCWaveFunctions/SPOSetSingle.h"
#include "QMCWaveFunctions/Fermion/determinant_update.h"
#include "Numerics/CUDA/cuda_inverse.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/Batching.h"

namespace qmcplusplus
{

class DiracDeterminantSingle : public DiracDeterminantBase,
                               public DiracDeterminantEval<Batching::SINGLE>
{
public:
  
  using SPOSetPtr = SPOSetSingle*;
  using Walker_t =  ParticleSet::Walker_t;
  SPOSetPtr Phi; //Out local Phi_

  DiracDeterminantSingle(SPOSetPtr const &spos, int first=0);
  //DiracDeterminantSingle(const DiracDeterminantSingle& s) = delete;

  
  virtual SPOSet* getPhi() { return dynamic_cast<SPOSet*>(Phi); }

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

  ValueType
  ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  void updateAfterSweep(ParticleSet& P,
			ParticleSet::ParticleGradient_t& G,
			ParticleSet::ParticleLaplacian_t& L);
  
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

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
  ValueType ratio(ParticleSet& P, int iat);
  
public:
  DiracDeterminantSingle* makeCopy(SPOSetPtr spo) const;
  DiracDeterminantBase* makeCopy(SPOSet* spo) const;

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  GradType evalGrad(ParticleSet& P, int iat);

  void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  GradType evalGradSource(ParticleSet& P, ParticleSet& source,
			  int iat);

  GradType evalGradSourcep(ParticleSet& P, ParticleSet& source,int iat,
			   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
			   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  GradType evalGradSource(ParticleSet& P, ParticleSet& source,int iat,
			  TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
			  TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi);
    
  void recompute(ParticleSet& P);
    
  void NLratios_CPU (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                     std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);
};

}
#endif // QMCPLUSPLUS_DIRAC_DETERMINANT_SINGLE_H
