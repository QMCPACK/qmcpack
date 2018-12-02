//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_BASE_H
#define QMCPLUSPLUS_DIRACDETERMINANT_BASE_H
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"

namespace qmcplusplus
{

class DiracDeterminantBase: public WaveFunctionComponent
{
public:
  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantBase(SPOSetPtr const &spos)
    : Phi(spos),
    UpdateTimer("DiracDeterminantBase::update",timer_level_fine),
    RatioTimer("DiracDeterminantBase::ratio",timer_level_fine),
    InverseTimer("DiracDeterminantBase::inverse",timer_level_fine),
    BufferTimer("DiracDeterminantBase::buffer",timer_level_fine),
    SPOVTimer("DiracDeterminantBase::spoval",timer_level_fine),
    SPOVGLTimer("DiracDeterminant::spovgl",timer_level_fine)
  {
    Optimizable = Phi->Optimizable;
    ClassName="DiracDeterminantBase";
    registerTimers();
  }

  ///default destructor
  virtual ~DiracDeterminantBase() { };

  // copy constructor and assign operator disabled
  DiracDeterminantBase(const DiracDeterminantBase& s) = delete;
  DiracDeterminantBase& operator=(const DiracDeterminantBase& s) = delete;

  // get the SPO pointer
  SPOSetPtr getPhi()
  {
    return Phi;
  };

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  virtual void set(int first, int nel, int delay=1) { };

  ///set BF pointers
  virtual void setBF(BackflowTransformation* BFTrans) {}

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

  // To be removed with AoS
  void resetTargetParticleSet(ParticleSet& P) final
  {
    Phi->resetTargetParticleSet(P);
    targetPtcl = &P;
  }

  inline void reportStatus(std::ostream& os) final { }

  using WaveFunctionComponent::evaluateDerivatives;
  using WaveFunctionComponent::evaluateLog;
  using WaveFunctionComponent::recompute;

  using WaveFunctionComponent::registerData;
  using WaveFunctionComponent::updateBuffer;
  using WaveFunctionComponent::copyFromBuffer;

  using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::ratioGrad;
  using WaveFunctionComponent::evalGrad;
  using WaveFunctionComponent::acceptMove;
  using WaveFunctionComponent::restore;
  using WaveFunctionComponent::completeUpdates;

  using WaveFunctionComponent::evalGradSource;
  using WaveFunctionComponent::evaluateRatios;
  using WaveFunctionComponent::evaluateRatiosAlltoOne;
  using WaveFunctionComponent::evaluateHessian;

  // Stop makeClone
  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const final
  {
    APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::makeClone");
    return 0;
  }

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  virtual DiracDeterminantBase* makeCopy(SPOSet* spo) const = 0;

#ifdef QMC_CUDA
  using WaveFunctionComponent::recompute;
  using WaveFunctionComponent::reserve;
  using WaveFunctionComponent::addLog;
  using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::addRatio;
  using WaveFunctionComponent::calcRatio;
  using WaveFunctionComponent::addGradient;
  using WaveFunctionComponent::calcGradient;
  using WaveFunctionComponent::update;
  using WaveFunctionComponent::gradLapl;
  using WaveFunctionComponent::NLratios;
#endif

  protected:
  /** flag to set the optimization mode */
  using WaveFunctionComponent::IsOptimizing;
  // Optimizable flag
  using WaveFunctionComponent::Optimizable;
  // Optimizable parameters
  using WaveFunctionComponent::myVars;
  /// Timers
  NewTimer UpdateTimer, RatioTimer, InverseTimer, BufferTimer, SPOVTimer, SPOVGLTimer;
  /// a set of single-particle orbitals used to fill in the  values of the matrix
  SPOSetPtr Phi;
  /// targetPtcl pointer. YE: to be removed.
  ParticleSet *targetPtcl;

  /// register all the timers
  void registerTimers()
  {
    UpdateTimer.reset();
    RatioTimer.reset();
    TimerManager.addTimer (&UpdateTimer);
    TimerManager.addTimer (&RatioTimer);
    TimerManager.addTimer (&InverseTimer);
    TimerManager.addTimer (&BufferTimer);
    TimerManager.addTimer (&SPOVTimer);
    TimerManager.addTimer (&SPOVGLTimer);
  }

};



}
#endif
