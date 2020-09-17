//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file DiracDeterminantBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_BASE_H
#define QMCPLUSPLUS_DIRACDETERMINANT_BASE_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/TimerManager.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"

namespace qmcplusplus
{
class DiracDeterminantBase : public WaveFunctionComponent
{
public:
  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantBase(SPOSetPtr const spos, int first = 0)
      : UpdateTimer(*timer_manager.createTimer("DiracDeterminantBase::update", timer_level_fine)),
        RatioTimer(*timer_manager.createTimer("DiracDeterminantBase::ratio", timer_level_fine)),
        InverseTimer(*timer_manager.createTimer("DiracDeterminantBase::inverse", timer_level_fine)),
        BufferTimer(*timer_manager.createTimer("DiracDeterminantBase::buffer", timer_level_fine)),
        SPOVTimer(*timer_manager.createTimer("DiracDeterminantBase::spoval", timer_level_fine)),
        SPOVGLTimer(*timer_manager.createTimer("DiracDeterminantBase::spovgl", timer_level_fine)),
        Phi(spos),
        FirstIndex(first),
        LastIndex(first + spos->size()),
        NumOrbitals(spos->size()),
        NumPtcls(spos->size())
  {
    Optimizable  = Phi->isOptimizable();
    is_fermionic = true;
    ClassName    = "DiracDeterminantBase";
    registerTimers();
  }

  ///default destructor
  virtual ~DiracDeterminantBase() {}

  // copy constructor and assign operator disabled
  DiracDeterminantBase(const DiracDeterminantBase& s) = delete;
  DiracDeterminantBase& operator=(const DiracDeterminantBase& s) = delete;

  // get the SPO pointer
  inline SPOSetPtr getPhi() const { return Phi; }

  // get FirstIndex, Last Index
  inline int getFirstIndex() const { return FirstIndex; }
  inline int getLastIndex() const { return LastIndex; }

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  virtual void set(int first, int nel, int delay = 1){};

  ///set BF pointers
  virtual void setBF(BackflowTransformation* BFTrans) {}

  ///optimizations  are disabled
  virtual inline void checkInVariables(opt_variables_type& active) override { Phi->checkInVariables(active); }

  virtual inline void checkOutVariables(const opt_variables_type& active) override { Phi->checkOutVariables(active); }

  virtual void resetParameters(const opt_variables_type& active) override { Phi->resetParameters(active); }

  // To be removed with AoS
  void resetTargetParticleSet(ParticleSet& P) override final
  {
    Phi->resetTargetParticleSet(P);
    targetPtcl = &P;
  }

  inline void reportStatus(std::ostream& os) override final {}

  // expose CPU interfaces
  using WaveFunctionComponent::evaluateDerivatives;
  using WaveFunctionComponent::evaluateLog;
  using WaveFunctionComponent::mw_evaluateLog;
  using WaveFunctionComponent::recompute;

  using WaveFunctionComponent::copyFromBuffer;
  using WaveFunctionComponent::registerData;
  using WaveFunctionComponent::updateBuffer;

  using WaveFunctionComponent::acceptMove;
  using WaveFunctionComponent::completeUpdates;
  using WaveFunctionComponent::evalGrad;
  using WaveFunctionComponent::mw_accept_rejectMove;
  using WaveFunctionComponent::mw_calcRatio;
  using WaveFunctionComponent::mw_completeUpdates;
  using WaveFunctionComponent::mw_evalGrad;
  using WaveFunctionComponent::mw_ratioGrad;
  using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::ratioGrad;
  using WaveFunctionComponent::restore;

  using WaveFunctionComponent::evalGradSource;
  using WaveFunctionComponent::evaluateHessian;
  using WaveFunctionComponent::evaluateRatios;
  using WaveFunctionComponent::evaluateRatiosAlltoOne;
  using WaveFunctionComponent::mw_evaluateRatios;

  // used by DiracDeterminantWithBackflow
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& active,
                                   int offset,
                                   Matrix<RealType>& dlogpsi,
                                   Array<GradType, 3>& dG,
                                   Matrix<RealType>& dL)
  {
    APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::evaluateDerivatives");
  }

  // Stop makeClone
  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const override final
  {
    APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::makeClone");
    return 0;
  }

  virtual PsiValueType ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad) override
  {
    APP_ABORT("  DiracDeterminantBase::ratioGradWithSpins():  Implementation required\n");
    return 0.0;
  }
  virtual GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override
  {
    APP_ABORT("  DiracDeterminantBase::evalGradWithSpins():  Implementation required\n");
    return GradType();
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
  // expose GPU interfaces
  //using WaveFunctionComponent::recompute;
  using WaveFunctionComponent::addLog;
  using WaveFunctionComponent::reserve;
  //using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::addGradient;
  using WaveFunctionComponent::addRatio;
  using WaveFunctionComponent::calcGradient;
  using WaveFunctionComponent::calcRatio;
  using WaveFunctionComponent::det_lookahead;
  using WaveFunctionComponent::gradLapl;
  using WaveFunctionComponent::NLratios;
  using WaveFunctionComponent::update;
#endif

protected:
  /// Timers
  NewTimer &UpdateTimer, &RatioTimer, &InverseTimer, &BufferTimer, &SPOVTimer, &SPOVGLTimer;
  /// a set of single-particle orbitals used to fill in the  values of the matrix
  SPOSetPtr const Phi;
  ///index of the first particle with respect to the particle set
  int FirstIndex;
  ///index of the last particle with respect to the particle set
  int LastIndex;
  ///number of single-particle orbitals which belong to this Dirac determinant
  int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  int NumPtcls;
  /// targetPtcl pointer. YE: to be removed.
  ParticleSet* targetPtcl;

  /// register all the timers
  void registerTimers()
  {
    UpdateTimer.reset();
    RatioTimer.reset();
  }
};

} // namespace qmcplusplus
#endif
