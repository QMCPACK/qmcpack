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
   *@param spos the single-particle orbital set.
   *  shared_ptr is intended neither for sharing between spin up and down electrons nor for sharing between clones.
   *  The sharing aspect is for the determinants used by the the multi-determinant slow implementation.
   *@param first index of the first particle
   */
  DiracDeterminantBase(const std::string& class_name, std::shared_ptr<SPOSet>&& spos, int first = 0)
      : WaveFunctionComponent(class_name),
        UpdateTimer(*timer_manager.createTimer(class_name + "::update", timer_level_fine)),
        RatioTimer(*timer_manager.createTimer(class_name + "::ratio", timer_level_fine)),
        InverseTimer(*timer_manager.createTimer(class_name + "::inverse", timer_level_fine)),
        BufferTimer(*timer_manager.createTimer(class_name + "::buffer", timer_level_fine)),
        SPOVTimer(*timer_manager.createTimer(class_name + "::spoval", timer_level_fine)),
        SPOVGLTimer(*timer_manager.createTimer(class_name + "::spovgl", timer_level_fine)),
        Phi(std::move(spos)),
        FirstIndex(first),
        LastIndex(first + Phi->size()),
        NumOrbitals(Phi->size()),
        NumPtcls(Phi->size())
  {
    Optimizable  = Phi->isOptimizable();
    is_fermionic = true;
    registerTimers();
  }

  ///default destructor
  virtual ~DiracDeterminantBase() {}

  // copy constructor and assign operator disabled
  DiracDeterminantBase(const DiracDeterminantBase& s) = delete;
  DiracDeterminantBase& operator=(const DiracDeterminantBase& s) = delete;

  // get the SPO pointer
  inline SPOSetPtr getPhi() const { return Phi.get(); }

  // get FirstIndex, Last Index
  inline int getFirstIndex() const { return FirstIndex; }
  inline int getLastIndex() const { return LastIndex; }

#ifndef NDEBUG
  virtual ValueMatrix_t& getPsiMinv() { return dummy_vmt; }
#endif

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

  inline void reportStatus(std::ostream& os) override final {}

  // expose CPU interfaces
  using WaveFunctionComponent::evaluateDerivatives;
  using WaveFunctionComponent::evaluateGL;
  using WaveFunctionComponent::evaluateLog;
  using WaveFunctionComponent::mw_evaluateGL;
  using WaveFunctionComponent::mw_evaluateLog;
  using WaveFunctionComponent::recompute;

  using WaveFunctionComponent::copyFromBuffer;
  using WaveFunctionComponent::registerData;
  using WaveFunctionComponent::updateBuffer;

  using WaveFunctionComponent::createResource;
  using WaveFunctionComponent::acquireResource;
  using WaveFunctionComponent::releaseResource;

  using WaveFunctionComponent::acceptMove;
  using WaveFunctionComponent::completeUpdates;
  using WaveFunctionComponent::evalGrad;
  using WaveFunctionComponent::mw_accept_rejectMove;
  using WaveFunctionComponent::mw_calcRatio;
  using WaveFunctionComponent::mw_completeUpdates;
  using WaveFunctionComponent::mw_evalGrad;
  using WaveFunctionComponent::mw_ratioGrad;
  using WaveFunctionComponent::mw_ratioGradAsync;
  using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::ratioGrad;
  using WaveFunctionComponent::ratioGradAsync;
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
  virtual DiracDeterminantBase* makeCopy(std::shared_ptr<SPOSet>&& spo) const = 0;

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
  /** a set of single-particle orbitals used to fill in the  values of the matrix
   *  shared_ptr is intended neither for sharing between spin up and down electrons nor for sharing between clones.
   *  The sharing aspect is for the determinants used by the the multi-determinant slow implementation.
   */
  const std::shared_ptr<SPOSet> Phi;
  ///index of the first particle with respect to the particle set
  int FirstIndex;
  ///index of the last particle with respect to the particle set
  int LastIndex;
  ///number of single-particle orbitals which belong to this Dirac determinant
  int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  int NumPtcls;

#ifndef NDEBUG
  ValueMatrix_t dummy_vmt;
#endif

  bool checkG(const GradType& g)
  {
    auto g_mag = std::abs(dot(g, g));
    if (std::isnan(g_mag))
      throw std::runtime_error("gradient of NaN");
    if (std::isinf(g_mag))
      throw std::runtime_error("gradient of inf");
    if (g_mag < std::abs(std::numeric_limits<RealType>::epsilon()))
    {
      std::cerr << "evalGrad gradient is " << g[0] << ' ' << g[1] << ' ' << g[2] << '\n';
      throw std::runtime_error("gradient of zero");
    }
    return true;
  }

  /// register all the timers
  void registerTimers()
  {
    UpdateTimer.reset();
    RatioTimer.reset();
  }
};

} // namespace qmcplusplus
#endif
