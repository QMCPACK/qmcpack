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
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include <strstream>
#include <stdexcept>

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/IteratorUtility.h"
#include "Concurrency/Info.hpp"

namespace qmcplusplus
{
typedef enum
{
  V_TIMER,
  VGL_TIMER,
  ACCEPT_TIMER,
  NL_TIMER,
  RECOMPUTE_TIMER,
  BUFFER_TIMER,
  DERIVS_TIMER,
  TIMER_SKIP
} TimerEnum;


TrialWaveFunction::TrialWaveFunction(Communicate* c)
    : MPIObjectBase(c),
      BufferCursor(0),
      BufferCursor_scalar(0),
      PhaseValue(0.0),
      PhaseDiff(0.0),
      LogValue(0.0),
      OneOverM(1.0)
{
  ClassName = "TrialWaveFunction";
  myName    = "psi0";
}

/** Destructor
*
*@warning Have not decided whether Z is cleaned up by TrialWaveFunction
*  or not. It will depend on I/O implementation.
*/
TrialWaveFunction::~TrialWaveFunction()
{
  delete_iter(Z.begin(), Z.end());
  //delete_iter(SPOSet.begin(),SPOSet.end());
  //delete_iter(myTimers.begin(),myTimers.end());
}

void TrialWaveFunction::resetTargetParticleSet(ParticleSet& P)
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->resetTargetParticleSet(P);
}

void TrialWaveFunction::startOptimization()
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->IsOptimizing = true;
}

void TrialWaveFunction::stopOptimization()
{
  for (int i = 0; i < Z.size(); i++)
  {
    Z[i]->finalizeOptimization();
    Z[i]->IsOptimizing = false;
  }
}

/** Takes owndership of aterm
 */
void TrialWaveFunction::addComponent(WaveFunctionComponent* aterm, std::string aname)
{
  Z.push_back(aterm);
  if (aterm->is_fermionic)
    app_log() << "  Added a fermionic WaveFunctionComponent " << aname << std::endl;

  std::vector<std::string> suffixes(7);
  suffixes[0] = "_V";
  suffixes[1] = "_VGL";
  suffixes[2] = "_accept";
  suffixes[3] = "_NLratio";
  suffixes[4] = "_recompute";
  suffixes[5] = "_buffer";
  suffixes[6] = "_derivs";
  for (int i = 0; i < suffixes.size(); i++)
  {
    std::string name = "WaveFunction::" + aname + suffixes[i];
    myTimers.push_back(TimerManager.createTimer(name));
  }
}


/** return log(|psi|)
*
* PhaseValue is the phase for the complex wave function
*/
TrialWaveFunction::RealType TrialWaveFunction::evaluateLog(ParticleSet& P)
{
  P.G = 0.0;
  P.L = 0.0;
  LogValueType logpsi(0.0);
  for (int i = 0, ii = RECOMPUTE_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    ScopedTimer local_timer(myTimers[ii]);
    logpsi += Z[i]->evaluateLog(P, P.G, P.L);
  }
  LogValue = std::real(logpsi);
  PhaseValue = std::imag(logpsi);
  return LogValue;
}

void TrialWaveFunction::flex_evaluateLog(const RefVector<TrialWaveFunction>& wf_list,
                                         const RefVector<ParticleSet>& p_list)
{
  if (wf_list.size() > 1)
  {
    constexpr RealType czero(0);
    const auto g_list(TrialWaveFunction::extractGRefList(wf_list));
    const auto l_list(TrialWaveFunction::extractLRefList(wf_list));

    int num_particles = p_list[0].get().getTotalNum();
    auto initGandL    = [num_particles, czero](TrialWaveFunction& twf, ParticleSet::ParticleGradient_t& grad,
                                            ParticleSet::ParticleLaplacian_t& lapl) {
      grad.resize(num_particles);
      lapl.resize(num_particles);
      grad           = czero;
      lapl           = czero;
      twf.LogValue   = czero;
      twf.PhaseValue = czero;
    };
    for (int iw = 0; iw < wf_list.size(); iw++)
      initGandL(wf_list[iw], g_list[iw], l_list[iw]);

    auto& wavefunction_components = wf_list[0].get().Z;
    const int num_wfc             = wf_list[0].get().Z.size();
    for (int i = 0, ii = RECOMPUTE_TIMER; i < num_wfc; ++i, ii += TIMER_SKIP)
    {
      ScopedTimer local_timer(wf_list[0].get().myTimers[ii]);
      const auto wfc_list(extractWFCRefList(wf_list, i));
      wavefunction_components[i]->mw_evaluateLog(convert_ref_to_ptr_list(wfc_list), convert_ref_to_ptr_list(p_list),
                                                 convert_ref_to_ptr_list(g_list), convert_ref_to_ptr_list(l_list));
      auto accumulateLogAndPhase = [](TrialWaveFunction& twf, WaveFunctionComponent& wfc) {
        twf.LogValue += std::real(wfc.LogValue);
        twf.PhaseValue += std::imag(wfc.LogValue);
      };
      for (int iw = 0; iw < wf_list.size(); iw++)
        accumulateLogAndPhase(wf_list[iw], wfc_list[iw]);
    }
    auto copyToP = [](ParticleSet& pset, TrialWaveFunction& twf) {
      pset.G = twf.G;
      pset.L = twf.L;
    };
    // Ye: temporal workaround to have P.G/L always defined.
    // remove when KineticEnergy use WF.G/L instead of P.G/L
    for (int iw = 0; iw < wf_list.size(); iw++)
      copyToP(p_list[iw], wf_list[iw]);
  }
  else if (wf_list.size() == 1)
  {
    wf_list[0].get().evaluateLog(p_list[0]);
    // Ye: temporal workaround to have WF.G/L always defined.
    // remove when KineticEnergy use WF.G/L instead of P.G/L
    wf_list[0].get().G = p_list[0].get().G;
    wf_list[0].get().L = p_list[0].get().L;
  }
}

void TrialWaveFunction::recompute(ParticleSet& P)
{
  std::vector<WaveFunctionComponent*>::iterator it(Z.begin());
  std::vector<WaveFunctionComponent*>::iterator it_end(Z.end());
  for (int ii = RECOMPUTE_TIMER; it != it_end; ++it, ii += TIMER_SKIP)
  {
    ScopedTimer local_timer(myTimers[ii]);
    (*it)->recompute(P);
  }
}

/** evaluate the log value of a many-body wave function
 * @param P input configuration containing N particles
 * @param recomputeall recompute all orbitals from scratch
 * @return the value of \f$ \log( \Pi_i \Psi_i) \f$  many-body wave function
 *
 * @if recomputeall == true
 *   all orbitals have "evaluateLog" called on them, including the non-optimized ones.
 * @else
 *   default value.  call evaluateLog only on optimizable orbitals.  OK if nonlocal pp's aren't used.
 *
 * To save time, logpsi, G, and L are only computed for orbitals that change over the course of the optimization.
 * It is assumed that the fixed components are stored elsewhere.  See evaluateDeltaLog(P,logpsi_fixed_r,logpsi_opt,fixedG,fixedL)
 * defined below.  Nonlocal pseudopotential evaluation requires temporary information like matrix inverses, so while
 * the logpsi, G, and L don't change, evaluateLog is called anyways to compute these auxiliary quantities from scratch.
 * logpsi, G, and L associated with these non-optimizable orbitals are discarded explicitly and with dummy variables. 
 */
TrialWaveFunction::RealType TrialWaveFunction::evaluateDeltaLog(ParticleSet& P, bool recomputeall)
{
  P.G = 0.0;
  P.L = 0.0;
  LogValueType logpsi(0.0);
  std::vector<WaveFunctionComponent*>::iterator it(Z.begin());
  std::vector<WaveFunctionComponent*>::iterator it_end(Z.end());
  int ii = RECOMPUTE_TIMER;
  for (; it != it_end; ++it, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    if ((*it)->Optimizable)
    {
      logpsi += (*it)->evaluateLog(P, P.G, P.L);
    }
    myTimers[ii]->stop();
  }
  LogValue = std::real(logpsi);
  PhaseValue = std::imag(logpsi);

  //In case we need to recompute orbitals, initialize dummy vectors for G and L.
  //evaluateLog dumps into these variables, and logPsi contribution is discarded.
  //Only called for non-optimizable orbitals.
  if (recomputeall)
  {
    ParticleSet::ParticleGradient_t dummyG(P.G);
    ParticleSet::ParticleLaplacian_t dummyL(P.L);

    it     = Z.begin();
    it_end = Z.end();

    for (; it != it_end; ++it)
    {
      if (!(*it)->Optimizable)
        (*it)->evaluateLog(P, dummyG,
                           dummyL); //update orbitals if its not flagged optimizable, AND recomputeall is true
    }
  }
  return LogValue;
}


/** evalaute the sum of log value of optimizable many-body wavefunctions
* @param P  input configuration containing N particles
* @param logpsi_fixed log(std::abs(psi)) of the invariant orbitals
* @param logpsi_opt log(std::abs(psi)) of the variable orbitals
* @param fixedG gradients of log(psi) of the fixed wave functions
* @param fixedL laplacians of log(psi) of the fixed wave functions
*
* This function is introduced for optimization only.
* fixedG and fixedL save the terms coming from the wave functions
* that are invarient during optimizations.
* It is expected that evaluateLog(P,false) is called later
* and the external object adds the varying G and L and the fixed terms.
*/
void TrialWaveFunction::evaluateDeltaLog(ParticleSet& P,
                                         RealType& logpsi_fixed_r,
                                         RealType& logpsi_opt_r,
                                         ParticleSet::ParticleGradient_t& fixedG,
                                         ParticleSet::ParticleLaplacian_t& fixedL)
{
  P.G    = 0.0;
  P.L    = 0.0;
  fixedL = 0.0;
  fixedG = 0.0;
  LogValueType logpsi_fixed(0.0);
  LogValueType logpsi_opt(0.0);
  std::vector<WaveFunctionComponent*>::iterator it(Z.begin());
  std::vector<WaveFunctionComponent*>::iterator it_end(Z.end());
  int ii = RECOMPUTE_TIMER;
  for (; it != it_end; ++it, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    if ((*it)->Optimizable)
      logpsi_opt += (*it)->evaluateLog(P, P.G, P.L);
    else
      logpsi_fixed += (*it)->evaluateLog(P, fixedG, fixedL);
    myTimers[ii]->stop();
  }
  P.G += fixedG;
  P.L += fixedL;
  convert(logpsi_fixed, logpsi_fixed_r);
  convert(logpsi_opt, logpsi_opt_r);
}

/*void TrialWaveFunction::evaluateHessian(ParticleSet & P, int iat, HessType& grad_grad_psi)
{
  std::vector<WaveFunctionComponent*>::iterator it(Z.begin());
  std::vector<WaveFunctionComponent*>::iterator it_end(Z.end());
  
  grad_grad_psi=0.0;
  
  for (; it!=it_end; ++it)
  {	
	  HessType tmp_hess;
	  (*it)->evaluateHessian(P, iat, tmp_hess);
	  grad_grad_psi+=tmp_hess;
  }
}*/

void TrialWaveFunction::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  std::vector<WaveFunctionComponent*>::iterator it(Z.begin());
  std::vector<WaveFunctionComponent*>::iterator it_end(Z.end());

  grad_grad_psi.resize(P.getTotalNum());

  for (int i = 0; i < Z.size(); i++)
  {
    HessVector_t tmp_hess(grad_grad_psi);
    tmp_hess = 0.0;
    Z[i]->evaluateHessian(P, tmp_hess);
    grad_grad_psi += tmp_hess;
    //  app_log()<<"TrialWavefunction::tmp_hess = "<<tmp_hess<< std::endl;
    //  app_log()<< std::endl<< std::endl;
  }
  // app_log()<<" TrialWavefunction::Hessian = "<<grad_grad_psi<< std::endl;
}

TrialWaveFunction::ValueType TrialWaveFunction::calcRatio(ParticleSet& P, int iat, ComputeType ct)
{
  PsiValueType r(1.0);
  for (int i = 0, ii = V_TIMER; i < Z.size(); i++, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    if (ct == ComputeType::ALL || (Z[i]->is_fermionic && ct == ComputeType::FERMIONIC) ||
        (!Z[i]->is_fermionic && ct == ComputeType::NONFERMIONIC))
      r *= Z[i]->ratio(P, iat);
    myTimers[ii]->stop();
  }
  return static_cast<ValueType>(r);
}

void TrialWaveFunction::flex_calcRatio(const RefVector<TrialWaveFunction>& wf_list,
                                       const RefVector<ParticleSet>& p_list,
                                       int iat,
                                       std::vector<PsiValueType>& ratios,
                                       ComputeType ct)
{
  const int num_wf = wf_list.size();
  ratios.resize(num_wf);
  std::fill(ratios.begin(), ratios.end(), PsiValueType(1));

  if (num_wf > 1)
  {
    const int num_wfc             = wf_list[0].get().Z.size();
    auto& wavefunction_components = wf_list[0].get().Z;

    std::vector<PsiValueType> ratios_z(num_wf);
    for (int i = 0, ii = V_TIMER; i < num_wfc; i++, ii += TIMER_SKIP)
    {
      if (ct == ComputeType::ALL || (wavefunction_components[i]->is_fermionic && ct == ComputeType::FERMIONIC) ||
          (!wavefunction_components[i]->is_fermionic && ct == ComputeType::NONFERMIONIC))
      {
        ScopedTimer local_timer(wf_list[0].get().get_timers()[ii]);
        const auto wfc_list(extractWFCRefList(wf_list, i));
        wavefunction_components[i]->mw_calcRatio(convert_ref_to_ptr_list(wfc_list), convert_ref_to_ptr_list(p_list),
                                                 iat, ratios_z);
        for (int iw = 0; iw < wf_list.size(); iw++)
          ratios[iw] *= ratios_z[iw];
      }
    }
    for (int iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].get().PhaseDiff = std::imag(std::arg(ratios[iw]));
  }
  else if (wf_list.size() == 1)
    ratios[0] = wf_list[0].get().calcRatio(p_list[0], iat);
}

TrialWaveFunction::GradType TrialWaveFunction::evalGrad(ParticleSet& P, int iat)
{
  GradType grad_iat;
  for (int i = 0, ii = VGL_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    grad_iat += Z[i]->evalGrad(P, iat);
    myTimers[ii]->stop();
  }
  return grad_iat;
}

void TrialWaveFunction::flex_evalGrad(const std::vector<std::reference_wrapper<TrialWaveFunction>>& wf_list,
                                      const std::vector<std::reference_wrapper<ParticleSet>>& p_list,
                                      int iat,
                                      std::vector<GradType>& grad_now)
{
  const int num_wf = wf_list.size();
  grad_now.resize(num_wf);
  std::fill(grad_now.begin(), grad_now.end(), GradType(0));


  if (num_wf > 1)
  {
    // Right now mw_evalGrad can only be called through an concrete instance of a wavefunctioncomponent
    const int num_wfc             = wf_list[0].get().Z.size();
    auto& wavefunction_components = wf_list[0].get().Z;

    std::vector<GradType> grad_now_z(num_wf);
    for (int i = 0, ii = VGL_TIMER; i < num_wfc; ++i, ii += TIMER_SKIP)
    {
      ScopedTimer localtimer(wf_list[0].get().get_timers()[ii]);
      const auto wfc_list(extractWFCRefList(wf_list, i));
      wavefunction_components[i]->mw_evalGrad(wfc_list, p_list, iat, grad_now_z);
      for (int iw = 0; iw < wf_list.size(); iw++)
        grad_now[iw] += grad_now_z[iw];
    }
  }
  else if (wf_list.size() == 1)
    grad_now[0] = wf_list[0].get().evalGrad(p_list[0], iat);
}


// Evaluates the gradient w.r.t. to the source of the Laplacian
// w.r.t. to the electrons of the wave function.
TrialWaveFunction::GradType TrialWaveFunction::evalGradSource(ParticleSet& P, ParticleSet& source, int iat)
{
  GradType grad_iat = GradType();
  for (int i = 0; i < Z.size(); ++i)
    grad_iat += Z[i]->evalGradSource(P, source, iat);
  return grad_iat;
}

TrialWaveFunction::GradType TrialWaveFunction::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat,
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
    TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad)
{
  GradType grad_iat = GradType();
  for (int dim = 0; dim < OHMMS_DIM; dim++)
    for (int i = 0; i < grad_grad[0].size(); i++)
    {
      grad_grad[dim][i] = GradType();
      lapl_grad[dim][i] = 0.0;
    }
  for (int i = 0; i < Z.size(); ++i)
    grad_iat += Z[i]->evalGradSource(P, source, iat, grad_grad, lapl_grad);
  return grad_iat;
}

TrialWaveFunction::ValueType TrialWaveFunction::calcRatioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  grad_iat = 0.0;
  PsiValueType r(1.0);
  for (int i = 0, ii = VGL_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    r *= Z[i]->ratioGrad(P, iat, grad_iat);
    myTimers[ii]->stop();
  }

  LogValueType logratio = convertValueToLog(r);
  PhaseDiff = std::imag(logratio);
  return static_cast<ValueType>(r);
}

void TrialWaveFunction::flex_calcRatioGrad(const RefVector<TrialWaveFunction>& wf_list,
                                       const RefVector<ParticleSet>& p_list,
                                       int iat,
                                       std::vector<PsiValueType>& ratios,
                                       std::vector<GradType>& grad_new)
{
  const int num_wf = wf_list.size();
  grad_new.resize(num_wf);
  std::fill(grad_new.begin(), grad_new.end(), GradType(0));
  grad_new.resize(num_wf);
  std::fill(ratios.begin(), ratios.end(), PsiValueType(1));

  if (wf_list.size() > 1)
  {
    const int num_wfc             = wf_list[0].get().Z.size();
    auto& wavefunction_components = wf_list[0].get().Z;

    std::vector<PsiValueType> ratios_z(wf_list.size());
    for (int i = 0, ii = VGL_TIMER; i < num_wfc; ++i, ii += TIMER_SKIP)
    {
      ScopedTimer localtimer(wf_list[0].get().get_timers()[ii]);
      //ScopedTimer local_timer(myTimers[ii]);
      const auto wfc_list(extractWFCRefList(wf_list, i));
      wavefunction_components[i]->mw_ratioGrad(convert_ref_to_ptr_list(wfc_list), convert_ref_to_ptr_list(p_list), iat,
                                               ratios_z, grad_new);
      for (int iw = 0; iw < wf_list.size(); iw++)
        ratios[iw] *= ratios_z[iw];
    }
    for (int iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].get().PhaseDiff = std::imag(std::arg(ratios[iw]));
  }
  else if (wf_list.size() == 1)
    ratios[0] = wf_list[0].get().calcRatioGrad(p_list[0], iat, grad_new[0]);
}

void TrialWaveFunction::printGL(ParticleSet::ParticleGradient_t& G,
                                ParticleSet::ParticleLaplacian_t& L,
                                std::string tag)
{
  std::ostringstream o;
  o << "---  reporting " << tag << std::endl << "  ---" << std::endl;
  for (int iat = 0; iat < L.size(); iat++)
    o << "index: " << std::fixed << iat << std::scientific << "   G: " << G[iat][0] << "  " << G[iat][1] << "  "
      << G[iat][2] << "   L: " << L[iat] << std::endl;
  o << "---  end  ---" << std::endl;
  std::cout << o.str();
}

/** restore to the original state
 * @param iat index of the particle with a trial move
 *
 * The proposed move of the iath particle is rejected.
 * All the temporary data should be restored to the state prior to the move.
 */
void TrialWaveFunction::rejectMove(int iat)
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->restore(iat);
  PhaseDiff = 0;
}

void TrialWaveFunction::flex_rejectMove(const RefVector<TrialWaveFunction>& wf_list, int iat)
{
  if (wf_list.size() > 1)
  {
    const int num_wfc             = wf_list[0].get().Z.size();
    auto& wavefunction_components = wf_list[0].get().Z;

    for (int i = 0; i < num_wfc; i++)
    {
      const auto wfc_list(extractWFCRefList(wf_list, i));
      wavefunction_components[i]->mw_restore(convert_ref_to_ptr_list(wfc_list), iat);
    }
  }
  else if (wf_list.size() == 1)
    wf_list[0].get().rejectMove(iat);
}

/** update the state with the new data
 * @param P ParticleSet
 * @param iat index of the particle with a trial move
 *
 * The proposed move of the iath particle is accepted.
 * All the temporary data should be incorporated so that the next move is valid.
 */
void TrialWaveFunction::acceptMove(ParticleSet& P, int iat)
{
  for (int i = 0, ii = ACCEPT_TIMER; i < Z.size(); i++, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->acceptMove(P, iat);
    myTimers[ii]->stop();
  }
  PhaseValue += PhaseDiff;
  PhaseDiff = 0.0;
  LogValue  = 0;
  for (int i = 0; i < Z.size(); i++)
    LogValue += std::real(Z[i]->LogValue);
}

void TrialWaveFunction::flex_acceptMove(const RefVector<TrialWaveFunction>& wf_list,
                                        const RefVector<ParticleSet>& p_list,
                                        int iat)
{
  if (wf_list.size() > 1)
  {
    const int num_wfc             = wf_list[0].get().Z.size();
    auto& wavefunction_components = wf_list[0].get().Z;

    for (int iw = 0; iw < wf_list.size(); iw++)
    {
      wf_list[iw].get().LogValue   = 0;
      wf_list[iw].get().PhaseValue = 0;
    }

    for (int i = 0, ii = ACCEPT_TIMER; i < num_wfc; i++, ii += TIMER_SKIP)
    {
      ScopedTimer localtimer(wf_list[0].get().get_timers()[ii]);
      const auto wfc_list(extractWFCRefList(wf_list, i));
      wavefunction_components[i]->mw_acceptMove(convert_ref_to_ptr_list(wfc_list), convert_ref_to_ptr_list(p_list),
                                                iat);
      for (int iw = 0; iw < wf_list.size(); iw++)
      {
        wf_list[iw].get().LogValue += std::real(wfc_list[iw].get().LogValue);
        wf_list[iw].get().PhaseValue += std::imag(wfc_list[iw].get().LogValue);
      }
    }
  }
  else if (wf_list.size() == 1)
    wf_list[0].get().acceptMove(p_list[0], iat);
}

void TrialWaveFunction::completeUpdates()
{
  for (int i = 0, ii = ACCEPT_TIMER; i < Z.size(); i++, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->completeUpdates();
    myTimers[ii]->stop();
  }
}

void TrialWaveFunction::flex_completeUpdates(const std::vector<TrialWaveFunction*>& WF_list) const
{
  if (WF_list.size() > 1)
  {
    for (int i = 0, ii = ACCEPT_TIMER; i < Z.size(); i++, ii += TIMER_SKIP)
    {
      ScopedTimer local_timer(myTimers[ii]);
      std::vector<WaveFunctionComponent*> WFC_list(extractWFCPtrList(WF_list, i));
      Z[i]->mw_completeUpdates(WFC_list);
    }
  }
  else if (WF_list.size() == 1)
    WF_list[0]->completeUpdates();
}

void TrialWaveFunction::checkInVariables(opt_variables_type& active)
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->checkInVariables(active);
}

void TrialWaveFunction::checkOutVariables(const opt_variables_type& active)
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->checkOutVariables(active);
}

void TrialWaveFunction::resetParameters(const opt_variables_type& active)
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->resetParameters(active);
}

void TrialWaveFunction::reportStatus(std::ostream& os)
{
  for (int i = 0; i < Z.size(); i++)
    Z[i]->reportStatus(os);
}

void TrialWaveFunction::getLogs(std::vector<RealType>& lvals)
{
  lvals.resize(Z.size(), 0);
  for (int i = 0; i < Z.size(); i++)
  {
    lvals[i] = std::real(Z[i]->LogValue);
  }
}

void TrialWaveFunction::getPhases(std::vector<RealType>& pvals)
{
  pvals.resize(Z.size(), 0);
  for (int i = 0; i < Z.size(); i++)
  {
    pvals[i] = std::imag(Z[i]->LogValue);
  }
}

void TrialWaveFunction::registerData(ParticleSet& P, WFBufferType& buf)
{
  //save the current position
  BufferCursor        = buf.current();
  BufferCursor_scalar = buf.current_scalar();
  for (int i = 0, ii = BUFFER_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->registerData(P, buf);
    myTimers[ii]->stop();
  }
  buf.add(PhaseValue);
  buf.add(LogValue);
}

void TrialWaveFunction::flex_registerData(const UPtrVector<TrialWaveFunction>& wf_list,
                                          const UPtrVector<ParticleSet>& P_list,
                                          const RefVector<WFBufferType>& buf_list)
{
  auto setBufferCursors = [](TrialWaveFunction& twf, WFBufferType& wb) {
    twf.BufferCursor        = wb.current();
    twf.BufferCursor_scalar = wb.current_scalar();
  };
  for (int iw = 0; iw < wf_list.size(); iw++)
  {
    setBufferCursors(*(wf_list[iw]), buf_list[iw]);
  }

  const int num_wfc             = wf_list[0]->Z.size();
  auto& wavefunction_components = wf_list[0]->Z;
  for (int i = 0, ii = BUFFER_TIMER; i < num_wfc; i++, ii += TIMER_SKIP)
  {
    ScopedTimer local_timer(wf_list[0]->myTimers[ii]);
    std::vector<WaveFunctionComponent*> wfc_list(extractWFCPtrList(wf_list, i));

    wavefunction_components[i]->mw_registerData(wfc_list, convertUPtrToPtrVector(P_list),
                                                convert_ref_to_ptr_list(buf_list));
  }

  auto addPhaseAndLog = [](WFBufferType& wfb, TrialWaveFunction& twf) {
    wfb.add(twf.PhaseValue);
    wfb.add(twf.LogValue);
  };
  for (int iw = 0; iw < wf_list.size(); iw++)
    addPhaseAndLog(buf_list[iw], *(wf_list[iw]));
}

void TrialWaveFunction::debugOnlyCheckBuffer(WFBufferType& buffer)
{
#ifndef NDEBUG
  if (buffer.size() < buffer.current() + buffer.current_scalar() * sizeof(FullPrecRealType))
  {
    std::strstream assert_message;
    assert_message << "On thread:" << Concurrency::getThreadId<>() << "  buf_list[iw].get().size():" << buffer.size()
                   << " < buf_list[iw].get().current():" << buffer.current()
                   << " + buf.current_scalar():" << buffer.current_scalar()
                   << " * sizeof(FullPrecRealType):" << sizeof(FullPrecRealType) << '\n';
    throw std::runtime_error(assert_message.str());
  }
#endif
}

TrialWaveFunction::RealType TrialWaveFunction::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  P.G = 0.0;
  P.L = 0.0;
  buf.rewind(BufferCursor, BufferCursor_scalar);
  LogValueType logpsi(0.0);
  for (int i = 0, ii = BUFFER_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    logpsi += Z[i]->updateBuffer(P, buf, fromscratch);
    myTimers[ii]->stop();
  }

  LogValue = std::real(logpsi);
  PhaseValue = std::imag(logpsi);
  //printGL(P.G,P.L);
  buf.put(PhaseValue);
  buf.put(LogValue);
  // Ye: temperal added check, to be removed
  debugOnlyCheckBuffer(buf);
  return LogValue;
}

/** updates "buffer" for multiple wavefunctions SIDEFFECT: reduces wfc values and phases to TWFs
 *
 *  NO UNIT TEST
 */
void TrialWaveFunction::flex_updateBuffer(const RefVector<TrialWaveFunction>& wf_list,
                                          const RefVector<ParticleSet>& p_list,
                                          const RefVector<WFBufferType>& buf_list,
                                          bool fromscratch)
{
  for (int iw = 0; iw < wf_list.size(); iw++)
  {
    constexpr RealType czero(0);

    p_list[iw].get().G           = czero; // Ye: remove when updateBuffer of all the WFC uses WF.G/L
    p_list[iw].get().L           = czero; // Ye: remove when updateBuffer of all the WFC uses WF.G/L
    wf_list[iw].get().LogValue   = czero;
    wf_list[iw].get().PhaseValue = czero;
    buf_list[iw].get().rewind(wf_list[iw].get().BufferCursor, wf_list[iw].get().BufferCursor_scalar);
  }

  auto& wavefunction_components = wf_list[0].get().Z;
  const int num_wfc             = wavefunction_components.size();

  for (int i = 0, ii = BUFFER_TIMER; i < num_wfc; ++i, ii += TIMER_SKIP)
  {
    ScopedTimer local_timer(wf_list[0].get().myTimers[ii]);
    const auto wfc_list(extractWFCRefList(wf_list, i));
    wavefunction_components[i]->mw_updateBuffer(wfc_list, p_list, buf_list, fromscratch);
    for (int iw = 0; iw < wf_list.size(); iw++)
    {
      wf_list[iw].get().LogValue += std::real(wfc_list[iw].get().LogValue);
      wf_list[iw].get().PhaseValue +=  std::imag(wfc_list[iw].get().LogValue);
    }
  }

  for (int iw = 0; iw < wf_list.size(); iw++)
  {
    buf_list[iw].get().put(wf_list[iw].get().PhaseValue);
    buf_list[iw].get().put(wf_list[iw].get().LogValue);
    debugOnlyCheckBuffer(buf_list[iw]);
  }
}


void TrialWaveFunction::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  buf.rewind(BufferCursor, BufferCursor_scalar);
  for (int i = 0, ii = BUFFER_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->copyFromBuffer(P, buf);
    myTimers[ii]->stop();
  }
  //get the gradients and laplacians from the buffer
  buf.get(PhaseValue);
  buf.get(LogValue);
  debugOnlyCheckBuffer(buf);
}

void TrialWaveFunction::flex_copyFromBuffer(const RefVector<TrialWaveFunction>& wf_list,
                                            const RefVector<ParticleSet>& p_list,
                                            const RefVector<WFBufferType>& buf_list) const
{
  auto rewind = [](WFBufferType& buf, TrialWaveFunction& twf) {
    buf.rewind(twf.BufferCursor, twf.BufferCursor_scalar);
  };
  for (int iw = 0; iw < wf_list.size(); iw++)
    rewind(buf_list[iw], wf_list[iw]);

  for (int i = 0, ii = BUFFER_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    ScopedTimer local_timer(myTimers[ii]);
    const auto wfc_list(extractWFCRefList(wf_list, i));
    Z[i]->mw_copyFromBuffer(wfc_list, p_list, buf_list);
  }

  auto bufGetTwfLog = [](WFBufferType& buf, TrialWaveFunction& twf) {
    buf.get(twf.PhaseValue);
    buf.get(twf.LogValue);
    debugOnlyCheckBuffer(buf);
  };
  for (int iw = 0; iw < wf_list.size(); iw++)
    bufGetTwfLog(buf_list[iw], wf_list[iw]);
}

void TrialWaveFunction::evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  assert(VP.getTotalNum() == ratios.size());
  std::vector<ValueType> t(ratios.size());
  std::fill(ratios.begin(), ratios.end(), 1.0);
  for (int i = 0, ii = NL_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->evaluateRatios(VP, t);
    for (int j = 0; j < ratios.size(); ++j)
      ratios[j] *= t[j];
    myTimers[ii]->stop();
  }
}

void TrialWaveFunction::evaluateDerivRatios(VirtualParticleSet& VP,
                                            const opt_variables_type& optvars,
                                            std::vector<ValueType>& ratios,
                                            Matrix<ValueType>& dratio)
{
#if defined(QMC_COMPLEX)
  APP_ABORT("TrialWaveFunction::evaluateDerivRatios not available for complex wavefunctions");
#else
  std::fill(ratios.begin(), ratios.end(), 1.0);
  std::vector<ValueType> t(ratios.size());
  for (int i = 0; i < Z.size(); ++i)
  {
    Z[i]->evaluateDerivRatios(VP, optvars, t, dratio);
    for (int j = 0; j < ratios.size(); ++j)
      ratios[j] *= t[j];
  }
#endif
}

bool TrialWaveFunction::put(xmlNodePtr cur) { return true; }

void TrialWaveFunction::reset() {}

TrialWaveFunction* TrialWaveFunction::makeClone(ParticleSet& tqp) const
{
  TrialWaveFunction* myclone   = new TrialWaveFunction(myComm);
  myclone->BufferCursor        = BufferCursor;
  myclone->BufferCursor_scalar = BufferCursor_scalar;
  for (int i = 0; i < Z.size(); ++i)
    myclone->addComponent(Z[i]->makeClone(tqp), Z[i]->ClassName);
  myclone->OneOverM = OneOverM;
  return myclone;
}

/** evaluate derivatives of KE wrt optimizable varibles
 *
 * @todo WaveFunctionComponent objects should take the mass into account.
 */
void TrialWaveFunction::evaluateDerivatives(ParticleSet& P,
                                            const opt_variables_type& optvars,
                                            std::vector<ValueType>& dlogpsi,
                                            std::vector<ValueType>& dhpsioverpsi,
                                            bool project)
{
  //     // First, zero out derivatives
  //  This should only be done for some variables.
  //     for (int j=0; j<dlogpsi.size(); j++)
  //       dlogpsi[j] = dhpsioverpsi[j] = 0.0;
  for (int i = 0; i < Z.size(); i++)
  {
    if (Z[i]->dPsi)
      (Z[i]->dPsi)->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
    else
      Z[i]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }
  //orbitals do not know about mass of particle.
  for (int i = 0; i < dhpsioverpsi.size(); i++)
    dhpsioverpsi[i] *= OneOverM;

  if (project)
  {
    for (int i = 0; i < Z.size(); i++)
    {
      if (Z[i]->dPsi)
        (Z[i]->dPsi)->multiplyDerivsByOrbR(dlogpsi);
      else
        Z[i]->multiplyDerivsByOrbR(dlogpsi);
    }
    RealType psiValue = std::exp(-LogValue) * std::cos(PhaseValue);
    for (int i = 0; i < dlogpsi.size(); i++)
      dlogpsi[i] *= psiValue;
  }
}

void TrialWaveFunction::evaluateDerivativesWF(ParticleSet& P,
                                              const opt_variables_type& optvars,
                                              std::vector<ValueType>& dlogpsi)
{
  for (int i = 0; i < Z.size(); i++)
  {
    if (Z[i]->dPsi)
      (Z[i]->dPsi)->evaluateDerivativesWF(P, optvars, dlogpsi);
    else
      Z[i]->evaluateDerivativesWF(P, optvars, dlogpsi);
  }
}

void TrialWaveFunction::evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
                                                std::vector<ValueType>& dgradlogpsi)
{
  for (int i = 0; i < Z.size(); i++)
  {
    Z[i]->evaluateGradDerivatives(G_in, dgradlogpsi);
  }
}

TrialWaveFunction::RealType TrialWaveFunction::KECorrection() const
{
  RealType sum = 0.0;
  for (int i = 0; i < Z.size(); ++i)
    sum += Z[i]->KECorrection();
  return sum;
}

void TrialWaveFunction::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  std::fill(ratios.begin(), ratios.end(), 1.0);
  std::vector<ValueType> t(ratios.size());
  for (int i = 0, ii = V_TIMER; i < Z.size(); ++i, ii += TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->evaluateRatiosAlltoOne(P, t);
    for (int j = 0; j < t.size(); ++j)
      ratios[j] *= t[j];
    myTimers[ii]->stop();
  }
}

std::vector<std::reference_wrapper<WaveFunctionComponent>> TrialWaveFunction::extractWFCRefList(
    const std::vector<std::reference_wrapper<TrialWaveFunction>>& wf_list,
    int id)
{
  std::vector<std::reference_wrapper<WaveFunctionComponent>> wfc_list;
  wfc_list.reserve(wf_list.size());
  for (auto wf : wf_list)
    wfc_list.push_back(*(wf.get().Z[id]));
  return wfc_list;
}

std::vector<WaveFunctionComponent*> TrialWaveFunction::extractWFCPtrList(const UPtrVector<TrialWaveFunction>& WF_list,
                                                                         int id)
{
  std::vector<WaveFunctionComponent*> WFC_list;
  WFC_list.reserve(WF_list.size());
  for (auto& WF : WF_list)
    WFC_list.push_back(WF->Z[id]);
  return WFC_list;
}

std::vector<WaveFunctionComponent*> TrialWaveFunction::extractWFCPtrList(const std::vector<TrialWaveFunction*>& WF_list,
                                                                         int id) const
{
  std::vector<WaveFunctionComponent*> WFC_list;
  WFC_list.reserve(WF_list.size());
  for (auto WF : WF_list)
    WFC_list.push_back(WF->Z[id]);
  return WFC_list;
}

std::vector<ParticleSet::ParticleGradient_t*> TrialWaveFunction::extractGPtrList(
    const std::vector<TrialWaveFunction*>& WF_list) const
{
  std::vector<ParticleSet::ParticleGradient_t*> G_list;
  for (auto WF : WF_list)
    G_list.push_back(&(WF->G));
  return G_list;
}

std::vector<ParticleSet::ParticleLaplacian_t*> TrialWaveFunction::extractLPtrList(
    const std::vector<TrialWaveFunction*>& WF_list) const
{
  std::vector<ParticleSet::ParticleLaplacian_t*> L_list;
  for (auto WF : WF_list)
    L_list.push_back(&(WF->L));
  return L_list;
}

RefVector<ParticleSet::ParticleGradient_t> TrialWaveFunction::extractGRefList(
    const RefVector<TrialWaveFunction>& wf_list)
{
  RefVector<ParticleSet::ParticleGradient_t> g_list;
  for (TrialWaveFunction& wf : wf_list)
    g_list.push_back(wf.G);
  return g_list;
}

RefVector<ParticleSet::ParticleLaplacian_t> TrialWaveFunction::extractLRefList(
    const RefVector<TrialWaveFunction>& wf_list)
{
  RefVector<ParticleSet::ParticleLaplacian_t> l_list;
  for (TrialWaveFunction& wf : wf_list)
    l_list.push_back(wf.L);
  return l_list;
}


} // namespace qmcplusplus
