//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{

typedef enum { V_TIMER, VGL_TIMER, ACCEPT_TIMER, NL_TIMER,
               RECOMPUTE_TIMER, DERIVS_TIMER, TIMER_SKIP
             } TimerEnum;

////////////////////////////////
// Vectorized member fuctions //
///////////////////////////////
void
TrialWaveFunction::freeGPUmem()
{
  for (int i=Z.size()-1; i>=0; i--)
    Z[i]->freeGPUmem();
}

void
TrialWaveFunction::recompute
(MCWalkerConfiguration &W, bool firstTime)
{
  for (int i=0,ii=RECOMPUTE_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->recompute(W, firstTime);
    myTimers[ii]->stop();
  }
}


void
TrialWaveFunction::reserve
(PointerPool<gpu::device_vector<CudaValueType> > &pool,
 bool onlyOptimizable)
{
  for(int i=0; i<Z.size(); i++)
    if (!onlyOptimizable || Z[i]->Optimizable)
      Z[i]->reserve(pool);
}

void
TrialWaveFunction::evaluateLog (MCWalkerConfiguration &W,
                                std::vector<RealType> &logPsi)
{
  for (int iw=0; iw<logPsi.size(); iw++)
    logPsi[iw] = RealType();
  for (int i=0; i<Z.size(); i++)
    Z[i]->addLog (W, logPsi);
}

void
TrialWaveFunction::calcGradient (MCWalkerConfiguration &W, int iat,
                                 std::vector<GradType> &grad)
{
  for (int iw=0; iw<grad.size(); iw++)
    grad[iw] = GradType();
  for(int i=0; i<Z.size(); i++)
    Z[i]->calcGradient(W, iat, grad);
}

void
TrialWaveFunction::addGradient (MCWalkerConfiguration &W, int iat,
                                std::vector<GradType> &grad)
{
  for(int i=0; i<Z.size()-1; i++)
    Z[i]->addGradient(W, iat, grad);
  Z[Z.size()-1]->addGradient(W,iat,grad);
}

void
TrialWaveFunction::getGradient (MCWalkerConfiguration &W, int iat,
                                std::vector<GradType> &grad)
{
  for (int iw=0; iw<grad.size(); iw++)
    grad[iw] = GradType();
  for(int i=0; i<Z.size(); i++)
    Z[i]->addGradient(W, iat, grad);
}

void
TrialWaveFunction::ratio (MCWalkerConfiguration &W, int iat,
                          std::vector<ValueType> &psi_ratios,
                          std::vector<GradType> &newG)
{
  for (int iw=0; iw<W.WalkerList.size(); iw++)
  {
    psi_ratios[iw] = 1.0;
    newG[iw] = GradType();
  }
  for (int i=0,ii=0; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->ratio (W, iat, psi_ratios, newG);
    myTimers[ii]->stop();
  }
}


void
TrialWaveFunction::ratio (MCWalkerConfiguration &W, int iat,
                          std::vector<ValueType> &psi_ratios,
                          std::vector<GradType> &newG, std::vector<ValueType> &newL)
{
  for (int iw=0; iw<W.WalkerList.size(); iw++)
  {
    psi_ratios[iw] = 1.0;
    newG[iw] = GradType();
    newL[iw] = ValueType();
  }
  for (int i=0,ii=1; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->ratio (W, iat, psi_ratios, newG, newL);
    myTimers[ii]->stop();
  }
}


void
TrialWaveFunction::calcRatio (MCWalkerConfiguration &W, int iat,
                              std::vector<ValueType> &psi_ratios,
                              std::vector<GradType> &newG, std::vector<ValueType> &newL)
{
  for (int iw=0; iw<W.WalkerList.size(); iw++)
  {
    psi_ratios[iw] = 1.0;
    newG[iw] = GradType();
    newL[iw] = ValueType();
  }
  for (int i=0,ii=1; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->calcRatio (W, iat, psi_ratios, newG, newL);
    myTimers[ii]->stop();
  }
}


void
TrialWaveFunction::addRatio (MCWalkerConfiguration &W, int iat,
                             std::vector<ValueType> &psi_ratios,
                             std::vector<GradType> &newG, std::vector<ValueType> &newL)
{
  for (int i=0,ii=1; i<Z.size()-1; i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->addRatio (W, iat, psi_ratios, newG, newL);
    myTimers[ii]->stop();
  }
  gpu::synchronize();
  myTimers[1+TIMER_SKIP*(Z.size()-1)]->start();
  Z[Z.size()-1]->addRatio (W, iat, psi_ratios, newG, newL);
  myTimers[1+TIMER_SKIP*(Z.size()-1)]->stop();
}

#if defined (QMC_COMPLEX)
void
TrialWaveFunction::convertRatiosFromComplexToReal (std::vector<ValueType> &psi_ratios,
                                                   std::vector<RealType> &psi_ratios_real)
{
  RealType PhaseValue;

  if (psi_ratios.size() != psi_ratios_real.size() ) {
    std::cerr << "Error: In " << __FILE__ << " , line " << __LINE__ << " , "
         << "input vector and output vector sizes unmatched." << std::endl; 
    abort();
  }

  for (int iw=0; iw<psi_ratios.size(); iw++)
  {
    psi_ratios_real[iw] = evaluateLogAndPhase(psi_ratios[iw], PhaseValue);
    psi_ratios_real[iw] = std::exp(psi_ratios_real[iw]);
  }
}
#endif

void
TrialWaveFunction::ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
                          std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
                          std::vector<GradType> &newG, std::vector<ValueType> &newL)
{
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] = 1.0;
    newG[iw] = GradType();
    newL[iw] = ValueType();
  }
  for (int i=0,ii=1; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->ratio (walkers, iatList, rNew, psi_ratios, newG, newL);
    myTimers[ii]->stop();
  }
}

void
TrialWaveFunction::ratio (MCWalkerConfiguration &W, int iat,
                          std::vector<ValueType> &psi_ratios)
{
  for (int iw=0; iw<W.WalkerList.size(); iw++)
    psi_ratios[iw] = 1.0;
  for (int i=0,ii=V_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->ratio (W, iat, psi_ratios);
    myTimers[ii]->stop();
  }
}

void
TrialWaveFunction::update (std::vector<Walker_t*> &walkers, int iat)
{
  for (int i=0,ii=ACCEPT_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->update(walkers, iat);
    myTimers[ii]->stop();
  }
}

void
TrialWaveFunction::update (const std::vector<Walker_t*> &walkers,
                           const std::vector<int> &iatList)
{
  for (int i=0,ii=ACCEPT_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->update(walkers, iatList);
    myTimers[ii]->stop();
  }
}


void
TrialWaveFunction::gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                             ValueMatrix_t &lapl)
{
  for (int i=0; i<grads.rows(); i++)
    for (int j=0; j<grads.cols(); j++)
    {
      grads(i,j) = GradType();
      lapl(i,j)  = ValueType();
    }
  for (int i=0,ii=VGL_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->gradLapl (W, grads, lapl);
    myTimers[ii]->stop();
  }
  for (int iw=0; iw<W.WalkerList.size(); iw++)
    for (int ptcl=0; ptcl<grads.cols(); ptcl++)
    {
      W[iw]->G[ptcl] = grads(iw, ptcl);
      W[iw]->L[ptcl]  = lapl(iw, ptcl);
    }
}

void
TrialWaveFunction::NLratios (MCWalkerConfiguration &W,
                             std::vector<NLjob> &jobList,
                             std::vector<PosType> &quadPoints,
                             std::vector<ValueType> &psi_ratios)
{
  for (int i=0; i<psi_ratios.size(); i++)
    psi_ratios[i] = 1.0;
  for (int i=0,ii=NL_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->NLratios(W, jobList, quadPoints, psi_ratios);
    myTimers[ii]->stop();
  }
}

void
TrialWaveFunction::NLratios (MCWalkerConfiguration &W,
                             gpu::device_vector<CUDA_PRECISION*> &Rlist,
                             gpu::device_vector<int*>            &ElecList,
                             gpu::device_vector<int>             &NumCoreElecs,
                             gpu::device_vector<CUDA_PRECISION*> &QuadPosList,
                             gpu::device_vector<CUDA_PRECISION*> &RatioList,
                             int numQuadPoints)
{
  for (int i=0,ii=NL_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->NLratios(W, Rlist, ElecList, NumCoreElecs,
                   QuadPosList, RatioList, numQuadPoints);
    myTimers[ii]->stop();
  }
}

void
TrialWaveFunction::evaluateDeltaLog(MCWalkerConfiguration &W,
                                    std::vector<RealType>& logpsi_opt)
{
  for (int iw=0; iw<logpsi_opt.size(); iw++)
    logpsi_opt[iw] = RealType();
  for (int i=0,ii=RECOMPUTE_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    if (Z[i]->Optimizable)
      Z[i]->addLog(W, logpsi_opt);
    myTimers[ii]->stop();
  }
}


void
TrialWaveFunction::evaluateDeltaLog (MCWalkerConfiguration &W,
                                     std::vector<RealType>& logpsi_fixed,
                                     std::vector<RealType>& logpsi_opt,
                                     GradMatrix_t&  fixedG,
                                     ValueMatrix_t& fixedL)
{
  for (int iw=0; iw<logpsi_fixed.size(); iw++)
  {
    logpsi_opt[iw] = RealType();
    logpsi_fixed[iw] = RealType();
  }
  fixedG = GradType();
  fixedL = RealType();
  // First, sum optimizable part, using fixedG and fixedL as temporaries
  for (int i=0,ii=RECOMPUTE_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    if (Z[i]->Optimizable)
    {
      Z[i]->addLog(W, logpsi_opt);
      Z[i]->gradLapl(W, fixedG, fixedL);
    }
  }
  for (int iw=0; iw<W.WalkerList.size(); iw++)
    for (int ptcl=0; ptcl<fixedG.cols(); ptcl++)
    {
      W[iw]->G[ptcl] = fixedG(iw, ptcl);
      W[iw]->L[ptcl]  = fixedL(iw, ptcl);
    }
  // Reset them, then accumulate the fixe part
  fixedG = GradType();
  fixedL = RealType();
  for (int i=0,ii=NL_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    if (!Z[i]->Optimizable)
    {
      Z[i]->addLog(W, logpsi_fixed);
      Z[i]->gradLapl(W, fixedG, fixedL);
    }
    myTimers[ii]->stop();
  }
  // Add on the fixed part to the total laplacian and gradient
  for (int iw=0; iw<W.WalkerList.size(); iw++)
    for (int ptcl=0; ptcl<fixedG.cols(); ptcl++)
    {
      W[iw]->G[ptcl] += fixedG(iw, ptcl);
      W[iw]->L[ptcl]  += fixedL(iw, ptcl);
    }
}

void
TrialWaveFunction::evaluateOptimizableLog (MCWalkerConfiguration &W,
    std::vector<RealType>& logpsi_opt,
    GradMatrix_t&  optG,
    ValueMatrix_t& optL)
{
  for (int iw=0; iw<W.getActiveWalkers(); iw++)
    logpsi_opt[iw] = RealType();
  optG = GradType();
  optL = RealType();
  // Sum optimizable part of log Psi
  for (int i=0,ii=RECOMPUTE_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    if (Z[i]->Optimizable)
    {
      Z[i]->addLog(W, logpsi_opt);
      Z[i]->gradLapl(W, optG, optL);
    }
    myTimers[ii]->stop();
  }
}



void
TrialWaveFunction::evaluateDerivatives (MCWalkerConfiguration &W,
                                        const opt_variables_type& optvars,
                                        RealMatrix_t &dlogpsi,
                                        RealMatrix_t &dhpsioverpsi)
{
  for (int i=0,ii=DERIVS_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    if (Z[i]->Optimizable)
      Z[i]->evaluateDerivatives(W, optvars, dlogpsi, dhpsioverpsi);
    myTimers[ii]->stop();
  }
}
}
