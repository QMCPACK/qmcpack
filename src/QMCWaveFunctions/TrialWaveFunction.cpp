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
    
    
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

typedef enum { V_TIMER, VGL_TIMER, ACCEPT_TIMER, NL_TIMER,
               RECOMPUTE_TIMER, DERIVS_TIMER, TIMER_SKIP
             } TimerEnum;


TrialWaveFunction::TrialWaveFunction(Communicate* c)
  : MPIObjectBase(c)
  , Ordered(true), RecomputeNeedsDistanceTable(false), NumPtcls(0), TotalDim(0), BufferCursor(0), BufferCursor_DP(0)
  , PhaseValue(0.0),LogValue(0.0),OneOverM(1.0), PhaseDiff(0.0), FermionWF(0)
{
  ClassName="TrialWaveFunction";
  myName="psi0";
}

///private and cannot be used
TrialWaveFunction::TrialWaveFunction()
  : MPIObjectBase(0)
  , Ordered(true), RecomputeNeedsDistanceTable(false), NumPtcls(0), TotalDim(0), BufferCursor(0), BufferCursor_DP(0)
  ,  PhaseValue(0.0),LogValue(0.0) ,OneOverM(1.0), PhaseDiff(0.0)
{
  ClassName="TrialWaveFunction";
  myName="psi0";
}

/** Destructor
*
*@warning Have not decided whether Z is cleaned up by TrialWaveFunction
*  or not. It will depend on I/O implementation.
*/
TrialWaveFunction::~TrialWaveFunction()
{
  delete_iter(Z.begin(),Z.end());
  //delete_iter(SPOSet.begin(),SPOSet.end());
  //delete_iter(myTimers.begin(),myTimers.end());
}

void
TrialWaveFunction::resetTargetParticleSet(ParticleSet& P)
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->resetTargetParticleSet(P);
}

void TrialWaveFunction::startOptimization()
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->IsOptimizing=true;
}

void TrialWaveFunction::stopOptimization()
{
  for (int i=0; i<Z.size(); i++)
  {
    Z[i]->finalizeOptimization();
    Z[i]->IsOptimizing=false;
  }
}

/** add an ObritalBase
 * @param aterm an OrbitalBase
 * @param aname  name of aterm
 * @param fermion if true, set aterm to FermionWF
 */
void
TrialWaveFunction::addOrbital(OrbitalBase* aterm, const std::string& aname, bool fermion)
{
  Z.push_back(aterm);
  aterm->IsFermionWF=fermion;
  if(fermion) 
  {
    app_log() << "  FermionWF=" << aname << std::endl;
    FermionWF=dynamic_cast<FermionBase*>(aterm);
  }

  if (aterm->RecomputeNeedsDistanceTable) RecomputeNeedsDistanceTable = true;

//#if defined(QMC_CUDA)
  char name1[64],name2[64],name3[64],name4[64], name5[64], name6[64];
  sprintf(name1,"WaveFunction::%s_V",aname.c_str());
  sprintf(name2,"WaveFunction::%s_VGL",aname.c_str());
  sprintf(name3,"WaveFunction::%s_accept",aname.c_str());
  sprintf(name4,"WaveFunction::%s_NLratio",aname.c_str());
  sprintf(name5,"WaveFunction::%s_recompute",aname.c_str());
  sprintf(name6,"WaveFunction::%s_derivs",aname.c_str());
  NewTimer *vtimer=new NewTimer(name1);
  NewTimer *vgltimer=new NewTimer(name2);
  NewTimer *accepttimer=new NewTimer(name3);
  NewTimer *NLtimer=new NewTimer(name4);
  NewTimer *recomputetimer=new NewTimer(name5);
  NewTimer *derivstimer=new NewTimer(name6);
  myTimers.push_back(vtimer);
  myTimers.push_back(vgltimer);
  myTimers.push_back(accepttimer);
  myTimers.push_back(NLtimer);
  myTimers.push_back(recomputetimer);
  myTimers.push_back(derivstimer);
  TimerManager.addTimer(vtimer);
  TimerManager.addTimer(vgltimer);
  TimerManager.addTimer(accepttimer);
  TimerManager.addTimer(NLtimer);
  TimerManager.addTimer(recomputetimer);
  TimerManager.addTimer(derivstimer);
//#endif
}


/** return log(|psi|)
*
* PhaseValue is the phase for the complex wave function
*/
TrialWaveFunction::RealType
TrialWaveFunction::evaluateLog(ParticleSet& P)
{
  //TAU_PROFILE("TrialWaveFunction::evaluateLog","ParticleSet& P", TAU_USER);
  P.G = 0.0;
  P.L = 0.0;
  ValueType logpsi(0.0);
  PhaseValue=0.0;
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  //WARNING: multiplication for PhaseValue is not correct, fix this!!
  for (; it!=it_end; ++it)
  {
    logpsi += (*it)->evaluateLog(P, P.G, P.L);
    PhaseValue += (*it)->PhaseValue;
  }
  convert(logpsi,LogValue);
  return LogValue;
  //return LogValue=real(logpsi);
}

void TrialWaveFunction::recompute(ParticleSet& P)
{
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (int ii=RECOMPUTE_TIMER; it!=it_end; ++it,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    (*it)->recompute(P);
    myTimers[ii]->stop();
  }
}

/** return log(|psi|)
*
* PhaseValue is the phase for the complex wave function
*/
TrialWaveFunction::RealType
TrialWaveFunction::evaluateLogOnly(ParticleSet& P)
{
  //TAU_PROFILE("TrialWaveFunction::evaluateLogOnly","ParticleSet& P", TAU_USER);
  tempP->R=P.R;
  tempP->L=0.0;
  tempP->G=0.0;
  ValueType logpsi(0.0);
  PhaseValue=0.0;
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  //WARNING: multiplication for PhaseValue is not correct, fix this!!
  for (; it!=it_end; ++it)
  {
    logpsi += (*it)->evaluateLog(*tempP, tempP->G, tempP->L);
    PhaseValue += (*it)->PhaseValue;
  }
  convert(logpsi,LogValue);
  return LogValue;
  //return LogValue=real(logpsi);
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
  ValueType logpsi(0.0);
  PhaseValue=0.0;
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  int ii=RECOMPUTE_TIMER;
  for (; it!=it_end; ++it,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    if ((*it)->Optimizable )
    {
      logpsi += (*it)->evaluateLog(P, P.G, P.L);
      PhaseValue += (*it)->PhaseValue;
    }
    myTimers[ii]->stop();
  }

  //In case we need to recompute orbitals, initialize dummy vectors for G and L.
  //evaluateLog dumps into these variables, and logPsi contribution is discarded.  
  //Only called for non-optimizable orbitals.   
  if (recomputeall)
  {
    ParticleSet::ParticleGradient_t dummyG(P.G);
    ParticleSet::ParticleLaplacian_t dummyL(P.L);
    
    it=Z.begin();
    it_end=Z.end();

    for( ; it!=it_end; ++it)
    {
      if (!(*it)->Optimizable) (*it)->evaluateLog(P,dummyG, dummyL); //update orbitals if its not flagged optimizable, AND recomputeall is true
    }
    
  }
  convert(logpsi,LogValue);
  return LogValue;
}

TrialWaveFunction::RealType TrialWaveFunction::evaluateDeltaLog(ParticleSet& P, PooledData<RealType>& buf)
{
  P.G = 0.0;
  P.L = 0.0;
  ValueType logpsi(0.0);
  PhaseValue=0.0;
  buf.rewind();
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  int ii=RECOMPUTE_TIMER;
  for (; it!=it_end; ++it,ii+=TIMER_SKIP)
  {
// mmorales: I don't remember if I did this, but eliminating the "if ((*it)->Optimizable)"
//           forces everything to be evaluated. This was probably done because for optm with the
//           nonlocal component in the cost function, the slater determinant might not be optimizable
//           but this must be called anyway to load the inverse. CHECK CHECK CHECK, FIX FIX FIX
    myTimers[ii]->start();
    if ((*it)->Optimizable)
    {
      logpsi += (*it)->evaluateLog(P, P.G, P.L,buf,false);
      PhaseValue += (*it)->PhaseValue;
    }
    else
//          ValueType x = (*it)->evaluateLog(P, P.G, P.L,buf,false);
      (*it)->copyFromDerivativeBuffer(P,buf);//keep buffer synched
    myTimers[ii]->stop();
  }
  convert(logpsi,LogValue);
  return LogValue;
  //return LogValue=real(logpsi);
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
* Additionally, dumpToBuffer and dumpFromBuffer is used to manage
* necessary data for ratio evaluations.
*/
void
TrialWaveFunction::evaluateDeltaLog(ParticleSet& P
                                    , RealType& logpsi_fixed_r, RealType& logpsi_opt_r
                                    , ParticleSet::ParticleGradient_t& fixedG, ParticleSet::ParticleLaplacian_t& fixedL)
{
  //TAU_PROFILE("TrialWaveFunction::evaluateDeltaLog","ParticleSet& P", TAU_USER);
  P.G = 0.0;
  P.L = 0.0;
  fixedL = 0.0;
  fixedG = 0.0;
  ValueType logpsi_fixed(0.0);
  ValueType logpsi_opt(0.0);
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  int ii=RECOMPUTE_TIMER;
  for (; it!=it_end; ++it,ii+=TIMER_SKIP)
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
  convert(logpsi_fixed,logpsi_fixed_r);
  convert(logpsi_opt,logpsi_opt_r);
  //logpsi_fixed_r = real(logpsi_fixed);
  //logpsi_opt_r = real(logpsi_opt);
}

void
TrialWaveFunction::evaluateDeltaLog(ParticleSet& P
                                    , RealType& logpsi_fixed_r, RealType& logpsi_opt_r
                                    , ParticleSet::ParticleGradient_t& fixedG
                                    , ParticleSet::ParticleLaplacian_t& fixedL
                                    , PooledData<RealType>& buf)
{
  //TAU_PROFILE("TrialWaveFunction::evaluateDeltaLog","ParticleSet& P", TAU_USER);
  P.G = 0.0;
  P.L = 0.0;
  fixedG = 0.0;
  fixedL = 0.0;
  ValueType logpsi_fixed(0.0);
  ValueType logpsi_opt(0.0);
  buf.rewind();
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  int ii=RECOMPUTE_TIMER;
  for (; it!=it_end; ++it,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    if ((*it)->Optimizable)
      logpsi_opt += (*it)->evaluateLog(P, P.G, P.L,buf,true);
    else
      logpsi_fixed += (*it)->evaluateLog(P, fixedG, fixedL,buf,true);
    myTimers[ii]->stop();
  }
  P.G += fixedG;
  P.L += fixedL;
  convert(logpsi_fixed,logpsi_fixed_r);
  convert(logpsi_opt,logpsi_opt_r);
  //logpsi_fixed_r = real(logpsi_fixed);
  //logpsi_opt_r = real(logpsi_opt);
}

/*void TrialWaveFunction::evaluateHessian(ParticleSet & P, int iat, HessType& grad_grad_psi)
{
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  
  grad_grad_psi=0.0;
  
  for (; it!=it_end; ++it)
  {	
	  HessType tmp_hess;
	  (*it)->evaluateHessian(P, iat, tmp_hess);
	  grad_grad_psi+=tmp_hess;
  }
}*/

void TrialWaveFunction::evaluateHessian(ParticleSet & P, HessVector_t& grad_grad_psi)
{
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  
  grad_grad_psi.resize(P.getTotalNum());
  
  for (int i=0; i<Z.size(); i++)
  {	
	  HessVector_t tmp_hess(grad_grad_psi);
	  tmp_hess=0.0;
	  Z[i]->evaluateHessian(P, tmp_hess);
	  grad_grad_psi+=tmp_hess;
	//  app_log()<<"TrialWavefunction::tmp_hess = "<<tmp_hess<< std::endl;
	//  app_log()<< std::endl<< std::endl;
  }
 // app_log()<<" TrialWavefunction::Hessian = "<<grad_grad_psi<< std::endl;
}


/** evaluate the value of a many-body wave function
*@param P input configuration containing N particles
*@return the value of many-body wave function
*
*Upon return, the gradient and laplacian operators are added by the components.
*/
TrialWaveFunction::ValueType TrialWaveFunction::evaluate(ParticleSet& P)
{
  //TAU_PROFILE("TrialWaveFunction::evaluate","ParticleSet& P", TAU_USER);
  P.G = 0.0;
  P.L = 0.0;
  ValueType psi(1.0);
  for (int i=0; i<Z.size(); i++)
  {
    psi *= Z[i]->evaluate(P, P.G, P.L);
  }
  //for(int iat=0; iat<P.getTotalNum(); iat++)
  // std::cout << P.G[iat] << " " << P.L[iat] << std::endl;
  LogValue = evaluateLogAndPhase(psi,PhaseValue);
  return psi;
}

TrialWaveFunction::RealType TrialWaveFunction::ratio(ParticleSet& P,int iat)
{
  //TAU_PROFILE("TrialWaveFunction::ratio","(ParticleSet& P,int iat)", TAU_USER);
  ValueType r(1.0);
  for (int i=0,ii=V_TIMER; i<Z.size(); ++i,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    r *= Z[i]->ratio(P,iat);
    myTimers[ii]->stop();
  }
#if defined(QMC_COMPLEX)
  //return std::exp(evaluateLogAndPhase(r,PhaseValue));
  RealType logr=evaluateLogAndPhase(r,PhaseDiff);
  return std::exp(logr);
#else
  if (r<0)
    PhaseDiff=M_PI;
  //     else PhaseDiff=0.0;
  return r;
#endif
}

TrialWaveFunction::ValueType TrialWaveFunction::full_ratio(ParticleSet& P,int iat)
{
  ValueType r(1.0);
  for(int i=0;i<Z.size();++i)
    r *= Z[i]->ratio(P,iat);
  return r;
}

TrialWaveFunction::RealType TrialWaveFunction::ratioVector(ParticleSet& P, int iat, std::vector<RealType>& ratios)
{
  //TAU_PROFILE("TrialWaveFunction::ratio","(ParticleSet& P,int iat)", TAU_USER);
  ratios.resize(Z.size(),0);
  ValueType r(1.0);
  for (int i=0,ii=V_TIMER; i<Z.size(); ++i,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    ValueType zr=Z[i]->ratio(P,iat);
    r *= zr;
#if defined(QMC_COMPLEX)
    ratios[i] = std::abs(zr);
#else
    ratios[i] = zr;
#endif
    myTimers[ii]->stop();
  }
#if defined(QMC_COMPLEX)
  //return std::exp(evaluateLogAndPhase(r,PhaseValue));
  RealType logr=evaluateLogAndPhase(r,PhaseDiff);
  return std::exp(logr);
#else
  if (r<0)
    PhaseDiff=M_PI;
  //     else PhaseDiff=0.0;
  return r;
#endif
}

TrialWaveFunction::RealType TrialWaveFunction::alternateRatio(ParticleSet& P)
{
  //TAU_PROFILE("TrialWaveFunction::ratio","(ParticleSet& P,int iat)", TAU_USER);
  ValueType r(1.0);
  for (int i=0,ii=0; i<Z.size(); ++i,ii+=2)
  {
    r *= Z[i]->alternateRatio(P);
  }
#if defined(QMC_COMPLEX)
  //return std::exp(evaluateLogAndPhase(r,PhaseValue));
  RealType logr=evaluateLogAndPhase(r,PhaseDiff);
  return std::exp(logr);
#else
  if (r<0)
    PhaseDiff=M_PI;
  return r;
#endif
}

TrialWaveFunction::GradType TrialWaveFunction::evalGrad(ParticleSet& P,int iat)
{
  //TAU_PROFILE("TrialWaveFunction::evalGrad","(ParticleSet& P,int iat)", TAU_USER);
  GradType grad_iat;
  for (int i=0; i<Z.size(); ++i)
    grad_iat += Z[i]->evalGrad(P,iat);
  return grad_iat;
}

TrialWaveFunction::GradType TrialWaveFunction::alternateEvalGrad(ParticleSet& P,int iat)
{
  //TAU_PROFILE("TrialWaveFunction::evalGrad","(ParticleSet& P,int iat)", TAU_USER);
  GradType grad_iat;
  for (int i=0; i<Z.size(); ++i)
    grad_iat += Z[i]->alternateEvalGrad(P,iat);
  return grad_iat;
}

// Evaluates the gradient w.r.t. to the source of the Laplacian
// w.r.t. to the electrons of the wave function.
TrialWaveFunction::GradType TrialWaveFunction::evalGradSource(ParticleSet& P
    , ParticleSet &source, int iat)
{
  GradType grad_iat = GradType();
  for (int i=0; i<Z.size(); ++i)
    grad_iat += Z[i]->evalGradSource(P, source, iat);
  return grad_iat;
}

TrialWaveFunction::GradType TrialWaveFunction::evalGradSource(ParticleSet& P
    , ParticleSet &source, int iat
    , TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad
    , TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  GradType grad_iat = GradType();
  for (int dim=0; dim<OHMMS_DIM; dim++)
    for (int i=0; i<grad_grad[0].size(); i++)
    {
      grad_grad[dim][i] = GradType();
      lapl_grad[dim][i] = 0.0;
    }
  for (int i=0; i<Z.size(); ++i)
    grad_iat += Z[i]->evalGradSource(P, source, iat,
                                     grad_grad, lapl_grad);
  return grad_iat;
}




TrialWaveFunction::RealType TrialWaveFunction::ratioGrad(ParticleSet& P
    ,int iat, GradType& grad_iat )
{
  //TAU_PROFILE("TrialWaveFunction::ratioGrad","(ParticleSet& P,int iat)", TAU_USER);
  grad_iat=0.0;
  ValueType r(1.0);
  for (int i=0,ii=VGL_TIMER; i<Z.size(); ++i,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    r *= Z[i]->ratioGrad(P,iat,grad_iat );
    myTimers[ii]->stop();
  }
#if defined(QMC_COMPLEX)
  //return std::exp(evaluateLogAndPhase(r,PhaseValue));
  RealType logr=evaluateLogAndPhase(r,PhaseValue);
  return std::exp(logr);
#else
  if (r<0)
    PhaseDiff=M_PI;
  //     else PhaseDiff=0.0;
  return r;
#endif
}

TrialWaveFunction::RealType TrialWaveFunction::alternateRatioGrad(ParticleSet& P
    ,int iat, GradType& grad_iat )
{
  //TAU_PROFILE("TrialWaveFunction::ratioGrad","(ParticleSet& P,int iat)", TAU_USER);
  grad_iat=0.0;
  ValueType r(1.0);
  for (int i=0,ii=VGL_TIMER; i<Z.size(); ++i,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    r *= Z[i]->alternateRatioGrad(P,iat,grad_iat);
    myTimers[ii]->stop();
  }
#if defined(QMC_COMPLEX)
  //return std::exp(evaluateLogAndPhase(r,PhaseValue));
  RealType logr=evaluateLogAndPhase(r,PhaseValue);
  return std::exp(logr);
#else
  if (r<0)
    PhaseDiff=M_PI;
  //     else PhaseDiff=0.0;
  return r;
#endif
}

void   TrialWaveFunction::update(ParticleSet& P,int iat)
{
  //ready to collect "changes" in the gradients and laplacians by the move
  delta_G=0.0;
  delta_L=0.0;
  for (int i=0, ii=ACCEPT_TIMER; i<Z.size(); i++,ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    Z[i]->update(P,delta_G,delta_L,iat);
    myTimers[ii]->stop();
  }
  P.G += delta_G;
  P.L += delta_L;
}


/** evaluate \f$ frac{\Psi({\bf R}_i^{'})}{\Psi({\bf R}_i)}\f$
 * @param P ParticleSet
 * @param iat index of the particle with a trial move
 * @param dG total differentcal gradients
 * @param dL total differential laplacians
 * @return ratio
 *
 * Each OrbitalBase object adds the differential gradients and lapacians.
 */
TrialWaveFunction::RealType TrialWaveFunction::ratio(ParticleSet& P, int iat
    , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL)
{
  //TAU_PROFILE("TrialWaveFunction::ratio","(P,iat,dG,dL)", TAU_USER);
  dG = 0.0;
  dL = 0.0;
  ValueType r(1.0);
  for (int i=0, ii=VGL_TIMER; i<Z.size(); ++i, ii+=TIMER_SKIP)
  {
    myTimers[ii]->start();
    r *= Z[i]->ratio(P,iat,dG,dL);
    myTimers[ii]->stop();
  }
#if defined(QMC_COMPLEX)
  return std::exp(evaluateLogAndPhase(r,PhaseDiff));
#else
  if (r<0)
    PhaseDiff=M_PI;
  //     else PhaseDiff=0.0;
  return r;
#endif
}

void TrialWaveFunction::printGL(ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L, std::string tag)
{
  std::ostringstream o;
  o << "---  reporting " << tag << std::endl << "  ---" << std::endl;
  for(int iat=0; iat<L.size(); iat++)
    o << "index: " << std::fixed << iat << std::scientific
      << "   G: " << G[iat][0] << "  " << G[iat][1] << "  " << G[iat][2]
      << "   L: " << L[iat] << std::endl;
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
  for (int i=0; i<Z.size(); i++)
  {
    Z[i]->restore(iat);
  }
  PhaseDiff=0;
}

/** update the state with the new data
 * @param P ParticleSet
 * @param iat index of the particle with a trial move
 *
 * The proposed move of the iath particle is accepted.
 * All the temporary data should be incorporated so that the next move is valid.
 */
void   TrialWaveFunction::acceptMove(ParticleSet& P,int iat)
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->acceptMove(P,iat);
  PhaseValue += PhaseDiff;
  PhaseDiff=0.0;
  LogValue=0;
  for (int i=0; i<Z.size(); i++)
    LogValue+= Z[i]->LogValue;
}

//void TrialWaveFunction::resizeByWalkers(int nwalkers){
//  for(int i=0; i<Z.size(); i++) Z[i]->resizeByWalkers(nwalkers);
//}

void TrialWaveFunction::checkInVariables(opt_variables_type& active)
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->checkInVariables(active);
}

void TrialWaveFunction::checkOutVariables(const opt_variables_type& active)
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->checkOutVariables(active);
}

void TrialWaveFunction::resetParameters(const opt_variables_type& active)
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->resetParameters(active);
}

void TrialWaveFunction::reportStatus(std::ostream& os)
{
  for (int i=0; i<Z.size(); i++)
    Z[i]->reportStatus(os);
}

void TrialWaveFunction::getLogs(std::vector<RealType>& lvals)
{
  lvals.resize(Z.size(),0);
  for (int i=0; i<Z.size(); i++)
  {
    lvals[i] = Z[i]->LogValue;
  }
}

void TrialWaveFunction::getPhases(std::vector<RealType>& pvals)
{
  pvals.resize(Z.size(),0);
  for (int i=0; i<Z.size(); i++)
  {
    pvals[i] = Z[i]->PhaseValue;
  }
}


TrialWaveFunction::RealType TrialWaveFunction::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  delta_G.resize(P.getTotalNum());
  delta_L.resize(P.getTotalNum());
  P.G = 0.0;
  P.L = 0.0;
  //save the current position
  BufferCursor=buf.current();
  BufferCursor_DP=buf.current_DP();
  ValueType logpsi(0.0);
  PhaseValue=0.0;
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (; it!=it_end; ++it)
  {
    logpsi += (*it)->registerData(P,buf);
    PhaseValue += (*it)->PhaseValue;
  }
  convert(logpsi,LogValue);
  //LogValue=real(logpsi);
//append current gradients and laplacians to the buffer
  NumPtcls = P.getTotalNum();
  TotalDim = PosType::Size*NumPtcls;
  buf.add(PhaseValue);
  buf.add(LogValue);
  buf.add(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
  buf.add(&(P.L[0]), &(P.L[P.getTotalNum()]));
  return LogValue;
}

TrialWaveFunction::RealType TrialWaveFunction::registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf, int storageType)
{
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (; it!=it_end; ++it)
    (*it)->registerDataForDerivatives(P,buf,storageType);
  return 1.0;
}

void TrialWaveFunction::memoryUsage_DataForDerivatives(ParticleSet& P,long& orbs_only,long& orbs, long& invs, long& dets)
{
  orbs_only=orbs=invs=dets=0;
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (; it!=it_end; ++it)
    (*it)->memoryUsage_DataForDerivatives(P,orbs_only,orbs,invs,dets);
}

TrialWaveFunction::RealType TrialWaveFunction::updateBuffer(ParticleSet& P
    , PooledData<RealType>& buf, bool fromscratch)
{
  //TAU_PROFILE("TrialWaveFunction::updateBuffer","(P,..)", TAU_USER);
  P.G = 0.0;
  P.L = 0.0;
  buf.rewind(BufferCursor,BufferCursor_DP);
  ValueType logpsi(0.0);
  PhaseValue=0.0;
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (; it!=it_end; ++it)
  {
    logpsi += (*it)->updateBuffer(P,buf,fromscratch);
    PhaseValue += (*it)->PhaseValue;
  }
  //printGL(P.G,P.L);
  convert(logpsi,LogValue);
  //LogValue=real(logpsi);
  buf.put(PhaseValue);
  buf.put(LogValue);
  buf.put(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
  buf.put(&(P.L[0]), &(P.L[0])+NumPtcls);
  return LogValue;
}

void TrialWaveFunction::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.rewind(BufferCursor,BufferCursor_DP);
  //TAU_PROFILE("TrialWaveFunction::copyFromBuffer","(P,..)", TAU_USER);
  for (int i=0; i<Z.size(); i++)
    Z[i]->copyFromBuffer(P,buf);
  //get the gradients and laplacians from the buffer
  buf.get(PhaseValue);
  buf.get(LogValue);
  buf.get(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
  buf.get(&(P.L[0]), &(P.L[0])+NumPtcls);
}

/** Dump data that are required to evaluate ratios to the buffer
* @param P active ParticleSet
* @param buf anonymous buffer to which the data will be dumped.
*
* This function lets the OrbitalBase objects store the minimal data
* that are required to evaluate the ratio, even though the components
* are invariant during the optimizations.
*/
void TrialWaveFunction::dumpToBuffer(ParticleSet& P, BufferType& buf)
{
  buf.rewind(BufferCursor,BufferCursor_DP);
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (; it!=it_end; ++it)
  {
    (*it)->dumpToBuffer(P,buf);
  }
}

/** copy data that are required to evaluate ratios from the buffer
* @param P active ParticleSet
* @param buf anonymous buffer from which the data will be copied.
*
* This function lets the OrbitalBase objects get the minimal data
* that are required to evaluate the ratio from the buffer.
* Only the data registered by dumToBuffer will be available.
*/
void TrialWaveFunction::dumpFromBuffer(ParticleSet& P, BufferType& buf)
{
  buf.rewind(BufferCursor,BufferCursor_DP);
  std::vector<OrbitalBase*>::iterator it(Z.begin());
  std::vector<OrbitalBase*>::iterator it_end(Z.end());
  for (; it!=it_end; ++it)
    (*it)->dumpFromBuffer(P,buf);
}

TrialWaveFunction::RealType
TrialWaveFunction::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.rewind(BufferCursor,BufferCursor_DP);
  LogValue=0.0;
  PhaseValue=0.0;
  for (int i=0; i<Z.size(); i++)
  {
    LogValue += Z[i]->evaluateLog(P,buf);
    PhaseValue += Z[i]->PhaseValue;
  }
  buf.put(PhaseValue);
  buf.put(LogValue);
  //buf.put(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
  //buf.put(&(P.L[0]), &(P.L[0])+NumPtcls);
  return LogValue;
}

void TrialWaveFunction::evaluateRatios(VirtualParticleSet& VP, std::vector<RealType>& ratios)
{
#if defined(QMC_COMPLEX)
  std::vector<ValueType> t(ratios.size()),r(ratios.size(),1.0);;
  for (int i=0; i<Z.size(); ++i)
  {
    Z[i]->evaluateRatios(VP,t);
    for (int j=0; j<ratios.size(); ++j)
      r[j]*=t[j];
  }
  RealType pdiff;
  for(int j=0; j<ratios.size(); ++j)
  {
    RealType logr=evaluateLogAndPhase(r[j],pdiff);
    ratios[j]=std::exp(logr)*std::cos(pdiff);
    //ratios[j]=std::abs(r)*std::cos(std::arg(r[j]));
  }
#else
  std::fill(ratios.begin(),ratios.end(),1.0);
  std::vector<ValueType> t(ratios.size());
  for (int i=0; i<Z.size(); ++i)
  {
    Z[i]->evaluateRatios(VP,t);
    for (int j=0; j<ratios.size(); ++j)
      ratios[j]*=t[j];
  }
#endif
}

void TrialWaveFunction::evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
    std::vector<RealType>& ratios, Matrix<RealType>& dratio)
{
#if defined(QMC_COMPLEX)
  APP_ABORT("TrialWaveFunction::evaluateDerivRatios not available for complex wavefunctions");
#else
  std::fill(ratios.begin(),ratios.end(),1.0);
  std::vector<ValueType> t(ratios.size());
  for (int i=0; i<Z.size(); ++i)
  {
    Z[i]->evaluateDerivRatios(VP,optvars,t,dratio);
    for (int j=0; j<ratios.size(); ++j)
      ratios[j]*=t[j];
  }
#endif
}

//TrialWaveFunction::RealType
//TrialWaveFunction::evaluate(ParticleSet& P, PooledData<RealType>& buf) {

//  ValueType psi(1.0);
//  for(int i=0; i<Z.size(); i++) psi *= Z[i]->evaluate(P,buf);
//  buf.put(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
//  buf.put(&(P.L[0]), &(P.L[0])+NumPtcls);
//  return real(psi);
//}
//bool TrialWaveFunction::hasSPOSet(const std::string& aname) {
//  return false;
//  //bool notfoundit=true;
//  //vector<OhmmsElementBase*>::iterator it(SPOSet.begin());
//  //vector<OhmmsElementBase*>::iterator it_end(SPOSet.end());
//  //while(notfoundit && it != it_end) {
//  //  if((*it)->getName() == aname) notfoundit=false;
//  //  ++it;
//  //}
//  //return !notfoundit;
//}

//OhmmsElementBase*
//TrialWaveFunction::getSPOSet(const std::string& aname) {
//  //bool notfoundit=true;
//  //vector<OhmmsElementBase*>::iterator it(SPOSet.begin());
//  //vector<OhmmsElementBase*>::iterator it_end(SPOSet.end());
//  //while(notfoundit && it != it_end) {
//  //  if((*it)->getName() == aname) return *it;
//  //  ++it;
//  //}
//  return 0;
//}

//void TrialWaveFunction::addSPOSet(OhmmsElementBase* spo) {
//  //SPOSet.push_back(spo);
//}

bool TrialWaveFunction::put(xmlNodePtr cur)
{
  return true;
}

void TrialWaveFunction::reset()
{
}

void TrialWaveFunction::reverse()
{
  Ordered=false;
  //vector<OrbitalBase*> zcopy(Z);
  //int n=Z.size()-1;
  //for(int i=0; i<Z.size(); ++i) Z[n-i]=zcopy[i];
}

TrialWaveFunction* TrialWaveFunction::makeClone(ParticleSet& tqp)  const
{
  TrialWaveFunction* myclone = new TrialWaveFunction(myComm);
  myclone->BufferCursor=BufferCursor;
  myclone->BufferCursor_DP=BufferCursor_DP;
  for (int i=0; i<Z.size(); ++i)
    myclone->addOrbital(Z[i]->makeClone(tqp),"dummy",Z[i]->IsFermionWF);
  myclone->OneOverM=OneOverM;
  return myclone;
}

/** evaluate derivatives of KE wrt optimizable varibles
 *
 * @todo OrbitalBase objects should take the mass into account.
 */
void TrialWaveFunction::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& optvars,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi,
    bool project)
{
  //     // First, zero out derivatives
  //  This should only be done for some variables.
  //     for (int j=0; j<dlogpsi.size(); j++)
  //       dlogpsi[j] = dhpsioverpsi[j] = 0.0;
  for (int i=0; i<Z.size(); i++)
  {
    if (Z[i]->dPsi)
      (Z[i]->dPsi)->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
    else
      Z[i]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }
  //orbitals do not know about mass of particle.
  for (int i=0; i<dhpsioverpsi.size(); i++)
    dhpsioverpsi[i]*=OneOverM;

  if (project)
  {
    for (int i=0; i<Z.size(); i++)
    {
      if (Z[i]->dPsi)
        (Z[i]->dPsi)->multiplyDerivsByOrbR(dlogpsi);
      else
        Z[i]->multiplyDerivsByOrbR(dlogpsi);
    }
    RealType psiValue=std::exp(-LogValue)*std::cos(PhaseValue);
    for (int i=0; i<dlogpsi.size(); i++)
      dlogpsi[i] *= psiValue;
  }
}

TrialWaveFunction::RealType
TrialWaveFunction::KECorrection() const
{
  RealType sum = 0.0;
  for (int i=0; i<Z.size(); ++i)
    sum += Z[i]->KECorrection();
  return sum;
}

void TrialWaveFunction::get_ratios(ParticleSet& P, std::vector<ValueType>& ratios)
{
  std::fill(ratios.begin(),ratios.end(),1.0);
  std::vector<ValueType> t(ratios.size());
  for (int i=0; i<Z.size(); ++i)
  {
    Z[i]->get_ratios(P,t);
    for (int j=0; j<t.size(); ++j)
      ratios[j]*=t[j];
  }
}

}
/***************************************************************************
* $RCSfile$   $Author$
* $Revision$   $Date$
* $Id$
***************************************************************************/
