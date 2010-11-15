// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
  {

  TrialWaveFunction::TrialWaveFunction(Communicate* c)
      : MPIObjectBase(c)
      , Ordered(true), NumPtcls(0), TotalDim(0), BufferCursor(0)
      , PhaseValue(0.0),LogValue(0.0),OneOverM(1.0), PhaseDiff(0.0)
  {
    ClassName="TrialWaveFunction";
    myName="psi0";
  }

  ///private and cannot be used
  TrialWaveFunction::TrialWaveFunction()
      : MPIObjectBase(0)
      , Ordered(true), NumPtcls(0), TotalDim(0), BufferCursor(0)
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
    for (int i=0; i<Z.size(); i++) Z[i]->resetTargetParticleSet(P);
  }

  void TrialWaveFunction::startOptimization()
  {
    for (int i=0; i<Z.size(); i++) Z[i]->IsOptimizing=true;
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
  *@param aterm a many-body wavefunction
  */
  void
  //TrialWaveFunction::addOrbital(OrbitalBase* aterm)
  TrialWaveFunction::addOrbital(OrbitalBase* aterm, const string& aname)
  {

    Z.push_back(aterm);
    //int n=Z.size();
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
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());

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
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());

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
   * @param needratio users request ratio evaluation
   * @param buf anonymous storage for the reusable data
   * @return the value of \f$ \log( \Pi_i \Psi_i) \f$  many-body wave function
   *
   * @if needratio == true
   *  need to update the data from buf, since external objects need to evaluate ratios, e.g., non-local pseudopotentials
   * @else
   *  evaluate the value only
   *
   * Upon return, the gradient and laplacian operators are added by the components.
   * Each OrbitalBase evaluates PhaseValue and LogValue = log(abs(psi_i))
   * Jastrow functions always have PhaseValue=1.
   */
  TrialWaveFunction::RealType TrialWaveFunction::evaluateDeltaLog(ParticleSet& P)
  {
    P.G = 0.0;
    P.L = 0.0;
    ValueType logpsi(0.0);
    PhaseValue=0.0;
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (; it!=it_end; ++it)
      {
        if ((*it)->Optimizable)
          {
            logpsi += (*it)->evaluateLog(P, P.G, P.L);
            PhaseValue += (*it)->PhaseValue;
          }
      }
    convert(logpsi,LogValue);
    return LogValue;
    //return LogValue=real(logpsi);
  }

  TrialWaveFunction::RealType TrialWaveFunction::evaluateDeltaLog(ParticleSet& P, PooledData<RealType>& buf)
  {
    P.G = 0.0;
    P.L = 0.0;
    ValueType logpsi(0.0);
    PhaseValue=0.0;
    buf.rewind();
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (; it!=it_end; ++it)
      {
        if ((*it)->Optimizable)
          {
            logpsi += (*it)->evaluateLog(P, P.G, P.L,buf,false);
            PhaseValue += (*it)->PhaseValue;
          }
      }
    convert(logpsi,LogValue);
    return LogValue;
    //return LogValue=real(logpsi);
  }


  /** evalaute the sum of log value of optimizable many-body wavefunctions
  * @param P  input configuration containing N particles
  * @param logpsi_fixed log(abs(psi)) of the invariant orbitals
  * @param logpsi_opt log(abs(psi)) of the variable orbitals
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
    fixedG = 0.0;
    fixedL = 0.0;
    ValueType logpsi_fixed(0.0);
    ValueType logpsi_opt(0.0);
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (; it!=it_end; ++it)
      {
        if ((*it)->Optimizable)
          logpsi_opt += (*it)->evaluateLog(P, P.G, P.L);
        else
          {
            logpsi_fixed += (*it)->evaluateLog(P, fixedG, fixedL);
          }
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
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (; it!=it_end; ++it)
      {
        if ((*it)->Optimizable)
          logpsi_opt += (*it)->evaluateLog(P, P.G, P.L,buf,true);
        else
          {
            logpsi_fixed += (*it)->evaluateLog(P, fixedG, fixedL);
          }
      }
    P.G += fixedG;
    P.L += fixedL;
    convert(logpsi_fixed,logpsi_fixed_r);
    convert(logpsi_opt,logpsi_opt_r);
    //logpsi_fixed_r = real(logpsi_fixed);
    //logpsi_opt_r = real(logpsi_opt);
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
    // cout << P.G[iat] << " " << P.L[iat] << endl;
    LogValue = evaluateLogAndPhase(psi,PhaseValue);
    return psi;
  }

  TrialWaveFunction::RealType TrialWaveFunction::ratio(ParticleSet& P,int iat)
  {
    //TAU_PROFILE("TrialWaveFunction::ratio","(ParticleSet& P,int iat)", TAU_USER);
    ValueType r(1.0);
    for (int i=0,ii=0; i<Z.size(); ++i,ii+=2)
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
    if (r<0) PhaseDiff=M_PI;
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
        myTimers[ii]->start();
        r *= Z[i]->alternateRatio(P);
        myTimers[ii]->stop();
      }
#if defined(QMC_COMPLEX)
    //return std::exp(evaluateLogAndPhase(r,PhaseValue));
    RealType logr=evaluateLogAndPhase(r,PhaseDiff);
    return std::exp(logr);
#else
    if (r<0) PhaseDiff=M_PI;
    return r;
#endif
  }

  TrialWaveFunction::GradType TrialWaveFunction::evalGrad(ParticleSet& P,int iat)
  {
    //TAU_PROFILE("TrialWaveFunction::evalGrad","(ParticleSet& P,int iat)", TAU_USER);
    GradType grad_iat;
    for (int i=0; i<Z.size(); ++i) grad_iat += Z[i]->evalGrad(P,iat);
    return grad_iat;
  }

TrialWaveFunction::GradType TrialWaveFunction::alternateEvalGrad(ParticleSet& P,int iat)
{
  //TAU_PROFILE("TrialWaveFunction::evalGrad","(ParticleSet& P,int iat)", TAU_USER);
  GradType grad_iat;
  for (int i=0; i<Z.size(); ++i) grad_iat += Z[i]->alternateEvalGrad(P,iat);
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
    for (int i=0,ii=1; i<Z.size(); ++i,ii+=2)
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
    if (r<0) PhaseDiff=M_PI;
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
  for (int i=0,ii=1; i<Z.size(); ++i,ii+=2)
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
  if (r<0) PhaseDiff=M_PI;
  //     else PhaseDiff=0.0;
  return r;
  #endif
}

  void   TrialWaveFunction::update(ParticleSet& P,int iat)
  {
    //ready to collect "changes" in the gradients and laplacians by the move
    delta_G=0.0;
    delta_L=0.0;
    for (int i=0; i<Z.size(); i++) Z[i]->update(P,delta_G,delta_L,iat);
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
    for (int i=0, ii=1; i<Z.size(); ++i, ii+=2)
      {
        myTimers[ii]->start();
        r *= Z[i]->ratio(P,iat,dG,dL);
        myTimers[ii]->stop();
      }

#if defined(QMC_COMPLEX)
    return std::exp(evaluateLogAndPhase(r,PhaseDiff));
#else
    if (r<0) PhaseDiff=M_PI;
    //     else PhaseDiff=0.0;
    return r;
#endif
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
    for (int i=0; i<Z.size(); i++) Z[i]->acceptMove(P,iat);
    PhaseValue += PhaseDiff;
    PhaseDiff=0.0;
  }

//void TrialWaveFunction::resizeByWalkers(int nwalkers){
//  for(int i=0; i<Z.size(); i++) Z[i]->resizeByWalkers(nwalkers);
//}

  void TrialWaveFunction::checkInVariables(opt_variables_type& active)
  {
    for (int i=0; i<Z.size(); i++) Z[i]->checkInVariables(active);
  }

  void TrialWaveFunction::checkOutVariables(const opt_variables_type& active)
  {
    for (int i=0; i<Z.size(); i++) Z[i]->checkOutVariables(active);
  }

  void TrialWaveFunction::resetParameters(const opt_variables_type& active)
  {
    for (int i=0; i<Z.size(); i++) Z[i]->resetParameters(active);
  }

  void TrialWaveFunction::reportStatus(ostream& os)
  {
    for (int i=0; i<Z.size(); i++) Z[i]->reportStatus(os);
  }

  TrialWaveFunction::RealType TrialWaveFunction::registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    delta_G.resize(P.getTotalNum());
    delta_L.resize(P.getTotalNum());

    P.G = 0.0;
    P.L = 0.0;

    //save the current position
    BufferCursor=buf.current();

    ValueType logpsi(0.0);
    PhaseValue=0.0;
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (;it!=it_end; ++it)
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

  TrialWaveFunction::RealType TrialWaveFunction::registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf)
  {
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (; it!=it_end; ++it)
      {
        if ((*it)->Optimizable)
          {
            (*it)->registerDataForDerivatives(P,buf);
          }
      }
    return 1.0;
  }

  TrialWaveFunction::RealType TrialWaveFunction::updateBuffer(ParticleSet& P
      , PooledData<RealType>& buf, bool fromscratch)
  {
    //TAU_PROFILE("TrialWaveFunction::updateBuffer","(P,..)", TAU_USER);
    P.G = 0.0;
    P.L = 0.0;

    buf.rewind(BufferCursor);

    ValueType logpsi(0.0);
    PhaseValue=0.0;
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (int ii=1; it!=it_end; ++it,ii+=2)
      {
        myTimers[ii]->start();
        logpsi += (*it)->updateBuffer(P,buf,fromscratch);
        PhaseValue += (*it)->PhaseValue;
        myTimers[ii]->stop();
      }
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

    buf.rewind(BufferCursor);

    //TAU_PROFILE("TrialWaveFunction::copyFromBuffer","(P,..)", TAU_USER);
    for (int i=0; i<Z.size(); i++) Z[i]->copyFromBuffer(P,buf);

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
    buf.rewind(BufferCursor);
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
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
    buf.rewind(BufferCursor);
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    for (; it!=it_end; ++it)
      (*it)->dumpFromBuffer(P,buf);
  }

  TrialWaveFunction::RealType
  TrialWaveFunction::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.rewind(BufferCursor);
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

//TrialWaveFunction::RealType
//TrialWaveFunction::evaluate(ParticleSet& P, PooledData<RealType>& buf) {

//  ValueType psi(1.0);
//  for(int i=0; i<Z.size(); i++) psi *= Z[i]->evaluate(P,buf);
//  buf.put(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
//  buf.put(&(P.L[0]), &(P.L[0])+NumPtcls);
//  return real(psi);
//}
//bool TrialWaveFunction::hasSPOSet(const string& aname) {
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
//TrialWaveFunction::getSPOSet(const string& aname) {
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
      for (int i=0; i<Z.size(); ++i)
        myclone->addOrbital(Z[i]->makeClone(tqp),"dummy");
      for (int i=0; i<myTimers.size(); i++)
        myclone->myTimers[i]->set_name(myTimers[i]->get_name());
      myclone->OneOverM=OneOverM;
      return myclone;
    }

  void TrialWaveFunction::evaluateDerivatives(ParticleSet& P,
      const opt_variables_type& optvars,
      vector<RealType>& dlogpsi,
      vector<RealType>& dhpsioverpsi)
  {
    //     // First, zero out derivatives
    //  This should only be done for some variables.
    //     for (int j=0; j<dlogpsi.size(); j++) 
    //       dlogpsi[j] = dhpsioverpsi[j] = 0.0;
      
    for (int i=0; i<Z.size(); i++) {
      if (Z[i]->dPsi)(Z[i]->dPsi)->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
      else Z[i]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
    }
    //orbitals do not know about mass of particle.
    for (int i=0;i<dhpsioverpsi.size();i++) 
      dhpsioverpsi[i]*=OneOverM;
  }

  TrialWaveFunction::RealType
  TrialWaveFunction::KECorrection() const
    {
      RealType sum = 0.0;
      for (int i=0; i<Z.size(); ++i)
        sum += Z[i]->KECorrection();
      return sum;
    }

  void TrialWaveFunction::get_ratios(ParticleSet& P, vector<ValueType>& ratios)
  {
    std::fill(ratios.begin(),ratios.end(),1.0);
    vector<ValueType> t(ratios.size());
    for (int i=0; i<Z.size(); ++i)
      {
        Z[i]->get_ratios(P,t);
        for (int j=0; j<t.size(); ++j) ratios[j]*=t[j];
      }
  }

}
/***************************************************************************
* $RCSfile$   $Author$
* $Revision$   $Date$
* $Id$
***************************************************************************/
