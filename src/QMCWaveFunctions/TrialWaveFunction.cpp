//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/Walker.h"
#include "Particle/WalkerSetRef.h"
#include "Utilities/OhmmsInfo.h"

namespace ohmmsqmc {

  TrialWaveFunction::TrialWaveFunction():SignValue(1.0),LogValue(0.0){ }

  /**@warning Have not decided whether Z is cleaned up by TrialWaveFunction 
   *  or not. It will depend on I/O implementation.
   */
  TrialWaveFunction::~TrialWaveFunction(){
    DEBUGMSG("TrialWaveFunction::~TrialWaveFunction")
      for(int i=0; i<Z.size(); i++) delete Z[i];
  }
  
  /**@param aterm a many-body wavefunction */
  void 
  TrialWaveFunction::add(OrbitalBase* aterm) {
    Z.push_back(aterm);
  }
  
  /** evaluate the log value of a many-body wave function
   *@param P input configuration containing N particles
   *@param all select the wave functions
   *@return the value of log( PI_i psi_i)  many-body wave function
   *
   * @if all == true
   *  all the wave functions evaluate the log value
   * @else
   *  only the wave functions whose Optimizable is set to zero.
   *
   * Upon return, the gradient and laplacian operators are added by the components.
   * Each OrbitalBase evaluates SignValue and LogValue = log(abs(psi_i))
   * Jastrow functions always have SignValue=1.
   */
  TrialWaveFunction::ValueType 
  TrialWaveFunction::evaluateLog(ParticleSet& P, bool all) {
    P.G = 0.0;
    P.L = 0.0;
    //ValueType logpsi=0.0;
    LogValue=0.0;
    SignValue=1.0;
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    while(it != it_end) {
      if(all || (*it)->Optimizable) {
        LogValue += (*it)->evaluateLog(P, P.G, P.L); 
        SignValue *= (*it)->SignValue;
      }
      ++it;
    }
    return LogValue;
  }

  TrialWaveFunction::ValueType 
  TrialWaveFunction::evaluateLog(ParticleSet& P) {
    P.G = 0.0;
    P.L = 0.0;
    //ValueType logpsi=0.0;
    LogValue=0.0;
    SignValue=1.0;
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    while(it != it_end) {
      LogValue += (*it)->evaluateLog(P, P.G, P.L); 
      SignValue *= (*it)->SignValue;
      ++it;
    }
    return LogValue;
  }
  /** evalaute the sum of log value of optimizable many-body wavefunctions
   * @param P  input configuration containing N particles
   * @param fixedG gradients of log(psi) of the fixed wave functions
   * @param fixedL laplacians of log(psi) of the fixed wave functions
   *
   * This function is introduced for optimization only.
   * fixedG and fixedL save the terms coming from the wave functions
   * that are invarient during optimizations.
   * It is expected that evaluateLog(P,false) is called later
   * and the external object adds the varying G and L and the fixed terms.
   */
  TrialWaveFunction::ValueType 
  TrialWaveFunction::evaluateLog(ParticleSet& P,
        ParticleSet::ParticleGradient_t& fixedG,
        ParticleSet::ParticleLaplacian_t& fixedL) {
    P.G = 0.0;
    P.L = 0.0;
    fixedG = 0.0;
    fixedL = 0.0;
    ValueType logpsi(0.0),t(0.0);
    vector<OrbitalBase*>::iterator it(Z.begin());
    vector<OrbitalBase*>::iterator it_end(Z.end());
    while(it != it_end) {
      if((*it)->Optimizable) 
        logpsi += (*it)->evaluateLog(P, P.G, P.L); 
      else
        t += (*it)->evaluateLog(P, fixedG, fixedL); 
      ++it;
    }
    P.G += fixedG; P.L += fixedL;
    return logpsi;
  }
  
  /** evaluate the value of a many-body wave function
   *@param P input configuration containing N particles
   *@return the value of many-body wave function
   *
   *Upon return, the gradient and laplacian operators are added by the components.
  */
  TrialWaveFunction::ValueType 
  TrialWaveFunction::evaluate(ParticleSet& P) {
    P.G = 0.0;
    P.L = 0.0;
    ValueType psi = 1.0;
    for(int i=0; i<Z.size(); i++) {
      psi *= Z[i]->evaluate(P, P.G, P.L);
    }
    //for(int iat=0; iat<P.getTotalNum(); iat++)
    // cout << P.G[iat] << " " << P.L[iat] << endl;
    return psi;
  }
  
  /**
     @param W the input set of walkers
     @param psi a array containing Nw wave function values of each walker
     @brief Upon return, the gradient and laplacian operators are summed 
     by the components. 
  */
  void 
  TrialWaveFunction::evaluate(WalkerSetRef& W, 
			      OrbitalBase::ValueVectorType& psi)
  {
    W.G = 0.0;
    W.L = 0.0;
    psi = 1.0;
    for(int i=0; i<Z.size(); i++) {
      Z[i]->evaluate(W,psi,W.G,W.L);
    }
    //for(int iw=0; iw<psi.size(); iw++) W.Properties(iw,Sign) = psi[iw];
  }
  
  TrialWaveFunction::ValueType
  TrialWaveFunction::ratio(ParticleSet& P,int iat) {
    RealType r=1.0;
    for(int i=0; i<Z.size(); i++) r *= Z[i]->ratio(P,iat);
    return r;
  }
  void   
  TrialWaveFunction::update(ParticleSet& P,int iat) {
    //ready to collect "changes" in the gradients and laplacians by the move
    delta_G=0.0; delta_L=0.0;
    for(int i=0; i<Z.size(); i++) Z[i]->update(P,delta_G,delta_L,iat);
    P.G += delta_G;
    P.L += delta_L;
  }

  
  TrialWaveFunction::ValueType 
  TrialWaveFunction::ratio(ParticleSet& P, int iat, 
			   ParticleSet::ParticleGradient_t& dG,
			   ParticleSet::ParticleLaplacian_t& dL) {
    dG = 0.0;
    dL = 0.0;
    RealType r=1.0;
    for(int i=0; i<Z.size(); i++) r *= Z[i]->ratio(P,iat,dG,dL);
    return r;
  }

  void 
  TrialWaveFunction::restore(int iat) {
    for(int i=0; i<Z.size(); i++) Z[i]->restore(iat);
  }

  void   
  TrialWaveFunction::update2(ParticleSet& P,int iat) {
    for(int i=0; i<Z.size(); i++) Z[i]->update(P,iat);
  }

  void TrialWaveFunction::resizeByWalkers(int nwalkers){
    for(int i=0; i<Z.size(); i++) Z[i]->resizeByWalkers(nwalkers);
  }
  
  void TrialWaveFunction::reset(){
    for(int i=0; i<Z.size(); i++) Z[i]->reset();
  }
  
  void TrialWaveFunction::registerData(ParticleSet& P, PooledData<RealType>& buf) {
    delta_G.resize(P.getTotalNum());
    delta_L.resize(P.getTotalNum());
    P.G = 0.0;
    P.L = 0.0;

    for(int i=0; i<Z.size(); i++) Z[i]->registerData(P,buf);

    //append current gradients and laplacians to the buffer
    TotalDim = OHMMS_DIM*P.getTotalNum();
    buf.add(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
    buf.add(P.L.begin(), P.L.end());

//     cout << "Registering gradients and laplacians " << endl;
//     for(int i=0; i<P.getLocalNum(); i++) {
//       cout << P.G[i] << " " << P.L[i] << endl;
//     }
  }

  void TrialWaveFunction::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

    for(int i=0; i<Z.size(); i++) Z[i]->copyFromBuffer(P,buf);

    //get the gradients and laplacians from the buffer
    buf.get(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
    buf.get(P.L.begin(), P.L.end());

//     cout << "Checking out gradients and laplacians " << endl;
//     for(int i=0; i<P.getLocalNum(); i++) {
//       cout << P.G[i] << " " << P.L[i] << endl;
//     }
  }

  TrialWaveFunction::ValueType
  TrialWaveFunction::evaluate(ParticleSet& P, PooledData<RealType>& buf) {

    ValueType psi = 1.0;
    for(int i=0; i<Z.size(); i++) psi *= Z[i]->evaluate(P,buf);
    buf.put(&(P.G[0][0]), &(P.G[0][0])+TotalDim);
    buf.put(P.L.begin(), P.L.end());

//     cout << "Checking in gradients and laplacians " << endl;
//     for(int i=0; i<P.getLocalNum(); i++) {
//       cout << P.G[i] << " " << P.L[i] << endl;
//     }
    return psi;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

