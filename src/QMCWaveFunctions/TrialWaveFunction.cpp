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

  TrialWaveFunction::TrialWaveFunction(){ }

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
      psi *= Z[i]->evaluate(P,P.G,P.L);
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
    cout << "TrialWaveFunction::ratio " << r << endl;
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

  void TrialWaveFunction::resizeByWalkers(int nwalkers){
    for(int i=0; i<Z.size(); i++) Z[i]->resizeByWalkers(nwalkers);
  }
  
  void TrialWaveFunction::reset(){
    for(int i=0; i<Z.size(); i++) Z[i]->reset();
  }
  
  void TrialWaveFunction::registerData(ParticleSet& P, PooledData<RealType>& buf) {
    delta_G.resize(P.getTotalNum());
    delta_L.resize(P.getTotalNum());
    for(int i=0; i<Z.size(); i++) Z[i]->registerData(P,buf);
  }

  void TrialWaveFunction::putData(ParticleSet& P, PooledData<RealType>& buf) {
    P.G=0.0; P.L=0.0;
    for(int i=0; i<Z.size(); i++) Z[i]->putData(P,buf);
  }

  TrialWaveFunction::ValueType
  TrialWaveFunction::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    ValueType psi = 1.0;
    for(int i=0; i<Z.size(); i++) psi *= Z[i]->evaluate(P,buf);
    return psi;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

