
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMC/VMC.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

//This macro is a temporary one to test the efficiency of all-walker moves. So far, not as much
//as to convince JK to implement it globally.
#define MOVE_ONE

  /// Constructor.
  VMC::VMC(MCWalkerConfiguration& w, 
	   TrialWaveFunction& psi, 
	   QMCHamiltonian& h, 
	   xmlNodePtr q): 
    QMCDriver(w,psi,h,q) { 
    RootName = "vmc";
    QMCType ="vmc";
  }
  
  /** Run the VMC algorithm.
   * \param nblocks number of blocks
   * \param nsteps number of steps
   * \param tau the timestep
   *
   * Advance the walkers nblocks*nsteps timesteps. 
   * For each timestep:
   * <ul>
   * <li> Move all the particles of a walker.
   * <li> Calculate the properties for the new walker configuration.
   * <li> Accept/reject the new configuration.
   * <li> Accumulate the estimators.
   * </ul>
   * For each block:
   * <ul>
   * <li> Flush the estimators and print to file.
   * <li> (Optional) Print the ensemble of walker configurations.
   * </ul>
   *
   * Default mode: Print the ensemble of walker configurations 
   * at the end of the run.
   */
  bool VMC::run() { 
    
#ifdef MOVE_ONE
    DistanceTable::create(1);
#else
    DistanceTable::create(W.getActiveWalkers());
    Psi.resizeByWalkers(W.getActiveWalkers());
#endif

    if(put(qmc_node)){
      
      //set the data members to start a new run
      getReady();
      
      //probably unnecessary
      for(MCWalkerConfiguration::iterator it = W.begin(); 
	  it != W.end(); ++it) {
	(*it)->Properties(Weight) = 1.0;
      }
      
      Estimators.reset();

      //create an output engine
      HDFWalkerOutput WO(RootName);
      
      IndexType block = 0;
      
      Pooma::Clock timer;
      
      double wh=0.0;
      IndexType accstep=0;
      IndexType nAcceptTot = 0;
      IndexType nRejectTot = 0;
      do {
	IndexType step = 0;
	timer.start();
	nAccept = 0; nReject=0;
	do {
#ifdef MOVE_ONE
	  advanceWalkerByWalker();
#else
	  advanceAllWalkers();
#endif
	  step++;accstep++;
	  Estimators.accumulate(W);
	} while(step<nSteps);
	
	timer.stop();
	nAcceptTot += nAccept;
	nRejectTot += nReject;
	
	Estimators.flush();
	Estimators.setColumn(AcceptIndex,
			     static_cast<double>(nAccept)/static_cast<double>(nAccept+nReject));
	Estimators.report(accstep);
	
	LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << endl;
	if(pStride) WO.get(W);
	nAccept = 0; nReject = 0;
	block++;
      } while(block<nBlocks);
      
      LogOut->getStream() 
	<< "Ratio = " 
	<< static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
	<< endl;
      
      if(!pStride) WO.get(W);
      
      Estimators.finalize();
      
      return true;
    } else {
      ERRORMSG("Error with Input")
      return false;
    }
  }


  bool 
  VMC::put(xmlNodePtr q){
    xmlNodePtr qsave=q;
    bool success = putQMCInfo(q);
    success = Estimators.put(qsave);
    return success;
  }

  /**  Advance all the walkers one timstep. 
   *
   Propose a move for each walker from its old 
   position \f${\bf R'}\f$ to a new position \f${\bf R}\f$
   \f[ 
   {\bf R'} + {\bf \chi} + 
   \tau {\bf v_{drift}}({\bf R'}) =  {\bf R},
   \f]
   where \f$ {\bf \chi} \f$ is a 3N-diminsional Gaussian 
   of mean zero and variance \f$ \tau \f$ and 
   \f$ {\bf v_{drift}} \f$ is the drift velocity
   \f[ 
   {\bf v_{drift}}({\bf R'}) = {\bf \nabla} 
   \ln |\Psi_T({\bf R'})| = \Psi_T({\bf R'})^{-1} 
   {\bf \nabla} \Psi_T({\bf R'}). 
   \f]
   Metropolis accept/reject with probability
   \f[
   P_{accept}(\mathbf{R'}\rightarrow\mathbf{R}) = 
   \min\left[1,\frac{G(\mathbf{R}\rightarrow\mathbf{R'})
   \Psi_T(\mathbf{R})^2}{G(\mathbf{R'}\rightarrow\mathbf{R})
   \Psi_T(\mathbf{R'})^2}\right],
   \f] 
   where \f$ G \f$ is the drift-diffusion Green's function 
   \f[
   G(\mathbf{R'} \rightarrow 
   \mathbf{R}) = (2\pi\tau)^{-3/2}\exp \left[ -
   (\mathbf{R}-\mathbf{R'}-\tau \mathbf{v_{drift}}
   (\mathbf{R'}))^2/2\tau \right].
   \f]
   If the move is accepted, update the walker configuration and
   properties.  For rejected moves, do not update except for the
   Age which needs to be incremented by one.
   *
   */

  /**@ingroup advanceWalkersVMC
   * \brief Loop through the walkers and advance one at a time.
   */
  
  void 
  VMC::advanceWalkerByWalker() {
    
    //Pooma::Clock timer;
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    MCWalkerConfiguration::PropertyContainer_t Properties;
    static ParticleSet::ParticlePos_t deltaR(W.getTotalNum());
    static ParticleSet::ParticlePos_t drift(W.getTotalNum());
    int nh = H.size()+1;
    
    for (MCWalkerConfiguration::iterator it = W.begin();
	 it != W.end(); ++it) {
      
      //copy the properties of the working walker
      Properties = (*it)->Properties;
      
      //save old local energy
      ValueType eold = Properties(LocalEnergy);
      
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      
      W.R = g*deltaR + (*it)->R + (*it)->Drift;
      
      //update the distance table associated with W
      DistanceTable::update(W);
      
      //evaluate wave function
      ValueType psi = Psi.evaluate(W);
      //update the properties
      Properties(PsiSq) = psi*psi;
      Properties(Sign) = psi;
      Properties(LocalEnergy) = H.evaluate(W);
      
      // deltaR = W.R - (*it)->R - (*it)->Drift;
      //  RealType forwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
      RealType forwardGF = exp(-0.5*Dot(deltaR,deltaR));
      
      //converting gradients to drifts, D = tau*G (reuse G)
      //W.G *= Tau;//original implementation with bare drift
      ValueType vsq = Dot(W.G,W.G);
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      drift = scale*W.G;
      //   W.G *= scale;
      
      deltaR = (*it)->R - W.R - drift;
      RealType backwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
      
      //forwardGF/backwardGF*Properties(PsiSq)/(*it)->Properties(PsiSq)
      //RealType prob   = min(Properties(PsiSq)/(*it)->Properties(PsiSq),1.0);
      //cout << "forward/backward " << forwardGF << " " << backwardGF << endl;
      if(Random() > 
	 backwardGF/forwardGF*Properties(PsiSq)/(*it)->Properties(PsiSq)) {
	(*it)->Properties(Age)++;     
	++nReject; 
      } else {
	Properties(Age) = 0;
	(*it)->R = W.R;
	(*it)->Drift = drift;
	(*it)->Properties = Properties;
	H.get((*it)->E);
	++nAccept;
      }
      
    }
  }
  
  /**  Advance all the walkers simultaneously. 
   */
  
  void VMC::advanceAllWalkers() {
    
    static ParticleSet::ParticlePos_t deltaR(W.getTotalNum());
    
    WalkerSetRef Wref(W);
    Wref.resize(W.getActiveWalkers(),W.getTotalNum());
    
    //Pooma::Clock timer;
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    MCWalkerConfiguration::PropertyContainer_t Properties;
    makeGaussRandom(Wref.R);
    
    Wref.R *= g;
    
    int nptcl = W.getTotalNum();
    int iw = 0;
    MCWalkerConfiguration::iterator it = W.begin();
    while(it !=  W.end()) {
      const ParticleSet::ParticlePos_t& r = (*it)->R;
      for(int jat=0; jat<nptcl; jat++) {
	Wref.R(iw,jat) += r(jat) + (*it)->Drift(jat);
      }
      iw++; it++;
    }
    
    DistanceTable::update(Wref);
    
    OrbitalBase::ValueVectorType   psi(iw), energy(iw);
    
    Psi.evaluate(Wref,psi);
    
    H.evaluate(Wref,energy);
    
    //multiply tau to convert gradient to drift term
    Wref.G *= Tau;
    
    iw = 0;
    it = W.begin();
    while(it !=  W.end()) {
      
      ValueType eold = Properties(LocalEnergy);
      
      for(int iat=0; iat<nptcl; iat++)
	deltaR(iat) = Wref.R(iw,iat) - (*it)->R(iat) - (*it)->Drift(iat);
      RealType forwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
      
      for(int iat=0; iat<nptcl; iat++)
	deltaR(iat) = (*it)->R(iat) - Wref.R(iw,iat) - Wref.G(iw,iat);
      
      RealType backwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
      
      ValueType psisq = psi(iw)*psi(iw);
      if(Random() > 
	 backwardGF/forwardGF*psisq/(*it)->Properties(PsiSq)) {
	++nReject; 
	(*it)->Properties(Age) += 1;
      } else {
	(*it)->Properties(Age) = 0;
	for(int iat=0; iat<nptcl; iat++) (*it)->R(iat) = Wref.R(iw,iat);
	for(int iat=0; iat<nptcl; iat++) (*it)->Drift(iat) = Wref.G(iw,iat);
	(*it)->Properties(Sign) = psi(iw);
	(*it)->Properties(PsiSq) = psisq;
	(*it)->Properties(LocalEnergy) = energy(iw);
	++nAccept;
      }
      iw++;it++;
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
