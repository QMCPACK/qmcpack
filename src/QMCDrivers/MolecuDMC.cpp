//////////////////////////////////////////////////////////////////
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
#include "QMCDrivers/MolecuDMC.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc {

  MolecuDMC::MolecuDMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h),BranchInfo("default"),branchEngine(0){ 
    RootName = "dmc";
    QMCType ="dmc";
  }

  MolecuDMC::~MolecuDMC() {
    if(branchEngine) delete branchEngine;
  }

  void MolecuDMC::setBranchInfo(const string& afile) {
    BranchInfo=afile;
  }
  
  bool MolecuDMC::put(xmlNodePtr cur){
    return true;
  }
  
  /** Advance the walkers nblocks*nsteps timesteps. 
   *
   * For each block:
   * <ul>
   *  <li> Advance walkers for nsteps
   *  For each timestep:
   *   <ul>
   *   <li> Move all the particles of a walker.
   *   <li> Calculate the properties for the new walker configuration.
   *   <li> Accept/reject the new configuration.
   *   <li> Accumulate the estimators.
   *   <li> Update the trial energy \f$ E_T \f$
   *   <li> Branch the population of walkers (birth/death algorithm).
   *   </ul>
   * <li> Flush the estimators and print to file.
   * <li> Update the estimate of the local energy.
   * <li> (Optional) Print the ensemble of walker configurations.
   * </ul>
   * Default mode: Print the ensemble of walker configurations 
   * at the end of the run.
   */
  bool MolecuDMC::run() { 

    //add columns
    IndexType PopIndex = Estimators->addColumn("Population");
    IndexType EtrialIndex = Estimators->addColumn("Etrial");
    //write the header
    Estimators->reportHeader();

    if(branchEngine == 0) {
      branchEngine=new BranchEngineType(Tau,W.getActiveWalkers());
      RealType e_ref = W.getLocalEnergy();
      branchEngine->setEguess(e_ref);
      branchEngine->put(qmcNode,LogOut);
    }

    branchEngine->flush(0);

    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    while(it != it_end) {
      (*it)->Weight= 1.0;
      (*it)->Multiplicity=1;
      ++it;
    }
    
    IndexType block = 0;
    Pooma::Clock timer;
    int Population = W.getActiveWalkers();
    int tPopulation = W.getActiveWalkers();
    RealType Eest = branchEngine->E_T;
    IndexType accstep=0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    do {
      IndexType step = 0;
      timer.start();
      IndexType pop_acc=0.0; 
      do {
        pop_acc += W.getActiveWalkers();
        advanceWalkerByWalker(*branchEngine);
        step++; accstep++;
        Estimators->accumulate(W);
        Eest = branchEngine->update(W.getActiveWalkers(), Eest);
        branchEngine->branch(accstep,W);
      } while(step<nSteps);
      timer.stop();
      
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      Estimators->flush();
      
      Estimators->setColumn(PopIndex,static_cast<RealType>(pop_acc/static_cast<RealType>(nSteps)));
      Estimators->setColumn(EtrialIndex,Eest);
      Estimators->setColumn(AcceptIndex,
      	            static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
      Estimators->report(accstep);
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time()
      		    << " " << Population << endl;
      Eest = Estimators->average(0);
      nAccept = 0; nReject = 0;
      block++;
      if(pStride) {
        //create an output engine: could accumulate the configurations
        HDFWalkerOutput WO(RootName);
        WO.get(W);
      }
      W.reset();
    } while(block<nBlocks);
    
    LogOut->getStream() 
      << "ratio = " << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << endl;
    
    if(!pStride) {
      //create an output engine: could accumulate the configurations
      HDFWalkerOutput WO(RootName);
      WO.get(W);
    }

    Estimators->finalize();
    return true;
  }

  /**  Advance all the walkers one timstep. 
   * @param Branch class that controls the trial energy and branching
   */
  template<class BRANCHER>
  void 
  MolecuDMC::advanceWalkerByWalker(BRANCHER& Branch) {
    

    //Pooma::Clock timer;
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);

    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    while(it != it_end) {
      
      Walker_t& thisWalker(**it);
      thisWalker.Weight= 1.0;
      thisWalker.Multiplicity=1;
      
      //save old local energy
      RealType eold    = thisWalker.Properties(LOCALENERGY);
      RealType signold = thisWalker.Properties(SIGN);
      RealType emixed  = eold;

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      
      W.R = g*deltaR + thisWalker.R + thisWalker.Drift;
      
      //update the distance table associated with W
      DistanceTable::update(W);
      
      //evaluate wave function
      ValueType logpsi(Psi.evaluateLog(W));

      bool accepted=false; 
      if(Branch(Psi.getSign(),thisWalker.Properties(SIGN))) {
        thisWalker.Age++;
        thisWalker.willDie();
      } else {
        RealType enew(H.evaluate(W));
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        ValueType vsq = Dot(W.G,W.G);
        //converting gradients to drifts, D = tau*G (reuse G)
        ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
        drift = scale*W.G;
        deltaR = (*it)->R - W.R - drift;
        RealType logGb = -oneover2tau*Dot(deltaR,deltaR);

        //set acceptance probability
        //RealType prob= std::min(exp(logGb-logGf +2.0*(W.Properties(LOGPSI)-thisWalker.Properties(LOGPSI))),1.0);
        RealType prob= std::min(exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
        if(Random() > prob){
          thisWalker.Age++;
          eold=enew;
        } else {
          accepted=true;  
          thisWalker.R = W.R;
          thisWalker.Drift = drift;
          thisWalker.resetProperty(logpsi,Psi.getSign(),enew);
          H.saveProperty(thisWalker.getPropertyBase());
          emixed = (emixed+enew)*0.5;
        }

        //calculate the weight and multiplicity
        ValueType M = Branch.branchGF(Tau,emixed,0.0); //1.0-prob);
        if(thisWalker.Age > 3) M = min(0.5,M);
        else if(thisWalker.Age > 0) M = min(1.0,M);
        thisWalker.Weight = M; 
        thisWalker.Multiplicity = M + Random();
      }

      Branch.accumulate(eold,thisWalker.Weight);

      /*

      RealType enew(H.evaluate(W));

      //deltaR = W.R - (*it)->R - (*it)->Drift;
      RealType logGf = -0.5*Dot(deltaR,deltaR);

      //scale the drift term to prevent persistent cofigurations
      ValueType vsq = Dot(W.G,W.G);

      //converting gradients to drifts, D = tau*G (reuse G)
      //   W.G *= Tau;//original implementation with bare drift
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      drift = scale*W.G;
      deltaR = (*it)->R - W.R - drift;

      //RealType backwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
      RealType logGb = -oneover2tau*Dot(deltaR,deltaR);

      //set acceptance probability
      //RealType prob= std::min(exp(logGb-logGf +2.0*(W.Properties(LOGPSI)-thisWalker.Properties(LOGPSI))),1.0);
      RealType prob= std::min(exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);

      if(Random() > prob){
        thisWalker.Age++;
      } else {
        accepted=true;  
        thisWalker.R = W.R;
        thisWalker.Drift = drift;
        thisWalker.resetProperty(logpsi,Psi.getSign(),enew);
        H.saveProperty(thisWalker.getPropertyBase());
        emixed = (emixed+enew)*0.5;
      }
      
      //calculate the weight and multiplicity
      ValueType M = Branch.branchGF(Tau,emixed,0.0); //1.0-prob);
      if(thisWalker.Age > 3) M = min(0.5,M);
      else if(thisWalker.Age > 0) M = min(1.0,M);
      thisWalker.Weight = M; 
      thisWalker.Multiplicity = M + Random();
      
      Branch.accumulate(emixed,M);

      //node-crossing: kill it for the time being
      if(Branch(signold,thisWalker.Properties(SIGN))) {
        accepted=false;     
        thisWalker.willDie();
      }
      */

      if(accepted) 
        ++nAccept;
      else 
        ++nReject;
      ++it;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
