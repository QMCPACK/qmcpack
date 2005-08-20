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
      do {
        Population = W.getActiveWalkers();
        advanceWalkerByWalker(*branchEngine);
        step++; accstep++;
        Estimators->accumulate(W);
        Eest = branchEngine->update(Population);
        branchEngine->branch(accstep,W);
      } while(step<nSteps);
      timer.stop();
      
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      Estimators->flush();
      
      Estimators->setColumn(PopIndex,static_cast<RealType>(Population));
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
   * 
   Propose a move for each walker from its old 
   position \f${\bf R'}\f$ to a new position \f${\bf R}\f$ 
   \f[ 
   {\bf R'} + {\bf \chi} + 
   \tau {\bf v_{drift}}({\bf R'}) =  {\bf R},
   \f]
   where \f$ {\bf \chi} \f$ is a 3N-diminsional 
   gaussian of mean zero and variance \f$ \tau \f$
   and \f$ {\bf v_{drift}} \f$ is the drift velocity
   \f[
   {\bf v_{drift}}({\bf R'}) = {\bf \nabla} 
   \ln |\Psi_T({\bf R'})| = \Psi_T({\bf R'})^{-1} 
   {\bf \nabla} \Psi_T({\bf R'}). 
   \f]
   For DMC it is necessary to check if the walker 
   crossed the nodal surface, if this is the case 
   then reject the move, otherwise Metropolis 
   accept/reject with probability
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
   Assign a weight and multiplicity for each walker
   \f[ weight = \exp \left[-\tau(E_L(\mathbf{R})+
   E_L(\mathbf{R})-2E_T)/2 \right]. \f]
   \f[ multiplicity = \exp \left[-\tau(E_L(\mathbf{R})+
   E_L(\mathbf{R})-2E_T)/2 \right] + \nu, \f]
   where \f$ \nu \f$ is a uniform random number.
   *
   Due to the fact that the drift velocity diverges on the nodal
   surface of the trial function \f$ \Psi_T \f$, it is possible
   for walkers close to the nodes to make excessively large proposed
   moves \f$ {\bf R'} \longrightarrow {\bf R} \f$.  With the
   accept/reject step this can lead to persistent configurations;
   a remedy is to impose a cutoff on the magnitude of the drift
   velocity.  We use the smooth cutoff proposed by Umrigar, 
   Nightingale and Runge [J. Chem. Phys.,  99, 2865, (1993)]
   \f[
   {\bf \bar{v}_{drift}} = \frac{-1+\sqrt{1+2 \tau v^2_{drift}}}
   {\tau v^2_{drift}}{\bf v_{drift}},
   \f]
   where \f$ {\bf v_{drift}} \f$ is evaluated at 
   \f$ {\bf R'} \f$ and the magnitude of the drift
   \f$ \tau {\bf v_{drift}} \f$ is unchanged for small
   \f$ \tau v^2_{drift} \f$ and is limited to \f$ \sqrt{2\tau} \f$
   for large \f$ \tau v^2_{drift} \f$. 
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
      //if(Branch(W.Properties(SIGN),(*it)->Properties(SIGN))) {
      if(Branch(signold,thisWalker.Properties(SIGN))) {
        accepted=false;     
        thisWalker.willDie();
      }

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
