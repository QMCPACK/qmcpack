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
#include "QMC/ReptationMC.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMC/MolecuFixedNodeBranch.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc {

  ReptationMC::ReptationMC(MCWalkerConfiguration& w, 
			   TrialWaveFunction& psi, 
			   QMCHamiltonian& h, 
			   xmlNodePtr q): 
    QMCDriver(w,psi,h,q), MoveHead(true),NumCut(1),UseBounce(true){ 
    RootName = "rmc";
    QMCType ="rmc";
  }
  
  /** Advance the walkers nblocks*nsteps timesteps. 
   * @param nblocks number of blocks
   * @param nsteps number of steps
   * @param tau the timestep
   *
   * For each timestep:
   * <ul>
   * <li> Move all the particles of a walker.
   * <li> Calculate the properties for the new walker configuration.
   * <li> Accept/reject the new configuration.
   * <li> Accumulate the estimators.
   * <li> Update the trial energy \f$ E_T. \f$
   * <li> Branch the population of walkers (birth/death algorithm).
   * </ul>
   * For each block:
   * <ul>
   * <li> Flush the estimators and print to file.
   * <li> Update the estimate of the local energy.
   * <li> (Optional) Print the ensemble of walker configurations.
   * </ul>
   * Default mode: Print the ensemble of walker configurations 
   * at the end of the run.
   */
  bool ReptationMC::run() { 

    //create a distance table for one walker
    DistanceTable::create(1);
    deltaR.resize(W.getTotalNum());
    drift.resize(W.getTotalNum());
    
    if(put(qmc_node)){
      
      //set the data members to start a new run
      //    getReady();
      int PopIndex, E_TIndex;
      Estimators.resetReportSettings(RootName);
      AcceptIndex = Estimators.addColumn("AcceptRatio");
      Estimators.reportHeader();
      
      for(MCWalkerConfiguration::iterator it = W.begin(); 
	  it != W.end(); ++it) {
	(*it)->Properties(WEIGHT) = 1.0;
	(*it)->Properties(MULTIPLICITY) = 1.0;
      }

      //construct a polymer
      Polymer.insert(Polymer.end(),W.begin(),W.end());      

      if(WalkerRepository.empty()) {
	for(int i=0; i<5; i++)//Maximum beads to be cut=5 (just checking)
	  WalkerRepository.push_back(new Walker_t(*W.begin()));
      }

      IndexType block = 0;
      Pooma::Clock timer;
      IndexType accstep=0;
      IndexType nAcceptTot = 0;
      IndexType nRejectTot = 0;
      
      //accumulate configuration: probably need to reorder
      HDFWalkerOutput WO(RootName);
      do {

	IndexType step = 0;
	timer.start();
	do {
          movePolymer();
	  step++; accstep++;
	  //how to accumulate
	  Estimators.accumulate(W);
	} while(step<nSteps);
	timer.stop();
	
	nAcceptTot += nAccept;
	nRejectTot += nReject;

	Estimators.flush();
	Estimators.setColumn(AcceptIndex,
			     static_cast<double>(nAccept)/static_cast<double>(nAccept+nReject));
	Estimators.report(accstep);


        //change NumBeads to make accstep ~ 50%

	LogOut->getStream() << "Block " << block << " " << timer.cpu_time()
			    << " " << Population << endl;
	
	nAccept = 0; nReject = 0;
	block++;

	if(pStride) WO.get(W);

      } while(block<nBlocks);
      
      LogOut->getStream() 
	<< "ratio = " 
	<< static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
	<< endl;

      Estimators.finalize();
      return true;
    } else 
      return false;
  }

  bool 
  ReptationMC::put(xmlNodePtr q){
    xmlNodePtr qsave=q;
    bool success = putQMCInfo(q);
    success = Estimators.put(qsave);
    return success;
  }
  
  void 
  ReptationMC::movePolymer(){
    
    //Pooma::Clock timer;
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    typedef MCWalkerConfiguration::PropertyContainer_t PropertyContainer_t;
    typedef MCWalkerConfiguration::Walker_t Walker_t;

    Walker_t* anchor=0;
    vector<Walker_t*> heads(WalkerRepository.begin(), WalkerRepository.begin()+NumBeads);
    vector<Walker_t*> tails(NumBeads+1);

    //copy NumBeads+1 to make tails
    if(MoveHead) {
      anchor=Polymer.front();
      for(int i=0,j=Polymer.size()-1; i<NumBeads+1; i++,j--) {
	tails[i] = Polymer[j];
      }
    } else {
      anchor=Polymer.back();
      for(int i=0; i<NumBeads+1; i++) {
	tails[i] = Polymer[i];
      }
    }
    
    RealType Wpolymer=0.0;

    //save the local energies of the anchor and tails
    RealType eloc_x=anchor->Properties(LOCALENERGY);
    RealType eloc_xp=tails[0]->Properties(LOCALENERGY);
    RealType eloc_yp=tails[1]->Properties(LOCALENERGY);

    for(int i=0; i<NumBeads; i++) {

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      W.R = g*deltaR + anchor->R + anchor->Drift;

      //update the distance table associated with W
      DistanceTable::update(W);
      //evaluate wave function
      ValueType psi = Psi.evaluate(W);

      //update the properties of the proposed head
      heads[i]->Properties(LOCALENERGY) = H.evaluate(W);
      H.copy(heads[i]->getEnergyBase());
      heads[i]->R = W.R;
      heads[i]->Drift = tau*W.G;

      //\f${x-y-\tau\nabla \ln \Psi_{T}(y))\f$
      //deltaR = anchor->R - W.R - heads[i]->Drift;
      //Gdrift *= exp(-oneover2tau*Dot(deltaR,deltaR));
      RealType eloc_y=heads[i]->Properties(LOCALENERGY);
      Wpolymer += eloc_x+eloc_y-eloc_xp-eloc_yp;

      //move the anchor and swap the local energies for Wpolyer
      anchor=heads[i];
      eloc_x  = eloc_y;
      eloc_xp = eloc_yp;
      eloc_yp = tails[i+1]->Properties(LOCALENERGY);
    }

    if(Random() < exp(-oneover2tau*Wpolyer)){//move accepted
      if(MoveHead){
	for(int i=0; i<NumBeads; i++) {
	  Polymer.push_front(heads[i]);
	  Polymer.pop_back();
	}
      }else {
	for(int i=0; i<NumBeads; i++) {
	  Polymer.push_back(heads[i]);
	  Polymer.pop_front();
	}
      }
      //copy NumBeads of the tails to the repository for next step
      for(int i=0; i<NumBeads; i++) {WalkerRepository[i] = tails[i];}
      ++nAccept;
    } else {
      ++nReject; 
    }

    RealType Bounce = UseBounce ? 1.0-Wpolymer: 0.5;
    if(Random()<Bounce) {
      MoveHead = ~MoveHead; //flip the direction
    }
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
