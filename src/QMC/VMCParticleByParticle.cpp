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
#include "QMC/VMCParticleByParticle.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"

namespace ohmmsqmc { 

  /// Constructor.
  VMCParticleByParticle::VMCParticleByParticle(MCWalkerConfiguration& w, 
					       TrialWaveFunction& psi, 
					       QMCHamiltonian& h, 
					       xmlNodePtr q): 
    QMCDriver(w,psi,h,q) { 
    RootName = "vmc";
    QMCType ="vmc";
  }
  
  bool VMCParticleByParticle::run() { 

    W.setUpdateMode(MCWalkerConfiguration::Update_Particle);
    //DistanceTable::create(1);

    if(put(qmc_node)){
      
      //set the data members to start a new run
      getReady();
      

      //going to add routines to calculate how much we need
      bool require_register =  W.createAuxDataSet();
      int iwalker=0;
      MCWalkerConfiguration::iterator it = W.begin();
      if(require_register) {
        while(it != W.end()) {
	  W.DataSet[iwalker]->rewind();
  	  W.registerData(**it,*(W.DataSet[iwalker]));
  	  Psi.registerData(W,*(W.DataSet[iwalker]));
	  it++;iwalker++;
        } 
      }      

      Estimators.reset();

      //create an output engine
      HDFWalkerOutput WO(RootName);
      
      IndexType block = 0;
      
      Pooma::Clock timer;

      RealType oneovertau = 1.0/Tau;
      RealType oneover2tau = 0.5*oneovertau;
      RealType g = sqrt(Tau);
      
      MCWalkerConfiguration::PropertyContainer_t Properties;
      ParticleSet::ParticlePos_t deltaR(W.getTotalNum());
      ParticleSet::ParticleGradient_t G(W.getTotalNum());
      ParticleSet::ParticleLaplacian_t L(W.getTotalNum());

      IndexType accstep=0;
      IndexType nAcceptTot = 0;
      IndexType nRejectTot = 0;
      ofstream fout("test.txt");
      do {
	IndexType step = 0;
	timer.start();
	nAccept = 0; nReject=0;
	IndexType nAllRejected = 0;
	do {
	  it = W.begin();	 
	  iwalker=0; 
	  while(it != W.end()) {

	    MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
	    W.R = (*it)->R;
	    w_buffer.rewind();
	    W.copyFromBuffer(w_buffer);
	    Psi.copyFromBuffer(W,w_buffer);

            ValueType psi_old = (*it)->Properties(PSI);
	    ValueType psi = psi_old;
	    //create a 3N-Dimensional Gaussian with variance=1
	    makeGaussRandom(deltaR);
	    bool moved = false;
	    for(int iat=0; iat<W.getTotalNum(); iat++) {
	      //PosType dr = g*deltaR[iat]+(*it)->Drift[iat];
	      PosType dr = g*deltaR[iat];
	      W.makeMove(iat,dr);
	      RealType ratio = Psi.ratio(W,iat);
	      RealType ratio2 = pow(ratio,2);
	      if(Random() < ratio2) {
		moved = true;
		++nAccept;
		W.acceptMove(iat);
		Psi.update(W,iat);
	      } else {
		++nReject; 
	      }
	    }

	    if(moved) {
 	      //(*it)->Data.rewind();
	      //W.copyToBuffer((*it)->Data);
 	      //psi = Psi.evaluate(W,(*it)->Data);
	      w_buffer.rewind();
	      W.copyToBuffer(w_buffer);
	      psi = Psi.evaluate(W,w_buffer);

	      (*it)->R = W.R;
	      (*it)->Properties(PSISQ) = psi*psi;
	      (*it)->Properties(PSI) = psi;
	      (*it)->Properties(LOCALENERGY) = H.evaluate(W);
	      H.copy((*it)->getEnergyBase());
	      //GreenFunction is not used here
	      //ValueType vsq = Dot(W.G,W.G);
	      //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
	      //(*it)->Drift = scale*W.G;
	    }
	    //Keep until everything is tested: debug routine
	    // 	    if(moved) {
	    // 	      //(*it)->Data.rewind();
	    // 	      //W.copyToBuffer((*it)->Data);
	    //   	      //psi = Psi.evaluate(W,(*it)->Data);
	    // 	      DataSet[iwalker]->rewind();
	    // 	      W.copyToBuffer(*DataSet[iwalker]);
	    // 	      psi = Psi.evaluate(W,*DataSet[iwalker]);
	    
	    // 	      cout << "What is the ratio " << psi/psi_old << endl;
	    // 	      // 		 std::cout << "Are they the same ratio=" << psi << endl;
	    
	    // 	      RealType please =  H.evaluate(W);

	    // 	      /**@note Debug the particle-by-particle move */
	    // 	      G = W.G;
	    // 	      L = W.L;
	    // 	      DistanceTable::update(W);
	    // 	      ValueType psinew = Psi.evaluate(W);
	    // 	      ValueType dsum=pow(psi-psinew,2);
	    // 	      std::cout<< "Diff of updated and calculated gradients/laplacians" << endl;
	    // 	      for(int iat=0; iat<W.getTotalNum(); iat++) {
	    // 		dsum += pow(G[iat][0]-W.G[iat][0],2)+pow(G[iat][1]-W.G[iat][1],2)
	    // 		  +pow(G[iat][2]-W.G[iat][2],2)+pow(L[iat]-W.L[iat],2);
	    // 		std::cout << G[iat]-W.G[iat] << " " << L[iat]-W.L[iat] <<endl;
	    // 	      }
	    // 	      std::cout << "Are they the same ratio=" << psi << " eval=" << psinew << " " << sqrt(dsum) << endl;
	    // 	      std::cout << "============================ " << endl;
	    // 	      for(int i=0; i<W.getLocalNum(); i++) {
	    // 		cout << W.G[i] << " " << W.L[i] << endl;
	    // 	      }
	    
	    // 	      (*it)->Properties(PSISQ) = psi*psi;
	    // 	      (*it)->Properties(PSI) = psi;
	    // 	      (*it)->Properties(LOCALENERGY) = H.evaluate(W);
	    // 	      std::cout << "Energy " << please << " " << (*it)->Properties(LOCALENERGY) 
	    // 			<< endl;
	    // 	      H.copy((*it)->getEnergyBase());
	    // 	      ValueType vsq = Dot(W.G,W.G);
	    // 	      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
	    // 	      (*it)->Drift = scale*W.G;
	    // 	      (*it)->R =W.R;
	    //	    } 
	    else {
	      //WARNMSG("All the particle moves are rejected.")
	      nAllRejected++;
	    }
	    it++; iwalker++;
	  }
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
	
	LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << " Fixed_configs " 
			    << static_cast<double>(nAllRejected)/static_cast<double>(step*W.getActiveWalkers()) << endl;
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
  VMCParticleByParticle::put(xmlNodePtr q){
    xmlNodePtr qsave=q;
    bool success = putQMCInfo(q);
    success = Estimators.put(qsave);
    return success;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
