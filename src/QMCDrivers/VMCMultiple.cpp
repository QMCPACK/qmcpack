//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
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
#include "QMCDrivers/VMCMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"
#include "Estimators/MultipleEnergyEstimator.h"

namespace ohmmsqmc { 

  /// Constructor.
  VMCMultiple::VMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), multiEstimator(0) { 
    RootName = "vmc";
    QMCType ="vmc-multi";

    //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
    add_H_and_Psi(&h,&psi);
  }

  /** allocate internal data here before run() is called
   * @author SIMONE
   *
   * See QMCDriver::process
   */
  bool VMCMultiple::put(xmlNodePtr q){


    nPsi=Psi1.size();	

    logpsi.resize(nPsi);
    sumratio.resize(nPsi);	
    invsumratio.resize(nPsi);	

    for(int ipsi=0; ipsi<nPsi; ipsi++) 
      H1[ipsi]->add2WalkerProperty(W);

    if(Estimators == 0) {
      Estimators = new ScalarEstimatorManager(H);
      multiEstimator = new MultipleEnergyEstimator(H,nPsi);
      Estimators->add(multiEstimator,"elocal");
    }

    LOGMSG("Number of H and Psi " << nPsi)

    H1[0]->setPrimary(true);
    for(int ipsi=1; ipsi<nPsi; ipsi++) {
      H1[ipsi]->setPrimary(false);
    }
    return true;
  }
  
  /** Run the VMCMultiple algorithm.
   *
   * Similar to VMC::run 
   */
  bool VMCMultiple::run() { 

    Estimators->reportHeader();

    bool require_register=false;

    //this is where the first values are evaulated
    multiEstimator->initialize(W,H1,Psi1,Tau,require_register);
    
    Estimators->reset();

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
        advanceWalkerByWalker();
        step++;accstep++;
        Estimators->accumulate(W);
      } while(step<nSteps);

      timer.stop();
      nAcceptTot += nAccept; nRejectTot += nReject;
      Estimators->flush();
      RealType TotalConfig=static_cast<RealType>(nAccept+nReject);		
      Estimators->setColumn(AcceptIndex,nAccept/TotalConfig);		

      Estimators->report(accstep);
      
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
    
    Estimators->finalize();

    return true;
  }

  /**  Advance all the walkers one timstep. 
   */
  void 
  VMCMultiple::advanceWalkerByWalker() {
    
    RealType oneovertau = 1.0/Tau;
    RealType oneover2tau = 0.5*oneovertau;
    RealType g = sqrt(Tau);
    
    //MCWalkerConfiguration::PropertyContainer_t Properties;
    
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    int iwlk(0); 
    int nPsi_minus_one(nPsi-1);

    while(it != it_end) {

      MCWalkerConfiguration::Walker_t &thisWalker(**it);

      //copy the properties of the working walker
      //Properties = thisWalker.Properties;   

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      
      W.R = g*deltaR + thisWalker.R + thisWalker.Drift;
      
      //update the distance table associated with W
      DistanceTable::update(W);
      
      //Evaluate Psi and graidients and laplacians
      //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
      for(int ipsi=0; ipsi< nPsi;ipsi++) {			  
	logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W); 		  
        Psi1[ipsi]->L=W.L; 
        Psi1[ipsi]->G=W.G; 
	sumratio[ipsi]=1.0;					  
      } 							  

      // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i	   
      for(int ipsi=0; ipsi< nPsi_minus_one; ipsi++) {			  
	for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++){     		  
	  RealType ratioij=exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));  
	  sumratio[ipsi] += ratioij;                              
	  sumratio[jpsi] += 1.0/ratioij;			  
	}                                                         
      }                                                           

      for(int ipsi=0; ipsi< nPsi; ipsi++)			  
	invsumratio[ipsi]=1.0/sumratio[ipsi];			  

      // Only these properties need to be updated                 
      // Using the sum of the ratio Psi^2[j]/Psi^2[iwref]	 
      // because these are number of order 1. Potentially	 
      // the sum of Psi^2[j] can get very big			 
      //Properties(LOGPSI) =logpsi[0];		   		 
      //Properties(SUMRATIO) = sumratio[0];			 
 
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      ValueType scale = Tau; // ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);	

      //accumulate the weighted drift
      drift = invsumratio[0]*Psi1[0]->G;
      for(int ipsi=1; ipsi< nPsi ;ipsi++) {               		
        drift += invsumratio[ipsi]*Psi1[ipsi]->G;  
      } 							    	
      drift *= scale; 						   	

      deltaR = thisWalker.R - W.R - drift;
      RealType logGb = -oneover2tau*Dot(deltaR,deltaR);
      
      //Original
      //RealType g = Properties(SUMRATIO)/thisWalker.Properties(SUMRATIO)*   		
      //	exp(logGb-logGf+2.0*(Properties(LOGPSI)-thisWalker.Properties(LOGPSI)));	
      //Reuse Multiplicity to store the sumratio[0]
      RealType g = sumratio[0]/thisWalker.Multiplicity*   		
       	exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));	

      if(Random() > g) {
	thisWalker.Age++;     
	++nReject; 
      } else {
	thisWalker.Age=0;
        thisWalker.Multiplicity=sumratio[0];
	thisWalker.R = W.R;
	thisWalker.Drift = drift;
	for(int ipsi=0; ipsi<nPsi; ipsi++){ 		            
          W.L=Psi1[ipsi]->L; 
          W.G=Psi1[ipsi]->G; 
	  RealType et = H1[ipsi]->evaluate(W);
          thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
          thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getSign();
          thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
          thisWalker.Properties(ipsi,LOCALENERGY)=et;
          //multiEstimator->updateSample(iwlk,ipsi,et,invsumratio[ipsi]); 
          H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
	}

	++nAccept;
      }
      ++it; 
      ++iwlk;
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
