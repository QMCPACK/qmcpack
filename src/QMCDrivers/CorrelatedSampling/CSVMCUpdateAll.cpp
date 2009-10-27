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
#include "QMCDrivers/CorrelatedSampling/CSVMCUpdateAll.h"
//#include "Utilities/OhmmsInfo.h"
//#include "Particle/MCWalkerConfiguration.h"
//#include "Particle/HDFWalkerIO.h"
//#include "ParticleBase/ParticleUtility.h"
//#include "ParticleBase/RandomSeqGenerator.h"
//#include "ParticleBase/ParticleAttribOps.h"
//#include "Message/Communicate.h"
//#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus { 

  /// Constructor.
  CSVMCUpdateAll::CSVMCUpdateAll(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg):
    CSUpdateBase(w,psi,h,rg) 
    { }

  /**  Advance all the walkers one timstep. 
   */
  void CSVMCUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
  {
    
    int iwlk(0); 
    int nPsi_minus_one(nPsi-1);

    while(it != it_end) 
    {

      MCWalkerConfiguration::Walker_t &thisWalker(**it);

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(deltaR,RandomGen);
      
      if(useDrift)
        W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
      else
        W.R = m_sqrttau*deltaR + thisWalker.R;
      
      //update the distance table associated with W
      //DistanceTable::update(W);
      W.update();
      
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
	  RealType ratioij=avgNorm[ipsi]/avgNorm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));  
	  sumratio[ipsi] += ratioij;                              
	  sumratio[jpsi] += 1.0/ratioij;			  
	}                                                         
      }                                                           

      for(int ipsi=0; ipsi<nPsi; ipsi++)
        invsumratio[ipsi]=1.0/sumratio[ipsi];			  

      RealType g = sumratio[0]/thisWalker.Multiplicity*   		
       	std::exp(2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));	

      if(useDrift)
      {
        //forward green function
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
        for(int ipsi=1; ipsi< nPsi ;ipsi++) {               		
          PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
        } 							    	
        setScaledDrift(Tau,drift);
        //backward green function
        deltaR = thisWalker.R - W.R - drift;
        RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
        g *= std::exp(logGb-logGf);
      }
      //Original
      //RealType g = Properties(SUMRATIO)/thisWalker.Properties(SUMRATIO)*   		
      //	exp(logGb-logGf+2.0*(Properties(LOGPSI)-thisWalker.Properties(LOGPSI)));	
      //Reuse Multiplicity to store the sumratio[0]
      //This is broken up into two pieces
      //RealType g = sumratio[0]/thisWalker.Multiplicity*   		
      // 	std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));	

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
          thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
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
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.cpp 1593 2007-01-04 23:23:27Z jnkim $ 
 ***************************************************************************/
